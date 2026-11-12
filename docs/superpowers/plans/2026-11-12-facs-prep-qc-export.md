# FACS prep / QC loop + final export Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A per-panel prep script saves processed FACS data and renders a QC
report on every iteration; a `final` switch exports that saved data to DuckDB,
dropping and recording the exclusions from the meta workbook.

**Architecture:** New `facs_prep_save()` / `facs_prep_load()` persist a prep
object (`data` + `info`) as `.rds` in `<panel folder>/processed/`.
`facs_match_exclusions()` resolves the meta `exclusions` sheet against that
data (used by the QC report and the export). `report_knit_facs_qc()` renders
a skeleton template into `<panel folder>/qc/`. `db_export_facs()` replaces a
root in one transaction, writes `import_log` + new `import_exclusions` rows.

**Tech Stack:** R, DuckDB (`DBI`/`duckdb`), `rmarkdown`, `testthat` 3e, `withr`.

**Spec:** `docs/superpowers/specs/2026-11-12-facs-prep-qc-export-design.md`
(parent: `2026-11-12-facs-pipeline-workflow-design.md`)

## Global Constraints

- Naming `domain_verb[_qualifier]()`, domains from CLAUDE.md only; files `domain_verb.R`.
- Base pipe `|>`; every external call namespaced `pkg::fn()`; deselect with `dplyr::select(!col)`.
- Non-ASCII in R strings as `\uXXXX`; roxygen text plain ASCII.
- Every export: `@param`, `@returns` (non-empty), `@export`; `@examples` in `\dontrun{}` (need files) or omitted.
- No new `Imports` (all needed packages already listed). `knitr` is only used inside the `.Rmd` template, not in package code.
- `devtools::document()` after roxygen changes; `devtools::check()` 0 errors / 0 warnings, only the NOTEs listed in CLAUDE.md.
- No self-referencing foreign keys (breaks `EXPORT DATABASE` snapshots).
- `source_file` / `qc_report` stored relative to the project's `data/` folder.
- Prep and report contain **all** samples; exclusions are applied only in `db_export_facs()`.
- An applicable exclusion row matching no file: **error** in the export, **warning** in the report.

## Review Focus

1. **Exclusions sheet as `meta_read()` returns it:** `janitor::remove_empty()` drops all-empty columns (e.g. no row uses `file_name` -> column absent), a header-only sheet comes back with 0 rows / 0 columns, and partly-empty columns are logical `NA`. Matching must treat absent columns as all-`NA` and compare as character. Tests in Task 2.
2. **Re-export of an existing root answered "no"** must leave data and log untouched; answered "yes" with the *same* sample keys must succeed (delete + reinsert in one DuckDB transaction). Tests in Task 5.
3. **Stale prep:** `.wsp` saved in FlowJo after the last prep run -> `facs_prep_load()` stops. Test in Task 1 (`Sys.setFileTime`).
4. **Every file excluded** (or none left for the panel) -> export errors, writes nothing. Test in Task 5.
5. **Exclusion row for another panel** is ignored (not "unmatched"); row with empty `panel` applies to all. Tests in Task 2.

---

### Task 1: Prep save / load + file-clash check

**Files:**
- Create: `R/facs_prep.R`
- Test: `tests/testthat/test-facs_prep.R`
- Create: `tests/testthat/helper-prep.R`

**Interfaces:**
- Produces:
  - `facs_check_file_clashes(data) -> tibble` with columns `file_name, mouse_ID, tissue, condition`: every file whose `mouse_ID x tissue x condition` is shared with another file; 0 rows if none.
  - `facs_prep_save(data, wsp_path, assay, panel_dir, condition_keyword = NA_character_) -> invisible(character)` path of the written `.rds`.
  - `facs_prep_load(wsp_path, panel_dir) -> list(data, info, paths)`.
  - Prep object written to disk: `list(data = <tibble>, info = list(wsp_path, wsp_root, wsp_version (Date), assay, condition_keyword, panel_dir, created_at (POSIXct), package_version (chr)))`. `facs_prep_load()` adds `paths = list(prep = <rds>, qc = <html>)`.
  - Internal `facs_prep_paths_(wsp_path, panel_dir) -> list(prep, qc)`: `file.path(panel_dir, "processed", paste0(base, "_prep.rds"))` and `file.path(panel_dir, "qc", paste0(base, "_qc.html"))`, `base = tools::file_path_sans_ext(basename(wsp_path))`.
  - Test helpers (in `helper-prep.R`): `local_prep_dirs(env)` -> list(`data_dir`, `panel_dir`, `wsp_path`) in a `withr::local_tempdir()`, with an empty file `<data_dir>/25-7_exp/raw/FACS/overview_25-7_20261101.wsp` and `panel_dir = <data_dir>/25-7_exp/raw/FACS/overview`; `fake_prep_data(files)` -> long tibble for a tibble `files` of `file_name, mouse_ID, tissue` (2 populations x metrics `count`, `pct_of_live_CD45+`, plus columns `population_full_path`, `population` (factor), `metric`, `value`).

- [ ] **Step 1: Write the failing tests**

```r
test_that("facs_prep_save writes rds with data + info and returns its path", {
  d <- local_prep_dirs()
  files <- tibble::tibble(file_name = c("a.fcs", "b.fcs"),
                          mouse_ID = c("25-7-1", "25-7-2"), tissue = "spleen")
  path <- facs_prep_save(fake_prep_data(files), d$wsp_path, "overview", d$panel_dir)
  expect_equal(path, file.path(d$panel_dir, "processed", "overview_25-7_20261101_prep.rds"))
  prep <- readRDS(path)
  expect_equal(prep$info$wsp_root, "overview_25-7")
  expect_equal(prep$info$wsp_version, as.Date("2026-11-01"))
  expect_equal(prep$info$assay, "overview")
  expect_equal(unique(prep$data$condition), "none")
  expect_setequal(unique(prep$data$file_name), c("a.fcs", "b.fcs"))
})

test_that("facs_prep_save drops files without mouse_ID or tissue", {
  # compensation / FMO files
  files <- tibble::tibble(file_name = c("a.fcs", "comp.fcs"),
                          mouse_ID = c("25-7-1", NA), tissue = c("spleen", NA))
  # ... expect only "a.fcs" in saved data
})

test_that("facs_prep_save takes condition from the keyword column", {
  # data with column restimulation = c("unstim", "MPO"); condition_keyword = "restimulation"
  # expect prep$data$condition == restimulation; NA in the keyword column -> error matching "restimulation"
})

test_that("facs_prep_save rejects a wsp name without _YYYYMMDD", {
  # wsp_path ".../overview_25-7.wsp" -> error matching "_YYYYMMDD"
})

test_that("facs_prep_save warns on files sharing mouse x tissue x condition", {
  # two files, same mouse_ID + tissue -> expect_warning(regexp = "share"); rds still written
})

test_that("facs_check_file_clashes lists every clashing file", {
  # 3 files: two share 25-7-1/spleen/none -> 2 rows; unique files -> 0 rows
})

test_that("facs_prep_load returns data, info and paths", {
  # save, then write an empty file at the qc path; load -> paths$qc == <panel_dir>/qc/overview_25-7_20261101_qc.html
})

test_that("facs_prep_load stops when the rds or report is missing", {
  # no rds -> error matching "final <- FALSE"; rds but no qc html -> same message
})

test_that("facs_prep_load stops when the wsp is newer than the rds", {
  # save + qc file, then Sys.setFileTime(d$wsp_path, Sys.time() + 60) -> error matching "newer"
})
```

- [ ] **Step 2: Run tests, verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_prep.R")'`
Expected: FAIL, `could not find function "facs_prep_save"`.

- [ ] **Step 3: Implement `R/facs_prep.R`**

Required columns for `data`: `file_name, population_full_path, population, metric, value, mouse_ID, tissue` (+ `condition_keyword` column when not `NA`); error naming missing ones. `wsp_root` = basename without `_YYYYMMDD.wsp` (regex `^(.+)_([0-9]{8})\.wsp$`), version via `as.Date(x, "%Y%m%d")`. `dir.create(recursive = TRUE)` for `processed/`. `package_version = as.character(utils::packageVersion("JanisHelpers"))`. Clash warning prints the clash tibble in its message and tells the user to fix keywords in FlowJo. Load messages: missing rds/report -> "Run the prep script with final <- FALSE first"; stale -> "<wsp> is newer than its prep result -- rerun with final <- FALSE". Roxygen with `\dontrun{}` examples.

- [ ] **Step 4: Run tests, verify they pass** (same command, all PASS)

- [ ] **Step 5: Commit**

```bash
git add R/facs_prep.R tests/testthat/test-facs_prep.R tests/testthat/helper-prep.R
git commit -m "feat(facs): save/load prep results and check file clashes"
```

---

### Task 2: Exclusion matching

**Files:**
- Create: `R/facs_exclude.R`
- Test: `tests/testthat/test-facs_exclude.R`

**Interfaces:**
- Consumes: prep `data` shape from Task 1 (needs `file_name, mouse_ID, tissue, condition`).
- Produces: `facs_match_exclusions(data, exclusions, panel) -> list(applicable, excluded_files, unmatched)`:
  - `applicable`: tibble, one row per exclusion row that applies to `panel`, columns `mouse_ID, tissue, panel, condition, file_name, reason, date` (all character except `date` Date) + `n_files_matched` (int).
  - `excluded_files`: tibble `file_name, reason` (distinct files; if several rows hit one file, reasons joined with `"; "`).
  - `unmatched`: the `applicable` rows with `n_files_matched == 0`.
  - `exclusions = NULL` or 0 rows -> all three empty tibbles with those columns.

- [ ] **Step 1: Write the failing tests**

```r
files <- tibble::tibble(
  file_name = c("m1_spl.fcs", "m1_lung.fcs", "m2_spl.fcs"),
  mouse_ID = c("25-7-1", "25-7-1", "25-7-2"),
  tissue = c("spleen", "lung", "spleen"), condition = "none")

test_that("whole-mouse row excludes all its files", {
  ex <- tibble::tibble(mouse_ID = "25-7-1", reason = "died early")  # other columns absent
  m <- facs_match_exclusions(files, ex, panel = "overview")
  expect_setequal(m$excluded_files$file_name, c("m1_spl.fcs", "m1_lung.fcs"))
  expect_equal(m$applicable$n_files_matched, 2L)
  expect_equal(nrow(m$unmatched), 0)
})

test_that("tissue and file_name narrow the match", {
  # mouse_ID+tissue="lung" -> "m1_lung.fcs"; mouse_ID+file_name="m2_spl.fcs" -> that file
})

test_that("rows for another panel are ignored, empty panel applies", {
  # panel = c("treg", NA) for 25-7-2 and 25-7-1 -> only 25-7-1 files; applicable has 1 row
})

test_that("unmatched applicable rows are returned, not an error", {
  # tissue = "Spleen" (case) -> unmatched has 1 row, excluded_files 0 rows
})

test_that("logical NA columns and NULL/0-row sheets work", {
  # ex with tissue = NA (logical), date as Date; NULL; tibble() with 0 cols -> empty results
})

test_that("missing mouse_ID/reason column or empty reason errors", {
  # expect_error(..., "reason")
})
```

- [ ] **Step 2: Run, verify fail** — `testthat::test_file("tests/testthat/test-facs_exclude.R")`, FAIL "could not find function".

- [ ] **Step 3: Implement `facs_match_exclusions(data, exclusions, panel)` in `R/facs_exclude.R`**

Add absent optional columns (`tissue, panel, condition, file_name, date`) as `NA`, coerce the five match columns + `reason` to character, `date` via `as.Date()`. Applicable = `is.na(panel) | panel == !!panel`. For each applicable row, a file matches when every non-`NA` match column equals the file's value (exact, no case folding). Error if `mouse_ID` or `reason` column is missing, or any applicable row has `NA`/empty `mouse_ID` or `reason`.

- [ ] **Step 4: Run, verify pass**

- [ ] **Step 5: Commit** — `git commit -m "feat(facs): match exclusion sheet rows to prep files"`

---

### Task 3: QC report template + `report_knit_facs_qc()`

**Files:**
- Create: `inst/rmarkdown/templates/facs-qc-template/template.yaml`
- Create: `inst/rmarkdown/templates/facs-qc-template/skeleton/skeleton.Rmd`
- Modify: `R/report_knit.R` (append)
- Test: `tests/testthat/test-report_knit.R` (append)

**Interfaces:**
- Consumes: `facs_prep_save()` output path, `facs_check_file_clashes()`, `facs_match_exclusions()`, `meta_read()`.
- Produces: `report_knit_facs_qc(prep_path, meta_path = NULL) -> character` absolute path of the rendered html at `<panel_dir>/qc/<base>_qc.html` (`base` = prep basename minus `_prep.rds`). Creates `qc/`.

- [ ] **Step 1: Write the failing test**

```r
test_that("report_knit_facs_qc renders next to the prep file and flags exclusions", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  files <- tibble::tibble(file_name = c("a.fcs", "b.fcs"),
                          mouse_ID = c("25-7-1", "25-7-2"), tissue = "spleen")
  prep <- facs_prep_save(fake_prep_data(files), d$wsp_path, "overview", d$panel_dir)
  meta <- file.path(d$data_dir, "meta.xlsx")
  # write meta.xlsx with sheet "exclusions": mouse_ID="25-7-2", reason="clots"
  # (use writexl if installed: skip_if_not_installed("writexl"))
  out <- report_knit_facs_qc(prep, meta)
  expect_equal(normalizePath(out),
    normalizePath(file.path(d$panel_dir, "qc", "overview_25-7_20261101_qc.html")))
  html <- paste(readLines(out), collapse = "\n")
  expect_match(html, "b.fcs")
  expect_match(html, "clots")
})
```

Also a test without `meta_path` (renders, no exclusion section error). If `writexl` isn't installed, the exclusion test is skipped; don't add it to `Suggests`.

- [ ] **Step 2: Run, verify fail**

- [ ] **Step 3: Implement**

Template `params: prep: ""` and `meta: ""`; sections, in order: header (wsp, assay, prep `created_at`, package version); **Files** table (one row per file: `file_name, mouse_ID, tissue, condition`, number of populations, `excluded` reason or blank) via `knitr::kable`; **Keyword clashes** (warning text + table, or "none"); **Exclusions** (applicable rows; unmatched rows shown under a visible "Warning: these rows match no file" heading); **QC checks** placeholder headings ("Gating", "Counts", "Outliers") each with a `<!-- add checks -->` note. All package calls `JanisHelpers::`. `template.yaml`: `name: FACS QC report`, `description`, `create_dir: false`.
`report_knit_facs_qc()`: `rmarkdown::render(system.file(...skeleton.Rmd), output_dir = qc_dir, output_file = basename(html), params = list(prep = normalizePath(prep_path), meta = meta_path or ""), intermediates_dir = tempdir(), envir = new.env(), quiet = TRUE)` (intermediates outside OneDrive). Return its value.

- [ ] **Step 4: Run, verify pass**

- [ ] **Step 5: Commit** — `git commit -m "feat(report): facs QC report skeleton rendered from prep result"`

---

### Task 4: `import_exclusions` table

**Files:**
- Modify: `R/db_connect.R` (after `import_log` block, ~line 191)
- Test: `tests/testthat/test-db_connect.R`

**Interfaces:**
- Produces table:
  ```sql
  import_exclusions (
    import_exclusion_id INTEGER PRIMARY KEY DEFAULT nextval('seq_import_exclusions'),
    import_id       INTEGER NOT NULL REFERENCES import_log(import_id),
    mouse_id        VARCHAR NOT NULL,
    tissue          VARCHAR,
    panel           VARCHAR,
    condition       VARCHAR,
    file_name       VARCHAR,
    reason          VARCHAR NOT NULL,
    decided_on      DATE,
    n_files_matched INTEGER NOT NULL
  )
  ```

- [ ] **Step 1: Failing tests** — extend "db_connect creates all tables" to expect `import_exclusions`; extend the EXPORT DATABASE test to insert one `import_log` row (via `JanisHelpers:::write_import_log_()`) + one `import_exclusions` row and expect `import_exclusions.parquet`.
- [ ] **Step 2: Run, verify fail** — `testthat::test_file("tests/testthat/test-db_connect.R")`
- [ ] **Step 3: Add sequence + `CREATE TABLE IF NOT EXISTS`** (comment: FK to another table, fine for snapshots); mention the table in `db_connect()` roxygen.
- [ ] **Step 4: Run, verify pass** (whole `test-db_connect.R`)
- [ ] **Step 5: Commit** — `git commit -m "feat(db): import_exclusions table"`

---

### Task 5: `db_export_facs()`

**Files:**
- Modify: `R/db_write_measurements.R` (extract internal insert from `db_write_facs()`)
- Create: `R/db_export.R`
- Test: `tests/testthat/test-db_export.R`

**Interfaces:**
- Consumes: `facs_prep_load()` result (Task 1), `facs_match_exclusions()`, `facs_check_file_clashes()`, `write_import_log_()`, `experiment_codes_()`, `check_facs_root_()`, `check_facs_root_keyword_()`, `lookup_sample_id_()`, `lookup_assay_id_()`, `db_write_samples()`, `db_write_assay()`.
- Produces:
  - Internal `insert_facs_stain_(con, sample_id, assay_id, data, wsp_root, source_file, condition, condition_keyword, stain_date = NA) -> integer` rows inserted; **opens no transaction**; `data` has snake-case `population_full_path, population, metric, value`. `db_write_facs()` keeps its signature/behaviour and calls it inside its own `dbWithTransaction()`.
  - `db_export_facs(con, prep, exclusions = NULL, data_dir = "data", overwrite = NULL) -> invisible(integer)` the new `import_id`, or `invisible(NULL)` when the user declines.

- [ ] **Step 1: Failing tests** (fixture: `local_test_db()`, experiment `25-7`, subjects `25-7-1`, `25-7-2`; prep from `local_prep_dirs()` + `fake_prep_data()` saved and an empty qc html written, then `facs_prep_load()`; `data_dir = d$data_dir`)

```r
test_that("exports all files, logs the import with relative paths", {
  id <- db_export_facs(con, prep, data_dir = d$data_dir)
  log <- DBI::dbGetQuery(con, "SELECT * FROM import_log WHERE import_id = ?", list(id))
  expect_equal(log$import_unit, "overview_25-7")
  expect_equal(log$source_file, "25-7_exp/raw/FACS/overview_25-7_20261101.wsp")
  expect_equal(log$qc_report, "25-7_exp/raw/FACS/overview/qc/overview_25-7_20261101_qc.html")
  expect_true(log$qc_passed)
  expect_equal(as.Date(log$source_version), as.Date("2026-11-01"))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM facs_stains")$n, 2)
  expect_equal(log$rows_written,
    DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM facs_measurements")$n)
})

test_that("applies exclusions and records them with the import", {
  ex <- tibble::tibble(mouse_ID = "25-7-2", reason = "clots")
  id <- db_export_facs(con, prep, ex, data_dir = d$data_dir)
  # 1 stain left (25-7-1); import_exclusions has 1 row, import_id = id,
  # reason "clots", n_files_matched 1, tissue/panel/file_name NA
})

test_that("unmatched exclusion row stops the export, nothing written", {
  # tissue = "Spleen" -> expect_error(regexp = "match no file"); 0 stains, 0 import_log rows
})

test_that("all files excluded stops the export", { ... "nothing left" })

test_that("keyword clash stops the export", { ... two files same mouse/tissue -> error "share" })

test_that("existing root: overwrite = FALSE aborts without changes", {
  # export once, then a second prep version (new wsp file _20261105, changed value) ->
  # db_export_facs(..., overwrite = FALSE) returns NULL; values + log unchanged (1 row)
})

test_that("existing root: overwrite = TRUE replaces data and links the log", {
  # second export -> measurements carry new source_file; import_log has 2 rows;
  # new row's supersedes_import_id == first id; first import's exclusions still present
})

test_that("warns when a newer dated wsp of the root exists on disk", {
  # file.create(<FACS dir>/overview_25-7_20261201.wsp) -> expect_warning(regexp = "newer")
})
```

- [ ] **Step 2: Run, verify fail** — `testthat::test_file("tests/testthat/test-db_export.R")`

- [ ] **Step 3: Refactor `db_write_facs()`** onto `insert_facs_stain_()`; run `test-db_write_measurements.R` — all PASS unchanged.

- [ ] **Step 4: Implement `db_export_facs()` in `R/db_export.R`**

Order (all checks before any write):
1. `prep` must have `paths` (came from `facs_prep_load()`), else error.
2. Clashes -> error listing them.
3. `facs_match_exclusions(prep$data, exclusions, prep$info$assay)`; unmatched -> error listing rows ("match no file"); drop `excluded_files`; 0 files left -> error ("nothing left to export").
4. Newer-version warning: files in `dirname(wsp_path)` matching `^<wsp_root>_[0-9]{8}\.wsp$` with a later date.
5. `source_file` / `qc_report` = path relative to `normalizePath(data_dir)`; error if outside it.
6. `db_write_samples()` (distinct `mouse_id`, `tissue`), `db_write_assay(assay, "facs")`.
7. Existing root? (`facs_stains.wsp_root`); previous id = max `import_id` in `import_log` for domain `facs` + `import_unit = wsp_root`. If existing: `overwrite` `NULL` -> `utils::askYesNo("Root '<root>' is already in the database (from <old source_file>). Overwrite?")`; not `TRUE` -> `message("Not exported -- nothing changed.")`, return `invisible(NULL)`.
8. One `DBI::dbWithTransaction()`: delete root's measurements, then stains; per file: `lookup_sample_id_()`, `check_facs_root_()`, `check_facs_root_keyword_()`, `insert_facs_stain_()`; `write_import_log_(con, "facs", experiment_codes_(...), wsp_root, source_file, n, source_version = wsp_version, qc_report, qc_passed = TRUE, supersedes_import_id = previous)`; insert `applicable` rows into `import_exclusions` (`decided_on` <- `date`).
9. `message()` stains / rows / excluded files; return `invisible(import_id)`.

Roxygen: `@param overwrite` "NULL asks interactively; TRUE/FALSE skips the prompt". Update `db_write_facs()` roxygen: "Stage 3 export" -> `[db_export_facs()]`.

- [ ] **Step 5: Run, verify pass** — `test-db_export.R` and `test-db_write_measurements.R`.

- [ ] **Step 6: Commit** — `git commit -m "feat(db): db_export_facs replaces a root with QC'd prep data and logs exclusions"`

---

### Task 6: Script template, docs, full check

**Files:**
- Create: `inst/scripts/facs_prep_template.R`
- Modify: `CLAUDE.md` (Data structures: prep/exclusions/export bullets; RMarkdown templates table: `facs-qc-template` -> `report_knit_facs_qc()`)
- Modify: `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md` (add `import_exclusions`)
- Test: `tests/testthat/test-facs_prep.R` (append)

- [ ] **Step 1: Failing test** — `expect_no_error(parse(system.file("scripts", "facs_prep_template.R", package = "JanisHelpers")))` and `expect_true(nzchar(...))`.
- [ ] **Step 2: Write the template.** Config block: `final <- FALSE`, `wsp`, `panel_dir` (folder name under `raw/FACS/`), `assay`, `condition_keyword <- NA`, `meta_path`, `data_dir <- "data"`. Then `if (!final)`: `facs_read_wsp(wsp_path, keywords = c("mouse_ID", "tissue"))$data` -> population selection vector (placeholder) -> `facs_calc_*()` (placeholder) -> `facs_prep_save()` -> `report_knit_facs_qc()`. `else`: `prep <- facs_prep_load()`, `exclusions <- meta_read(meta_path)$exclusions` (`NULL` if no sheet), `con <- db_open_...` placeholder, `db_export_facs()`, snapshot placeholder. Comments explain the loop in 3-4 lines.
- [ ] **Step 3: Docs** — CLAUDE.md bullets + templates row; schema spec section for `import_exclusions` (columns as Task 4, why kept on overwrite).
- [ ] **Step 4: Verify** — `Rscript -e 'devtools::document(); devtools::test()'` all PASS; `Rscript -e 'devtools::check()'` -> 0 errors, 0 warnings, only known NOTEs.
- [ ] **Step 5: Commit** — `git commit -m "docs: prep script template, CLAUDE.md and schema spec for prep/QC/export"`

---

### Task 7: Rework the 26-1 project scripts (outside this repo)

Project dir: `~/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Shared Documents/MPO_mRNA/data/26-1_BMT-mRNA-V4/src/`. Don't commit there (ask the user). Git/IO on OneDrive is slow; use long timeouts.

- [ ] **Step 1:** `db_import_26-1_common.R`: delete `db_import_facs_wsp()`; keep config, `db_open_26_1()`, `db_snapshot_close()`; add `meta_path`.
- [ ] **Step 2:** Create `facs_prep_26-1_overview.R` from the template, carrying over the user's population selection and calculations from `db_import_26-1_overview.R`: `wsp_path` built from `wsp` (was undefined `ov_path`), `dplyr::mutate`, meta for `facs_calc_count_per_g()` from `meta_clean(meta_read(meta_path))` — confirm the sheet/column names with the user if `meta$calculations` doesn't exist — and drop the `exclude` filter. `panel_dir <- "overview"`, `assay <- "overview"`.
- [ ] **Step 3:** Same for `treg` (`panel_dir <- "Treg"`) and `mpo-restim` (`panel_dir <- "MPO-Restim"`, `condition_keyword <- "restimulation"`), with empty selection/calculation placeholders for the user. Remove the old `db_import_26-1_<panel>.R` files.
- [ ] **Step 4:** Tell the user to add the `exclusions` sheet to `meta_26-1.xlsx` in Excel (header row `mouse_ID | tissue | panel | condition | file_name | reason | date`) — not written programmatically, to avoid damaging the workbook's formatting.
- [ ] **Step 5:** Dry run of overview with `final <- FALSE` against the real `.wsp` (needs JanisHelpers installed in the project's renv — report if missing rather than installing).
