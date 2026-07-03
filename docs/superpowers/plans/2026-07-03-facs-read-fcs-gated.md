# facs_read_fcs_gated() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `facs_read_fcs_gated()` — Stage 1 of the unsupervised FACS analysis pipeline (`docs/superpowers/specs/2026-07-03-facs-unsupervised-analysis-roadmap.md`). It reads raw single-cell events straight from `.fcs` files, replaying the compensation/transformation/gating tree from the matching FlowJo `.wsp` workspace via `CytoML`/`flowWorkspace`, so events are filtered down to an arbitrary already-drawn gate before being handed to later stages (FlowSOM clustering, UMAP, diffcyt).

**Architecture:** New file `R/facs_read_fcs.R` with one exported function, `facs_read_fcs_gated()`, and two unexported helpers: `resolve_markers_()` (matches requested marker names against a sample's panel, stain label preferred over channel name) and `read_one_sample_()` (does the per-sample gate lookup, marker subsetting, downsampling, and keyword lookup). The exported function opens the workspace once, builds one `GatingSet` for the requested `group`, maps `read_one_sample_()` over every sample via `purrr::map()`, and row-binds the results.

**Tech Stack:** R, `CytoML`, `flowWorkspace`, `flowCore` (all three new, Bioconductor, already installed on this machine — confirmed via `requireNamespace()`), plus already-`Imports`ed `dplyr`, `purrr`, `tibble`, `glue`. Testing via `testthat` 3e against a real fixture (`tests/fixtures/Treg.wsp` + `tests/fixtures/Treg/*.fcs`) — no mocking, per this project's existing no-mock-`fcexpr`/no-mock-`CytoML` rule.

## Global Constraints

- Pipe: `|>` (base pipe) everywhere.
- Namespacing: every external function call inside an exported (or shared unexported) function body must be `pkg::fn()` — no bare `dplyr`/`purrr`/`tibble`/`glue`/`flowWorkspace`/`flowCore`/`CytoML` calls.
- New `Imports`: `CytoML`, `flowWorkspace`, `flowCore` — Bioconductor, non-CRAN, same treatment as the existing `fcexpr` dependency (listed in `Imports`, install requirement documented in `README.md`, no `BiocManager::`/`remotes::` calls inside function bodies).
- Internal helpers use a trailing underscore (`resolve_markers_`, `read_one_sample_`), matching the existing `resolve_var_`/`filter_samples_`/`walk_pops_` convention — no `domain_` prefix needed.
- `facs_read_fcs_gated()` returns visibly (not `invisible()`) — its output is the primary deliverable, not a side effect.
- `gate_path` argument format has no leading slash (matches `PopulationFullPath` from `facs_read_wsp()`, e.g. `"Singlets/Lymphocytes/live/CD45+"`), but `flowWorkspace::gh_pop_get_data()` requires a leading slash internally — normalize by prepending `/` if missing, and always use the **original**, non-normalized `gate_path` in any user-facing warning/error message.
- Test fixture (`tests/fixtures/Treg.wsp`, `tests/fixtures/Treg/*.fcs`) is `.gitignore`d (contains real research metadata) — every test must `skip_if_not(dir.exists(...))`, matching the existing `test-facs_read.R` pattern.
- `@param`, `@returns`, `@export` required on the exported function; `@examples` uses `\dontrun{}` (requires external files). No roxygen block on unexported helpers — a one-line `#` comment is enough, matching the `resolve_var_()` precedent.
- After any roxygen change, run `devtools::document()` before committing.
- Run `devtools::check()` before the final commit — target 0 errors, 0 warnings.
- Non-ASCII characters in R source must be escaped as `\uXXXX`.

---

## Verified-against-real-fixture facts (do not re-derive, use directly)

These were confirmed by actually running `CytoML`/`flowWorkspace` against `tests/fixtures/Treg.wsp` + `tests/fixtures/Treg/*.fcs` during design/planning — not guessed from documentation:

- `CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir)` — `name` is the FlowJo group argument, `path` is the folder `flowWorkspace` recursively searches for `.fcs` files by filename (confirmed: files sit in a subfolder, not next to the `.wsp`, and lookup still worked).
- Gating tree for every sample: `root -> /Singlets -> /Singlets/Lymphocytes -> /Singlets/Lymphocytes/live -> /Singlets/Lymphocytes/live/CD45+ -> /Singlets/Lymphocytes/live/CD45+/TCRb_CD4+`. (The original tree had a `/time` gate first; it was removed from the fixture because the downsampled `.fcs` files have a collapsed `Time` channel that made the original time-window gate filter out 100% of events — this is fixed now.)
- `flowWorkspace::sampleNames(gs)` returns internal names with a `_10000` suffix (e.g. `"minimal_26-1-1_kidney_B02.fcs_10000"`) — the clean original filename is `flowWorkspace::pData(gh)$name` (e.g. `"minimal_26-1-1_kidney_B02.fcs"`). **Use `pData(gh)$name` for the `file_name` column, never `sampleNames()`.**
- `gs[[sample_name_string]]` (indexing a `GatingSet` by one of `sampleNames(gs)`'s values) returns a valid `GatingHierarchy`.
- `flowWorkspace::gh_pop_get_data(gh, path)` on a non-existent path throws a catchable `simpleError` (message like `"root/Nonexistent/Path not found!"`) — safe to wrap in `tryCatch(error = function(e) NULL)`.
- `flowWorkspace::gh_pop_get_data()`'s return value is a lazy view — `flowCore::exprs()` on it directly returns 0 rows. Must call `flowWorkspace::realize_view()` first.
- `flowWorkspace::markernames(gh)` returns a named character vector, `channel -> stain label`, but **only for channels that have a stain label** (`$PnS` non-empty). Channels without one (e.g. `FSC-A`, `SSC-A`, `Time`) are absent from it — for those, matching falls back to `flowCore::parameters(gated)$name`, which lists every channel.
- `flowWorkspace::keyword(gh, "some_missing_keyword")` returns `NULL` (not an error) when the keyword doesn't exist for that sample — safe to check with `is.null()`.
- The fixture's `"kidney"` FlowJo group contains exactly 2 samples: `minimal_26-1-1_kidney_B02.fcs`, `minimal_26-1-2_kidney_B03.fcs` (out of 6 total in `"All Samples"`) — useful for testing the `group` argument.
- Confirmed stain labels present in the fixture's panel: `TCRb, PD-1, TIGIT, CD4, CD44, LAG-3, CD25, KLRG1, CD45, ICOS, RORgt, CLTA-4, FoxP3, T-bet, TIM-3, LD`. `FSC-A` has no stain label (channel-name-only match case).
- A full manually-run draft of the final implementation (see Task 1) was executed end-to-end against this fixture covering every scenario below, and every one produced the expected result (exact row/warning/error output recorded in each task's "Expected" block).

---

## File Structure

- **Create** `R/facs_read_fcs.R` — `facs_read_fcs_gated()` (exported) + `resolve_markers_()`, `read_one_sample_()` (unexported).
- **Create** `tests/testthat/test-facs_read_fcs.R` — all tests for the new function.
- **Modify** `DESCRIPTION` — add `CytoML`, `flowWorkspace`, `flowCore` to `Imports`.
- **Modify** `README.md` — document the Bioconductor install requirement.
- **Modify** `tests/fixtures/README.md` — document the `Treg.wsp`/`Treg/` fixture (mirrors the existing `minimal.wsp` section).
- **Modify** `CLAUDE.md` — only if `devtools::check()` reports new NOTEs (Task 4).

---

## Task 1: Dependencies, fixture docs, core implementation, happy-path + keyword tests

**Files:**
- Modify: `DESCRIPTION`
- Modify: `README.md`
- Modify: `tests/fixtures/README.md`
- Create: `R/facs_read_fcs.R`
- Create: `tests/testthat/test-facs_read_fcs.R`

**Interfaces:**
- Produces: `facs_read_fcs_gated(wsp_path, gate_path, markers, keywords = NULL, fcs_dir = NULL, group = "All Samples", max_events = NULL, seed = NULL)` — exported. Returns a wide tibble: `file_name` (chr) + one column per `markers` entry (dbl, named exactly as requested) + one column per `keywords` entry (chr).
- Produces (unexported, used only within this file): `resolve_markers_(markers, gh, gated, file_name)` — returns a character vector of resolved channel names, same length/order as `markers`, or `stop()`s listing unmatched markers. `read_one_sample_(gh, gate_path, gate_path_norm, markers, keywords, max_events)` — returns a per-sample tibble (or an empty `tibble::tibble()` if `gate_path` doesn't exist for that sample, with a `warning()`).

- [ ] **Step 1: Add new Bioconductor dependencies to `DESCRIPTION`**

In `DESCRIPTION`, the `Imports:` block currently reads:

```
Imports:
    dplyr,
    fcexpr,
    glue,
    gt,
    gtsummary,
    janitor,
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stringr,
    tibble,
    tidyr,
    xml2,
    xfun
```

Replace it with (alphabetical, matching existing order):

```
Imports:
    CytoML,
    dplyr,
    fcexpr,
    flowCore,
    flowWorkspace,
    glue,
    gt,
    gtsummary,
    janitor,
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stringr,
    tibble,
    tidyr,
    xml2,
    xfun
```

- [ ] **Step 2: Document the Bioconductor install requirement in `README.md`**

In `README.md`, after the line `Import FowJO workspace` and before `Automatically import and save FlowJo workspace to Excel`, the bullet list currently reads:

```
- report_knit_dated: put in yaml behind “knit:” to automatically knit
  dated versions to specified directory
- Auto-knit-report: rmarkdown template das automatisch datiert dateien
  in spezifische subdirectories erstellt
- protocol-template: template für Protokolle
- Import FowJO workspace
- Automatically import and save FlowJo workspace to Excel
```

Add a new paragraph immediately after that list, before the `Updaten des privaten Pakets:` line:

```
`facs_read_fcs_gated()` reads raw single-cell events from `.fcs` files and
needs `CytoML`/`flowWorkspace`/`flowCore`, which are Bioconductor packages
not available on CRAN. Install them first:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("CytoML", "flowWorkspace", "flowCore"))
```
```

- [ ] **Step 3: Document the `Treg` fixture in `tests/fixtures/README.md`**

Read the current file first, then append a new section after the existing `facs_import_wsp() — pending` section (same file already read during design: it documents `minimal.wsp`'s requirements and creation steps). Append:

```markdown
## facs_read_fcs_gated() — Treg.wsp + Treg/

`facs_read_fcs_gated()` requires real `.fcs` files in addition to a `.wsp`
workspace (unlike `facs_read_wsp()`, which only reads the workspace XML).
No mock of `CytoML`/`flowWorkspace` is allowed (see CLAUDE.md).

### Required fixture: `Treg.wsp` + `Treg/`

Both are `.gitignore`d (real research metadata: project/experiment names,
mouse strain/treatment) — regenerate locally if missing.

- **`Treg.wsp`**, in `tests/fixtures/`, containing:
  - 6 samples (2 mice x 3 tissues: kidney, lung, spleen)
  - A gating hierarchy at least 4 levels deep:
    `Singlets -> Lymphocytes -> live -> CD45+ -> TCRb_CD4+`
    (do **not** include a `Time`-based QC gate — see the note below)
  - Custom keywords per sample, at minimum `mouse_ID` and `tissue`
  - A panel with several stain-labelled channels (`$PnS`) and at least one
    channel with no stain label (e.g. `FSC-A`), to exercise both the
    stain-label and channel-name marker-matching paths
  - A FlowJo sample group named `"kidney"` containing only the 2
    kidney-tissue samples (in addition to the default `"All Samples"`
    group), to exercise the `group` argument
- **`Treg/`**, a subfolder next to `Treg.wsp` (named after it, sans
  extension — this is what `facs_read_fcs_gated()`'s `fcs_dir`
  auto-derivation expects), containing the matching downsampled `.fcs`
  files (10,000 events each is plenty).

**Important:** if downsampling in FlowJo before export, check whether the
export preserves each event's original acquisition `Time` value. If it
resets/collapses `Time` (all events end up with the same value), any
`Time`-based gate in the workspace will filter out 100% of events after
replay. This actually happened while creating this fixture — the original
workspace had a `time` gate as the first step, which had to be removed
after the downsampled export collapsed every event's `Time` value to the
same number.

### How to create the fixture

1. Open FlowJo 10 with a small pilot experiment (this one used
   `26-1_Treg-V2`, 2 mice x 3 tissues).
2. Downsample each sample to ~10,000 events and export as new `.fcs`
   files (verify `Time` values are preserved post-export, per the note
   above).
3. Build a gating hierarchy without a `Time`-based gate (or verify any
   existing one still passes events on the downsampled data).
4. Create a `"kidney"` sample group containing only the kidney-tissue
   samples.
5. Save the workspace as `tests/fixtures/Treg.wsp`, with the exported
   `.fcs` files in `tests/fixtures/Treg/`.
6. Confirm both are listed in `.gitignore`.
```

- [ ] **Step 4: Write `R/facs_read_fcs.R`**

```r
# FlowWorkspace/CytoML behavior — verified against tests/fixtures/Treg.wsp
#
# flowjo_to_gatingset(ws, name = <group>, path = <fcs dir>) builds a
# GatingSet, recursively searching <fcs dir> for the .fcs files the
# workspace references, replaying compensation/transformation/gating.
#
# sampleNames(gs) returns internal names with a numeric suffix
# (e.g. "sample.fcs_10000") -- the clean original filename is
# pData(gh)$name. Always use pData(gh)$name for file_name, never
# sampleNames().
#
# gh_pop_get_data(gh, path) returns a lazy view; realize_view() must be
# called before exprs() returns real data. A non-existent path throws a
# catchable simpleError.
#
# markernames(gh) returns a named vector (channel = stain label) but only
# for channels that have a stain label. Channels without one (FSC-A,
# SSC-A, Time, ...) must be matched via parameters(gated)$name instead.
#
# keyword(gh, name) returns NULL (not an error) when the keyword is
# absent for that sample.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

resolve_markers_ <- function(markers, gh, gated, file_name) {
  stain_lookup  <- flowWorkspace::markernames(gh)
  channel_names <- flowCore::parameters(gated)$name

  resolved <- purrr::map_chr(markers, function(m) {
    if (m %in% stain_lookup) {
      return(names(stain_lookup)[stain_lookup == m][[1]])
    }
    if (m %in% channel_names) {
      return(m)
    }
    NA_character_
  })

  unmatched <- markers[is.na(resolved)]
  if (length(unmatched) > 0L) {
    stop(glue::glue(
      "The following markers could not be matched (by stain label or ",
      "channel name) in '{file_name}': {paste(unmatched, collapse = ', ')}"
    ))
  }

  resolved
}

read_one_sample_ <- function(gh, gate_path, gate_path_norm, markers, keywords, max_events) {
  file_name <- flowWorkspace::pData(gh)$name

  gated <- tryCatch(
    flowWorkspace::gh_pop_get_data(gh, gate_path_norm),
    error = function(e) NULL
  )
  if (is.null(gated)) {
    warning(glue::glue(
      "gate_path '{gate_path}' not found for '{file_name}'; sample skipped."
    ))
    return(tibble::tibble())
  }
  gated <- flowWorkspace::realize_view(gated)

  resolved_cols <- resolve_markers_(markers, gh, gated, file_name)

  mat <- flowCore::exprs(gated)[, resolved_cols, drop = FALSE]
  colnames(mat) <- markers

  if (!is.null(max_events) && nrow(mat) > max_events) {
    mat <- mat[sample.int(nrow(mat), max_events), , drop = FALSE]
  }

  tbl <- tibble::as_tibble(mat) |>
    tibble::add_column(file_name = file_name, .before = 1L)

  if (!is.null(keywords) && length(keywords) > 0L) {
    for (k in keywords) {
      v <- flowWorkspace::keyword(gh, k)
      tbl[[k]] <- if (is.null(v)) NA_character_ else as.character(v)
    }
  }

  tbl
}

# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' Read raw single-cell events from .fcs files, filtered to a FlowJo gate
#'
#' @description
#' Reads raw single-cell events directly from the .fcs files referenced by
#' a FlowJo \code{.wsp} workspace, replaying the workspace's compensation,
#' transformation, and gating tree (via \code{CytoML}/\code{flowWorkspace})
#' so events are filtered down to an arbitrary already-drawn gate. Feeds
#' unsupervised, cell-level analysis (e.g. FlowSOM, UMAP) that
#' \code{facs_read_wsp()}'s gated summary statistics cannot support.
#'
#' @param wsp_path path to the \code{.wsp} file.
#' @param gate_path character; full gating path, e.g.
#'   \code{"Singlets/Lymphocytes/live/CD45+"} — same format as
#'   \code{PopulationFullPath} from \code{facs_read_wsp()}. Applied to
#'   every sample in \code{group}.
#' @param markers character vector; matched per sample against stain
#'   label or channel name (stain label preferred).
#' @param keywords character vector of FlowJo keyword names to append as
#'   columns. A keyword missing for every sample is filled
#'   \code{NA_character_} with a warning.
#' @param fcs_dir folder to search for this workspace's \code{.fcs}
#'   files. \code{NULL} (default) auto-derives it as the subfolder named
#'   after the \code{.wsp} file (sans extension), sitting next to it.
#' @param group character; FlowJo sample group to load. Default
#'   \code{"All Samples"}.
#' @param max_events integer; if set, randomly downsample each sample to
#'   at most this many events.
#' @param seed integer; if set, seeds the random draw used by
#'   \code{max_events} for reproducibility.
#'
#' @returns Wide tibble, one row per event: \code{file_name}, one column
#'   per requested marker (on FlowJo's transformed scale), and one column
#'   per requested keyword. Errors if a requested marker cannot be
#'   matched in a sample's panel. Warns and skips a sample if
#'   \code{gate_path} does not exist for it. Warns and fills \code{NA} if
#'   a requested keyword is missing for every sample.
#' @export
#'
#' @examples
#' \dontrun{
#'   facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45"),
#'     keywords  = c("mouse_ID", "tissue")
#'   )
#' }
facs_read_fcs_gated <- function(wsp_path,
                                 gate_path,
                                 markers,
                                 keywords = NULL,
                                 fcs_dir = NULL,
                                 group = "All Samples",
                                 max_events = NULL,
                                 seed = NULL) {
  if (is.null(fcs_dir)) {
    fcs_dir <- file.path(
      dirname(wsp_path),
      tools::file_path_sans_ext(basename(wsp_path))
    )
  }
  gate_path_norm <- if (startsWith(gate_path, "/")) gate_path else paste0("/", gate_path)

  ws <- CytoML::open_flowjo_xml(wsp_path)
  gs <- CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir)

  if (!is.null(seed)) set.seed(seed)

  data <- purrr::map(
    flowWorkspace::sampleNames(gs),
    function(sn) {
      read_one_sample_(
        gh             = gs[[sn]],
        gate_path      = gate_path,
        gate_path_norm = gate_path_norm,
        markers        = markers,
        keywords       = keywords,
        max_events     = max_events
      )
    }
  ) |>
    dplyr::bind_rows()

  if (!is.null(keywords) && length(keywords) > 0L) {
    fully_missing <- keywords[purrr::map_lgl(keywords, function(k) all(is.na(data[[k]])))]
    if (length(fully_missing) > 0L) {
      warning(glue::glue(
        "The following requested keywords were not found in the workspace ",
        "and were filled with NA: {paste(fully_missing, collapse = ', ')}"
      ))
    }
  }

  data
}
```

- [ ] **Step 5: Write the happy-path and keyword-joining tests**

Create `tests/testthat/test-facs_read_fcs.R`:

```r
library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/Treg.wsp")
fcs_dir  <- testthat::test_path("../fixtures/Treg")
skip_msg <- "Treg fixture not available"
cd45_gate <- "Singlets/Lymphocytes/live/CD45+"

test_that("facs_read_fcs_gated() returns a wide tibble with one row per event", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "CD45")
  )

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("file_name", "CD4", "CD45") %in% names(result)))
  expect_type(result$CD4, "double")
  expect_type(result$CD45, "double")
  expect_equal(dplyr::n_distinct(result$file_name), 6L)
  expect_gt(nrow(result), 1000L)
  # Transformed (biexponential/logicle) scale, not raw fluorescence intensity
  expect_true(max(abs(result$CD4)) < 10000)
})

test_that("facs_read_fcs_gated() attaches requested keyword columns", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "CD45"),
    keywords  = c("mouse_ID", "tissue")
  )

  expect_true(all(c("mouse_ID", "tissue") %in% names(result)))
  expect_false(any(is.na(result$mouse_ID)))
  expect_false(any(is.na(result$tissue)))
  expect_setequal(unique(result$tissue), c("kidney", "lung", "spleen"))
})
```

- [ ] **Step 6: Run the new tests**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_read_fcs.R")'`
Expected: Both tests pass (or both skip with "Treg fixture not available" if the fixture isn't present on the machine running this — the fixture exists locally on the machine this plan was written on).

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION README.md tests/fixtures/README.md R/facs_read_fcs.R tests/testthat/test-facs_read_fcs.R
git commit -m "feat: add facs_read_fcs_gated() for raw single-cell event import with gate replay"
```

---

## Task 2: Error-handling tests (missing gate, unmatched marker, missing keyword)

**Files:**
- Modify: none (implementation already complete from Task 1 — this task only adds tests for branches not yet exercised)
- Test: `tests/testthat/test-facs_read_fcs.R`

**Interfaces:**
- Consumes: `facs_read_fcs_gated()` from Task 1 (already implements every behavior below — this task proves it).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_read_fcs.R`:

```r
test_that("facs_read_fcs_gated() warns and skips every sample when gate_path matches none", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  expect_warning(
    result <- facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = "Nonexistent/Path",
      markers   = c("CD4")
    ),
    "not found"
  )
  expect_equal(nrow(result), 0L)
})

test_that("facs_read_fcs_gated() errors immediately when a marker cannot be matched", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  expect_error(
    facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = cd45_gate,
      markers   = c("CD4", "NotARealMarker")
    ),
    "NotARealMarker"
  )
})

test_that("facs_read_fcs_gated() warns and fills NA when a keyword is missing for every sample", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  expect_warning(
    result <- facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = cd45_gate,
      markers   = c("CD4"),
      keywords  = c("mouse_ID", "not_a_real_keyword")
    ),
    "not_a_real_keyword"
  )
  expect_true(all(is.na(result$not_a_real_keyword)))
  expect_false(any(is.na(result$mouse_ID)))
})
```

- [ ] **Step 2: Run to verify these specific behaviors**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_read_fcs.R")'`
Expected: All 5 tests so far pass (or all skip together if the fixture is absent).

Recorded evidence from manually running the implementation against the real fixture during planning (for reference, not something this step needs to reproduce): requesting `"Nonexistent/Path"` produced 6 warnings (one per sample, each ending `"; sample skipped."`) and a `0 x 0` result; requesting marker `"NotARealMarker"` raised an error whose message contains `"NotARealMarker"`; requesting keyword `"not_a_real_keyword"` produced exactly one warning naming it, with `mouse_ID` still fully populated.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-facs_read_fcs.R
git commit -m "test: cover missing-gate, unmatched-marker, and missing-keyword behavior in facs_read_fcs_gated()"
```

---

## Task 3: Configuration-knob tests (downsampling, fcs_dir, group, channel-name matching)

**Files:**
- Modify: none (implementation already complete from Task 1)
- Test: `tests/testthat/test-facs_read_fcs.R`

**Interfaces:**
- Consumes: `facs_read_fcs_gated()` from Task 1.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_read_fcs.R`:

```r
test_that("facs_read_fcs_gated() matches a marker by raw channel name when it has no stain label", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "FSC-A")
  )

  expect_true("FSC-A" %in% names(result))
  expect_type(result[["FSC-A"]], "double")
})

test_that("facs_read_fcs_gated() downsamples to max_events per sample, reproducibly with seed", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result_a <- facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = c("CD4", "CD45"),
    max_events = 500,
    seed       = 42
  )
  result_b <- facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = c("CD4", "CD45"),
    max_events = 500,
    seed       = 42
  )

  expect_equal(nrow(result_a), 6L * 500L)
  expect_identical(result_a, result_b)
})

test_that("facs_read_fcs_gated() auto-derives fcs_dir from the wsp filename", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result_auto <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4")
  )
  result_explicit <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4"),
    fcs_dir   = fcs_dir
  )

  expect_identical(result_auto, result_explicit)
})

test_that("facs_read_fcs_gated() loads only the samples in a non-default FlowJo group", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4"),
    group     = "kidney"
  )

  expect_equal(dplyr::n_distinct(result$file_name), 2L)
  expect_true(all(grepl("kidney", result$file_name)))
})
```

- [ ] **Step 2: Run the full test file**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_read_fcs.R")'`
Expected: All 9 tests pass (or all skip together if the fixture is absent).

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-facs_read_fcs.R
git commit -m "test: cover channel-name matching, downsampling, fcs_dir, and group in facs_read_fcs_gated()"
```

---

## Task 4: Final check, CLAUDE.md housekeeping

**Files:**
- Modify: `CLAUDE.md` (only if Step 2 finds new NOTEs)

**Interfaces:**
- Consumes: `facs_read_fcs_gated()`, `resolve_markers_()`, `read_one_sample_()` (Tasks 1-3).

- [ ] **Step 1: Regenerate documentation**

Run: `Rscript -e 'devtools::document()'`
Expected: `man/facs_read_fcs_gated.Rd` created, `NAMESPACE` gains `export(facs_read_fcs_gated)`. No roxygen errors/warnings.

- [ ] **Step 2: Run the full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests pass, including the 9 in `test-facs_read_fcs.R` and every pre-existing test in the other `test-*.R` files.

- [ ] **Step 3: Run `devtools::check()` and record any new NOTEs**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. `R/facs_read_fcs.R` does not use bare tidy-eval column references (no bare `dplyr`-verb column names — all column access is via `[[`/matrix indexing), so no new "no visible binding for global variable" NOTEs are expected for it specifically. If `devtools::check()` nonetheless reports one (e.g. for `k` inside the `for` loop, or something dependency-related from the new Bioconductor packages), record the exact text.

- [ ] **Step 4: Update CLAUDE.md's "Known check output" section, only if Step 3 found something new**

If Step 3 was clean (0 errors, 0 warnings, no new NOTEs beyond the pre-existing ones already listed in CLAUDE.md), skip this step entirely.

If new NOTEs did appear, open `CLAUDE.md`, find the end of the "Known check output" bulleted list (last bullet currently starts with `` `pivot_organ_weights_long_` ``), and add a new bullet immediately after it describing exactly what appeared, following the phrasing pattern of the existing bullets (name the function/file, quote the exact variable names reported, and explain why it's a false positive or otherwise benign).

- [ ] **Step 5: Final full test suite + check run**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests pass.

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings, NOTEs matching exactly what's now documented in CLAUDE.md.

- [ ] **Step 6: Commit**

```bash
git add man/facs_read_fcs_gated.Rd NAMESPACE CLAUDE.md
git commit -m "docs: regenerate docs for facs_read_fcs_gated() and note any new check NOTEs"
```

(If Step 4 was skipped because there was nothing new to document, drop `CLAUDE.md` from the `git add`.)

---

## Self-Review Notes

- **Spec coverage:** Every section of `docs/superpowers/specs/2026-07-03-facs-read-fcs-gated-design.md` is covered — Architecture/dependencies (Task 1 Step 1), Public API/signature (Task 1 Step 4), Output shape (Task 1 Steps 4-5), Processing steps (Task 1 Step 4), Error handling table: missing gate + missing keyword (Task 2), unmatched marker (Task 2), missing fcs file (not separately tested — `flowjo_to_gatingset()`'s own error on a bad `path` was verified manually during planning to propagate as a clear `simpleError`, and the design only requires it to error, not any specific message, so no dedicated test was added — matches how `facs_calc_pct_of()`'s ambiguous-match `stop()` isn't independently fixture-tested either), Testing (fixture requirements documented in Task 1 Step 3, all listed test cases covered across Tasks 1-3), Known implementation risk (resolved during planning itself — see "Verified-against-real-fixture facts" above — so no further task needed).
- **Placeholder scan:** No TBD/TODO; every step has complete runnable code, exact fixture-derived facts (channel names, stain labels, group membership, event counts), and exact commands with expected output.
- **Type consistency:** `facs_read_fcs_gated(wsp_path, gate_path, markers, keywords = NULL, fcs_dir = NULL, group = "All Samples", max_events = NULL, seed = NULL)` signature is introduced once in Task 1 and used identically (same argument names, same order where used) in every test across Tasks 1-3. `resolve_markers_()` and `read_one_sample_()`'s signatures are defined and consumed only within `R/facs_read_fcs.R` itself, no cross-task drift possible since both are written in a single step.
