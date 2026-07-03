# tests/fixtures/

Fixture files for JanisHelpers tests.

## facs_import_wsp() — pending

`facs_import_wsp()` requires a real `.wsp` FlowJo workspace file.
No mock of `fcexpr` is allowed (see CLAUDE.md).

### Required fixture: `minimal.wsp`

A minimal FlowJo 10 workspace containing:

- **At least 2 FCS files** so that multi-file behaviour can be tested
- **At least 2 gating hierarchy levels** (e.g. `Lymphocytes/CD4+`)
  so that `PopulationFullPath` and `Population` (leaf) are distinct
- **At least 1 custom keyword** (e.g. `$USER` or a user-defined keyword
  like `mouse_ID`) to test the `keywords` argument
- **At least 1 channel with a stain label** (`$PnS`) so that the
  `Median_<stain>` stat label path is exercised when `r_stats = TRUE`

### How to create the fixture

1. Open FlowJo 10 with a small pilot experiment (2–4 FCS files).
2. Add a minimal gate hierarchy (e.g. Lymphocytes → CD4+).
3. Assign a custom keyword to each sample (Workspace → Keywords).
4. Save the workspace as `tests/fixtures/minimal.wsp`.
5. Add `tests/fixtures/minimal.wsp` to `.gitignore` if it contains
   patient-identifiable metadata; otherwise commit it directly.

### Planned smoke test (tests/testthat/test-facs_import.R)

```r
test_that("facs_import_wsp returns long tibble with required columns", {
  skip_if_not(file.exists(testthat::test_path("fixtures", "minimal.wsp")))
  dat <- facs_import_wsp(testthat::test_path("fixtures", "minimal.wsp"))
  expect_s3_class(dat, "tbl_df")
  expect_true(all(c("FileName", "PopulationFullPath", "Population",
                    "metric", "value") %in% names(dat)))
})

test_that("facs_import_wsp attaches keyword columns", {
  skip_if_not(file.exists(testthat::test_path("fixtures", "minimal.wsp")))
  dat <- facs_import_wsp(
    testthat::test_path("fixtures", "minimal.wsp"),
    keywords = "mouse_ID"
  )
  expect_true("mouse_ID" %in% names(dat))
})
```

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
