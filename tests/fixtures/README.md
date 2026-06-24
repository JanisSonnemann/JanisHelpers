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
