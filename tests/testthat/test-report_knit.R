# Helper: build a realistic nested path without creating the file
exp_input <- function() {
  root <- file.path(tempdir(), paste0("JHtest_", Sys.getpid()))
  file.path(root, "exp1", "scripts", "report.Rmd")
}

# ── report_knit_exp ───────────────────────────────────────────────────────────

test_that("report_knit_exp passes correct output_dir without rstudioapi", {
  input <- exp_input()

  captured <- list()
  local_mocked_bindings(
    render = function(...) { captured <<- list(...); invisible(NULL) },
    .package = "rmarkdown"
  )

  report_knit_exp(input = input, format = "html_document")

  expected_dir <- paste0(
    dirname(dirname(normalizePath(input, mustWork = FALSE))),
    "/analysis/output/",
    xfun::sans_ext(basename(input))
  )
  expect_equal(captured$output_dir, expected_dir)
})

test_that("report_knit_exp passes format argument through to render", {
  input <- exp_input()

  captured <- list()
  local_mocked_bindings(
    render = function(...) { captured <<- list(...); invisible(NULL) },
    .package = "rmarkdown"
  )

  report_knit_exp(input = input, format = "html_document")

  expect_true(length(captured) > 0)
  expect_equal(captured$output_format, "html_document")
})

# ── report_knit_wide ──────────────────────────────────────────────────────────

test_that("report_knit_wide passes correct output_dir without rstudioapi", {
  input <- exp_input()

  captured <- list()
  local_mocked_bindings(
    render = function(...) { captured <<- list(...); invisible(NULL) },
    .package = "rmarkdown"
  )

  report_knit_wide(input = input)

  expected_dir <- file.path(
    dirname(dirname(normalizePath(input, mustWork = FALSE))),
    "analysis/output",
    xfun::sans_ext(basename(input))
  )
  expect_equal(captured$output_dir, expected_dir)
})

test_that("report_knit_wide injects an html_document format object", {
  input <- exp_input()

  captured <- list()
  local_mocked_bindings(
    render = function(...) { captured <<- list(...); invisible(NULL) },
    .package = "rmarkdown"
  )

  report_knit_wide(input = input)

  expect_true(inherits(captured$output_format, "rmarkdown_output_format"))
})

# -- report_knit_facs_qc ------------------------------------------------------

test_that("report_knit_facs_qc renders into qc/ and flags exclusions", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  files <- tibble::tibble(file_name = c("a.fcs", "b.fcs"),
                          mouse_ID = c("25-7-1", "25-7-2"), tissue = "spleen")
  prep <- suppressMessages(
    facs_prep_save(fake_prep_data(files), d$wsp_path, "overview", d$panel_dir)
  )
  ex <- tibble::tibble(
    mouse_ID = c("25-7-2", "25-7-9"), reason = c("clots", "typo-mouse")
  )
  out <- report_knit_facs_qc(prep, exclusions = ex)
  expect_equal(
    normalizePath(out),
    normalizePath(file.path(d$panel_dir, "qc", "overview_25-7_20261101_qc.html"))
  )
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "b.fcs")
  expect_match(html, "clots")
  expect_match(html, "match no file")
  expect_match(html, "typo-mouse")
})

test_that("report_knit_facs_qc lists exclusion rows for other panels", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  files <- tibble::tibble(file_name = "a.fcs", mouse_ID = "25-7-1", tissue = "spleen")
  prep <- suppressMessages(
    facs_prep_save(fake_prep_data(files), d$wsp_path, "overview", d$panel_dir)
  )
  ex <- tibble::tibble(mouse_ID = "25-7-1", panel = "Overvew", reason = "typo-panel")
  html <- paste(readLines(report_knit_facs_qc(prep, ex), warn = FALSE), collapse = "\n")
  expect_match(html, "other panels")
  expect_match(html, "Overvew")
})

render_qc_html_ <- function(prep_path, exclusions = NULL) {
  paste(readLines(report_knit_facs_qc(prep_path, exclusions), warn = FALSE),
        collapse = "\n")
}

test_that("report_knit_facs_qc renders QC checks with flags", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  qd <- fake_qc_data()
  prep <- suppressMessages(facs_prep_save(
    qd, d$wsp_path, "overview", d$panel_dir,
    group_col = "grp", batch_col = "day",
    expected_mice = c(unique(qd$mouse_ID), "m-99"),
    qc_mfi_populations = "T"
  ))
  html <- render_qc_html_(prep, tibble::tibble(mouse_ID = "m-5", reason = "clots"))
  for (section in c("Completeness", "Fraction of parent", "Raw counts",
                    "Batch effects", "MFI",
                    "Flagged samples", "Batch")) {
    expect_match(html, section)
  }
  expect_match(html, "m-3")
  expect_match(html, "m-99")
  expect_match(html, "<img")
})

test_that("report_knit_facs_qc reports empty values, missing stats and keywords", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  qd <- fake_qc_data()
  qd$tissue[qd$mouse_ID %in% c("m-1", "m-2", "m-9", "m-10")] <- "liver"
  per_g <- dplyr::filter(qd, metric == "count") |>
    dplyr::mutate(metric = "count_per_g",
                  value = dplyr::if_else(tissue == "spleen", NA_real_, value))
  qd <- dplyr::bind_rows(qd, per_g) |>
    dplyr::filter(!(mouse_ID == "m-4" & metric == "median_CD25"))
  qd$grp[qd$mouse_ID == "m-5"] <- NA
  comp <- dplyr::mutate(qd[1, ], file_name = "comp.fcs", mouse_ID = NA,
                        tissue = NA)
  prep <- suppressMessages(facs_prep_save(
    dplyr::bind_rows(qd, comp), d$wsp_path, "overview", d$panel_dir,
    group_col = "grp", batch_col = "day"
  ))
  html <- render_qc_html_(prep)
  expect_match(html, "Empty values and missing stats")
  expect_match(html, "count_per_g \\(spleen\\)")
  expect_match(html, "all empty")
  expect_match(html, "missing in: m-4")
  expect_match(html, "comp.fcs")
  expect_match(html, "empty group / batch values")
})

test_that("report_knit_facs_qc renders an old prep without qc settings", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  path <- suppressMessages(
    facs_prep_save(fake_qc_data(), d$wsp_path, "overview", d$panel_dir)
  )
  prep <- readRDS(path)
  prep$info$qc <- NULL
  saveRDS(prep, path)
  html <- render_qc_html_(path)
  expect_match(html, "rerun the prep script")
  expect_match(html, "Expected mice not set")
})

test_that("report_knit_facs_qc notes when no MFI populations are selected", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  prep <- suppressMessages(facs_prep_save(
    fake_qc_data(), d$wsp_path, "overview", d$panel_dir,
    group_col = "grp", batch_col = "day"
  ))
  expect_match(render_qc_html_(prep), "No MFI populations selected")
})

test_that("report_knit_facs_qc renders the QC plots without warnings", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  prep <- suppressMessages(facs_prep_save(
    fake_qc_data(), d$wsp_path, "overview", d$panel_dir,
    group_col = "grp", batch_col = "day", qc_mfi_populations = "T"
  ))
  expect_no_warning(report_knit_facs_qc(prep))
})

test_that("report_knit_facs_qc renders without exclusions and shows clashes", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_prep_dirs()
  files <- tibble::tibble(file_name = c("a.fcs", "b.fcs"),
                          mouse_ID = "25-7-1", tissue = "spleen")
  prep <- suppressWarnings(suppressMessages(
    facs_prep_save(fake_prep_data(files), d$wsp_path, "overview", d$panel_dir)
  ))
  out <- report_knit_facs_qc(prep)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "share one mouse_ID")
  expect_match(html, "QC checks")
})

# -- report_knit_elisa_qc -----------------------------------------------------

test_that("report_knit_elisa_qc renders into qc/ and names the sections", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- withr::local_tempdir()
  path <- file.path(d, "elisa_curated.xlsx")
  file.copy(test_path("..", "fixtures", "elisa_curated.xlsx"), path)
  data <- suppressMessages(elisa_read_curated(path, tissue = "spleen")) |>
    dplyr::mutate(grp = ifelse(mouse_ID == "M1", "ctrl", "treated"))
  prep <- suppressMessages(elisa_prep_save(data, path, "restim",
                                           group_col = "grp",
                                           expected_mice = c("M1", "M2", "M9")))
  ex <- tibble::tibble(mouse_ID = c("M2", "M7"), reason = c("sick", "typo"))
  out <- expect_no_warning(report_knit_elisa_qc(prep, exclusions = ex))
  expect_equal(normalizePath(out),
               normalizePath(file.path(d, "qc", "elisa_curated_qc.html")))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  for (s in c("Overview", "Completeness", "Exclusions", "Empty values", "Range",
              "Replicate agreement", "Outliers", "Flagged samples", "sick",
              "M7", "M9")) {
    expect_match(html, s, fixed = TRUE)
  }
})

# -- report_knit_histo_qc -----------------------------------------------------

test_that("report_knit_histo_qc renders into qc/ and names the sections", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- withr::local_tempdir()
  path <- file.path(d, "histo_curated.xlsx")
  file.copy(test_path("..", "fixtures", "histo_curated.xlsx"), path)
  data <- suppressMessages(histo_read_curated(path, tissue = "kidney")) |>
    dplyr::mutate(grp = ifelse(mouse_ID == "M1", "ctrl", "treated"))
  prep <- suppressMessages(histo_prep_save(data, path, "PAS", group_col = "grp",
                                           rules = "crescents <= gloms"))
  ex <- tibble::tibble(mouse_ID = c("M2", "M7"), reason = c("sick", "typo"))
  out <- expect_no_warning(report_knit_histo_qc(prep, exclusions = ex))
  expect_equal(normalizePath(out),
               normalizePath(file.path(d, "qc", "histo_curated_qc.html")))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  for (s in c("Overview", "Completeness", "Exclusions", "Empty values",
              "Sanity rules", "Outliers", "Flagged samples", "sick", "M7", "M4")) {
    expect_match(html, s, fixed = TRUE)
  }
})

test_that("report_knit_histo_qc says so when no rules are set", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- withr::local_tempdir()
  path <- file.path(d, "histo_curated.xlsx")
  file.copy(test_path("..", "fixtures", "histo_curated.xlsx"), path)
  prep <- suppressMessages(histo_prep_save(
    suppressMessages(histo_read_curated(path, tissue = "kidney")), path, "PAS"))
  out <- report_knit_histo_qc(prep)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "No rules set", fixed = TRUE)
})
