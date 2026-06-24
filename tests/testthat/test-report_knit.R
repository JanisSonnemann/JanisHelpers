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

test_that("report_knit_exp passes input path and format through to render", {
  input <- exp_input()

  captured <- list()
  local_mocked_bindings(
    render = function(...) { captured <<- list(...); invisible(NULL) },
    .package = "rmarkdown"
  )

  report_knit_exp(input = input, format = "html_document")

  expect_equal(captured[[1]], input)               # positional arg: input
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
