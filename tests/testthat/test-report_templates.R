purl_template_ <- function(name) {
  path <- system.file("rmarkdown", "templates", name, "skeleton", "skeleton.Rmd",
                      package = "JanisHelpers")
  expect_true(nzchar(path) && file.exists(path))
  out <- tempfile(fileext = ".R")
  knitr::purl(path, output = out, quiet = TRUE, documentation = 0)
  code <- readLines(out)
  expect_no_error(parse(text = code))
  paste(code, collapse = "\n")
}

test_that("analysis-db template parses, loads read-only and builds dat + pops", {
  code <- purl_template_("analysis-db")
  expect_match(code, "dat <- ", fixed = TRUE)
  expect_match(code, "pops <- ", fixed = TRUE)
  expect_match(code, "read_only = TRUE", fixed = TRUE)
  expect_match(code, "db_load_experiment(", fixed = TRUE)
  expect_match(code, "JanisHelpers.group_colors", fixed = TRUE)
  expect_no_match(code, "rstudioapi", fixed = TRUE)
})

test_that("analysis-wsp template parses and builds the same dat + pops shape", {
  code <- purl_template_("analysis-wsp")
  expect_match(code, "dat <- ", fixed = TRUE)
  expect_match(code, "pops <- ", fixed = TRUE)
  expect_match(code, "facs_read_wsp(", fixed = TRUE)
  expect_match(code, "mouse_id = mouse_ID", fixed = TRUE)
  expect_match(code, "JanisHelpers.group_colors", fixed = TRUE)
  expect_no_match(code, "rstudioapi", fixed = TRUE)
})
