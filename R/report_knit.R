#' Function for automatic knitting of dated files in subdirectory
#'
#' @param input unchanged
#' @param path change to name of desired subdirectory
#' @param format specify desired files to output in quotes: "html_document", "pdf_document", "all"
#' @param ... further arguments
#'
#' @return
#' @export
#'
#'
#'
#' @examples
report_knit_dated <- function(input, path, format, ...) {
  rmarkdown::render(
    input,
    output_dir = paste0("analysis/0-output/", path),
    output_file = paste0(
      xfun::sans_ext(input), '-', Sys.Date()
    ),
    output_format = format,
    envir = globalenv()
  )
}


#' New function for automatic knitting of dated files in subdirectory based on new experiment based subdirectories
#'
#' @param input unchanged
#' @param format specify desired files to output in quotes: "html_document", "pdf_document", "all"
#' @param ... further arguments
#'
#' @returns
#' @export
#'
#' @examples
report_knit_exp <- function(input, format, ...) {
  rmarkdown::render(
    input,
    output_dir = paste0(dirname(dirname(normalizePath(input, mustWork = FALSE))), "/analysis/output/", xfun::sans_ext(basename(input))),
    output_file = paste0(
      xfun::sans_ext(input), '-', Sys.Date()
    ),
    output_format = format,
    envir = globalenv()
  )
}


#' function to generate wide html
#'
#' @param input file to be knitted
#' @param ... further arguments
#'
#' @returns
#' @export
#'
#' @examples
report_knit_wide <- function(input, ...) {
  css_path <- system.file(
    "resources", "wide-output.css",
    package = "JanisHelpers"
  )

  # Inject the CSS into the format object
  format <- rmarkdown::html_document(
    css = css_path,
    toc = TRUE,
    toc_float = TRUE,
    code_folding = "hide"
    )

  # Continue with rendering
  rmarkdown::render(
    input,
    output_dir = file.path(dirname(dirname(normalizePath(input, mustWork = FALSE))), "analysis/output", xfun::sans_ext(basename(input))),
    output_file = paste0(xfun::sans_ext(input), "-", Sys.Date()),
    output_format = format,
    envir = globalenv(),
    ...
  )
}
