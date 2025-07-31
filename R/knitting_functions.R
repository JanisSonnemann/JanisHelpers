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
knit_multiple_dated <- function(input, path, format, ...) {
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
knit_exp_structure <- function(input, format, ...) {
  rmarkdown::render(
    input,
    output_dir = paste0(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)), "/analysis/output/", xfun::sans_ext(basename(input))),
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
#' @param format defaults to wide html_document
#' @param ... further arguments
#'
#' @returns
#' @export
#'
#' @examples
knit_wide_html <- function(input, ...) {
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
    output_dir = file.path(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)), "analysis/output", xfun::sans_ext(basename(input))),
    output_file = paste0(xfun::sans_ext(input), "-", Sys.Date()),
    output_format = format,
    envir = globalenv(),
    ...
  )
}


#' Panel knitting
#' @description create panel
#'
#'
#' @param input unchanged
#' @param path path of experiment
#' @param format pdf_document
#' @param ... further parameters
#'
#' @return
#' @export
#'
#' @examples
knit_panel <- function(input, path, format, ...) {
  rmarkdown::render(
    input,
    output_dir = paste0("data/", path),
    output_file = paste0(
      xfun::sans_ext(input), '-', Sys.Date()
    ),
    output_format = format,
    envir = globalenv()
  )
}
