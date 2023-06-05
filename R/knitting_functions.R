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
    output_dir = paste0("doc/output/", path),
    output_file = paste0(
      xfun::sans_ext(input), '-', Sys.Date()
    ),
    output_format = format,
    envir = globalenv()
  )
}
