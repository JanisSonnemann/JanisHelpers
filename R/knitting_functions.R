#' Function for automatic knitting of dated files in subdirectory
#'
#' @param input unchanged
#' @param path change to name of desired subdirectory
#' @param ... further arguments
#'
#' @return
#' @export
#'
#' @examples
knit_multiple_dated <- function(input, path, ...) {
  rmarkdown::render(
    input,
    output_dir = paste0("doc/output/", path),
    output_file = paste0(
      xfun::sans_ext(input), '-', Sys.Date()
    ),
    output_format = "all",
    envir = globalenv()
  )
}
