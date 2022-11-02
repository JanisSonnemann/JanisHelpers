#' Knits multiple files into subdirectory and adds date to filename
#'
#' @export
knit_with_date <- function(input, ...) {
  rmarkdown::render(
    input,
    output_file = paste0(
      xfun::sans_ext(input), '-', Sys.Date(), '.'
    ),
    output_dir = "output",
    output_format = "all",
    envir = globalenv()
  )
}
