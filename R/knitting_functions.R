#' Custom Knit function for RStudio
#'
#' @export
knit_with_date <- function(input, ...) {
  rmarkdown::render(
    input,
    output_file = paste0(
      xfun::sans_ext(input), '-', Sys.Date(), '.'
    ),
    output_dir = "output",
    envir = globalenv()
  )
}
