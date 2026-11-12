#' Function for automatic knitting of dated files in subdirectory
#'
#' @param input unchanged
#' @param path change to name of desired subdirectory
#' @param format specify desired files to output in quotes: "html_document", "pdf_document", "all"
#' @param ... further arguments
#'
#' @returns path to the rendered output file, returned by `rmarkdown::render()`
#' @export
#'
#' @examples
#' \dontrun{
#'   report_knit_dated(
#'     input  = "analysis/flow-report.Rmd",
#'     path   = "flow-experiment-1",
#'     format = "html_document"
#'   )
#' }
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
#' @returns path to the rendered output file, returned by `rmarkdown::render()`
#' @export
#'
#' @examples
#' \dontrun{
#'   report_knit_exp(
#'     input  = "scripts/flow-report.Rmd",
#'     format = "html_document"
#'   )
#' }
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
#' @returns path to the rendered output file, returned by `rmarkdown::render()`
#' @export
#'
#' @examples
#' \dontrun{
#'   report_knit_wide(input = "scripts/flow-report.Rmd")
#' }
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


#' Render the FACS QC report for one prep result
#'
#' Renders the `facs-qc-template` from a prep result saved by
#' [facs_prep_save()] into the `qc/` folder beside its `processed/` folder:
#' `<panel_dir>/qc/<wsp name>_qc.html` (one per workspace version, overwritten
#' on rerun). The report lists every file, flags files hit by the exclusions
#' sheet, warns about keyword clashes and exclusion rows that match no file,
#' and has placeholder sections for QC checks.
#'
#' @param prep_path Path to the prep `.rds` returned by [facs_prep_save()].
#' @param exclusions Optional tibble -- the `exclusions` sheet of the meta
#'   workbook (e.g. `meta_read(path)$exclusions`), see
#'   [facs_match_exclusions()]. `NULL` = no exclusions.
#' @returns The path of the rendered `.html` file.
#' @examples
#' \dontrun{
#' prep_path <- facs_prep_save(data, wsp, "overview", panel_dir)
#' report_knit_facs_qc(prep_path, meta_read(meta_path)$exclusions)
#' }
#' @export
report_knit_facs_qc <- function(prep_path, exclusions = NULL) {
  prep_path <- normalizePath(prep_path, mustWork = TRUE)
  qc_dir <- file.path(dirname(dirname(prep_path)), "qc")
  dir.create(qc_dir, showWarnings = FALSE)
  # Intermediates go to a temp dir: knitting on OneDrive is slow and flaky.
  work_dir <- tempfile("facs_qc_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  rmarkdown::render(
    system.file("rmarkdown", "templates", "facs-qc-template", "skeleton",
                "skeleton.Rmd", package = "JanisHelpers"),
    output_dir = qc_dir,
    output_file = sub("_prep\\.rds$", "_qc.html", basename(prep_path)),
    params = list(prep = prep_path, exclusions = exclusions),
    intermediates_dir = work_dir,
    knit_root_dir = work_dir,
    envir = new.env(),
    quiet = TRUE
  )
}
