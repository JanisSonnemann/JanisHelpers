#' Plot a UMAP embedding
#'
#' @description
#' Renders a \code{ggplot2} scatter plot of a UMAP embedding (e.g.
#' \code{facs_reduce_umap()}'s output), colored by a chosen column and
#' optionally faceted. The color scale is chosen automatically: a
#' continuous viridis scale for a \code{double}-typed \code{color_by}
#' column (e.g. a marker's expression), a discrete viridis scale otherwise
#' (e.g. \code{metacluster}, \code{group}). Both are colorblind-safe and
#' perceptually uniform, and ship inside \code{ggplot2} -- no additional
#' dependency beyond \code{ggplot2} itself.
#'
#' @param data tibble shaped like \code{facs_reduce_umap()}'s output: must
#'   contain \code{UMAP1} and \code{UMAP2}.
#' @param color_by character; column in \code{data} to color points by.
#'   Default \code{"metacluster"}.
#' @param facet_by character or \code{NULL} (default); column in
#'   \code{data} to facet panels by via \code{ggplot2::facet_wrap()}.
#'   \code{NULL} produces a single panel.
#'
#' @returns A \code{ggplot} object (not printed or saved). Errors if
#'   \code{data} does not contain both \code{UMAP1} and \code{UMAP2}, if
#'   \code{color_by} is not a column in \code{data}, or if \code{facet_by}
#'   is supplied and is not a column in \code{data}.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_reduce_umap(facs_cluster_flowsom(facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   )))
#'   facs_plot_umap(dat, color_by = "metacluster", facet_by = "group")
#' }
facs_plot_umap <- function(data, color_by = "metacluster", facet_by = NULL) {
  if (!all(c("UMAP1", "UMAP2") %in% names(data))) {
    stop("`data` must contain both `UMAP1` and `UMAP2` columns (see facs_reduce_umap()).")
  }
  if (!color_by %in% names(data)) {
    stop(glue::glue("`color_by` ('{color_by}') not found in `data`."))
  }
  if (!is.null(facet_by) && !facet_by %in% names(data)) {
    stop(glue::glue("`facet_by` ('{facet_by}') not found in `data`."))
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = UMAP1, y = UMAP2, color = .data[[color_by]])) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2", color = color_by) +
    ggplot2::theme_minimal()

  p <- if (is.double(data[[color_by]])) {
    p + ggplot2::scale_color_viridis_c()
  } else {
    p + ggplot2::scale_color_viridis_d()
  }

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", facet_by)))
  }

  p
}
