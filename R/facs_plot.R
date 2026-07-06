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

#' Plot per-sample cluster abundance by group
#'
#' @description
#' Renders a \code{ggplot2} boxplot with jittered points of
#' \code{facs_calc_cluster_freq()}'s per-sample \code{fraction}, grouped by
#' \code{group_col} and faceted by \code{cluster_col}. Optionally restricts
#' the facets shown to clusters flagged significant by
#' \code{facs_test_cluster_abundance()}.
#'
#' @param freq_data tibble shaped like \code{facs_calc_cluster_freq()}'s
#'   output: must contain \code{fraction}, \code{group_col}, and
#'   \code{cluster_col}.
#' @param test_result optional tibble from \code{facs_test_cluster_abundance()};
#'   used only when \code{significant_only = TRUE}, to determine which
#'   \code{cluster_col} values have any \code{p_adj <= p_adj_threshold}.
#' @param group_col character; column in \code{freq_data} to plot on the
#'   x-axis.
#' @param cluster_col character; column in \code{freq_data} (and, if
#'   supplied, \code{test_result}) to facet by. Default \code{"metacluster"}.
#' @param significant_only logical; if \code{TRUE}, restrict facets to
#'   clusters with any \code{p_adj <= p_adj_threshold} in \code{test_result}.
#'   Default \code{FALSE}.
#' @param p_adj_threshold numeric; adjusted p-value cutoff used when
#'   \code{significant_only = TRUE}. Default \code{0.05}.
#'
#' @returns A \code{ggplot} object (not printed or saved). Errors if
#'   \code{group_col} or \code{cluster_col} are not columns in
#'   \code{freq_data}.
#' @export
#'
#' @examples
#' \dontrun{
#'   freq <- facs_calc_cluster_freq(facs_cluster_flowsom(facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   )))
#'   facs_plot_cluster_abundance(freq, group_col = "group")
#' }
facs_plot_cluster_abundance <- function(freq_data,
                                         test_result = NULL,
                                         group_col,
                                         cluster_col = "metacluster",
                                         significant_only = FALSE,
                                         p_adj_threshold = 0.05) {
  if (!group_col %in% names(freq_data)) {
    stop(glue::glue("`group_col` ('{group_col}') not found in `freq_data`."))
  }
  if (!cluster_col %in% names(freq_data)) {
    stop(glue::glue("`cluster_col` ('{cluster_col}') not found in `freq_data`."))
  }

  if (isTRUE(significant_only) && !is.null(test_result)) {
    significant_clusters <- test_result |>
      dplyr::filter(p_adj <= p_adj_threshold) |>
      dplyr::pull(dplyr::all_of(cluster_col)) |>
      unique()

    freq_data <- dplyr::filter(freq_data, .data[[cluster_col]] %in% significant_clusters)
  }

  ggplot2::ggplot(freq_data, ggplot2::aes(x = .data[[group_col]], y = fraction)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.6) +
    ggplot2::facet_wrap(stats::as.formula(paste0("~", cluster_col))) +
    ggplot2::labs(x = group_col, y = "fraction") +
    ggplot2::theme_minimal()
}
