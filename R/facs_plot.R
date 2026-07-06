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

#' Plot a cluster marker-expression heatmap
#'
#' @description
#' Renders a \code{ggplot2} tile heatmap of \code{facs_calc_cluster_marker_
#' medians()}'s per-cluster median marker expression -- the standard
#' cluster-annotation view for naming a metacluster (e.g. "Treg-like") by
#' its marker profile before interpreting a
#' \code{facs_test_cluster_abundance()} hit. Uses a diverging red/blue fill
#' scale (not the viridis scales used elsewhere in this package) since a
#' scale centered at a meaningful midpoint -- 0, whether that's a z-score's
#' mean or the transformed marker scale's own background/negative boundary
#' -- is the field-standard way to read this kind of heatmap.
#'
#' @param marker_medians tibble shaped like \code{facs_calc_cluster_marker_
#'   medians()}'s output: must contain \code{cluster_col}, \code{marker},
#'   and \code{median}.
#' @param cluster_col character; column in \code{marker_medians} to plot on
#'   the y-axis. Default \code{"metacluster"}.
#' @param scale character; one of \code{"zscore"} (default, z-scores
#'   \code{median} per \code{marker} across clusters before plotting) or
#'   \code{"raw"} (plots \code{median} unscaled).
#'
#' @returns A \code{ggplot} object (not printed or saved). Errors if
#'   \code{cluster_col}, \code{"marker"}, or \code{"median"} are not columns
#'   in \code{marker_medians}.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_cluster_flowsom(facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   ), seed = 1)
#'   medians <- facs_calc_cluster_marker_medians(dat)
#'   facs_plot_cluster_heatmap(medians)
#' }
facs_plot_cluster_heatmap <- function(marker_medians,
                                      cluster_col = "metacluster",
                                      scale = c("zscore", "raw")) {
  scale <- match.arg(scale)

  if (!cluster_col %in% names(marker_medians)) {
    stop(glue::glue("`cluster_col` ('{cluster_col}') not found in `marker_medians`."))
  }
  if (!"marker" %in% names(marker_medians)) {
    stop("`marker` column not found in `marker_medians`.")
  }
  if (!"median" %in% names(marker_medians)) {
    stop("`median` column not found in `marker_medians`.")
  }

  # `scale` (the argument) is a character string, not a function -- R's
  # call-position lookup still resolves `base::scale()` correctly here (see
  # the plan's "Verified mechanics" section), but it's spelled out
  # explicitly below for a human reader given the name collision.
  plot_data <- marker_medians |>
    dplyr::group_by(marker) |>
    dplyr::mutate(fill_value = if (scale == "zscore") {
      as.numeric(base::scale(median))
    } else {
      median
    }) |>
    dplyr::ungroup()

  ggplot2::ggplot(plot_data, ggplot2::aes(x = marker, y = .data[[cluster_col]], fill = fill_value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::labs(x = "marker", y = cluster_col, fill = if (scale == "zscore") "z-score" else "median") +
    ggplot2::theme_minimal()
}
