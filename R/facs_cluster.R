# FlowSOM behavior -- verified against tests/fixtures/Treg.wsp
# (installed: FlowSOM 2.16.0, ConsensusClusterPlus 1.72.0, current
# Bioconductor 3.21 release of both as of this writing)
#
# FlowSOM::FlowSOM()'s `xdim`/`ydim` are not named formals -- per its own
# roxygen docs they are absorbed by `...` and forwarded to the internal
# SOM-building step, so passing them by name still reaches it. It returns
# the fsom list directly (not the two-element list its own docs describe),
# with fsom$metaclustering already attached -- GetClusters()/
# GetMetaclusters() take that returned object as-is.
#
# Three reproducible upstream bugs constrain which grid/marker/metacluster
# combinations actually work, independent of anything in this file:
#
# 1. `n_metaclusters = 2` always throws ("argument of length 0" /
#    "Argument der Laenge 0") from ConsensusClusterPlus's internal
#    clusterTrackingPlot(), for *any* data -- reproduces directly via
#    ConsensusClusterPlus::ConsensusClusterPlus(matrix(rnorm(40), 10),
#    maxK = 2, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
#    plot = "pdf", clusterAlg = "hc", distance = "euclidean"). Its
#    consensus-color matrix only ever accumulates one row when maxK = 2,
#    and R's default `[` silently drops that 1-row matrix to a bare
#    vector, so nrow() on it is NULL. n_metaclusters must be >= 3 (and
#    < grid_xdim * grid_ydim, per bug 3 below and the check below) until
#    upstream fixes this.
# 2. Selecting exactly one marker column always throws ("argument must be
#    coercible to non-negative integer") from FlowSOM::BuildSOM()'s
#    `fsom$data[, colsToUse]`, which has no `drop = FALSE` -- a
#    single-column `colsToUse` collapses to a bare vector, and the
#    subsequent seq_len(ncol(data)) inside SOM() fails on its NULL
#    ncol(). `markers` (explicit or defaulted) must resolve to >= 2
#    columns until upstream fixes this.
# 3. `n_metaclusters == grid_xdim * grid_ydim` (exact equality with the
#    SOM's node count) always throws from ConsensusClusterPlus's internal
#    cutree() call, for *any* data. ConsensusClusterPlus's default
#    `pItem = 0.9` means every resampling iteration works on a subset of
#    the N SOM nodes, never all N of them, so its resampled hierarchical
#    tree has fewer than N leaves. Asking cutree() for exactly N clusters
#    then fails because that many clusters don't exist in the resampled
#    tree. n_metaclusters must be strictly less than grid_xdim * grid_ydim
#    (not just `<=`, per the check below) until upstream
#    fixes this.
#
# All three bugs are now guarded against, below: this function rejects
# n_metaclusters < 3 (bug 1), a resolved `markers` of length < 2 (bug 2),
# and n_metaclusters >= grid_xdim * grid_ydim (bug 3) outright, so a caller
# who would otherwise hit any of these opaque FlowSOM/ConsensusClusterPlus
# crashes instead gets a clear error from this package's own validation.

#' Cluster single-cell events with FlowSOM
#'
#' @description
#' Runs FlowSOM (self-organizing map + consensus metaclustering) on a
#' per-event tibble (e.g. \code{facs_read_fcs_gated()}'s output), appending
#' a raw SOM node assignment and a consensus metacluster assignment to each
#' event. Clustering is performed directly on the input's numeric scale (no
#' additional z-score normalization), matching the transformed
#' (logicle/biexponential) scale \code{facs_read_fcs_gated()} already
#' returns.
#'
#' @param data tibble shaped like \code{facs_read_fcs_gated()}'s output:
#'   one row per event, \code{dbl} marker columns, and optionally
#'   \code{chr} keyword columns.
#' @param markers character vector of column names in \code{data} to
#'   cluster on. \code{NULL} (default) uses every \code{dbl}-typed column
#'   in \code{data}. Must resolve to at least 2 columns, and, when
#'   explicitly supplied, every named column must be \code{double}-typed
#'   in \code{data}.
#' @param grid_xdim,grid_ydim integer; SOM grid dimensions. Default
#'   \code{10}/\code{10} (100 nodes).
#' @param n_metaclusters integer; target consensus metacluster count.
#'   Default \code{10}. Must be at least 3 and strictly less than
#'   \code{grid_xdim * grid_ydim}.
#' @param seed integer; if set, seeds SOM training and consensus
#'   metaclustering for reproducible assignments.
#'
#' @returns \code{data} with two columns appended: \code{cluster} (integer,
#'   raw SOM node, \code{1:(grid_xdim * grid_ydim)}) and \code{metacluster}
#'   (factor, consensus grouping, \code{1:n_metaclusters} levels). Errors if
#'   any \code{markers} name is absent from \code{data}, if an explicitly
#'   supplied \code{markers} column is not \code{double}-typed in
#'   \code{data} (listing the offending column(s)), if the resolved
#'   \code{markers} has fewer than 2 columns, if a selected marker column
#'   contains \code{NA}, if \code{n_metaclusters} is less than 3, or if
#'   \code{n_metaclusters} is not strictly less than the grid's node count.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   )
#'   facs_cluster_flowsom(dat, seed = 1)
#' }
facs_cluster_flowsom <- function(data,
                                  markers = NULL,
                                  grid_xdim = 10,
                                  grid_ydim = 10,
                                  n_metaclusters = 10,
                                  seed = NULL) {
  if (is.null(markers)) {
    is_dbl  <- purrr::map_lgl(data, is.double)
    markers <- names(data)[is_dbl]
  } else {
    missing_markers <- setdiff(markers, names(data))
    if (length(missing_markers) > 0L) {
      stop(glue::glue(
        "The following `markers` were not found in `data`: ",
        "{paste(missing_markers, collapse = ', ')}"
      ))
    }

    non_double_markers <- markers[!purrr::map_lgl(markers, function(m) is.double(data[[m]]))]
    if (length(non_double_markers) > 0L) {
      stop(glue::glue(
        "The following `markers` are not double-typed columns in `data`: ",
        "{paste(non_double_markers, collapse = ', ')}"
      ))
    }
  }

  if (length(markers) < 2L) {
    stop(glue::glue(
      "`markers` must resolve to at least 2 columns (found ",
      "{length(markers)}) -- FlowSOM::BuildSOM() crashes opaquely when ",
      "exactly one marker column is selected."
    ))
  }

  na_cols <- markers[purrr::map_lgl(markers, function(m) anyNA(data[[m]]))]
  if (length(na_cols) > 0L) {
    stop(glue::glue(
      "The following marker column(s) contain NA and cannot be clustered: ",
      "{paste(na_cols, collapse = ', ')}"
    ))
  }

  if (n_metaclusters < 3L) {
    stop(glue::glue(
      "`n_metaclusters` ({n_metaclusters}) must be at least 3 -- ",
      "ConsensusClusterPlus crashes opaquely whenever n_metaclusters < 3."
    ))
  }

  if (n_metaclusters >= grid_xdim * grid_ydim) {
    stop(glue::glue(
      "`n_metaclusters` ({n_metaclusters}) must be strictly less than the ",
      "SOM grid's node count (grid_xdim * grid_ydim = ",
      "{grid_xdim * grid_ydim})."
    ))
  }

  input <- flowCore::flowFrame(as.matrix(data[markers]))

  # See the "FlowSOM behavior" header comment at the top of this file for
  # the `xdim`/`ydim`/return-shape notes and the upstream bugs guarded
  # against above.
  fsom <- FlowSOM::FlowSOM(
    input     = input,
    colsToUse = markers,
    xdim      = grid_xdim,
    ydim      = grid_ydim,
    nClus     = n_metaclusters,
    scale     = FALSE,
    seed      = seed
  )

  data$cluster     <- as.integer(FlowSOM::GetClusters(fsom))
  data$metacluster <- factor(as.integer(FlowSOM::GetMetaclusters(fsom)))

  data
}

#' Compute per-sample cluster frequencies from FlowSOM assignments
#'
#' @description
#' Aggregates \code{facs_cluster_flowsom()}'s per-event cluster/metacluster
#' assignments to per-sample counts and fractions, one row per
#' \code{file_name} x \code{cluster_col} value. Every \code{file_name} x
#' \code{cluster_col} combination seen anywhere in \code{data} is
#' represented (zero-filled where a sample had no events in that cluster),
#' so the result is ready for count-based differential abundance testing
#' without further completion.
#'
#' @param data tibble shaped like \code{facs_cluster_flowsom()}'s output:
#'   must contain \code{file_name} and \code{cluster_col}.
#' @param cluster_col character; column in \code{data} to aggregate on.
#'   Default \code{"metacluster"}; pass \code{"cluster"} for
#'   raw-SOM-node frequencies instead.
#'
#' @returns Long tibble: \code{file_name}, \code{cluster_col} (name
#'   reused), \code{n} (event count), \code{fraction} (\code{n} divided by
#'   that sample's total event count), and any column from \code{data}
#'   that is constant within \code{file_name} (e.g. keyword columns),
#'   carried through automatically. Errors if \code{cluster_col} is not a
#'   column in \code{data}.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_cluster_flowsom(facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   ), seed = 1)
#'   facs_calc_cluster_freq(dat)
#' }
facs_calc_cluster_freq <- function(data, cluster_col = "metacluster") {
  if (!cluster_col %in% names(data)) {
    stop(glue::glue("`cluster_col` ('{cluster_col}') not found in `data`."))
  }

  totals <- data |>
    dplyr::count(file_name, name = "total")

  passthrough_candidates <- setdiff(names(data), c("file_name", cluster_col))
  is_constant_ <- function(col) {
    n_distinct_per_file <- data |>
      dplyr::group_by(file_name) |>
      dplyr::summarise(n_distinct = dplyr::n_distinct(.data[[col]]), .groups = "drop") |>
      dplyr::pull(n_distinct)
    all(n_distinct_per_file == 1L)
  }
  passthrough_cols <- passthrough_candidates[purrr::map_lgl(passthrough_candidates, is_constant_)]

  passthrough <- data |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c("file_name", passthrough_cols))))

  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("file_name", cluster_col)))) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    tidyr::complete(file_name, .data[[cluster_col]], fill = list(n = 0L)) |>
    dplyr::left_join(totals, by = "file_name") |>
    dplyr::mutate(fraction = n / total) |>
    dplyr::select(!total) |>
    dplyr::left_join(passthrough, by = "file_name") |>
    dplyr::arrange(file_name, .data[[cluster_col]])
}

#' Compute per-cluster median marker expression
#'
#' @description
#' Summarizes \code{facs_cluster_flowsom()}'s per-event marker columns to a
#' per-cluster median, one row per \code{cluster_col} x marker. Intended as
#' a cluster-annotation step -- pairs with \code{facs_plot_cluster_heatmap()}
#' to let a metacluster be named by its marker expression profile before
#' interpreting a \code{facs_test_cluster_abundance()} hit.
#'
#' @param data tibble shaped like \code{facs_cluster_flowsom()}'s output:
#'   one row per event, \code{dbl} marker columns, and a \code{cluster_col}
#'   column.
#' @param markers character vector of column names in \code{data} to
#'   summarize. \code{NULL} (default) uses every \code{dbl}-typed column in
#'   \code{data}, identical convention to \code{facs_cluster_flowsom()}'s
#'   own \code{markers} argument. When explicitly supplied, every named
#'   column must be \code{double}-typed in \code{data}.
#' @param cluster_col character; column in \code{data} to group by. Default
#'   \code{"metacluster"}, matching \code{facs_cluster_flowsom()}'s and
#'   \code{facs_calc_cluster_freq()}'s own default.
#'
#' @returns Long tibble: \code{{cluster_col name}} (type unchanged from
#'   \code{data}), \code{marker} (chr), \code{median} (dbl). One row per
#'   cluster x marker. Errors if \code{cluster_col} is not a column in
#'   \code{data}, if an explicitly supplied \code{markers} name is absent
#'   from \code{data}, if an explicitly supplied \code{markers} column is
#'   not \code{double}-typed in \code{data} (listing the offending
#'   column(s)), or if the resolved \code{markers} is empty.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_cluster_flowsom(facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   ), seed = 1)
#'   facs_calc_cluster_marker_medians(dat)
#' }
facs_calc_cluster_marker_medians <- function(data, markers = NULL, cluster_col = "metacluster") {
  if (!cluster_col %in% names(data)) {
    stop(glue::glue("`cluster_col` ('{cluster_col}') not found in `data`."))
  }

  if (is.null(markers)) {
    is_dbl  <- purrr::map_lgl(data, is.double)
    markers <- names(data)[is_dbl]
  } else {
    missing_markers <- setdiff(markers, names(data))
    if (length(missing_markers) > 0L) {
      stop(glue::glue(
        "The following `markers` were not found in `data`: ",
        "{paste(missing_markers, collapse = ', ')}"
      ))
    }

    non_double_markers <- markers[!purrr::map_lgl(markers, function(m) is.double(data[[m]]))]
    if (length(non_double_markers) > 0L) {
      stop(glue::glue(
        "The following `markers` are not double-typed columns in `data`: ",
        "{paste(non_double_markers, collapse = ', ')}"
      ))
    }
  }

  if (length(markers) == 0L) {
    stop("Resolved `markers` is empty -- no double-typed columns found in `data`.")
  }

  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(cluster_col))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(markers), stats::median), .groups = "drop") |>
    tidyr::pivot_longer(cols = dplyr::all_of(markers), names_to = "marker", values_to = "median")
}
