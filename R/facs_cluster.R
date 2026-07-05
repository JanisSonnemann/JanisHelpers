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
# Two reproducible upstream bugs constrain which grid/marker/metacluster
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
#    <= grid_xdim * grid_ydim, per the check below) until upstream fixes
#    this.
# 2. Selecting exactly one marker column always throws ("argument must be
#    coercible to non-negative integer") from FlowSOM::BuildSOM()'s
#    `fsom$data[, colsToUse]`, which has no `drop = FALSE` -- a
#    single-column `colsToUse` collapses to a bare vector, and the
#    subsequent seq_len(ncol(data)) inside SOM() fails on its NULL
#    ncol(). `markers` (explicit or defaulted) must resolve to >= 2
#    columns until upstream fixes this.
#
# Neither is guarded against here (matching the "don't silently paper
# over a mismatch" instruction) -- both surface as FlowSOM/
# ConsensusClusterPlus errors, not JanisHelpers ones, if triggered.

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
#'   in \code{data}.
#' @param grid_xdim,grid_ydim integer; SOM grid dimensions. Default
#'   \code{10}/\code{10} (100 nodes).
#' @param n_metaclusters integer; target consensus metacluster count.
#'   Default \code{10}. Must not exceed \code{grid_xdim * grid_ydim}.
#' @param seed integer; if set, seeds SOM training and consensus
#'   metaclustering for reproducible assignments.
#'
#' @returns \code{data} with two columns appended: \code{cluster} (integer,
#'   raw SOM node, \code{1:(grid_xdim * grid_ydim)}) and \code{metacluster}
#'   (factor, consensus grouping, \code{1:n_metaclusters} levels). Errors if
#'   any \code{markers} name is absent from \code{data}, if a selected
#'   marker column contains \code{NA}, or if \code{n_metaclusters} exceeds
#'   the grid's node count.
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
  }

  na_cols <- markers[purrr::map_lgl(markers, function(m) anyNA(data[[m]]))]
  if (length(na_cols) > 0L) {
    stop(glue::glue(
      "The following marker column(s) contain NA and cannot be clustered: ",
      "{paste(na_cols, collapse = ', ')}"
    ))
  }

  if (n_metaclusters > grid_xdim * grid_ydim) {
    stop(glue::glue(
      "`n_metaclusters` ({n_metaclusters}) cannot exceed the SOM grid's ",
      "node count (grid_xdim * grid_ydim = {grid_xdim * grid_ydim})."
    ))
  }

  input <- flowCore::flowFrame(as.matrix(data[markers]))

  # See the "FlowSOM behavior" header comment at the top of this file for
  # the `xdim`/`ydim`/return-shape notes and the two upstream bugs this
  # call is subject to.
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
