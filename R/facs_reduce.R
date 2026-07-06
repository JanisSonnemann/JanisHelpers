# uwot behavior -- verified against tests/fixtures/Treg.wsp
#
# uwot::umap()'s neighbor search and SGD optimization are parallelized by
# default (n_threads). Per uwot's own documentation, exact run-to-run
# reproducibility under a fixed set.seed() is only guaranteed with
# single-threaded execution, so facs_reduce_umap() forces n_threads = 1
# whenever `seed` is supplied (trading speed for honoring the `seed`
# contract), leaving uwot's own default threading in place otherwise.
# uwot::umap() also requires enough input rows relative to n_neighbors
# (roughly n_neighbors < nrow(data)) -- callers working with small or
# heavily downsampled inputs must pass a correspondingly small
# n_neighbors.

# Stratified-per-file_name downsampling for facs_reduce_umap(). Keeps
# every sample visible in the UMAP embedding regardless of how many
# events it happened to contribute -- a pure random pool would let a
# high-yield sample dominate the embedding independent of biology. A
# sample with fewer rows than its computed share contributes all of its
# rows; the shortfall is not redistributed to other samples, so the
# returned row count can come in slightly under `max_events`.
downsample_stratified_ <- function(data, max_events, seed) {
  if (!is.null(seed)) set.seed(seed)

  file_names <- sort(unique(data$file_name))
  n_samples  <- length(file_names)
  base_share <- max_events %/% n_samples
  remainder  <- max_events %% n_samples

  shares <- rep(base_share, n_samples)
  if (remainder > 0L) {
    shares[seq_len(remainder)] <- shares[seq_len(remainder)] + 1L
  }
  names(shares) <- file_names

  purrr::map(file_names, function(fn) {
    rows  <- data[data$file_name == fn, , drop = FALSE]
    share <- shares[[fn]]
    if (share <= 0L || nrow(rows) == 0L) {
      return(rows[0L, , drop = FALSE])
    }
    if (nrow(rows) <= share) {
      return(rows)
    }
    rows[sample.int(nrow(rows), share), , drop = FALSE]
  }) |>
    dplyr::bind_rows()
}

#' Compute a UMAP embedding of gated single-cell events
#'
#' @description
#' Runs UMAP (via \code{uwot::umap()}) on a per-event tibble (e.g.
#' \code{facs_read_fcs_gated()}'s output), appending a 2D embedding for
#' visualization. Embedding is performed directly on the input's numeric
#' scale (no additional z-score normalization), matching the transformed
#' (logicle/biexponential) scale \code{facs_read_fcs_gated()} already
#' returns.
#'
#' @param data tibble shaped like \code{facs_read_fcs_gated()}'s output:
#'   one row per event, \code{dbl} marker columns, \code{file_name}, and
#'   optionally \code{chr} keyword columns.
#' @param markers character vector of column names in \code{data} to embed
#'   on. \code{NULL} (default) uses every \code{dbl}-typed column in
#'   \code{data}. Must resolve to at least 2 columns, and, when explicitly
#'   supplied, every named column must be \code{double}-typed in
#'   \code{data}.
#' @param max_events integer or \code{NULL} (default). If set, downsamples
#'   the combined tibble to (approximately) this many total rows before
#'   running UMAP, stratified per \code{file_name}: each sample
#'   contributes an equal share of \code{max_events} (remainder
#'   distributed one extra row each to the first samples by
#'   \code{file_name} sort order), and a sample with fewer rows than its
#'   share contributes all of its rows without redistributing the
#'   shortfall elsewhere. No effect if \code{max_events} is at least
#'   \code{nrow(data)}.
#' @param n_neighbors,min_dist numeric; passed through to
#'   \code{uwot::umap()} unchanged. Defaults \code{15}/\code{0.1}
#'   (\code{uwot}'s own defaults).
#' @param seed integer; if set, seeds both the stratified downsampling
#'   draw and \code{uwot::umap()} (forcing single-threaded execution so
#'   the embedding itself is reproducible).
#'
#' @returns \code{data} (or its downsampled subset, if \code{max_events}
#'   triggered downsampling) with two columns appended: \code{UMAP1} and
#'   \code{UMAP2} (both \code{dbl}). Every returned row has a real
#'   embedding -- rows excluded by downsampling are dropped, not kept with
#'   \code{NA}. Errors if any \code{markers} name is absent from
#'   \code{data}, if an explicitly supplied \code{markers} column is not
#'   \code{double}-typed in \code{data} (listing the offending column(s)),
#'   if the resolved \code{markers} has fewer than 2 columns, if a
#'   selected marker column contains \code{NA}, or if \code{max_events} is
#'   supplied and is not a single positive integer.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   )
#'   facs_reduce_umap(dat, seed = 1)
#' }
facs_reduce_umap <- function(data,
                              markers = NULL,
                              max_events = NULL,
                              n_neighbors = 15,
                              min_dist = 0.1,
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
      "`markers` must resolve to at least 2 columns (found {length(markers)})."
    ))
  }

  na_cols <- markers[purrr::map_lgl(markers, function(m) anyNA(data[[m]]))]
  if (length(na_cols) > 0L) {
    stop(glue::glue(
      "The following marker column(s) contain NA and cannot be embedded: ",
      "{paste(na_cols, collapse = ', ')}"
    ))
  }

  if (!is.null(max_events)) {
    if (!is.numeric(max_events) || length(max_events) != 1L ||
        max_events != as.integer(max_events) || max_events <= 0L) {
      stop("`max_events` must be a single positive integer.")
    }
    max_events <- as.integer(max_events)
    if (max_events < nrow(data)) {
      data <- downsample_stratified_(data, max_events, seed)
    }
  }

  if (!is.null(seed)) set.seed(seed)

  umap_args <- list(
    X           = as.matrix(data[markers]),
    n_neighbors = n_neighbors,
    min_dist    = min_dist
  )
  if (!is.null(seed)) {
    umap_args$n_threads <- 1L
  }
  embedding <- do.call(uwot::umap, umap_args)

  data$UMAP1 <- embedding[, 1]
  data$UMAP2 <- embedding[, 2]

  data
}
