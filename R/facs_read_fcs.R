# FlowWorkspace/CytoML behavior -- verified against tests/fixtures/Treg.wsp
#
# flowjo_to_gatingset(ws, name = <group>, path = <fcs dir>) builds a
# GatingSet, recursively searching <fcs dir> for the .fcs files the
# workspace references, replaying compensation/transformation/gating.
#
# sampleNames(gs) returns internal names with a numeric suffix
# (e.g. "sample.fcs_10000") -- the clean original filename is
# pData(gh)$name. Always use pData(gh)$name for file_name, never
# sampleNames().
#
# gh_pop_get_data(gh, path) returns a lazy view; realize_view() must be
# called before exprs() returns real data. A non-existent path throws a
# catchable simpleError.
#
# markernames(gh) returns a named vector (channel = stain label) but only
# for channels that have a stain label. Channels without one (FSC-A,
# SSC-A, Time, ...) must be matched via parameters(gated)$name instead.
#
# keyword(gh, name) returns NULL (not an error) when the keyword is
# absent for that sample. It only sees keywords physically embedded in
# the raw .fcs file's own TEXT segment -- a keyword typed into FlowJo's
# UI but never written back into the .fcs file (a "workspace-only"
# keyword) is invisible to it, even though FlowJo's own "Inspect" view
# shows it. Keywords are therefore read from the .wsp XML directly (via
# parse_keywords_(), the same helper facs_read_wsp() uses), which is the
# authoritative superset of both raw-.fcs-embedded and workspace-only
# keywords.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

resolve_markers_ <- function(markers, gh, gated, file_name) {
  stain_lookup  <- flowWorkspace::markernames(gh)
  channel_names <- flowCore::parameters(gated)$name

  resolved <- purrr::map_chr(markers, function(m) {
    if (m %in% stain_lookup) {
      return(names(stain_lookup)[stain_lookup == m][[1]])
    }
    if (m %in% channel_names) {
      return(m)
    }
    NA_character_
  })

  unmatched <- markers[is.na(resolved)]
  if (length(unmatched) > 0L) {
    stop(glue::glue(
      "The following markers could not be matched (by stain label or ",
      "channel name) in '{file_name}': {paste(unmatched, collapse = ', ')}"
    ))
  }

  resolved
}

# Resolves the sample names belonging to a FlowJo group directly from the
# workspace XML metadata -- no GatingSet build required, so this can run
# once up front before any per-sample work is dispatched.
resolve_group_samples_ <- function(ws, group) {
  groups  <- CytoML::fj_ws_get_sample_groups(ws)
  samples <- CytoML::fj_ws_get_samples(ws)
  ids <- groups$sampleID[groups$groupName == group]
  names <- samples$name[samples$sampleID %in% ids]
  if (length(names) == 0L) {
    stop(glue::glue("No samples found for FlowJo group '{group}'."))
  }
  names
}

read_one_sample_solo_ <- function(wsp_path, fcs_dir, group, sample_name, sample_index,
                                   gate_path, gate_path_norm, markers,
                                   max_events, seed) {
  warnings_ <- character()
  add_warning_ <- function(w) {
    warnings_ <<- c(warnings_, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

  data <- withCallingHandlers(
    {
      # CytoML's workspace/GatingSet objects are C++ external pointers and
      # cannot be passed across a process boundary (forked or socket-based)
      # -- each call must open its own copy from the plain wsp_path string.
      ws <- CytoML::open_flowjo_xml(wsp_path)

      # CytoML::flowjo_to_gatingset()'s `subset` argument is resolved via
      # eval(substitute(subset)) internally -- this only sees a bare variable
      # name (not its value) when the argument is forwarded through a wrapper
      # function like this one, so it must be called via do.call() to pass
      # sample_name as an already-evaluated value rather than a promise.
      gs <- do.call(
        CytoML::flowjo_to_gatingset,
        list(ws = ws, name = group, path = fcs_dir, subset = sample_name)
      )
      gh <- gs[[1]]
      file_name <- flowWorkspace::pData(gh)$name

      gated <- tryCatch(
        flowWorkspace::gh_pop_get_data(gh, gate_path_norm),
        error = function(e) NULL
      )

      if (is.null(gated)) {
        warning(glue::glue(
          "gate_path '{gate_path}' not found for '{file_name}'; sample skipped."
        ))
        tibble::tibble()
      } else {
        gated <- flowWorkspace::realize_view(gated)

        resolved_cols <- resolve_markers_(markers, gh, gated, file_name)

        mat <- flowCore::exprs(gated)[, resolved_cols, drop = FALSE]
        colnames(mat) <- markers

        if (!is.null(max_events) && nrow(mat) > max_events) {
          if (!is.null(seed)) set.seed(seed + sample_index)
          mat <- mat[sample.int(nrow(mat), max_events), , drop = FALSE]
        }

        tibble::as_tibble(mat) |>
          tibble::add_column(file_name = file_name, .before = 1L)
      }
    },
    warning = add_warning_
  )

  list(data = data, warnings = warnings_)
}

# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' Read raw single-cell events from .fcs files, filtered to a FlowJo gate
#'
#' @description
#' Reads raw single-cell events directly from the .fcs files referenced by
#' a FlowJo \code{.wsp} workspace, replaying the workspace's compensation,
#' transformation, and gating tree (via \code{CytoML}/\code{flowWorkspace})
#' so events are filtered down to an arbitrary already-drawn gate. Feeds
#' unsupervised, cell-level analysis (e.g. FlowSOM, UMAP) that
#' \code{facs_read_wsp()}'s gated summary statistics cannot support.
#'
#' @param wsp_path path to the \code{.wsp} file.
#' @param gate_path character; full gating path, e.g.
#'   \code{"Singlets/Lymphocytes/live/CD45+"} -- same format as
#'   \code{PopulationFullPath} from \code{facs_read_wsp()}. Applied to
#'   every sample in \code{group}.
#' @param markers character vector; matched per sample against stain
#'   label or channel name (stain label preferred).
#' @param keywords character vector of FlowJo keyword names to append as
#'   columns. Read from the \code{.wsp} workspace XML (same source as
#'   \code{facs_read_wsp()}), not from the raw \code{.fcs} files -- this
#'   also picks up keywords typed into FlowJo's UI that were never written
#'   back into the \code{.fcs} file itself. A keyword missing for every
#'   sample is filled \code{NA_character_} with a warning.
#' @param fcs_dir folder to search for this workspace's \code{.fcs}
#'   files. \code{NULL} (default) auto-derives it as the subfolder named
#'   after the \code{.wsp} file (sans extension), sitting next to it.
#' @param group character; FlowJo sample group to load. Default
#'   \code{"All Samples"}.
#' @param max_events integer; if set, randomly downsample each sample to
#'   at most this many events.
#' @param seed integer; if set, seeds the random draw used by
#'   \code{max_events} for reproducibility.
#'
#' @returns Wide tibble, one row per event: \code{file_name}, one column
#'   per requested marker (on FlowJo's transformed scale), and one column
#'   per requested keyword. Errors if a requested marker cannot be
#'   matched in a sample's panel. Warns and skips a sample if
#'   \code{gate_path} does not exist for it. Warns and fills \code{NA} if
#'   a requested keyword is missing for every sample.
#' @export
#'
#' @examples
#' \dontrun{
#'   facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45"),
#'     keywords  = c("mouse_ID", "tissue")
#'   )
#' }
facs_read_fcs_gated <- function(wsp_path,
                                 gate_path,
                                 markers,
                                 keywords = NULL,
                                 fcs_dir = NULL,
                                 group = "All Samples",
                                 max_events = NULL,
                                 seed = NULL,
                                 workers = 1L) {
  if (is.null(fcs_dir)) {
    fcs_dir <- file.path(
      dirname(wsp_path),
      tools::file_path_sans_ext(basename(wsp_path))
    )
  }
  gate_path_norm <- if (startsWith(gate_path, "/")) gate_path else paste0("/", gate_path)

  ws <- CytoML::open_flowjo_xml(wsp_path)
  sample_names <- resolve_group_samples_(ws, group)
  rm(ws)

  # read_one_sample_solo_()'s own environment() attribute is the package
  # namespace. Its internal (non-::-qualified) call to resolve_markers_()
  # would otherwise resolve via that namespace when reconstructed on a PSOCK
  # worker -- which loads whatever's actually *installed* there, not the
  # code currently in memory here (e.g. under devtools::load_all()).
  # Re-homing it onto a plain environment before use fixes this; verified
  # via canary-function testing against the real fixture.
  export_env <- new.env()
  export_env$resolve_markers_ <- resolve_markers_
  read_one_sample_solo_export_ <- read_one_sample_solo_
  environment(read_one_sample_solo_export_) <- export_env

  worker_fn <- function(i) {
    read_one_sample_solo_export_(
      wsp_path       = wsp_path,
      fcs_dir        = fcs_dir,
      group          = group,
      sample_name    = sample_names[i],
      sample_index   = i,
      gate_path      = gate_path,
      gate_path_norm = gate_path_norm,
      markers        = markers,
      max_events     = max_events,
      seed           = seed
    )
  }

  results <- if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, seq_along(sample_names), worker_fn)
  } else {
    purrr::map(seq_along(sample_names), worker_fn)
  }

  # parallel::parLapply() already raises a per-worker error immediately as
  # a real error in this process (verified: the propagated message is
  # prefixed "one node produced an error: <original message>", so existing
  # regexp-based expect_error() checks still match on the original text) --
  # no try-error detection/re-throw needed here, unlike the fork-based
  # mclapply() approach this replaced.

  for (res in results) {
    for (w in res$warnings) warning(w, call. = FALSE)
  }

  data <- purrr::map(results, "data") |>
    dplyr::bind_rows()

  if (!is.null(keywords) && length(keywords) > 0L) {
    doc <- xml2::read_xml(wsp_path)
    user_kws <- parse_keywords_(doc) |>
      dplyr::filter(key %in% keywords) |>
      tidyr::pivot_wider(names_from = key, values_from = value)

    missing_kws <- setdiff(keywords, names(user_kws))
    if (length(missing_kws) > 0L) {
      warning(
        "The following requested keywords were not found in the workspace ",
        "and were filled with NA: ",
        paste(missing_kws, collapse = ", ")
      )
      user_kws[missing_kws] <- NA_character_
    }

    data <- dplyr::left_join(data, user_kws, by = "file_name")
  }

  data
}
