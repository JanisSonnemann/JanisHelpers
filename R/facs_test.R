#' Test differential cluster abundance between groups
#'
#' @description
#' Tests whether \code{facs_calc_cluster_freq()}'s per-sample cluster
#' frequencies differ between levels of a grouping column, using
#' Bioconductor's \code{diffcyt} package. One \code{diffcyt} test is run per
#' non-reference level of \code{fixed} (relative to \code{ref_level}),
#' producing one contrast per comparison. Builds the
#' \code{SummarizedExperiment} \code{diffcyt} expects directly from
#' \code{freq_data} (bypassing \code{diffcyt::prepareData()}/
#' \code{calcCounts()}, since the per-sample-per-cluster counts are already
#' computed).
#'
#' @param freq_data tibble shaped like \code{facs_calc_cluster_freq()}'s
#'   output: must contain \code{file_name}, \code{cluster_col}, \code{n},
#'   and the \code{fixed}/\code{random} columns (either already present, or
#'   supplied via \code{meta}).
#' @param meta optional tibble to left-join onto \code{freq_data} via
#'   \code{meta_annotate()} before testing. \code{NULL} (default) assumes
#'   \code{fixed}/\code{random} are already columns in \code{freq_data}.
#' @param fixed character; name of the fixed-effect column to test (e.g.
#'   \code{"group"}). Must resolve to a column with at least 2 levels.
#' @param random character vector or \code{NULL} (default); random-effect
#'   column(s) (e.g. \code{"mouse_ID"}) for pairing/blocking. Only
#'   supported with \code{method = "glmm"}.
#' @param by character vector; join key(s) forwarded to \code{meta_annotate()}
#'   when \code{meta} is supplied. Default \code{"mouse_ID"}. Ignored if
#'   \code{meta} is \code{NULL}.
#' @param cluster_col character; column in \code{freq_data} identifying the
#'   cluster. Default \code{"metacluster"}, matching
#'   \code{facs_calc_cluster_freq()}'s own default.
#' @param ref_level character or \code{NULL} (default); reference level of
#'   \code{fixed}. \code{NULL} uses \code{fixed}'s first factor level.
#' @param method character; one of \code{"glmm"} (default), \code{"edgeR"},
#'   \code{"voom"} -- which \code{diffcyt} test function to dispatch to.
#'
#' @returns Tibble, one row per \code{cluster_col} value x contrast:
#'   \code{{cluster_col name}}, \code{contrast} (chr,
#'   \code{"{level}_vs_{ref_level}"}), \code{p_val}, \code{p_adj}, plus
#'   method-specific columns (\code{method = "edgeR"} adds \code{logFC},
#'   \code{logCPM}, \code{LR}; \code{method = "voom"} adds \code{logFC},
#'   \code{AveExpr}, \code{t}, \code{B}; \code{method = "glmm"} adds no
#'   further columns). Errors if \code{fixed}, \code{random}, or
#'   \code{cluster_col} are not columns in \code{freq_data} (after any
#'   \code{meta} join), if \code{fixed} has fewer than 2 levels, if
#'   \code{ref_level} is not among \code{fixed}'s levels, if \code{random}
#'   is supplied with \code{method != "glmm"}, or if \code{fixed}/\code{random}
#'   contain \code{NA} after an unmatched \code{meta} join.
#' @export
#'
#' @examples
#' \dontrun{
#'   freq <- facs_calc_cluster_freq(facs_cluster_flowsom(facs_read_fcs_gated(
#'     wsp_path  = "experiment.wsp",
#'     gate_path = "Singlets/Lymphocytes/live/CD45+",
#'     markers   = c("CD4", "CD45")
#'   )))
#'   facs_test_cluster_abundance(freq, fixed = "group", method = "glmm")
#' }
facs_test_cluster_abundance <- function(freq_data,
                                         meta = NULL,
                                         fixed,
                                         random = NULL,
                                         by = "mouse_ID",
                                         cluster_col = "metacluster",
                                         ref_level = NULL,
                                         method = c("glmm", "edgeR", "voom")) {
  method <- match.arg(method)

  if (!is.null(meta)) {
    freq_data <- meta_annotate(freq_data, meta, by = by)
  }

  required_cols <- c("file_name", cluster_col, "n", fixed, random)
  missing_cols <- setdiff(required_cols, names(freq_data))
  if (length(missing_cols) > 0L) {
    stop(glue::glue(
      "The following column(s) were not found in `freq_data`: ",
      "{paste(missing_cols, collapse = ', ')}."
    ))
  }

  if (method != "glmm" && !is.null(random)) {
    stop(glue::glue(
      "random effects are only supported with method = 'glmm'; pass ",
      "random = NULL, or fold this column into `fixed`, for edgeR/voom."
    ))
  }

  fixed_vec <- freq_data[[fixed]]
  if (!is.factor(fixed_vec)) fixed_vec <- factor(fixed_vec)
  if (nlevels(fixed_vec) < 2L) {
    stop(glue::glue(
      "`fixed` column '{fixed}' must have at least 2 levels ",
      "(found {nlevels(fixed_vec)})."
    ))
  }
  if (is.null(ref_level)) ref_level <- levels(fixed_vec)[1]
  if (!ref_level %in% levels(fixed_vec)) {
    stop(glue::glue(
      "`ref_level` ('{ref_level}') is not among `fixed`'s levels: ",
      "{paste(levels(fixed_vec), collapse = ', ')}."
    ))
  }

  na_check_cols <- c(fixed, random)
  na_present <- na_check_cols[purrr::map_lgl(na_check_cols, function(col) anyNA(freq_data[[col]]))]
  if (length(na_present) > 0L) {
    stop(glue::glue(
      "The following column(s) contain NA (e.g. from an unmatched `meta` ",
      "join key): {paste(na_present, collapse = ', ')}."
    ))
  }

  experiment_info <- freq_data |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c("file_name", fixed, random)))) |>
    dplyr::rename(sample_id = file_name) |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(fixed),
      ~ stats::relevel(factor(.x), ref = ref_level)
    ))

  counts_wide <- freq_data |>
    dplyr::select(dplyr::all_of(c("file_name", cluster_col, "n"))) |>
    tidyr::pivot_wider(names_from = "file_name", values_from = "n")

  cluster_ids <- counts_wide[[cluster_col]]
  count_matrix <- as.matrix(counts_wide[, setdiff(names(counts_wide), cluster_col)])
  count_matrix <- count_matrix[, as.character(experiment_info$sample_id), drop = FALSE]

  d_counts <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = count_matrix),
    rowData = data.frame(cluster_id = cluster_ids, n_cells = rowSums(count_matrix)),
    colData = experiment_info
  )

  design_mm <- stats::model.matrix(stats::as.formula(paste0("~", fixed)), data = experiment_info)
  non_ref_cols <- setdiff(colnames(design_mm), "(Intercept)")

  results <- purrr::map_dfr(seq_along(non_ref_cols), function(i) {
    contrast_vec <- rep(0, ncol(design_mm))
    contrast_vec[i + 1] <- 1
    contrast <- diffcyt::createContrast(contrast_vec)
    level_name <- stringr::str_remove(non_ref_cols[i], paste0("^", fixed))

    res <- if (method == "glmm") {
      formula_obj <- diffcyt::createFormula(
        experiment_info,
        cols_fixed  = fixed,
        cols_random = c("sample_id", random)
      )
      diffcyt::testDA_GLMM(d_counts, formula_obj, contrast)
    } else {
      design <- diffcyt::createDesignMatrix(experiment_info, cols_design = fixed)
      if (method == "edgeR") {
        diffcyt::testDA_edgeR(d_counts, design, contrast)
      } else {
        diffcyt::testDA_voom(d_counts, design, contrast)
      }
    }

    tibble::as_tibble(as.data.frame(SummarizedExperiment::rowData(res))) |>
      dplyr::mutate(
        contrast = paste0(level_name, "_vs_", ref_level),
        .after = "cluster_id"
      )
  })

  dplyr::rename(results, !!cluster_col := cluster_id)
}
