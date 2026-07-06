# Stage 4: facs_test_cluster_abundance() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add differential cluster abundance testing (Stage 4 of the unsupervised FACS pipeline) via `diffcyt`: a new `facs_test_cluster_abundance()` in `R/facs_test.R`, plus a companion `facs_plot_cluster_abundance()` added to the existing `R/facs_plot.R`.

**Architecture:** `facs_test_cluster_abundance()` takes `facs_calc_cluster_freq()`'s output (a `file_name` x `cluster_col` tibble of `n`/`fraction`), builds a `SummarizedExperiment` by hand (bypassing `diffcyt::prepareData()`/`calcCounts()`, since the counts are already computed), and dispatches to `diffcyt::testDA_GLMM()`/`testDA_edgeR()`/`testDA_voom()` -- once per non-reference level of the `fixed` column relative to `ref_level` -- row-binding results into one tibble. `facs_plot_cluster_abundance()` is a `ggplot2` boxplot/jitter of `fraction` by group, faceted by cluster, in the style of the existing `facs_plot_umap()`.

**Tech Stack:** R, `diffcyt`/`SummarizedExperiment` (Bioconductor, both new), `dplyr`/`tidyr`/`purrr`/`glue`/`tibble`/`stringr`/`stats` (existing), `ggplot2` (existing since Stage 3), `testthat` edition 3.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-07-06-facs-test-cluster-abundance-design.md` -- read it before starting; this plan implements it exactly except for the Step-2 "recompute total" simplification noted below.
- **Correction to the approved spec:** the design's data-flow Step 2 ("recompute each sample's total event count as `sum(n)`") is unneeded. Verified empirically (see below) -- `testDA_GLMM()`/`testDA_edgeR()`/`testDA_voom()` derive each sample's total directly from the count matrix's own column sums, since `facs_calc_cluster_freq()`'s zero-filled complete grid guarantees every cluster is present per sample. Do not implement that step.
- Namespacing: every external call is `pkg::fn()` (no bare `map()`/`filter()`/etc.) per `CLAUDE.md`.
- Pipe: `|>` (base pipe) everywhere.
- `@param`/`@returns`/`@export` required on every exported function; `@returns` must be non-empty; `\dontrun{}` for examples requiring external files (per `CLAUDE.md`'s Roxygen conventions).
- No `invisible()` on `facs_test_cluster_abundance()`'s return -- the tibble is the primary output, not a side effect.
- Deselection uses `dplyr::select(!col)`, not `-col`.
- `devtools::check()` target: 0 errors, 0 warnings (some pre-existing NOTEs are expected and already documented in `CLAUDE.md`; add new NOTEs there rather than trying to eliminate tidy-eval false positives).
- Tests use `testthat` edition 3, synthetic data only (no `.wsp`/`.fcs` fixture dependency) -- mirrors Stage 3's precedent, since both new functions operate on plain tibbles.
- `diffcyt` (Bioconductor 1.28.0, confirmed installed and working against Bioconductor 3.21 during this plan's own verification) and `SummarizedExperiment` are new `Imports`.

---

## Verified diffcyt mechanics (read before Task 1)

This plan's code was written against a real, working prototype run against installed `diffcyt` 1.28.0 -- not assumed from documentation. Key facts, confirmed by direct experimentation:

- `calcCounts()`'s output shape (which `testDA_*()` expects as `d_counts`) is: `assay` = a plain numeric matrix, rows = clusters, columns = samples; `rowData` = a data frame with `cluster_id` (factor) and `n_cells` (numeric, `rowSums` of the assay); `colData` = exactly the `experiment_info` data frame (one row per sample, in the same order as the assay's columns). Building this by hand with `SummarizedExperiment::SummarizedExperiment(assays = list(counts = ...), rowData = ..., colData = ...)` works identically to going through `diffcyt::prepareData()`/`calcCounts()` -- confirmed by running `testDA_GLMM()`/`testDA_edgeR()`/`testDA_voom()` on a hand-built object and getting sensible results.
- `diffcyt::createFormula(experiment_info, cols_fixed = fixed, cols_random = c("sample_id", random))` -- `"sample_id"` is *always* included in `cols_random`, even when the caller's `random` is `NULL`, because it supplies the per-sample "observation-level random effect" (OLRE) that `diffcyt`'s GLMM method relies on for overdispersion control (confirmed via `createFormula.Rd`'s own worked example, which always includes `sample_id`). This produces a formula like `y ~ group_id + (1 | sample_id) + (1 | mouse_id)`.
- Contrast vectors correspond to `stats::model.matrix(~ fixed, data = experiment_info)`'s columns: `"(Intercept)"` followed by one column per non-reference level, in factor-level order once `fixed` is releveled via `stats::relevel(factor(fixed_col), ref = ref_level)`. For the i-th non-reference level, the contrast vector is all zeros except a `1` at position `i + 1`. Confirmed for both a 2-level and a 3-level factor.
- Output columns via `SummarizedExperiment::rowData()` on the result, confirmed by direct inspection:
  - `method = "glmm"`: `cluster_id`, `p_val`, `p_adj` only -- **no effect-size/direction column**.
  - `method = "edgeR"`: `cluster_id`, `logFC`, `logCPM`, `LR`, `p_val`, `p_adj`.
  - `method = "voom"`: `cluster_id`, `logFC`, `AveExpr`, `t`, `p_val`, `p_adj`, `B`.
- `"boundary (singular) fit"` / convergence warnings from `lme4` are expected with small synthetic sample sizes and are not suppressed -- they surface naturally to the caller, consistent with not silencing meaningful warnings elsewhere in the package.

---

### Task 1: Core statistical function -- happy path across all three methods

**Files:**
- Create: `R/facs_test.R`
- Modify: `DESCRIPTION` (add `diffcyt`, `SummarizedExperiment` to `Imports`, alphabetically)
- Test: `tests/testthat/test-facs_test.R`

**Interfaces:**
- Produces: `facs_test_cluster_abundance(freq_data, meta = NULL, fixed, random = NULL, by = "mouse_ID", cluster_col = "metacluster", ref_level = NULL, method = c("glmm", "edgeR", "voom"))` -- returns a tibble with columns `{cluster_col name}`, `contrast` (chr, `"{level}_vs_{ref_level}"`), `p_val`, `p_adj`, plus method-specific extra columns (see above). This task implements the full computational core assuming well-formed input (`fixed`/`random`/`cluster_col` present, `fixed` has >= 2 levels, no `NA`s) and `method %in% c("glmm","edgeR","voom")` used correctly (no `random` + non-GLMM misuse yet -- that guard is Task 2). `meta`/`by`-based joining is also Task 2; this task assumes `fixed`/`random` columns are already present in `freq_data`.

- [ ] **Step 1: Add new dependencies to `DESCRIPTION`**

Edit the `Imports:` block in `DESCRIPTION` to read (alphabetical order preserved):

```
Imports:
    CytoML,
    diffcyt,
    dplyr,
    fcexpr,
    flowCore,
    FlowSOM,
    flowWorkspace,
    ggplot2,
    glue,
    gt,
    gtsummary,
    janitor,
    parallel,
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stats,
    stringr,
    SummarizedExperiment,
    tibble,
    tidyr,
    tools,
    uwot,
    xml2,
    xfun
```

- [ ] **Step 2: Confirm diffcyt/SummarizedExperiment are installed**

Run: `Rscript -e 'requireNamespace("diffcyt", quietly = TRUE) && requireNamespace("SummarizedExperiment", quietly = TRUE)'`
Expected: `[1] TRUE`

If `FALSE`, install with:
```r
options(repos = BiocManager::repositories())
BiocManager::install(c("diffcyt", "SummarizedExperiment"), update = FALSE, ask = FALSE)
```

- [ ] **Step 3: Write the failing test file**

Create `tests/testthat/test-facs_test.R`:

```r
library(testthat)
library(JanisHelpers)

# Synthetic tibble shaped like facs_calc_cluster_freq()'s output: file_name x
# metacluster x n/fraction, plus mouse_ID/group passthrough columns. Cluster
# "1" gets a deliberate abundance shift for the non-reference group(s), so
# tests can assert it has the smallest p_val without pinning exact values.
freq_input <- function(groups = c("control", "control", "control", "control",
                                   "treated", "treated", "treated", "treated")) {
  mice  <- paste0("mouse", seq_along(groups))
  files <- paste0(mice, ".fcs")
  clusters <- factor(1:4)

  base <- tidyr::expand_grid(file_name = files, metacluster = clusters)
  set.seed(7)
  base$n <- stats::rpois(nrow(base), lambda = 50)

  shifted_files <- files[groups != groups[1]]
  base$n[base$metacluster == "1" & base$file_name %in% shifted_files] <-
    base$n[base$metacluster == "1" & base$file_name %in% shifted_files] + 40L

  lookup <- tibble::tibble(file_name = files, mouse_ID = mice, group = factor(groups))
  base <- dplyr::left_join(base, lookup, by = "file_name")

  totals <- base |> dplyr::group_by(file_name) |> dplyr::summarise(total = sum(n), .groups = "drop")
  base |>
    dplyr::left_join(totals, by = "file_name") |>
    dplyr::mutate(fraction = n / total) |>
    dplyr::select(!total)
}

test_that("facs_test_cluster_abundance() flags the shifted cluster with method='glmm'", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "glmm")

  expect_true(all(c("metacluster", "contrast", "p_val", "p_adj") %in% names(result)))
  expect_equal(nrow(result), 4L)
  expect_equal(unique(result$contrast), "treated_vs_control")

  shifted_p <- result$p_val[result$metacluster == "1"]
  expect_true(shifted_p == min(result$p_val))
})

test_that("facs_test_cluster_abundance() supports random effects with method='glmm'", {
  result <- facs_test_cluster_abundance(
    freq_input(), fixed = "group", random = "mouse_ID", method = "glmm"
  )
  expect_true(all(c("metacluster", "contrast", "p_val", "p_adj") %in% names(result)))
  expect_equal(nrow(result), 4L)
})

test_that("facs_test_cluster_abundance() flags the shifted cluster with method='edgeR'", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "edgeR")

  expect_true(all(c("metacluster", "contrast", "logFC", "logCPM", "p_val", "p_adj") %in% names(result)))
  shifted_p <- result$p_val[result$metacluster == "1"]
  expect_true(shifted_p == min(result$p_val))
})

test_that("facs_test_cluster_abundance() flags the shifted cluster with method='voom'", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "voom")

  expect_true(all(c("metacluster", "contrast", "logFC", "t", "p_val", "p_adj") %in% names(result)))
  shifted_p <- result$p_val[result$metacluster == "1"]
  expect_true(shifted_p == min(result$p_val))
})

test_that("facs_test_cluster_abundance() tests one contrast per non-reference level", {
  groups <- rep(c("low", "mid", "high"), length.out = 8)
  result <- facs_test_cluster_abundance(
    freq_input(groups), fixed = "group", ref_level = "mid", method = "glmm"
  )

  expect_equal(nrow(result), 8L)
  expect_setequal(unique(result$contrast), c("low_vs_mid", "high_vs_mid"))
})

test_that("facs_test_cluster_abundance() defaults ref_level to the first factor level", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "glmm")
  expect_equal(unique(result$contrast), "treated_vs_control")
})
```

- [ ] **Step 4: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_test.R")'`
Expected: FAIL with `could not find function "facs_test_cluster_abundance"`

- [ ] **Step 5: Implement `facs_test_cluster_abundance()`**

Create `R/facs_test.R`:

```r
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

  fixed_vec <- freq_data[[fixed]]
  if (!is.factor(fixed_vec)) fixed_vec <- factor(fixed_vec)
  if (is.null(ref_level)) ref_level <- levels(fixed_vec)[1]

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
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_test.R")'`
Expected: PASS (6/6). `lme4` convergence warnings ("boundary (singular) fit") are expected and not failures.

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION R/facs_test.R tests/testthat/test-facs_test.R
git commit -m "feat: add facs_test_cluster_abundance() core (glmm/edgeR/voom)"
```

---

### Task 2: Validation, `random`/method restriction, and `meta` join

**Files:**
- Modify: `R/facs_test.R`
- Test: `tests/testthat/test-facs_test.R`

**Interfaces:**
- Consumes: `facs_test_cluster_abundance()` from Task 1 (same signature); `meta_annotate(data, meta, by = "mouse_ID")` from `R/meta_wrangle.R` (existing, unmodified).
- Produces: the same function, now validating `fixed`/`random`/`cluster_col` presence, `fixed`'s level count, `ref_level` validity, the `random`+non-GLMM restriction, `NA` handling after a `meta` join, and honoring `meta`/`by` when `meta` is supplied.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_test.R`:

```r
test_that("facs_test_cluster_abundance() errors when fixed is missing from freq_data", {
  expect_error(
    facs_test_cluster_abundance(freq_input(), fixed = "NotAColumn", method = "glmm"),
    "NotAColumn"
  )
})

test_that("facs_test_cluster_abundance() errors when random is missing from freq_data", {
  expect_error(
    facs_test_cluster_abundance(freq_input(), fixed = "group", random = "NotAColumn", method = "glmm"),
    "NotAColumn"
  )
})

test_that("facs_test_cluster_abundance() errors when cluster_col is missing from freq_data", {
  expect_error(
    facs_test_cluster_abundance(freq_input(), fixed = "group", cluster_col = "NotAColumn", method = "glmm"),
    "NotAColumn"
  )
})

test_that("facs_test_cluster_abundance() errors when fixed has fewer than 2 levels", {
  dat <- freq_input() |> dplyr::mutate(group = factor("control"))
  expect_error(
    facs_test_cluster_abundance(dat, fixed = "group", method = "glmm"),
    "at least 2 levels"
  )
})

test_that("facs_test_cluster_abundance() errors when ref_level is not among fixed's levels", {
  expect_error(
    facs_test_cluster_abundance(freq_input(), fixed = "group", ref_level = "NotALevel", method = "glmm"),
    "NotALevel"
  )
})

test_that("facs_test_cluster_abundance() errors when random is set with method != 'glmm'", {
  expect_error(
    facs_test_cluster_abundance(freq_input(), fixed = "group", random = "mouse_ID", method = "edgeR"),
    "only supported with method"
  )
})

test_that("facs_test_cluster_abundance() joins meta and errors on NA after an unmatched key", {
  dat  <- freq_input() |> dplyr::select(!group)
  meta <- tibble::tibble(mouse_ID = paste0("mouse", 1:3), group = factor(c("control", "control", "treated")))

  expect_error(
    suppressWarnings(
      facs_test_cluster_abundance(dat, meta = meta, fixed = "group", by = "mouse_ID", method = "glmm")
    ),
    "NA"
  )
})

test_that("facs_test_cluster_abundance() via meta join matches the passthrough-column path", {
  dat_with_group <- freq_input()
  dat_without_group <- dplyr::select(dat_with_group, !group)
  meta <- dplyr::distinct(dat_with_group, mouse_ID, group)

  result_passthrough <- facs_test_cluster_abundance(dat_with_group, fixed = "group", method = "glmm")
  result_via_meta <- facs_test_cluster_abundance(
    dat_without_group, meta = meta, fixed = "group", by = "mouse_ID", method = "glmm"
  )

  expect_equal(result_passthrough$p_val, result_via_meta$p_val)
})
```

- [ ] **Step 2: Run the tests to verify the new ones fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_test.R")'`
Expected: the 8 new tests FAIL (missing-column tests fail because the function does not yet validate; the `meta`-join tests fail with "argument \"meta\" ... not used" style errors since no join is wired in the way the test expects, or with unrelated errors -- either way, not passing for the intended reason).

- [ ] **Step 3: Add validation and the meta join to `facs_test_cluster_abundance()`**

In `R/facs_test.R`, replace everything from the function's opening line
(`facs_test_cluster_abundance <- function(...`) through Task 1's
`if (is.null(ref_level)) ref_level <- levels(fixed_vec)[1]` line
inclusive -- i.e. every line up to but NOT including
`experiment_info <- freq_data |> ...` -- with the following expanded block.
Everything from `experiment_info <- freq_data |> ...` onward is unchanged
from Task 1:

```r
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

```

(The rest of the function -- `experiment_info <- freq_data |> ...` onward -- stays as written in Task 1.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_test.R")'`
Expected: PASS (14/14 total). The `meta`-join NA test expects a `warning()` from `meta_annotate()` (about the unmatched key) followed by a `stop()` from this function -- `suppressWarnings()` around the `expect_error()` call is why Step 1's test wraps it that way.

- [ ] **Step 5: Commit**

```bash
git add R/facs_test.R tests/testthat/test-facs_test.R
git commit -m "feat: add validation and meta join to facs_test_cluster_abundance()"
```

---

### Task 3: `facs_plot_cluster_abundance()`

**Files:**
- Modify: `R/facs_plot.R` (add the new function alongside `facs_plot_umap()`)
- Test: `tests/testthat/test-facs_plot.R`

**Interfaces:**
- Consumes: nothing from Tasks 1-2 except test-time convenience (the "significant_only" test builds a synthetic `test_result` tibble shaped like `facs_test_cluster_abundance()`'s output: `{cluster_col}`, `contrast`, `p_val`, `p_adj`).
- Produces: `facs_plot_cluster_abundance(freq_data, test_result = NULL, group_col, cluster_col = "metacluster", significant_only = FALSE, p_adj_threshold = 0.05)` -- returns a `ggplot` object.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_plot.R`:

```r
abundance_input <- function() {
  tibble::tibble(
    file_name   = rep(c("a.fcs", "b.fcs", "c.fcs", "d.fcs"), each = 3),
    metacluster = factor(rep(1:3, 4)),
    fraction    = runif(12),
    group       = rep(c("control", "treated"), each = 6)
  )
}

test_that("facs_plot_cluster_abundance() returns a ggplot object", {
  p <- facs_plot_cluster_abundance(abundance_input(), group_col = "group")
  expect_s3_class(p, "ggplot")
})

test_that("facs_plot_cluster_abundance() errors when group_col is not a column", {
  expect_error(
    facs_plot_cluster_abundance(abundance_input(), group_col = "NotAColumn"),
    "NotAColumn"
  )
})

test_that("facs_plot_cluster_abundance() errors when cluster_col is not a column", {
  expect_error(
    facs_plot_cluster_abundance(abundance_input(), group_col = "group", cluster_col = "NotAColumn"),
    "NotAColumn"
  )
})

test_that("facs_plot_cluster_abundance() facets by cluster_col", {
  p <- facs_plot_cluster_abundance(abundance_input(), group_col = "group")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 3L)
})

test_that("facs_plot_cluster_abundance() filters to significant clusters when significant_only = TRUE", {
  test_result <- tibble::tibble(
    metacluster = factor(1:3),
    contrast    = "treated_vs_control",
    p_val       = c(0.001, 0.2, 0.5),
    p_adj       = c(0.003, 0.3, 0.5)
  )

  p <- facs_plot_cluster_abundance(
    abundance_input(),
    test_result      = test_result,
    group_col        = "group",
    significant_only = TRUE,
    p_adj_threshold  = 0.05
  )
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 1L)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_plot.R")'`
Expected: FAIL with `could not find function "facs_plot_cluster_abundance"`

- [ ] **Step 3: Implement `facs_plot_cluster_abundance()`**

Append to `R/facs_plot.R`:

```r
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_plot.R")'`
Expected: PASS (all `facs_plot_umap()` tests plus the 5 new ones).

- [ ] **Step 5: Commit**

```bash
git add R/facs_plot.R tests/testthat/test-facs_plot.R
git commit -m "feat: add facs_plot_cluster_abundance()"
```

---

### Task 4: Documentation, full check, and CLAUDE.md updates

**Files:**
- Modify: `man/` (regenerated by roxygen -- do not hand-edit)
- Modify: `NAMESPACE` (regenerated by roxygen -- do not hand-edit)
- Modify: `CLAUDE.md` (append new `Known check output` bullets, update Stage 4 status if referenced)

**Interfaces:**
- Consumes: both functions from Tasks 1-3, already fully implemented and tested.
- Produces: nothing new -- this task verifies the whole branch and documents remaining `devtools::check()` NOTEs.

- [ ] **Step 1: Regenerate documentation**

Run: `Rscript -e 'devtools::document()'`
Expected: exits 0; `man/facs_test_cluster_abundance.Rd` and updated `man/facs_plot_cluster_abundance` entry (or new file) created; `NAMESPACE` gains `export(facs_test_cluster_abundance)` and `export(facs_plot_cluster_abundance)`.

- [ ] **Step 2: Run the full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: 0 failing (all of `test-facs_test.R` and the augmented `test-facs_plot.R` passing, alongside every pre-existing test file).

- [ ] **Step 3: Run devtools::check() and record new NOTEs**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Read the NOTE list; for any new NOTE caused by `facs_test.R`/`facs_plot.R` bare-column tidy-eval references (e.g. `cluster_id`, `p_adj`, `fraction`, `group_col`/`cluster_col`-as-variable-holding-a-name), confirm it is the same false-positive pattern documented for other files (bare column name inside a `dplyr`/`tidyr` verb, or a formal argument referenced bare and resolved by codetools via scope) -- do not attempt to "fix" these away.

- [ ] **Step 4: Append a `CLAUDE.md` bullet documenting the new NOTEs**

In `CLAUDE.md`'s "Known check output" section, add (adjusting exact NOTE wording to match what Step 3 actually produced):

```markdown
  - All `facs_test_cluster_abundance` variable-binding notes (`cluster_id`) -- bare column name inside `dplyr::rename()`'s left-hand side reference and `dplyr::mutate()`, same false-positive pattern as elsewhere. `facs_plot_cluster_abundance` variable-binding notes (`fraction`, `p_adj`) -- bare column names inside `ggplot2::aes()`/`dplyr::filter()`, same pattern as `facs_plot_umap`.
```

- [ ] **Step 5: Commit**

```bash
git add man/ NAMESPACE CLAUDE.md
git commit -m "docs: regenerate roxygen docs for Stage 4, update CLAUDE.md check notes"
```
