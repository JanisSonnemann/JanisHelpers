# UMAP Dimensionality Reduction (Stage 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `facs_reduce_umap()` (UMAP embedding of Stage 1's per-event tibble) and `facs_plot_umap()` (ggplot2 visualization of that embedding) to the package, per `docs/superpowers/specs/2026-07-06-facs-reduce-umap-design.md`.

**Architecture:** New file `R/facs_reduce.R` with `facs_reduce_umap()` (runs `uwot::umap()` on the input tibble's marker columns, with optional per-`file_name`-stratified downsampling, appending `UMAP1`/`UMAP2`). New file `R/facs_plot.R` with `facs_plot_umap()` (a `ggplot2` scatter of that embedding, colored/faceted by caller-chosen columns, using `ggplot2`'s built-in viridis scales).

**Tech Stack:** R, `uwot` (new CRAN dependency), `ggplot2` (new CRAN dependency), `dplyr`/`purrr`/`glue`/`tibble` (already dependencies), `testthat` (edition 3).

## Global Constraints

- Namespace every external call as `pkg::fn()` — no bare `dplyr`/`purrr`/`ggplot2`/`uwot` calls (CLAUDE.md).
- Use `|>` (base pipe) everywhere in new code.
- Every exported function needs `@param`, `@returns`, `@export`; `@examples` must use `\dontrun{}` since both functions need real fixture data (or Stage 3's own output).
- After any roxygen change, run `devtools::document()` before committing.
- Target `devtools::check()`: 0 errors, 0 warnings. New NOTEs are expected and acceptable — bare-column tidy-eval false positives (`UMAP1`, `UMAP2`, `file_name` in `R/facs_reduce.R`; `.data`, `UMAP1`, `UMAP2` in `R/facs_plot.R`), matching the pattern already documented in CLAUDE.md's "Known check output" section. Task 3 adds a bullet there for these two files.
- No mocking `uwot`/`ggplot2` — tests run the real packages, same no-mock rule as every other `facs_` domain function. `facs_reduce_umap()`'s primary tests chain off the real `tests/fixtures/Treg.wsp` + `tests/fixtures/Treg/*.fcs` fixture via `facs_read_fcs_gated()` (Stage 1, already merged to `master`); validation-only edge cases use small synthetic tibbles (same pattern several `facs_cluster_flowsom()` tests already use).
- **Stage 2 status:** `facs_cluster_flowsom()`/`facs_calc_cluster_freq()` (Stage 2, PR #1) merged to `master` on 2026-07-06, so this worktree (rebased onto the post-merge `master`) already has `R/facs_cluster.R`, `FlowSOM`, and `parallel` in `DESCRIPTION`. This plan still keeps its tests decoupled from Stage 2 rather than reworking them to depend on it: `facs_reduce_umap()` selects marker columns via `is.double()`, and Stage 2's `cluster` (integer) / `metacluster` (factor) columns are not `double`, so they're automatically excluded from the default marker set with no special-casing required. `facs_plot_umap()`'s default `color_by = "metacluster"` is just a column-name convention — its tests build a small synthetic tibble containing a `metacluster` column directly (`tibble::tibble(..., metacluster = factor(...))`) rather than calling `facs_cluster_flowsom()`, keeping those tests fast and independent of FlowSOM's runtime.
- **Implementation risk to verify against the real fixture** (same "don't assume, verify" treatment Stage 1 gave `flowjo_to_gatingset()` and Stage 2 gave `FlowSOM`/`ConsensusClusterPlus`): `uwot::umap()` requires enough rows relative to `n_neighbors` (roughly `n_neighbors < nrow(data)`) — this plan's tests deliberately pass small `n_neighbors` values on small/downsampled inputs. Also, per `uwot`'s own documentation, exact run-to-run reproducibility under a fixed `set.seed()` is only guaranteed with single-threaded execution — `facs_reduce_umap()` therefore forces `n_threads = 1` whenever `seed` is supplied. If the installed `uwot` version's behavior differs from either assumption, adjust the implementation/tests and note the deviation in a code comment (same treatment as the header comments in `R/facs_read_fcs.R` and `R/facs_cluster.R`).
- If `uwot`/`ggplot2` aren't already installed locally, install them (`install.packages(c("uwot", "ggplot2"))`) before running tests.

---

### Task 1: `facs_reduce_umap()`

**Files:**
- Create: `R/facs_reduce.R`
- Modify: `DESCRIPTION` (insert `uwot,` between `tools,` and `xml2,`)
- Test: `tests/testthat/test-facs_reduce.R`

**Interfaces:**
- Consumes: `facs_read_fcs_gated()`'s output shape (`R/facs_read_fcs.R`) — tibble with `file_name` (chr), marker columns (`dbl`), optional keyword columns (`chr`). Not called directly by this task, but test fixtures are built by chaining it.
- Produces: `facs_reduce_umap(data, markers = NULL, max_events = NULL, n_neighbors = 15, min_dist = 0.1, seed = NULL)` returning `data` (or its downsampled subset) with `UMAP1` (dbl) and `UMAP2` (dbl) columns appended. Task 2 consumes this exact output shape (specifically, the `UMAP1`/`UMAP2` column names).

- [ ] **Step 1: Add `uwot` to `DESCRIPTION`**

Edit `DESCRIPTION` so the `Imports:` block reads:

```
Imports:
    CytoML,
    dplyr,
    fcexpr,
    flowCore,
    FlowSOM,
    flowWorkspace,
    glue,
    gt,
    gtsummary,
    janitor,
    parallel,
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stringr,
    tibble,
    tidyr,
    tools,
    uwot,
    xml2,
    xfun
```

- [ ] **Step 2: Write the failing tests**

Create `tests/testthat/test-facs_reduce.R`:

```r
library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/Treg.wsp")
fcs_dir  <- testthat::test_path("../fixtures/Treg")
skip_msg <- "Treg fixture not available"
cd45_gate <- "Singlets/Lymphocytes/live/CD45+"

reduce_input <- function(markers = c("CD4", "CD45", "TCRb"), max_events = 50, seed = 1) {
  facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = markers,
    keywords   = c("mouse_ID", "tissue"),
    max_events = max_events,
    seed       = seed
  )
}

test_that("facs_reduce_umap() appends UMAP1 and UMAP2 columns", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input()
  result <- facs_reduce_umap(dat, n_neighbors = 5, seed = 1)

  expect_true(all(c("UMAP1", "UMAP2") %in% names(result)))
  expect_type(result$UMAP1, "double")
  expect_type(result$UMAP2, "double")
  expect_equal(nrow(result), nrow(dat))
})

test_that("facs_reduce_umap() defaults to embedding on every double column", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input()
  result <- facs_reduce_umap(dat, n_neighbors = 5, seed = 1)
  expect_true(all(c("CD4", "CD45", "TCRb", "UMAP1", "UMAP2") %in% names(result)))
})

test_that("facs_reduce_umap() embeds on an explicit markers override", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input()
  result <- facs_reduce_umap(dat, markers = c("CD4", "CD45"), n_neighbors = 5, seed = 1)
  expect_true(all(c("UMAP1", "UMAP2") %in% names(result)))
})

test_that("facs_reduce_umap() errors when a marker override is not found in data", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    CD45      = c(3.3, 4.4)
  )
  expect_error(
    facs_reduce_umap(dat, markers = c("CD4", "NotAColumn")),
    "NotAColumn"
  )
})

test_that("facs_reduce_umap() errors when an explicit marker is not double-typed", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    mouse_ID  = c("M1", "M2")
  )
  expect_error(
    facs_reduce_umap(dat, markers = c("CD4", "mouse_ID")),
    "mouse_ID"
  )
})

test_that("facs_reduce_umap() errors when fewer than 2 markers resolve", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2)
  )
  expect_error(facs_reduce_umap(dat), "at least 2")
})

test_that("facs_reduce_umap() errors when a marker column contains NA", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "a.fcs", "b.fcs"),
    CD4       = c(1.1, NA, 2.2),
    CD45      = c(3.3, 4.4, 5.5)
  )
  expect_error(facs_reduce_umap(dat, n_neighbors = 2), "CD4")
})

test_that("facs_reduce_umap() errors when max_events is not a positive integer", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    CD45      = c(3.3, 4.4)
  )
  expect_error(facs_reduce_umap(dat, max_events = -1), "positive integer")
  expect_error(facs_reduce_umap(dat, max_events = 1.5), "positive integer")
})

test_that("facs_reduce_umap() stratifies downsampling per file_name", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input(max_events = 50)
  n_samples <- length(unique(dat$file_name))
  result <- facs_reduce_umap(dat, max_events = 12, n_neighbors = 5, seed = 1)

  expect_equal(nrow(result), 12L)
  counts <- table(result$file_name)
  expect_equal(length(counts), n_samples)
  expect_true(all(counts == 2L))
})

test_that("facs_reduce_umap() gives short samples all their rows without redistributing the shortfall", {
  dat <- tibble::tibble(
    file_name = c(rep("a.fcs", 2), rep("b.fcs", 10), rep("c.fcs", 10)),
    CD4       = rnorm(22),
    CD45      = rnorm(22)
  )
  result <- facs_reduce_umap(dat, max_events = 15, n_neighbors = 3, seed = 1)

  counts <- table(result$file_name)
  expect_equal(unname(counts[["a.fcs"]]), 2L)
  expect_lt(nrow(result), 15L)
})

test_that("facs_reduce_umap() drops samples entirely once max_events is smaller than n_samples", {
  file_names <- sprintf("s%02d.fcs", 1:20)
  dat <- tibble::tibble(
    file_name = file_names,
    CD4       = rnorm(20),
    CD45      = rnorm(20)
  )
  result <- facs_reduce_umap(dat, max_events = 16, n_neighbors = 5, seed = 1)

  expect_equal(nrow(result), 16L)
  expect_setequal(unique(result$file_name), file_names[1:16])
})

test_that("facs_reduce_umap() skips downsampling silently when max_events exceeds nrow", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input(max_events = 20)
  result <- facs_reduce_umap(dat, max_events = 10000L, n_neighbors = 5, seed = 1)
  expect_equal(nrow(result), nrow(dat))
})

test_that("facs_reduce_umap() is reproducible with the same seed", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input(max_events = 50)
  r1 <- facs_reduce_umap(dat, max_events = 12, n_neighbors = 5, seed = 42)
  r2 <- facs_reduce_umap(dat, max_events = 12, n_neighbors = 5, seed = 42)

  expect_equal(r1$file_name, r2$file_name)
  expect_equal(r1$UMAP1, r2$UMAP1)
  expect_equal(r1$UMAP2, r2$UMAP2)
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_reduce.R")'`
Expected: FAIL with `could not find function "facs_reduce_umap"`

- [ ] **Step 4: Implement `facs_reduce_umap()`**

Create `R/facs_reduce.R`:

```r
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
```

- [ ] **Step 5: Regenerate docs**

Run: `Rscript -e 'devtools::document()'`
Expected: `Writing NAMESPACE`, `Writing facs_reduce_umap.Rd`, no errors.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_reduce.R")'`
Expected: All `facs_reduce_umap()` tests PASS. If the installed `uwot` version's `umap()` signature, minimum-row requirements, or threading/reproducibility behavior differ from what's assumed above, adjust the implementation or the affected test's `n_neighbors`/`n_threads` handling and add a one-line comment noting the deviation (same treatment `R/facs_read_fcs.R`'s and `R/facs_cluster.R`'s header comments give their own upstream quirks) — don't silently paper over a mismatch.

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION R/facs_reduce.R tests/testthat/test-facs_reduce.R man/facs_reduce_umap.Rd NAMESPACE
git commit -m "feat: add facs_reduce_umap() for UMAP dimensionality reduction"
```

---

### Task 2: `facs_plot_umap()`

**Files:**
- Create: `R/facs_plot.R`
- Modify: `DESCRIPTION` (insert `ggplot2,` between `flowWorkspace,` and `glue,`; insert `stats,` between `rstatix,` and `stringr,`)
- Test: `tests/testthat/test-facs_plot.R`

**Interfaces:**
- Consumes: `facs_reduce_umap()`'s output shape from Task 1 — any tibble containing `UMAP1`/`UMAP2` (`dbl`) columns, plus whatever other columns exist as `color_by`/`facet_by` targets. Tests build a small synthetic tibble matching this shape directly (see the Stage 2 branch-dependency note in Global Constraints) rather than calling `facs_reduce_umap()` for every test, to keep this task's tests fast and independent of UMAP's actual numerics.
- Produces: `facs_plot_umap(data, color_by = "metacluster", facet_by = NULL)` returning a `ggplot` object. No later task consumes this.

- [ ] **Step 1: Add `ggplot2` and `stats` to `DESCRIPTION`**

Edit `DESCRIPTION` so the `Imports:` block reads:

```
Imports:
    CytoML,
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
    tibble,
    tidyr,
    tools,
    uwot,
    xml2,
    xfun
```

- [ ] **Step 2: Write the failing tests**

Create `tests/testthat/test-facs_plot.R`:

```r
library(testthat)
library(JanisHelpers)

umap_input <- function() {
  tibble::tibble(
    file_name   = rep(c("a.fcs", "b.fcs"), each = 5),
    UMAP1       = rnorm(10),
    UMAP2       = rnorm(10),
    CD4         = rnorm(10),
    metacluster = factor(rep(1:2, 5)),
    group       = rep(c("control", "treated"), 5)
  )
}

test_that("facs_plot_umap() returns a ggplot object", {
  p <- facs_plot_umap(umap_input())
  expect_s3_class(p, "ggplot")
})

test_that("facs_plot_umap() errors when UMAP1/UMAP2 are missing", {
  dat <- umap_input() |> dplyr::select(!UMAP1)
  expect_error(facs_plot_umap(dat), "UMAP1")
})

test_that("facs_plot_umap() errors when color_by is not a column", {
  expect_error(facs_plot_umap(umap_input(), color_by = "NotAColumn"), "NotAColumn")
})

test_that("facs_plot_umap() errors when facet_by is not a column", {
  expect_error(facs_plot_umap(umap_input(), facet_by = "NotAColumn"), "NotAColumn")
})

test_that("facs_plot_umap() uses a continuous viridis scale for a double color_by", {
  p <- facs_plot_umap(umap_input(), color_by = "CD4")
  expect_s3_class(p$scales$get_scales("colour"), "ScaleContinuous")
})

test_that("facs_plot_umap() uses a discrete viridis scale for the default metacluster color_by", {
  p <- facs_plot_umap(umap_input())
  expect_s3_class(p$scales$get_scales("colour"), "ScaleDiscrete")
})

test_that("facs_plot_umap() facets by facet_by", {
  p <- facs_plot_umap(umap_input(), facet_by = "group")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 2L)
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_plot.R")'`
Expected: FAIL with `could not find function "facs_plot_umap"`

- [ ] **Step 4: Implement `facs_plot_umap()`**

Create `R/facs_plot.R`:

```r
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
```

- [ ] **Step 5: Regenerate docs**

Run: `Rscript -e 'devtools::document()'`
Expected: `Writing NAMESPACE`, `Writing facs_plot_umap.Rd`, no errors.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_plot.R")'`
Expected: All `facs_plot_umap()` tests PASS. If the installed `ggplot2` version's internal scale class names (`ScaleContinuous`/`ScaleDiscrete`) or `ggplot_build()` layout structure differ from what's assumed above, adjust the two affected tests and add a one-line comment noting the deviation.

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION R/facs_plot.R tests/testthat/test-facs_plot.R man/facs_plot_umap.Rd NAMESPACE
git commit -m "feat: add facs_plot_umap() for UMAP visualization"
```

---

### Task 3: Full-suite verification and CLAUDE.md update

**Files:**
- Modify: `CLAUDE.md` (append a bullet to the "Known check output" section)

**Interfaces:**
- Consumes: Both functions from Tasks 1-2, plus the full existing test suite. No new interfaces produced.

- [ ] **Step 1: Run the full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests PASS (including the pre-existing suite), 0 failures.

- [ ] **Step 2: Run `devtools::check()`**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Compare the NOTEs against CLAUDE.md's existing "Known check output" list — any genuinely new NOTE should only be the expected bare-column tidy-eval false positives for `UMAP1`/`UMAP2`/`file_name` (in `facs_reduce_umap`/`downsample_stratified_`) and `.data`/`UMAP1`/`UMAP2`/`color_by` (in `facs_plot_umap`). If an unexpected error, warning, or NOTE appears, investigate before proceeding — don't paper over it.

- [ ] **Step 3: Update CLAUDE.md**

Add a new bullet to the "Known check output" section (after the `facs_calc_cluster_freq` bullet, the current last entry — Stage 2's merge already added it after `facs_read_fcs_gated`'s), matching the existing bullet-list style:

```markdown
  - `facs_reduce_umap`/`downsample_stratified_` variable-binding notes (`UMAP1`, `UMAP2`, `file_name`) — bare column names inside base-R subsetting and `dplyr::bind_rows()`, same false-positive pattern as the rest of the `facs_` domain.
  - `facs_plot_umap` variable-binding notes (`UMAP1`, `UMAP2`, `.data`) — bare column names and the rlang tidy-eval pronoun inside `ggplot2::aes()`, same false-positive pattern as `.data` usage elsewhere in the package (see the `analysis_stats.R` bullet above).
```

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: document facs_reduce_umap()/facs_plot_umap() check NOTEs in CLAUDE.md"
```
