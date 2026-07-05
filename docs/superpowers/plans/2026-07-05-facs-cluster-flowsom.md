# FlowSOM Clustering (Stage 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `facs_cluster_flowsom()` (FlowSOM clustering of Stage 1's per-event tibble) and `facs_calc_cluster_freq()` (per-sample cluster frequency table) to the package, per `docs/superpowers/specs/2026-07-05-facs-cluster-flowsom-design.md`.

**Architecture:** New file `R/facs_cluster.R` with the two exported functions. `facs_cluster_flowsom()` builds a `flowCore::flowFrame` from the input tibble's marker columns and runs `FlowSOM::FlowSOM()`, appending `cluster`/`metacluster` columns. `facs_calc_cluster_freq()` aggregates those assignments to a zero-filled long tibble of per-sample frequencies, auto-detecting keyword columns to carry through via a constant-within-`file_name` check (same pattern already used by `calc_restim_proportions_()` in `R/facs_calc.R`).

**Tech Stack:** R, `FlowSOM` (new Bioconductor dependency), `flowCore` (already a dependency), `dplyr`/`tidyr`/`purrr`/`glue`/`tibble`, `testthat` (edition 3).

## Global Constraints

- Namespace every external call as `pkg::fn()` — no bare `dplyr`/`purrr`/`tidyr` calls (CLAUDE.md).
- Use `|>` (base pipe) everywhere in new code.
- Every exported function needs `@param`, `@returns`, `@export`; `@examples` must use `\dontrun{}` since both functions need real fixture data.
- After any roxygen change, run `devtools::document()` before committing.
- Target `devtools::check()`: 0 errors, 0 warnings (new NOTEs for bare-column tidy-eval false positives are expected and acceptable, matching the pattern already documented in CLAUDE.md's "Known check output" section for `facs_calc.R`).
- No mocking `FlowSOM`/`flowCore` — tests chain off the real `tests/fixtures/Treg.wsp` + `tests/fixtures/Treg/*.fcs` fixture via `facs_read_fcs_gated()`, same no-mock rule as every other `facs_` domain function. Use small SOM grids (e.g. 2x2/3x3) and `max_events` downsampling in tests to keep `R CMD check` runtime reasonable.
- `FlowSOM::FlowSOM()`'s exact call signature/return structure should be verified against the real fixture during implementation (flagged as a known risk in the design doc) — if the installed version's API differs from what's written below, adjust the extraction calls and note the deviation in a code comment.

---

### Task 1: `facs_cluster_flowsom()`

**Files:**
- Create: `R/facs_cluster.R`
- Modify: `DESCRIPTION:16-18` (insert `FlowSOM,` between `flowCore,` and `flowWorkspace,`, alphabetical order)
- Test: `tests/testthat/test-facs_cluster.R`

**Interfaces:**
- Consumes: `facs_read_fcs_gated()`'s output shape (`R/facs_read_fcs.R`) — tibble with `file_name` (chr), marker columns (`dbl`), optional keyword columns (`chr`). Not called directly by this task, but test fixtures are built by chaining it.
- Produces: `facs_cluster_flowsom(data, markers = NULL, grid_xdim = 10, grid_ydim = 10, n_metaclusters = 10, seed = NULL)` returning `data` with `cluster` (integer) and `metacluster` (factor) columns appended. Task 2 consumes this exact output shape.

- [ ] **Step 1: Add `FlowSOM` to `DESCRIPTION`**

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
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stringr,
    tibble,
    tidyr,
    tools,
    xml2,
    xfun
```

- [ ] **Step 2: Write the failing tests**

Create `tests/testthat/test-facs_cluster.R`:

```r
library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/Treg.wsp")
fcs_dir  <- testthat::test_path("../fixtures/Treg")
skip_msg <- "Treg fixture not available"
cd45_gate <- "Singlets/Lymphocytes/live/CD45+"

cluster_input <- function(markers = c("CD4", "CD45"), max_events = 200, seed = 1) {
  facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = markers,
    keywords   = c("mouse_ID", "tissue"),
    max_events = max_events,
    seed       = seed
  )
}

test_that("facs_cluster_flowsom() appends cluster and metacluster columns", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  result <- facs_cluster_flowsom(
    dat,
    grid_xdim = 3, grid_ydim = 3, n_metaclusters = 3, seed = 1
  )

  expect_true(all(c("cluster", "metacluster") %in% names(result)))
  expect_type(result$cluster, "integer")
  expect_s3_class(result$metacluster, "factor")
  expect_true(all(result$cluster >= 1L & result$cluster <= 9L))
  expect_equal(nrow(result), nrow(dat))
})

test_that("facs_cluster_flowsom() defaults to clustering on every double column", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  result <- facs_cluster_flowsom(
    dat,
    grid_xdim = 2, grid_ydim = 2, n_metaclusters = 2, seed = 1
  )
  expect_true(all(c("CD4", "CD45", "cluster", "metacluster") %in% names(result)))
})

test_that("facs_cluster_flowsom() clusters on an explicit markers override", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  result <- facs_cluster_flowsom(
    dat, markers = "CD4",
    grid_xdim = 2, grid_ydim = 2, n_metaclusters = 2, seed = 1
  )
  expect_true(all(c("cluster", "metacluster") %in% names(result)))
})

test_that("facs_cluster_flowsom() errors when a marker override is not found in data", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  expect_error(
    facs_cluster_flowsom(dat, markers = c("CD4", "NotAColumn")),
    "NotAColumn"
  )
})

test_that("facs_cluster_flowsom() errors when a marker column contains NA", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "a.fcs", "b.fcs"),
    CD4       = c(1.1, NA, 2.2),
    CD45      = c(3.3, 4.4, 5.5)
  )
  expect_error(
    facs_cluster_flowsom(dat, grid_xdim = 2, grid_ydim = 2, n_metaclusters = 2),
    "CD4"
  )
})

test_that("facs_cluster_flowsom() errors when n_metaclusters exceeds the grid's node count", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    CD45      = c(3.3, 4.4)
  )
  expect_error(
    facs_cluster_flowsom(dat, grid_xdim = 1, grid_ydim = 1, n_metaclusters = 5),
    "n_metaclusters"
  )
})

test_that("facs_cluster_flowsom() is reproducible with the same seed", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  result_a <- facs_cluster_flowsom(dat, grid_xdim = 3, grid_ydim = 3, n_metaclusters = 3, seed = 42)
  result_b <- facs_cluster_flowsom(dat, grid_xdim = 3, grid_ydim = 3, n_metaclusters = 3, seed = 42)

  expect_identical(result_a$cluster, result_b$cluster)
  expect_identical(result_a$metacluster, result_b$metacluster)
})
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_cluster.R")'`
Expected: FAIL with `could not find function "facs_cluster_flowsom"`

- [ ] **Step 4: Implement `facs_cluster_flowsom()`**

Create `R/facs_cluster.R`:

```r
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
```

- [ ] **Step 5: Regenerate docs**

Run: `Rscript -e 'devtools::document()'`
Expected: `Writing NAMESPACE`, `Writing facs_cluster_flowsom.Rd`, no errors.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_cluster.R")'`
Expected: All `facs_cluster_flowsom()` tests PASS. If `FlowSOM::FlowSOM()`'s installed-version signature or `GetClusters()`/`GetMetaclusters()` differ from what's written above, adjust the implementation to match the installed API and add a one-line comment noting the deviation (same treatment `R/facs_read_fcs.R`'s header comment gives `flowjo_to_gatingset()` quirks) — don't silently paper over a mismatch.

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION R/facs_cluster.R tests/testthat/test-facs_cluster.R
git commit -m "feat: add facs_cluster_flowsom() for FlowSOM clustering"
```

---

### Task 2: `facs_calc_cluster_freq()`

**Files:**
- Modify: `R/facs_cluster.R` (append second function)
- Test: `tests/testthat/test-facs_cluster.R` (append tests)

**Interfaces:**
- Consumes: `facs_cluster_flowsom()`'s output from Task 1 — tibble with `file_name`, marker columns, `cluster` (int), `metacluster` (factor), and optional keyword columns.
- Produces: `facs_calc_cluster_freq(data, cluster_col = "metacluster")` returning a long tibble: `file_name`, `cluster_col` (name reused), `n`, `fraction`, plus any column constant within `file_name` carried through. No later task consumes this — it is Stage 2's second deliverable, feeding Stage 4 (not yet planned).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_cluster.R`:

```r
test_that("facs_calc_cluster_freq() returns one row per file_name x metacluster, zero-filled", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input() |>
    facs_cluster_flowsom(grid_xdim = 2, grid_ydim = 2, n_metaclusters = 4, seed = 1)
  freq <- facs_calc_cluster_freq(dat)

  expect_true(all(c("file_name", "metacluster", "n", "fraction") %in% names(freq)))
  expect_equal(nrow(freq), dplyr::n_distinct(dat$file_name) * dplyr::n_distinct(dat$metacluster))
  expect_true(any(freq$n == 0L))
})

test_that("facs_calc_cluster_freq() fractions sum to 1 per file_name", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input() |>
    facs_cluster_flowsom(grid_xdim = 2, grid_ydim = 2, n_metaclusters = 4, seed = 1)
  freq <- facs_calc_cluster_freq(dat)

  totals <- freq |>
    dplyr::group_by(file_name) |>
    dplyr::summarise(total_fraction = sum(fraction), .groups = "drop")
  expect_true(all(abs(totals$total_fraction - 1) < 1e-9))
})

test_that("facs_calc_cluster_freq() carries through keyword columns constant within file_name", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input() |>
    facs_cluster_flowsom(grid_xdim = 2, grid_ydim = 2, n_metaclusters = 4, seed = 1)
  freq <- facs_calc_cluster_freq(dat)

  expect_true(all(c("mouse_ID", "tissue") %in% names(freq)))
  expect_false(any(is.na(freq$mouse_ID)))
})

test_that("facs_calc_cluster_freq() supports cluster_col override", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input() |>
    facs_cluster_flowsom(grid_xdim = 2, grid_ydim = 2, n_metaclusters = 4, seed = 1)
  freq <- facs_calc_cluster_freq(dat, cluster_col = "cluster")

  expect_true("cluster" %in% names(freq))
  expect_true(all(freq$cluster %in% 1:4))
})

test_that("facs_calc_cluster_freq() errors when cluster_col is not found in data", {
  dat <- tibble::tibble(file_name = "a.fcs", metacluster = factor(1))
  expect_error(
    facs_calc_cluster_freq(dat, cluster_col = "not_a_col"),
    "not_a_col"
  )
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_cluster.R")'`
Expected: New tests FAIL with `could not find function "facs_calc_cluster_freq"`; Task 1's tests still PASS.

- [ ] **Step 3: Implement `facs_calc_cluster_freq()`**

Append to `R/facs_cluster.R`:

```r
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
    tidyr::complete(file_name, dplyr::all_of(cluster_col), fill = list(n = 0L)) |>
    dplyr::left_join(totals, by = "file_name") |>
    dplyr::mutate(fraction = n / total) |>
    dplyr::select(!total) |>
    dplyr::left_join(passthrough, by = "file_name") |>
    dplyr::arrange(file_name, .data[[cluster_col]])
}
```

- [ ] **Step 4: Regenerate docs**

Run: `Rscript -e 'devtools::document()'`
Expected: `Writing facs_calc_cluster_freq.Rd`, no errors.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_cluster.R")'`
Expected: All tests in the file PASS.

- [ ] **Step 6: Commit**

```bash
git add R/facs_cluster.R tests/testthat/test-facs_cluster.R man/facs_calc_cluster_freq.Rd NAMESPACE
git commit -m "feat: add facs_calc_cluster_freq() for per-sample cluster frequencies"
```

---

### Task 3: Full package check and CLAUDE.md update

**Files:**
- Modify: `CLAUDE.md` (dependency list + "Known check output" section)

**Interfaces:**
- Consumes: Both functions from Tasks 1–2. No new interfaces produced — this task verifies the whole package still checks cleanly and documents any new NOTEs.

- [ ] **Step 1: Run the full check**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Note any new NOTEs (expected: bare-column-in-dplyr-verb false positives for `facs_cluster_flowsom`/`facs_calc_cluster_freq`, e.g. `markers`, `file_name`, `n`, `total`, `n_distinct` — same category already documented for every other `facs_` function).

- [ ] **Step 2: Update `CLAUDE.md`**

In the "Dependency philosophy" section, add `FlowSOM` to the `Current Imports` list (alphabetical, between `flowCore` and `flowWorkspace`).

In the "Known check output" section, add a bullet describing the new NOTEs from Step 1, following the exact style of the existing `facs_calc_pct_of`/`facs_calc_count_per_g` bullet — name the actual flagged identifiers observed in Step 1's output (fill in from the real `devtools::check()` output, not guessed).

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: record facs_cluster_flowsom()/facs_calc_cluster_freq() in CLAUDE.md"
```
