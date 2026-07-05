# Parallelize `facs_read_fcs_gated()` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut the wall-clock time of `facs_read_fcs_gated()` by restructuring it to build one single-sample `GatingSet` per file (instead of one shared multi-sample `GatingSet`) and dispatching those per-sample builds across CPU cores via a new opt-in `workers` argument.

**Architecture:** `R/facs_read_fcs.R` is modified in place — no new files. `read_one_sample_()` is replaced by `read_one_sample_solo_()`, which builds its own single-sample `GatingSet` via `CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir, subset = sample_name)` instead of receiving a `GatingHierarchy` carved out of a shared group-wide `GatingSet`. A new helper, `resolve_group_samples_()`, resolves the sample-name list for a `group` directly from the workspace XML metadata (`CytoML::fj_ws_get_sample_groups()` + `CytoML::fj_ws_get_samples()`) without building any `GatingSet` — this is what lets dispatch happen before any of the expensive per-sample work starts. `facs_read_fcs_gated()` gains a `workers = 1L` argument: `workers == 1` dispatches via `purrr::map()` (today's sequential behavior, byte-identical output); `workers > 1` dispatches via a `parallel::makeCluster()` PSOCK cluster + `parallel::parLapply()` (base R, socket-based separate processes, works cross-platform including Windows).

**Revision note (discovered during Task 2 implementation, not known when this plan was first written):** the original design called for `parallel::mclapply()` (fork-based). This was tried first and **segfaults** — `CytoML`/`flowWorkspace`'s C++ objects are not fork-safe (confirmed: even having each forked child re-open its own workspace from scratch segfaults, indicating corrupted state inherited from the pre-fork process, a known class of issue with multithreaded C++ libraries and `fork()`). PSOCK clusters (`parallel::makeCluster()`/`parLapply()`) were verified working end-to-end against the real fixture instead — these spawn genuinely separate fresh R processes rather than forking, sidestepping the fork-safety problem entirely. This requires two structural changes from the original design: (1) `read_one_sample_solo_()` takes `wsp_path` (a plain string) instead of a pre-opened `ws` object, and opens its own copy internally — a C++ external pointer cannot be passed across a process boundary at all, forked or socket-based; (2) `read_one_sample_solo_()`/`resolve_markers_()` must be shipped to each PSOCK worker explicitly via `parallel::clusterExport()` rather than relying on the worker re-loading the package by name, because a fresh worker process may resolve a stale *installed* copy of the package that doesn't match the code currently running in-memory (e.g. under `devtools::load_all()` during development/testing).

**Tech Stack:** R, `CytoML`, `flowWorkspace`, `flowCore` (already `Imports`), `parallel` (ships with every R install, needs to be added to `Imports:` in `DESCRIPTION` since `R CMD check` requires every `pkg::fn()` call to be declared there — this is a correction to the design doc, which said "no dependency change"; the *practical* meaning of that claim still holds, since `parallel` requires no separate install). Testing via `testthat` 3e against the real `tests/fixtures/Treg.wsp` + `tests/fixtures/Treg/*.fcs` fixture (gitignored — copy it into your worktree from the main checkout before starting; see Task 0).

## Global Constraints

- Pipe: `|>` (base pipe) everywhere.
- Namespacing: every external function call inside an exported (or shared unexported) function body must be `pkg::fn()` — no bare `dplyr`/`purrr`/`tibble`/`glue`/`flowWorkspace`/`flowCore`/`CytoML`/`parallel` calls.
- `parallel` must be added to `DESCRIPTION`'s `Imports:` (alphabetical, matching existing order — falls between `janitor` and `purrr`).
- Internal helpers use a trailing underscore (`resolve_group_samples_`, `read_one_sample_solo_`), matching the existing convention.
- `facs_read_fcs_gated()` returns visibly (not `invisible()`) — unchanged.
- Test fixture (`tests/fixtures/Treg.wsp`, `tests/fixtures/Treg/*.fcs`) is `.gitignore`d — every test must `skip_if_not(dir.exists(...))`, matching the existing pattern in `tests/testthat/test-facs_read_fcs.R`.
- `@param`, `@returns`, `@export` required on the exported function; update `@param seed` and add `@param workers`. No roxygen block on unexported helpers — a one-line `#` comment is enough.
- After any roxygen change, run `devtools::document()` before committing.
- Run `devtools::check()` before the final commit — target 0 errors, 0 warnings.
- Non-ASCII characters in R source must be escaped as `\uXXXX`.

---

## Task 0: Restore the gitignored fixture into your worktree

**Files:** none modified — this is a one-time local setup step.

The `Treg` fixture is gitignored (real research metadata) and is not carried into a fresh `git worktree` checkout. Before writing any test, confirm it's present:

- [ ] **Step 1: Check whether the fixture is present**

Run: `ls tests/fixtures/Treg.wsp tests/fixtures/Treg/ 2>&1`

Expected: both exist. If you get "No such file or directory", copy them from the main checkout (the non-worktree `JanisHelpers` directory) into `tests/fixtures/` in your current worktree — same for `tests/fixtures/minimal.wsp` and `tests/fixtures/meta_minimal.xlsx` if those are also missing (they are untracked-but-present in the main checkout, not yet committed).

- [ ] **Step 2: Confirm the baseline test suite is green**

Run: `Rscript -e 'devtools::test(filter = "facs_read_fcs")'`
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 27 ]`

If this doesn't pass before you start, stop and investigate — you'd otherwise have no way to tell your changes apart from pre-existing breakage.

---

## Task 1: Restructure to one solo `GatingSet` per sample (regression-safe refactor, no new behavior)

**Files:**
- Modify: `R/facs_read_fcs.R`

**Interfaces:**
- Produces (unexported, internal to this file): `resolve_group_samples_(ws, group)` — returns a character vector of clean sample file names (matching `pData(gh)$name`) belonging to `group`, resolved from workspace XML metadata alone (no `GatingSet` build).
- Produces (unexported, replaces `read_one_sample_()`): `read_one_sample_solo_(ws, fcs_dir, group, sample_name, sample_index, gate_path, gate_path_norm, markers, max_events, seed)` — builds its own single-sample `GatingSet`, returns a plain tibble (same shape `read_one_sample_()` used to return). `sample_index` is that sample's 1-based position in the list returned by `resolve_group_samples_()` — used to derive a per-sample downsampling seed.
- `facs_read_fcs_gated()`'s public signature is **unchanged** in this task — this task is purely internal restructuring, verified by the existing test suite staying green throughout.

This task deliberately does **not** touch parallel dispatch yet — it only changes *how* each sample's `GatingSet` is built (solo vs. shared-group), keeping dispatch itself sequential (`purrr::map()`, same as today) so the refactor's correctness can be verified against the existing suite in isolation from the parallelism work in Task 2.

- [ ] **Step 1: Replace `read_one_sample_()` and add `resolve_group_samples_()`**

In `R/facs_read_fcs.R`, replace this block:

```r
read_one_sample_ <- function(gh, gate_path, gate_path_norm, markers, max_events) {
  file_name <- flowWorkspace::pData(gh)$name

  gated <- tryCatch(
    flowWorkspace::gh_pop_get_data(gh, gate_path_norm),
    error = function(e) NULL
  )
  if (is.null(gated)) {
    warning(glue::glue(
      "gate_path '{gate_path}' not found for '{file_name}'; sample skipped."
    ))
    return(tibble::tibble())
  }
  gated <- flowWorkspace::realize_view(gated)

  resolved_cols <- resolve_markers_(markers, gh, gated, file_name)

  mat <- flowCore::exprs(gated)[, resolved_cols, drop = FALSE]
  colnames(mat) <- markers

  if (!is.null(max_events) && nrow(mat) > max_events) {
    mat <- mat[sample.int(nrow(mat), max_events), , drop = FALSE]
  }

  tibble::as_tibble(mat) |>
    tibble::add_column(file_name = file_name, .before = 1L)
}
```

with:

```r
# Resolves the sample names belonging to a FlowJo group directly from the
# workspace XML metadata -- no GatingSet build required, so this can run
# once up front before any per-sample work is dispatched.
resolve_group_samples_ <- function(ws, group) {
  groups  <- CytoML::fj_ws_get_sample_groups(ws)
  samples <- CytoML::fj_ws_get_samples(ws)
  ids <- groups$sampleID[groups$groupName == group]
  samples$name[samples$sampleID %in% ids]
}

# Builds its own single-sample GatingSet (verified byte-identical to
# carving the same sample out of a shared multi-sample GatingSet) so the
# per-sample work this function does is independent and parallelizable.
read_one_sample_solo_ <- function(ws, fcs_dir, group, sample_name, sample_index,
                                   gate_path, gate_path_norm, markers,
                                   max_events, seed) {
  gs <- CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir, subset = sample_name)
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
    return(tibble::tibble())
  }
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
```

- [ ] **Step 2: Rewire `facs_read_fcs_gated()`'s body to use the new helpers**

Replace this block (inside `facs_read_fcs_gated()`):

```r
  ws <- CytoML::open_flowjo_xml(wsp_path)
  gs <- CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir)

  if (!is.null(seed)) set.seed(seed)

  data <- purrr::map(
    flowWorkspace::sampleNames(gs),
    function(sn) {
      read_one_sample_(
        gh             = gs[[sn]],
        gate_path      = gate_path,
        gate_path_norm = gate_path_norm,
        markers        = markers,
        max_events     = max_events
      )
    }
  ) |>
    dplyr::bind_rows()
```

with:

```r
  ws <- CytoML::open_flowjo_xml(wsp_path)
  sample_names <- resolve_group_samples_(ws, group)

  data <- purrr::map(
    seq_along(sample_names),
    function(i) {
      read_one_sample_solo_(
        ws             = ws,
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
  ) |>
    dplyr::bind_rows()
```

- [ ] **Step 3: Run the full existing suite, confirm no regressions**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); devtools::test(filter = "facs_read_fcs")'`
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 27 ]` — identical to the Task 0 baseline. If anything fails, the solo-build refactor changed observable behavior; do not proceed to Task 2 until this is green.

- [ ] **Step 4: Commit**

```bash
git add R/facs_read_fcs.R
git commit -m "refactor: build one solo GatingSet per sample in facs_read_fcs_gated()

Replaces the shared multi-sample GatingSet build with independent
per-sample GatingSet construction (verified byte-identical output),
which is what makes the per-sample work parallelizable in the next
commit. No observable behavior change -- full existing suite stays
green."
```

---

## Task 2: Add `workers` argument and parallel dispatch

**Files:**
- Modify: `R/facs_read_fcs.R`
- Modify: `tests/testthat/test-facs_read_fcs.R`

**Interfaces:**
- Produces: `facs_read_fcs_gated(..., workers = 1L)` — new argument. `workers == 1` behaves exactly as after Task 1 (plus the `wsp_path`-per-call change below). `workers > 1` dispatches per-sample work across a PSOCK cluster via `parallel::makeCluster()` + `parallel::parLapply()`.
- Modifies: `read_one_sample_solo_()`'s signature changes from `(ws, fcs_dir, ...)` to `(wsp_path, fcs_dir, ...)` — it now opens its own `ws <- CytoML::open_flowjo_xml(wsp_path)` internally, since a C++ external pointer (what `ws` is) cannot be passed across a process boundary, whether forked or socket-based; every call (sequential or parallel) must independently open the workspace from the plain path string.
- Modifies: `read_one_sample_solo_()` now returns `list(data = tibble, warnings = character())` instead of a bare tibble — necessary because neither fork-based nor PSOCK-based parallel workers propagate `warning()` calls back to the parent process as real conditions; they must be captured in the worker and re-issued in the parent.

**Note on this task's history:** the first implementation attempt used `parallel::mclapply()` (fork-based) per the plan's original wording and hit a hard segfault — `CytoML`/`flowWorkspace`'s C++ objects are not fork-safe, and this could not be worked around (confirmed by testing: even having each forked child independently re-open its own workspace still segfaulted, indicating the corruption is in inherited process/threading state from before the fork, not in passing the object itself). The steps below reflect the corrected, verified-working PSOCK-based design — follow them as written; do not reach for `mclapply()`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/testthat/test-facs_read_fcs.R` (after the existing `max_events`/`seed` test):

```r
test_that("facs_read_fcs_gated() with workers = 2 matches workers = 1 output", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result_seq <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "CD45"),
    keywords  = c("mouse_ID", "tissue"),
    workers   = 1
  )
  result_par <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "CD45"),
    keywords  = c("mouse_ID", "tissue"),
    workers   = 2
  )

  result_seq <- dplyr::arrange(result_seq, file_name, CD4, CD45)
  result_par <- dplyr::arrange(result_par, file_name, CD4, CD45)
  expect_equal(result_seq, result_par)
})

test_that("facs_read_fcs_gated() still errors on an unmatched marker under workers = 2", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  expect_error(
    facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = cd45_gate,
      markers   = c("CD4", "NotARealMarker"),
      workers   = 2
    ),
    "NotARealMarker"
  )
})

test_that("facs_read_fcs_gated() still warns on a missing gate_path under workers = 2", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  warns <- testthat::capture_warnings(
    result <- facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = "Nonexistent/Path",
      markers   = c("CD4"),
      workers   = 2
    )
  )
  expect_true(all(grepl("not found", warns)))
  expect_length(warns, 6L)
  expect_equal(nrow(result), 0L)
})
```

- [ ] **Step 2: Run the new tests, confirm they fail**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); devtools::test(filter = "facs_read_fcs")'`
Expected: FAIL — `unused argument (workers = ...)`, since `facs_read_fcs_gated()` doesn't accept `workers` yet.

- [ ] **Step 3: Change `read_one_sample_solo_()` to open its own workspace from `wsp_path` and capture warnings**

In `R/facs_read_fcs.R`, replace the full body of `read_one_sample_solo_()` (the whole function, signature included) with:

```r
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
```

Two things to note:
1. The `if (is.null(gated)) { ... } else { ... }` restructuring (replacing the old early `return(tibble::tibble())`): a bare `return()` inside the `withCallingHandlers({...})` block would exit `read_one_sample_solo_()` immediately, skipping the final `list(data =, warnings =)` wrap — the if/else must produce the tibble as the block's value instead.
2. This function now takes `wsp_path` instead of `ws` — every caller of `read_one_sample_solo_()` must be updated to pass `wsp_path` instead of a pre-opened `ws` object (done in Step 4 below).

- [ ] **Step 4: Add `workers` argument and dispatch branching to `facs_read_fcs_gated()`**

Replace this block (the one written in Task 1 Step 2):

```r
  ws <- CytoML::open_flowjo_xml(wsp_path)
  sample_names <- resolve_group_samples_(ws, group)

  data <- purrr::map(
    seq_along(sample_names),
    function(i) {
      read_one_sample_solo_(
        ws             = ws,
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
  ) |>
    dplyr::bind_rows()
```

with:

```r
  ws <- CytoML::open_flowjo_xml(wsp_path)
  sample_names <- resolve_group_samples_(ws, group)

  worker_fn <- function(i) {
    read_one_sample_solo_(
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
    # read_one_sample_solo_()/resolve_markers_() must be shipped to each
    # worker by value (not re-resolved by name on the worker), since a
    # PSOCK worker is a fresh R process that may see a different
    # installed copy of this package than the one currently running here
    # (e.g. under devtools::load_all() during development).
    parallel::clusterExport(
      cl,
      varlist = c("read_one_sample_solo_", "resolve_markers_"),
      envir   = environment(read_one_sample_solo_)
    )
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
```

Note: `ws` is still opened once here in the parent process — it's only used for `resolve_group_samples_()` (cheap XML metadata lookup), never passed to `read_one_sample_solo_()` anymore. Each call to `read_one_sample_solo_()` (whether sequential or parallel) re-opens its own `ws` from `wsp_path` internally, adding a small per-sample overhead (~0.04s, per Task 1's profiling) even when `workers == 1` — an intentional, required trade-off: `ws` is a C++ external pointer that cannot cross a process boundary, so keeping one code path that works identically in both dispatch modes is preferable to special-casing the sequential path to reuse a shared `ws`.

Then update the function signature. Replace:

```r
facs_read_fcs_gated <- function(wsp_path,
                                 gate_path,
                                 markers,
                                 keywords = NULL,
                                 fcs_dir = NULL,
                                 group = "All Samples",
                                 max_events = NULL,
                                 seed = NULL) {
```

with:

```r
facs_read_fcs_gated <- function(wsp_path,
                                 gate_path,
                                 markers,
                                 keywords = NULL,
                                 fcs_dir = NULL,
                                 group = "All Samples",
                                 max_events = NULL,
                                 seed = NULL,
                                 workers = 1L) {
```

- [ ] **Step 5: Run the tests, confirm all pass**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); devtools::test(filter = "facs_read_fcs")'`
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 30 ]` (27 existing + 3 new).

- [ ] **Step 6: Commit**

```bash
git add R/facs_read_fcs.R tests/testthat/test-facs_read_fcs.R
git commit -m "feat: parallelize facs_read_fcs_gated() across samples via workers=

Adds an opt-in workers argument (default 1, sequential, unchanged
behavior). workers > 1 dispatches per-sample GatingSet construction
across a PSOCK cluster via parallel::makeCluster()/parLapply(), since
profiling showed GatingSet construction -- not the post-build
extraction step -- dominates runtime and scales linearly per sample.
mclapply() (fork-based) was tried first and segfaults, since CytoML's
C++ objects are not fork-safe; PSOCK workers are separate fresh
processes and sidestep this. Each worker re-opens the workspace from
wsp_path independently (a C++ external pointer can't cross a process
boundary) and read_one_sample_solo_()/resolve_markers_() are shipped
to workers via clusterExport() rather than relying on the worker
re-loading the package by name, so a stale installed copy can't shadow
the code actually running. Warnings are captured per-worker and
re-issued in the parent process (neither fork nor PSOCK workers
propagate warning() conditions back on their own), and per-sample
downsampling uses a seed derived from each sample's position rather
than one shared sequential RNG stream, so results no longer depend on
execution order."
```

---

## Task 3: Documentation and final check

**Files:**
- Modify: `R/facs_read_fcs.R` (roxygen only)
- Modify: `DESCRIPTION`
- Modify: `man/facs_read_fcs_gated.Rd` (regenerated, not hand-edited)
- Modify: `CLAUDE.md` (only if `devtools::check()` reports a new NOTE)

- [ ] **Step 1: Add `parallel` to `DESCRIPTION`'s `Imports:`**

In `DESCRIPTION`, the `Imports:` block currently reads (alphabetical):

```
Imports:
    CytoML,
    dplyr,
    fcexpr,
    flowCore,
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

Add `parallel` between `janitor` and `purrr`:

```
Imports:
    CytoML,
    dplyr,
    fcexpr,
    flowCore,
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
    xml2,
    xfun
```

- [ ] **Step 2: Update the roxygen block on `facs_read_fcs_gated()`**

In `R/facs_read_fcs.R`, find the existing `@param seed` line:

```r
#' @param seed integer; if set, seeds the random draw used by
#'   \code{max_events} for reproducibility.
```

Replace it with:

```r
#' @param seed integer; if set, seeds the random draw used by
#'   \code{max_events} for reproducibility. Each sample's draw is seeded
#'   from \code{seed} offset by that sample's position in its FlowJo
#'   group, so results are reproducible across repeated calls with the
#'   same inputs regardless of \code{workers} -- but the exact per-sample
#'   offset is an internal implementation detail, not a public contract.
#' @param workers integer; number of samples to process in parallel via a
#'   \code{parallel::makeCluster()} PSOCK cluster. Default \code{1}
#'   processes samples sequentially. Works cross-platform (including
#'   Windows) since each worker is a separate process, not a fork --
#'   necessary because \code{CytoML}'s workspace/GatingSet objects are
#'   C++ external pointers that cannot be shared across a process
#'   boundary; each worker re-opens the workspace independently.
```

- [ ] **Step 3: Regenerate documentation**

Run: `Rscript -e 'devtools::document()'`
Expected: no errors; `man/facs_read_fcs_gated.Rd` is updated with the new `@param` entries.

- [ ] **Step 4: Run `devtools::check()`**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Compare the NOTE list against the ones already documented in `CLAUDE.md`'s "Known check output" section. If a genuinely new NOTE appears (e.g. a variable-binding note from `resolve_group_samples_()` or `read_one_sample_solo_()`), add one bullet describing it to that section, following the existing bullet style — do not skip this if a new NOTE appears, but do not add speculative entries if none do.

- [ ] **Step 5: Run the full test suite one more time**

Run: `Rscript -e 'devtools::test()'`
Expected: all suites pass, no regressions outside `test-facs_read_fcs.R`.

- [ ] **Step 6: Commit**

```bash
git add R/facs_read_fcs.R DESCRIPTION man/facs_read_fcs_gated.Rd CLAUDE.md
git commit -m "docs: document workers= argument and add parallel to Imports

devtools::check() requires every pkg::fn() call to be declared in
Imports even for base-distributed packages like parallel."
```

(Drop `CLAUDE.md` from the `git add` if Step 4 found no new NOTEs to document.)
