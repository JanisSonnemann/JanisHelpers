# Temporary project layout for prep / export tests:
#   <data_dir>/25-7_exp/raw/FACS/overview_25-7_20261101.wsp  (empty file)
#   <data_dir>/25-7_exp/raw/FACS/overview/                   (panel folder)
local_prep_dirs <- function(env = parent.frame()) {
  data_dir <- withr::local_tempdir(.local_envir = env)
  facs_dir <- file.path(data_dir, "25-7_exp", "raw", "FACS")
  panel_dir <- file.path(facs_dir, "overview")
  dir.create(panel_dir, recursive = TRUE)
  wsp_path <- file.path(facs_dir, "overview_25-7_20261101.wsp")
  file.create(wsp_path)
  # Make the wsp clearly older than anything written afterwards.
  Sys.setFileTime(wsp_path, Sys.time() - 3600)
  list(data_dir = data_dir, facs_dir = facs_dir, panel_dir = panel_dir,
       wsp_path = wsp_path)
}

# Long prep data (facs_read_wsp() shape after calculations) for `files`, a
# tibble of file_name, mouse_ID, tissue (+ any keyword columns). Two
# populations x two metrics per file.
fake_prep_data <- function(files, value = 100) {
  pops <- tibble::tibble(
    population_full_path = c("root/live_CD45+", "root/live_CD45+/CD4+"),
    population = factor(c("live_CD45+", "CD4+"), levels = c("live_CD45+", "CD4+"))
  )
  metrics <- tibble::tibble(metric = c("count", "pct_of_live_CD45+"))
  files |>
    dplyr::cross_join(pops) |>
    dplyr::cross_join(metrics) |>
    dplyr::mutate(value = value) |>
    dplyr::relocate(file_name, population_full_path, population, metric, value)
}

# Long QC test data: 2 groups (grp A / B) x n_per_group mice, one spleen file
# each, two batch days (alternating). Populations: live (~1e5 events),
# live/T (~3e4) and a rare live/T/rare (10 events), each with `count` and
# `fraction_of_parent`; plus `median_CD25` for T. Deterministic scatter of
# +-5 %; mouse `outlier` has a T fraction of 0.05 instead of ~0.3.
fake_qc_data <- function(n_per_group = 8, outlier = "m-3") {
  n <- 2 * n_per_group
  mice <- tibble::tibble(
    mouse_ID = paste0("m-", seq_len(n)),
    file_name = paste0("m-", seq_len(n), "_spleen.fcs"),
    tissue = "spleen",
    condition = "none",
    grp = rep(c("A", "B"), each = n_per_group),
    day = as.Date(c("2026-06-15", "2026-06-16"))[(seq_len(n) %% 2) + 1],
    s = 1 + 0.05 * sin(seq_len(n) * 1.7)
  )
  mice |>
    dplyr::mutate(
      live_n = 1e5 * s,
      t_frac = dplyr::if_else(mouse_ID == outlier, 0.05, 0.3 * (2 - s)),
      t_n = live_n * t_frac
    ) |>
    dplyr::reframe(
      population_full_path = c("root/live", "root/live", "root/live/T",
                               "root/live/T", "root/live/T/rare",
                               "root/live/T/rare", "root/live/T"),
      population = factor(c("live", "live", "T", "T", "rare", "rare", "T"),
                          levels = c("live", "T", "rare")),
      metric = c("count", "fraction_of_parent", "count", "fraction_of_parent",
                 "count", "fraction_of_parent", "median_CD25"),
      value = c(live_n, 0.6 * s, t_n, t_frac, 10, 10 / t_n, 900 * s),
      .by = c(file_name, mouse_ID, tissue, condition, grp, day)
    ) |>
    dplyr::relocate(file_name, population_full_path, population, metric, value)
}

# Resolved QC settings for fake_qc_data().
fake_qc_settings <- function(data, ...) {
  JanisHelpers:::qc_resolve_(data, utils::modifyList(
    list(group_col = "grp", batch_col = "day", qc_mfi_populations = "T"),
    list(...)
  ))
}

# Save a prep result for `files` and create an (empty) QC report beside it,
# then load it the way the final export does.
local_loaded_prep <- function(d, files, value = 100, wsp_path = d$wsp_path,
                              condition_keyword = NA_character_) {
  path <- facs_prep_save(
    fake_prep_data(files, value), wsp_path, "overview", d$panel_dir,
    condition_keyword = condition_keyword
  )
  qc <- sub("_prep\\.rds$", "_qc.html",
            file.path(d$panel_dir, "qc", basename(path)))
  dir.create(dirname(qc), showWarnings = FALSE)
  file.create(qc)
  facs_prep_load(wsp_path, d$panel_dir)
}
