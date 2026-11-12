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
