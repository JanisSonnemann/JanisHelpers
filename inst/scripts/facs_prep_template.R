# FACS prep script for one panel -- copy into the project and adapt.
# Template: system.file("scripts", "facs_prep_template.R", package = "JanisHelpers")
#
# The loop (facs_prep_run()):
#   final <- FALSE  read the .wsp, select populations, calculate, save the
#                   prep result and render the QC report (all samples shown;
#                   exclusions only flagged). Adjust gates in FlowJo (save a
#                   new dated .wsp for major changes) and/or edit the meta
#                   workbook's `exclusions` sheet, then rerun.
#   final <- TRUE   QC passed: export the saved prep result to the database,
#                   dropping and recording the exclusions, then snapshot.
#                   Asks before overwriting a root that is already stored.
#
# Run interactively from the project root (the folder holding data/).
# With several panels per experiment, move the shared part (paths, meta,
# db_open()) into one sourced script and keep only the panel part here.

library(JanisHelpers)

final <- FALSE

# ---- experiment (shared across panels) --------------------------------------

data_dir  <- "data"
exp_dir   <- "26-1_BMT-mRNA-V4"                 # inside data_dir
db_path   <- "~/databases/MPO_mRNA/MPO_mRNA.duckdb"
snapshot_root <- file.path(data_dir, "z_snapshots")

facs_dir   <- file.path(data_dir, exp_dir, "raw", "FACS")
meta       <- meta_read(file.path(data_dir, exp_dir, "meta", "meta_26-1.xlsx"))
mice       <- dplyr::filter(meta$meta, !is.na(mouse_ID))
exclusions <- meta$exclusions                   # NULL if the sheet doesn't exist

# Open the database and register the experiment + subjects (the column
# mapping depends on the meta workbook -- adapt).
db_open <- function() {
  con <- db_connect(path.expand(db_path))
  db_write_experiment(con, "26-1", "BMT-mRNA-V4", project = "MPO_mRNA")
  db_write_subjects(con, dplyr::transmute(
    mice, mouse_id = mouse_ID, experiment_code = "26-1",
    treatment_group = as.character(group)
  ))
  con
}

# ---- panel --------------------------------------------------------------------

wsp       <- "overview-26-1_20260806.wsp"       # in raw/FACS/, <root>_YYYYMMDD.wsp
panel_dir <- "overview"                         # panel folder in raw/FACS/ (holds the .fcs)
assay     <- "overview"                         # panel's assay name in the database

# Populations to keep, in display order, and this panel's calculations.
selection <- c("live_CD45+")

prepare <- function(data) {
  data |>
    dplyr::filter(population %in% selection) |>
    dplyr::mutate(population = factor(population, levels = selection)) |>
    facs_calc_pct_of(ref_pop = "live_CD45+")
}

con <- if (final) db_open()
facs_prep_run(
  final, file.path(facs_dir, wsp), file.path(facs_dir, panel_dir), assay,
  prepare = prepare, mice = mice, exclusions = exclusions, con = con,
  condition_keyword = NA_character_,            # e.g. "restimulation" for multi-well panels
  # QC report: outliers scored within groups, batch plots per staining day.
  group_col = "group", batch_col = "death_date",
  # Mice this panel should contain: all minus dropouts (`early` -- adapt).
  expected_mice = mice$mouse_ID[is.na(mice$early)],
  qc_mfi_populations = character(),             # e.g. "Treg_FoxP3+_CD25hi"
  data_dir = data_dir, snapshot_root = snapshot_root
)
if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
