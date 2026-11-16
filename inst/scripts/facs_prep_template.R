# FACS prep script for one panel -- copy into the project and adapt.
# Template: system.file("scripts", "facs_prep_template.R", package = "JanisHelpers")
#
# One self-contained script per panel; the pipeline logic lives in
# JanisHelpers (facs_prep_run(), db_write_meta()).
#
#   final <- FALSE  read the .wsp, select populations, calculate, save the
#                   prep result and render the QC report (all samples shown;
#                   exclusions only flagged). Adjust gates in FlowJo (save a
#                   new dated .wsp for major changes) and/or edit the meta
#                   workbook's `exclusions` sheet, then rerun.
#   final <- TRUE   QC passed: register experiment + mice, export the saved
#                   prep result to the database (dropping and recording the
#                   exclusions), then snapshot. Asks before overwriting a
#                   root that is already stored.
#
# Run interactively from the project root (the folder holding data/).
# Keep the live .duckdb local; only snapshots go on OneDrive.

library(JanisHelpers)

final <- FALSE

# ---- experiment ---------------------------------------------------------------

exp_path  <- file.path("data", "26-1_BMT-mRNA-V4")
meta      <- meta_read(file.path(exp_path, "meta", "meta_26-1.xlsx"))  # standard template
mice      <- dplyr::filter(meta$meta, !is.na(mouse_ID))
db_path   <- "~/databases/MPO_mRNA/MPO_mRNA.duckdb"

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

# ---- run ------------------------------------------------------------------------

con <- if (final) {
  db_connect(path.expand(db_path)) |>
    db_write_meta(meta, "26-1", "BMT-mRNA-V4", project = "MPO_mRNA")
}
facs_prep_run(
  final,
  wsp_path   = file.path(exp_path, "raw", "FACS", wsp),
  panel_dir  = file.path(exp_path, "raw", "FACS", panel_dir),
  assay      = assay,
  prepare    = prepare,
  mice       = mice,
  exclusions = meta$exclusions,                 # NULL if the sheet doesn't exist
  con        = con,
  condition_keyword = NA_character_,            # e.g. "restimulation" for multi-well panels
  # QC report: outliers scored within groups, batch plots per staining day,
  # completeness against all mice minus dropouts (`early`).
  group_col     = "group",
  batch_col     = "death_date",
  expected_mice = mice$mouse_ID[is.na(mice$early)],
  qc_mfi_populations = character(),             # e.g. "Treg_FoxP3+_CD25hi"
  snapshot_root = file.path("data", "z_snapshots")
)
if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
