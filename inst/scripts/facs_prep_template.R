# FACS prep script for one panel -- copy into the project and adapt.
# Copy a clean version of this template to the current working directory (then rename it):
# file.copy(system.file("scripts", "facs_prep_template.R", package = "JanisHelpers"), ".")
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

# ==============================================================================
# INPUT -- fill in everything in this section
# ==============================================================================

# ---- run mode ----------------------------------------------------------------

final <- FALSE                                  # TRUE = export to the database

# ---- experiment --------------------------------------------------------------

exp_code      <- "26-1"                         # experiment code in the database
exp_name      <- "BMT-mRNA-V4"                  # experiment name in the database
project       <- "MPO_mRNA"                     # project name in the database
exp_dir       <- "26-1_BMT-mRNA-V4"             # experiment folder in data/
meta_file     <- "meta_26-1.xlsx"               # meta workbook in <exp_dir>/meta/
db_path       <- "~/databases/MPO_mRNA/MPO_mRNA.duckdb"
snapshot_root <- file.path("data", "z_snapshots")

# ---- panel -------------------------------------------------------------------

wsp       <- "overview-26-1_20260806.wsp"       # in raw/FACS/, <root>_YYYYMMDD.wsp
panel_dir <- "overview"                         # panel folder in raw/FACS/ (holds the .fcs)
assay     <- "overview"                         # panel's assay name in the database
wsp_group <- NULL                               # FlowJo group to read, e.g. "Samples"; NULL = all samples
condition_keyword <- NA_character_              # e.g. "restimulation" for multi-well panels

# Populations to keep, in display order (workspace spelling -- a name that
# isn't in the workspace stops with the closest match).
# Display renames: facs_select_populations(labels = c(old = "new")).
selection <- c("live_CD45+")

# This panel's calculations. `calc_meta` (volumes + organ weights from the
# meta workbook, see below) is available inside this function.
prepare <- function(data) {
  data <- data |>
    facs_select_populations(selection) |>
    facs_calc_pct_of(ref_pop = "live_CD45+")
  # Counts per g organ (column names = the meta workbook's facs_volumes /
  # organ_weights headers -- check names(calc_meta)):
  # data <- facs_calc_count_per_g(
  #   data,
  #   meta = calc_meta,
  #   vol_total = "total_vol",
  #   vol_stained = "ov_vol",
  #   vol_resuspended = "resuspend_vol",
  #   vol_measured = "measured_vol",
  #   organ_piece_weight = "facs_weight"
  # )
  data
}

# ---- QC report (defaults are usually fine) -------------------------------------

group_col          <- "group"                   # outliers scored within these groups
batch_col          <- "death_date"              # batch plots per staining day
expected_mice      <- NULL                      # NULL = all mice in the meta sheet minus dropouts (`early`)
expected_tissues   <- NULL                      # NULL = tissues in the data
# Which populations / MFIs the QC checks cover (labels after prepare()).
qc_populations     <- NULL                      # hierarchy + count checks; NULL = all selected
qc_mfi_populations <- character()               # MFI check, e.g. "Treg_FoxP3+_CD25hi"; empty = none
qc_metrics         <- NULL                      # e.g. "median_CD25"; NULL = all median_*
min_events         <- 50                        # points outside the whiskers are labelled unless under this many events

# ==============================================================================
# PIPELINE -- do not change anything below this line
# ==============================================================================

exp_path  <- file.path("data", exp_dir)
meta_path <- file.path(exp_path, "meta", meta_file)
meta      <- meta_read(meta_path)
mice      <- dplyr::filter(meta$meta, !is.na(mouse_ID))
calc_meta <- meta_clean(meta)
if (is.null(expected_mice)) {
  dropped <- if ("early" %in% names(mice)) !is.na(mice$early) else FALSE
  expected_mice <- mice$mouse_ID[!dropped]
}

con <- if (final) {
  db_connect(path.expand(db_path)) |>
    db_write_meta(meta, exp_code, exp_name, project = project)
}
facs_prep_run(
  final,
  wsp_path           = file.path(exp_path, "raw", "FACS", wsp),
  panel_dir          = file.path(exp_path, "raw", "FACS", panel_dir),
  assay              = assay,
  prepare            = prepare,
  mice               = mice,
  exclusions         = meta$exclusions,         # NULL if the sheet doesn't exist
  con                = con,
  meta_path          = meta_path,               # final = TRUE also records the exclusions sheet
  wsp_group          = wsp_group,
  condition_keyword  = condition_keyword,
  group_col          = group_col,
  batch_col          = batch_col,
  expected_mice      = expected_mice,
  expected_tissues   = expected_tissues,
  qc_populations     = qc_populations,
  qc_mfi_populations = qc_mfi_populations,
  qc_metrics         = qc_metrics,
  min_events         = min_events,
  snapshot_root      = snapshot_root
)
if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
