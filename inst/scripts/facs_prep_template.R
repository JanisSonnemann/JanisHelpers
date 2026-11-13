# FACS prep script for one panel -- copy into the project and adapt.
# Template: system.file("scripts", "facs_prep_template.R", package = "JanisHelpers")
#
# The loop:
#   final <- FALSE  read the .wsp, select populations, calculate, save the
#                   prep result and render the QC report (all samples shown;
#                   exclusions only flagged). Adjust gates in FlowJo (save a
#                   new dated .wsp for major changes) and/or edit the meta
#                   workbook's `exclusions` sheet, then rerun.
#   final <- TRUE   QC passed: export the saved prep result to the database,
#                   dropping and recording the exclusions. Asks before
#                   overwriting a root that is already stored.
#
# Run interactively from the project root (the folder holding data/).

library(JanisHelpers)

# ---- config -----------------------------------------------------------------

final <- FALSE

data_dir  <- "data"
exp_dir   <- "26-1_BMT-mRNA-V4"                 # inside data_dir
wsp       <- "overview-26-1_20260806.wsp"       # in raw/FACS/, <root>_YYYYMMDD.wsp
panel_dir <- "overview"                         # panel folder in raw/FACS/ (holds the .fcs)
assay     <- "overview"                         # panel's assay name in the database
condition_keyword <- NA_character_              # e.g. "restimulation" for multi-well panels
meta_path <- file.path(data_dir, exp_dir, "meta", "meta_26-1.xlsx")
db_path   <- "~/databases/MPO_mRNA/MPO_mRNA.duckdb"

# QC report settings (columns joined from the meta workbook's `meta` sheet).
group_col  <- "group"                           # treatment group: outliers scored within groups
batch_col  <- "death_date"                      # staining / acquisition day for the batch plots
qc_populations     <- NULL                      # NULL = all selected populations
qc_mfi_populations <- character()               # e.g. "Treg_FoxP3+_CD25hi"; empty = no MFI check
qc_metrics         <- NULL                      # NULL = all median_* MFIs

facs_dir   <- file.path(data_dir, exp_dir, "raw", "FACS")
wsp_path   <- file.path(facs_dir, wsp)
panel_path <- file.path(facs_dir, panel_dir)
meta       <- meta_read(meta_path)
exclusions <- meta$exclusions                   # NULL if the sheet doesn't exist

if (!final) {
  # ---- prep -------------------------------------------------------------------
  data <- facs_read_wsp(
    wsp_path,
    keywords = c("mouse_ID", "tissue", if (!is.na(condition_keyword)) condition_keyword)
  )$data

  # Populations to keep, in display order.
  selection <- c("live_CD45+")

  data <- data |>
    dplyr::filter(population %in% selection) |>
    dplyr::mutate(population = factor(population, levels = selection))

  # Calculations for this panel (adapt).
  data <- data |>
    facs_calc_pct_of(ref_pop = "live_CD45+")

  # Group and batch columns for the QC report.
  mice <- meta$meta |> dplyr::filter(!is.na(mouse_ID))
  data <- meta_annotate(
    data, dplyr::select(mice, mouse_ID, dplyr::all_of(c(group_col, batch_col))),
    by = "mouse_ID"
  )
  # Mice this panel should contain: all mice minus dropouts. The dropout
  # column (`early` here) is experiment-specific -- adapt.
  expected_mice <- mice$mouse_ID[is.na(mice$early)]

  prep_path <- facs_prep_save(
    data, wsp_path, assay, panel_path,
    condition_keyword = condition_keyword,
    group_col = group_col, batch_col = batch_col,
    expected_mice = expected_mice,
    qc_populations = qc_populations,
    qc_mfi_populations = qc_mfi_populations,
    qc_metrics = qc_metrics
  )
  report_knit_facs_qc(prep_path, exclusions)
} else {
  # ---- final export -------------------------------------------------------------
  prep <- facs_prep_load(wsp_path, panel_path)
  con <- db_connect(path.expand(db_path))
  # Experiment and subjects must be registered first (db_write_experiment(),
  # db_write_subjects()).
  db_export_facs(con, prep, exclusions, data_dir = data_dir)
  # Snapshot: DBI::dbExecute(con, "EXPORT DATABASE '<dir>' (FORMAT PARQUET)")
  DBI::dbDisconnect(con, shutdown = TRUE)
}
