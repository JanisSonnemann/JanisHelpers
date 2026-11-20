# histology prep script for one curated file -- copy into the project and adapt.
# Copy a clean version of this template and a blank curated sheet to the current working directory (then rename them):
# file.copy(system.file(c("scripts/histo_prep_template.R", "templates/histo_curated_template.xlsx"), package = "JanisHelpers"), ".")
#
# One self-contained script per curated .xlsx; the pipeline logic lives in
# JanisHelpers (histo_prep_run(), db_write_meta()).
#
#   final <- FALSE  read the curated sheet, save the prep result and render
#                   the QC report (all samples shown; exclusions only
#                   flagged). Fix the .xlsx and/or the meta workbook's
#                   `exclusions` sheet, then rerun.
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

# ---- assay -------------------------------------------------------------------

histo_file <- "histo_26-1.xlsx"                 # curated workbook in <exp_dir>/raw/histology/
assay      <- "PAS"                             # assay name in the database
tissue     <- "kidney"                          # tissue of every sample in the file

# Changes to the long table before saving: derived ratio metrics, dropping a metric.
# New score columns in the sheet are picked up as metrics automatically; a
# `comment` column is free text shown in the QC report (not exported).
# Ratios: one named entry per new metric (new name = numerator), all over
# `denominator` (or give one denominator per ratio). Every metric named here
# must exist in the sheet.
prepare <- function(data) {
  histo_calc_ratio(
    data,
    numerator = c(pct_cresc    = "crescents",
                  pct_necrosis = "necrosis",
                  pct_abnormal = "abnormal_gloms"),
    denominator = "gloms"
  )
}

# Plausibility rules for the QC report (one TRUE / FALSE per sample, metrics
# as variable names; FALSE = violation, NA = not checked).
rules <- c("crescents <= gloms", "necrosis <= gloms", "abnormal_gloms <= gloms")

# ---- QC report (defaults are usually fine) -------------------------------------

group_col     <- "group"                        # outliers scored within these groups
batch_col     <- "death_date"                   # batch column shown in the overview
expected_mice <- NULL                           # NULL = all mice in the meta sheet minus dropouts (`early`)

# ==============================================================================
# PIPELINE -- do not change anything below this line
# ==============================================================================

exp_path  <- file.path("data", exp_dir)
meta_path <- file.path(exp_path, "meta", meta_file)
meta      <- meta_read(meta_path)
mice      <- dplyr::filter(meta$meta, !is.na(mouse_ID))
if (is.null(expected_mice)) {
  dropped <- if ("early" %in% names(mice)) !is.na(mice$early) else FALSE
  expected_mice <- mice$mouse_ID[!dropped]
}

con <- if (final) {
  db_connect(path.expand(db_path)) |>
    db_write_meta(meta, exp_code, exp_name, project = project)
}
histo_prep_run(
  final,
  path          = file.path(exp_path, "raw", "histology", histo_file),
  assay         = assay,
  tissue        = tissue,
  prepare       = prepare,
  mice          = mice,
  exclusions    = meta$exclusions,              # NULL if the sheet doesn't exist
  con           = con,
  meta_path     = meta_path,                    # final = TRUE also records the exclusions sheet
  group_col     = group_col,
  batch_col     = batch_col,
  expected_mice = expected_mice,
  rules         = rules,
  snapshot_root = snapshot_root
)
if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
