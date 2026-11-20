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

final <- FALSE

# ---- experiment ---------------------------------------------------------------

exp_path  <- file.path("data", "26-1_BMT-mRNA-V4")
meta_path <- file.path(exp_path, "meta", "meta_26-1.xlsx")
meta      <- meta_read(meta_path)  # standard template
mice      <- dplyr::filter(meta$meta, !is.na(mouse_ID))
db_path   <- "~/databases/MPO_mRNA/MPO_mRNA.duckdb"

# ---- assay --------------------------------------------------------------------

path      <- file.path(exp_path, "raw", "histology", "histo_26-1.xlsx")  # curated workbook
assay     <- "PAS"                              # assay name in the database
tissue    <- "kidney"                           # tissue of every sample in the file

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

# ---- run ------------------------------------------------------------------------

con <- if (final) {
  db_connect(path.expand(db_path)) |>
    db_write_meta(meta, "26-1", "BMT-mRNA-V4", project = "MPO_mRNA")
}
histo_prep_run(
  final,
  path       = path,
  assay      = assay,
  tissue     = tissue,
  prepare    = prepare,
  mice       = mice,
  exclusions = meta$exclusions,                 # NULL if the sheet doesn't exist
  con        = con,
  meta_path  = meta_path,                       # final = TRUE also records the exclusions sheet
  # QC report: outliers scored within groups (batch column shown in the
  # overview), completeness against all mice minus dropouts (`early`).
  group_col     = "group",
  batch_col     = "death_date",
  expected_mice = mice$mouse_ID[is.na(mice$early)],
  rules         = rules,
  snapshot_root = file.path("data", "z_snapshots")
)
if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
