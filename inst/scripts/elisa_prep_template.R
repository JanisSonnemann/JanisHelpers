# ELISA prep script for one curated file -- copy into the project and adapt.
# Copy this script and a blank curated sheet into the project (run once from the project root, adapt the targets):
# file.copy(system.file(c("scripts/elisa_prep_template.R", "templates/elisa_curated_template.xlsx"), package = "JanisHelpers"), c("data/26-1_BMT-mRNA-V4/src/elisa_prep_26-1_spleen-restim-48h.R", "data/26-1_BMT-mRNA-V4/raw/ELISA/spleen-restim-48h_26-1.xlsx"))
#
# One self-contained script per curated .xlsx; the pipeline logic lives in
# JanisHelpers (elisa_prep_run(), db_write_meta()).
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
meta      <- meta_read(file.path(exp_path, "meta", "meta_26-1.xlsx"))  # standard template
mice      <- dplyr::filter(meta$meta, !is.na(mouse_ID))
db_path   <- "~/databases/MPO_mRNA/MPO_mRNA.duckdb"

# ---- assay --------------------------------------------------------------------

path      <- file.path(exp_path, "raw", "ELISA", "spleen-restim-48h_26-1.xlsx")  # curated workbook
assay     <- "restim"                           # assay name in the database
tissue    <- "spleen"                           # tissue of every sample in the file

# Optional changes to the long table before saving (e.g. drop an analyte).
prepare <- function(data) {
  data
}

# ---- run ------------------------------------------------------------------------

con <- if (final) {
  db_connect(path.expand(db_path)) |>
    db_write_meta(meta, "26-1", "BMT-mRNA-V4", project = "MPO_mRNA")
}
elisa_prep_run(
  final,
  path       = path,
  assay      = assay,
  tissue     = tissue,
  prepare    = prepare,
  mice       = mice,
  exclusions = meta$exclusions,                 # NULL if the sheet doesn't exist
  con        = con,
  # QC report: outliers scored within groups, batch plots per day,
  # completeness against all mice minus dropouts (`early`).
  group_col     = "group",
  batch_col     = "death_date",
  expected_mice = mice$mouse_ID[is.na(mice$early)],
  max_cv        = 20,                           # replicate CV limit in percent
  log_scale     = FALSE,                        # TRUE = score / plot on a log axis
  snapshot_root = file.path("data", "z_snapshots")
)
if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
