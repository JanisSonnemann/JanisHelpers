# ELISA prep script for one curated file -- copy into the project and adapt.
# Copy a clean version of this template and a blank curated sheet to the current working directory (then rename them):
# file.copy(system.file(c("scripts/elisa_prep_template.R", "templates/elisa_curated_template.xlsx"), package = "JanisHelpers"), ".")
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

elisa_file   <- "spleen-restim-48h_26-1.xlsx"   # curated workbook in <exp_dir>/raw/ELISA/
assay        <- "restim"                        # assay name in the database
tissue       <- "spleen"                        # tissue of every sample in the file

# Wide sheet (one row per mouse, one column per date, `final` = end-point)?
# Set wide <- TRUE: it is reshaped to <stem>_long.xlsx first (e.g.
# elisa_file <- "anti-mpo-titer_26-1_wide.xlsx", wide_analyte <- "anti-MPO-titer").
wide         <- FALSE
wide_analyte <- NULL                            # analyte name if the sheet has no `analyte` column

# Optional changes to the long table before saving (e.g. drop an analyte).
prepare <- function(data) {
  data
}

# ---- QC report (defaults are usually fine) -------------------------------------

group_col     <- "group"                        # outliers scored within these groups
batch_col     <- "death_date"                   # batch plots per day
expected_mice <- NULL                           # NULL = all mice in the meta sheet minus dropouts (`early`)
max_cv        <- 20                             # replicate CV limit in percent
log_scale     <- FALSE                          # TRUE = score / plot on a log axis

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

path <- file.path(exp_path, "raw", "ELISA", elisa_file)
if (wide) path <- elisa_reshape_wide(path, analyte = wide_analyte)

con <- if (final) {
  db_connect(path.expand(db_path)) |>
    db_write_meta(meta, exp_code, exp_name, project = project)
}
elisa_prep_run(
  final,
  path          = path,
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
  max_cv        = max_cv,
  log_scale     = log_scale,
  snapshot_root = snapshot_root
)
if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
