# Builds the blank curated-sheet templates shipped in inst/templates/.
# Run from the package root:  Rscript data-raw/make_curated_templates.R
# Needs writexl (Suggests only, not Imports).

readme <- function(lines) data.frame(instructions = lines)

elisa_values <- data.frame(
  mouse_ID = character(), tissue = character(), date = as.Date(character()),
  TNFa = character(), `IL-2` = character(), check.names = FALSE
)
elisa_range <- data.frame(
  label = c("LLOQ", "ULOQ", "unit"),
  TNFa = c(NA, NA, "pg/ml"), `IL-2` = c(NA, NA, "pg/ml"), check.names = FALSE
)
elisa_readme <- readme(c(
  "Curated ELISA sheet -- read by JanisHelpers::elisa_read_curated(); one file per measurement set.",
  "Name the file <root>.xlsx (e.g. spleen-restim-48h_26-1.xlsx). The root identifies the import: re-importing the same file replaces it in the database.",
  "",
  "Sheet 'values' -- one row per well:",
  "  mouse_ID  required; as in the meta workbook (e.g. 26-1-3).",
  "  tissue    optional column and cell; empty = the prep script's tissue argument (e.g. spleen, serum).",
  "  date      optional column and cell; collection date. Empty = the mouse's end (death) date.",
  "  one column per analyte -- rename/add/remove TNFa, IL-2 as needed. Headers must be unique and non-empty.",
  "  Values: numbers, OOR< or OOR>. Any other text (n.d., <0.5, ...) stops the import.",
  "  Repeated rows for the same mouse x tissue x date are technical replicates (numbered in row order).",
  "",
  "Sheet 'range' (optional) -- column 1 holds the labels LLOQ, ULOQ, unit; one column per analyte, same headers as in 'values'.",
  "  Values outside LLOQ/ULOQ are stored as they are and counted in the QC report.",
  "",
  "This README sheet is ignored by the reader. Do not add other columns to 'values'."
))

histo_values <- data.frame(
  mouse_ID = character(), ID = character(), tissue = character(),
  date = as.Date(character()), gloms = numeric(), crescents = numeric(),
  necrosis = numeric()
)
histo_readme <- readme(c(
  "Curated histology sheet -- read by JanisHelpers::histo_read_curated(); one file per scoring set.",
  "Name the file <root>.xlsx (e.g. histo_26-1.xlsx). The root identifies the import: re-importing the same file replaces it in the database.",
  "",
  "Sheet 'values' -- one row per mouse (x tissue x date):",
  "  mouse_ID  required; as in the meta workbook (e.g. 26-1-3).",
  "  ID        optional slide label; ignored (as is any column ending in _ID).",
  "  tissue    optional column and cell; empty = the prep script's tissue argument (e.g. kidney).",
  "  date      optional column and cell; collection date. Empty = the mouse's end (death) date.",
  "  one numeric column per score -- rename/add/remove gloms, crescents, necrosis as needed. Headers must be unique and non-empty.",
  "  Leave a cell empty if not scored; a mouse with all cells empty is kept (shows as missing in the QC report).",
  "  Any text in a score cell stops the import. Other text columns (notes, ...) are not allowed.",
  "  Derived values (e.g. % crescents) are calculated in the prep script (histo_calc_ratio()), not in the sheet.",
  "",
  "This README sheet is ignored by the reader."
))

writexl::write_xlsx(
  list(values = elisa_values, range = elisa_range, README = elisa_readme),
  "inst/templates/elisa_curated_template.xlsx"
)
writexl::write_xlsx(
  list(values = histo_values, README = histo_readme),
  "inst/templates/histo_curated_template.xlsx"
)
