# Generator for elisa_curated.xlsx and histo_curated.xlsx.
# writexl is in DESCRIPTION Suggests only (not Imports); install it to regenerate.
# Run from the package root: Rscript tests/fixtures/make_curated_fixtures.R
elisa_values <- data.frame(
  mouse_ID = c("M1", "M1", "M2", "M2", "M1"),
  date     = as.Date(c(NA, NA, NA, NA, "2026-05-06")),
  TNFa     = c("4.2", "3.9", "6.0", "6.1", "5.0"),
  `IL-2`   = c("OOR<", "10.1", "20.0", "20.5", "12.0"),
  check.names = FALSE
)
# numeric cells stay numeric except OOR flags
elisa_values <- lapply(elisa_values, function(col) col)
elisa_values$TNFa <- as.numeric(elisa_values$TNFa)
il2 <- elisa_values[["IL-2"]]
elisa_values[["IL-2"]] <- NULL
elisa_values <- as.data.frame(elisa_values, check.names = FALSE)
elisa_values[["IL-2"]] <- il2   # character column mixing numbers and "OOR<"

elisa_range <- data.frame(
  param = c("LLOQ", "ULOQ", "unit"),
  TNFa  = c("7.81", "31.3", "pg/ml"),
  `IL-2` = c("15.6", "250", "pg/ml"),
  check.names = FALSE
)

histo_values <- data.frame(
  mouse_ID  = c("M1", "M2", "M3", "M4"),
  ID        = c("s1", "s2", "s3", "s4"),
  gloms     = c(85, 68, NA, 70),
  crescents = c(22, 3, NA, 80),
  necrosis  = c(15, 1, NA, 2)
)

writexl::write_xlsx(list(values = elisa_values, range = elisa_range),
                    "tests/fixtures/elisa_curated.xlsx")
writexl::write_xlsx(list(values = histo_values),
                    "tests/fixtures/histo_curated.xlsx")
