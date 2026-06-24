# ── Fixtures ──────────────────────────────────────────────────────────────────

# Two-group data, no tissue column
df_two <- data.frame(
  group   = rep(c("ctrl", "trt"), each = 10),
  CD4_pct = c(rep(30, 10), rep(45, 10)),
  CD8_pct = c(rep(20, 10), rep(28, 10))
)

# Two-group data with tissue column (5 obs per group per tissue)
df_tissue <- data.frame(
  group   = rep(c("ctrl", "trt"), each = 10),
  tissue  = rep(c("spleen", "blood"), times = 10),
  CD4_pct = c(rep(30, 10), rep(45, 10))
)

# Three-group data: maximally distinct groups → KW always significant
df_sig <- data.frame(
  group   = rep(c("ctrl", "trt_A", "trt_B"), each = 10),
  CD4_pct = c(rep(10, 10), rep(50, 10), rep(90, 10))
)

# Three-group data: identical values across groups → KW never significant (p = 1)
df_ns <- data.frame(
  group   = rep(c("ctrl", "trt_A", "trt_B"), each = 5),
  CD4_pct = c(30, 31, 29, 30, 31,
              30, 29, 31, 30, 29,
              31, 30, 30, 29, 31)
)

# ── analysis_summary_table ────────────────────────────────────────────────────

test_that("analysis_summary_table returns a gt_tbl when tissue_col = NULL", {
  result <- analysis_summary_table(df_two, group_col = "group", tissue_col = NULL)
  expect_s3_class(result, "gt_tbl")
})

test_that("analysis_summary_table returns named list of gt_tbl when tissue_col is set", {
  result <- analysis_summary_table(df_tissue, group_col = "group", tissue_col = "tissue")
  expect_type(result, "list")
  expect_named(result, c("spleen", "blood"), ignore.order = TRUE)
  expect_s3_class(result[["spleen"]], "gt_tbl")
  expect_s3_class(result[["blood"]], "gt_tbl")
})

# ── analysis_posthoc_tables ───────────────────────────────────────────────────

test_that("analysis_posthoc_tables returns named list keyed by variable name", {
  result <- suppressWarnings(
    analysis_posthoc_tables(df_sig, group_col = "group", tissue_col = NULL)
  )
  expect_type(result, "list")
  expect_named(result, "CD4_pct")
})

test_that("analysis_posthoc_tables returns gt_tbl when KW is significant", {
  result <- suppressWarnings(
    analysis_posthoc_tables(df_sig, group_col = "group", tissue_col = NULL)
  )
  expect_s3_class(result[["CD4_pct"]], "gt_tbl")
})

test_that("analysis_posthoc_tables returns character message when KW is not significant", {
  result <- suppressWarnings(
    analysis_posthoc_tables(df_ns, group_col = "group", tissue_col = NULL)
  )
  expect_type(result[["CD4_pct"]], "character")
  expect_match(result[["CD4_pct"]], "No post-hoc testing performed")
})

test_that("analysis_posthoc_tables returns nested list when tissue_col is set", {
  df <- data.frame(
    group   = rep(c("ctrl", "trt_A", "trt_B"), each = 10),
    tissue  = rep(c("spleen", "blood"), times = 15),
    CD4_pct = c(rep(10, 10), rep(50, 10), rep(90, 10))
  )
  result <- suppressWarnings(
    analysis_posthoc_tables(df, group_col = "group", tissue_col = "tissue")
  )
  expect_type(result, "list")
  expect_named(result, c("spleen", "blood"), ignore.order = TRUE)
  expect_type(result[["spleen"]], "list")
  expect_named(result[["spleen"]], "CD4_pct")
})
