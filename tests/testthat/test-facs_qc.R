qc_files <- tibble::tibble(
  file_name = c("a.fcs", "b.fcs"),
  mouse_ID = c("25-7-1", "25-7-2"),
  tissue = "spleen"
)

# -- settings -------------------------------------------------------------------

test_that("qc_settings_ fills defaults for a prep saved without qc settings", {
  prep <- list(data = fake_prep_data(qc_files), info = list())
  s <- qc_settings_(prep)
  expect_null(s$group_col)
  expect_null(s$batch_col)
  expect_null(s$expected_mice)
  expect_equal(s$expected_tissues, "spleen")
  expect_equal(s$qc_populations, c("live_CD45+", "CD4+"))
  expect_equal(s$qc_mfi_populations, character())
  expect_equal(s$qc_metrics, character())
})

test_that("qc_resolve_ errors on an unknown MFI metric", {
  expect_error(qc_resolve_(fake_qc_data(), list(qc_metrics = "median_X")), "median_X")
})

test_that("qc_settings_ returns stored settings when present", {
  stored <- list(group_col = "grp", batch_col = NULL, expected_mice = NULL,
                 expected_tissues = "spleen", qc_populations = "CD4+",
                 qc_mfi_populations = character(), qc_metrics = character())
  prep <- list(data = fake_prep_data(qc_files), info = list(qc = stored))
  expect_equal(qc_settings_(prep), stored)
})

# -- scoring --------------------------------------------------------------------

test_that("qc_score_ flags the planted outlier only", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d))
  expect_equal(unique(r$mouse_ID[r$flag]), "m-3")
  expect_setequal(unique(r$check), c("hierarchy", "counts", "mfi"))
})

test_that("qc_score_ does not score populations under 50 events", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d))
  rare <- r[r$population == "rare", ]
  expect_gt(nrow(rare), 0)
  expect_true(all(is.na(rare$z)))
  expect_false(any(rare$flag))
})

test_that("qc_score_ gives finite z when a box has zero spread", {
  d <- fake_qc_data()
  d$value[d$population == "live" & d$metric == "count"] <- 1e5
  r <- qc_score_(d, fake_qc_settings(d))
  z <- r$z[r$population == "live" & r$metric == "count"]
  expect_true(all(is.finite(z)))
})

test_that("qc_score_ keeps fraction 0 / 1 and count 0 finite", {
  d <- fake_qc_data()
  d$value[d$mouse_ID == "m-1" & d$population == "live" &
            d$metric == "fraction_of_parent"] <- 1
  d$value[d$mouse_ID == "m-2" & d$population == "T" &
            d$metric %in% c("count", "fraction_of_parent")] <- 0
  r <- qc_score_(d, fake_qc_settings(d))
  expect_false(any(is.infinite(r$t) | is.nan(r$t)))
  expect_false(any(is.infinite(r$z) | is.nan(r$z)))
  expect_false(anyNA(r$flag))
})

test_that("qc_score_ skips boxes under 4 samples", {
  d <- fake_qc_data(n_per_group = 3)
  r <- qc_score_(d, fake_qc_settings(d))
  expect_true(all(is.na(r$z)))
})

test_that("qc_score_ leaves excluded files out of the reference but scores them", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d), excluded_files = "m-3_spleen.fcs")
  out <- r[r$mouse_ID == "m-3", ]
  expect_true(all(out$excluded))
  expect_false(any(r$excluded[r$mouse_ID != "m-3"]))
  expect_true(any(out$flag))
  expect_false(any(r$flag[r$mouse_ID != "m-3"]))
})

test_that("qc_score_ excluded files do not shift the box median", {
  d <- fake_qc_data()
  # a second outlier in group A; excluding it must not change m-1's z
  d$value[d$mouse_ID == "m-5" & d$population == "T" &
            d$metric == "fraction_of_parent"] <- 0.9
  s <- fake_qc_settings(d)
  z_of <- function(r) {
    r$z[r$mouse_ID == "m-1" & r$population == "T" & r$metric == "fraction_of_parent"]
  }
  r_ex <- qc_score_(d, s, excluded_files = "m-5_spleen.fcs")
  r_drop <- qc_score_(d[d$mouse_ID != "m-5", ], s)
  expect_equal(z_of(r_ex), z_of(r_drop))
})

test_that("qc_score_ ignores mice with NA group", {
  d <- fake_qc_data()
  d$grp[d$mouse_ID == "m-1"] <- NA
  r <- qc_score_(d, fake_qc_settings(d))
  expect_true(all(is.na(r$z[r$mouse_ID == "m-1"])))
  expect_false(any(r$flag[r$mouse_ID == "m-1"]))
})

test_that("qc_score_ does not score without a group column", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d, group_col = NULL))
  expect_gt(nrow(r), 0)
  expect_true(all(is.na(r$z)))
  expect_false(any(r$flag))
})

test_that("qc_completeness_ lists missing mouse x tissue", {
  d <- fake_qc_data()
  d <- d[d$mouse_ID == "m-1", ]
  res <- qc_completeness_(d, fake_qc_settings(d, expected_mice = c("m-1", "m-99")))
  expect_equal(res$missing$mouse_ID, "m-99")
  expect_equal(res$missing$tissue, "spleen")
  expect_equal(res$n_expected, 2L)
  expect_equal(res$n_present, 1L)
})

test_that("qc_completeness_ counts a missing condition well as missing", {
  d <- dplyr::bind_rows(
    dplyr::mutate(fake_qc_data(), condition = "unstim"),
    dplyr::mutate(fake_qc_data(), condition = "MPO",
                  file_name = sub(".fcs", "_mpo.fcs", file_name, fixed = TRUE))
  )
  d <- d[!(d$mouse_ID == "m-1" & d$condition == "MPO"), ]
  res <- qc_completeness_(d, fake_qc_settings(d, expected_mice = c("m-1", "m-2")))
  expect_equal(nrow(res$missing), 1L)
  expect_equal(res$missing$mouse_ID, "m-1")
  expect_equal(res$missing$condition, "MPO")
  expect_equal(res$n_expected, 4L)
})

test_that("qc_completeness_ returns NULL without expected_mice", {
  d <- fake_qc_data()
  expect_null(qc_completeness_(d, fake_qc_settings(d)))
})

fake_scored <- function(day = as.POSIXct("2026-06-15", tz = "UTC")) {
  tibble::tibble(
    mouse_ID = c(rep("m-1", 6), "m-2", "m-2", "m-3"),
    tissue = c(rep("spleen", 7), "lung", "spleen"),
    condition = "none",
    grp = "A",
    day = day,
    population = c("a", "b", "c", "a", "b", "c", "a", "a", "a"),
    metric = c(rep("fraction_of_parent", 3), rep("count", 3),
               "fraction_of_parent", "fraction_of_parent", "count"),
    check = c(rep("hierarchy", 3), rep("counts", 3), "hierarchy", "hierarchy",
              "counts"),
    value = 1,
    z = c(4, -5, 6, 4, 4, 4, 9, -4, 1),
    flag = c(rep(TRUE, 8), FALSE),
    excluded = c(rep(FALSE, 6), TRUE, TRUE, FALSE)
  )
}

test_that("qc_flag_table_ counts a cascading gate once per sample and check", {
  s <- list(group_col = "grp", batch_col = "day")
  tbl <- qc_flag_table_(fake_scored(), s)
  expect_equal(nrow(tbl), 8L)
  expect_equal(unique(tbl$score[tbl$mouse_ID == "m-1"]), 2L)
  expect_equal(unique(tbl$score[tbl$mouse_ID == "m-2"]), 2L)
  expect_false("m-3" %in% tbl$mouse_ID)
  expect_equal(names(tbl), c("mouse_ID", "score", "tissue", "condition", "group",
                             "batch", "population", "metric", "value", "z",
                             "excluded"))
  # m-1 before m-2 (same score, mouse_ID order); within a mouse by |z| desc
  expect_equal(tbl$mouse_ID[1], "m-1")
  expect_equal(tbl$z[1], 6)
  expect_true(all(tbl$excluded[tbl$mouse_ID == "m-2"]))
})

test_that("qc_flag_table_ formats a POSIXct batch as a date", {
  tbl <- qc_flag_table_(fake_scored(), list(group_col = "grp", batch_col = "day"))
  expect_equal(unique(tbl$batch), "2026-06-15")
})

test_that("qc_flag_table_ returns zero rows when nothing is flagged", {
  sc <- fake_scored()
  sc$flag <- FALSE
  tbl <- qc_flag_table_(sc, list(group_col = "grp", batch_col = NULL))
  expect_equal(nrow(tbl), 0L)
})

test_that("qc_score_ only scores MFI for qc_mfi_populations and positive values", {
  d <- fake_qc_data()
  r0 <- qc_score_(d, fake_qc_settings(d, qc_mfi_populations = character()))
  expect_false("mfi" %in% r0$check)
  d$value[d$mouse_ID == "m-1" & d$metric == "median_CD25"] <- -5
  r <- qc_score_(d, fake_qc_settings(d))
  neg <- r[r$mouse_ID == "m-1" & r$metric == "median_CD25", ]
  expect_true(is.na(neg$t))
  expect_true(is.na(neg$z))
})
