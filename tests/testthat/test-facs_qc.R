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
  expect_equal(s$min_events, 50)
})

test_that("qc_resolve_ validates min_events", {
  d <- fake_qc_data()
  expect_error(qc_resolve_(d, list(min_events = -5)), "min_events")
  expect_error(qc_resolve_(d, list(min_events = "50")), "min_events")
  expect_equal(qc_resolve_(d, list(min_events = 0))$min_events, 0)
})

test_that("qc_resolve_ errors on an unknown MFI metric", {
  expect_error(qc_resolve_(fake_qc_data(), list(qc_metrics = "median_X")), "median_X")
})

test_that("qc_settings_ returns stored settings when present", {
  stored <- list(group_col = "grp", batch_col = NULL, expected_mice = NULL,
                 expected_tissues = "spleen", qc_populations = "CD4+",
                 qc_mfi_populations = character(), qc_metrics = character())
  prep <- list(data = fake_prep_data(qc_files), info = list(qc = stored))
  # older stored settings get the min_events default
  expect_equal(qc_settings_(prep), c(stored, list(min_events = 50)))
  stored$min_events <- 20
  prep$info$qc <- stored
  expect_equal(qc_settings_(prep), stored)
})

# -- scoring --------------------------------------------------------------------

# Outliers by ggplot2's own boxplot statistics, for comparison.
ggplot_outliers_ <- function(y) {
  st <- ggplot2::StatBoxplot$compute_group(data.frame(x = 1, y = y),
                                           scales = list(x = NULL))
  st$outliers[[1]]
}

test_that("qc_score_ flags the planted outlier only", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d))
  expect_equal(unique(r$mouse_ID[r$flag]), "m-3")
  expect_setequal(unique(r$check), c("hierarchy", "counts", "mfi"))
  t_frac <- r[r$mouse_ID == "m-3" & r$population == "T" &
                r$metric == "fraction_of_parent", ]
  expect_true(t_frac$flag)
  expect_equal(t_frac$side, "low")
})

test_that("qc_score_ runs without warnings", {
  d <- fake_qc_data()
  expect_no_warning(qc_score_(d, fake_qc_settings(d)))
  expect_no_warning(qc_score_(d, fake_qc_settings(d), excluded_files = "m-3_spleen.fcs"))
})

test_that("qc_score_ flags exactly the points ggplot2 draws outside the whiskers", {
  d <- fake_qc_data()
  # spread group A's T fractions so the box has a real IQR, plus a high value
  a_t <- d$grp == "A" & d$population == "T" & d$metric == "fraction_of_parent"
  d$value[a_t] <- c(0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.60)
  r <- qc_score_(d, fake_qc_settings(d))
  box <- r[r$grp == "A" & r$population == "T" & r$metric == "fraction_of_parent", ]
  expect_equal(sort(box$value[box$flag]), sort(ggplot_outliers_(box$value)))
  expect_equal(box$side[box$value == 0.60], "high")
  expect_equal(box$box_median, rep(stats::median(box$value), nrow(box)))
  # counts on the log axis: same rule on log10(count + 1)
  cnt <- r[r$grp == "A" & r$population == "T" & r$metric == "count", ]
  expect_equal(sort(log10(cnt$value[cnt$flag] + 1)),
               sort(ggplot_outliers_(log10(cnt$value + 1))))
})

test_that("qc_score_ uses the fraction itself, not a transformed scale", {
  d <- fake_qc_data()
  # 26-1-13-like box: seven mice at 3-8 %, one at 18 %
  a_t <- d$grp == "A" & d$population == "T" & d$metric == "fraction_of_parent"
  d$value[a_t] <- c(0.030, 0.036, 0.038, 0.039, 0.061, 0.077, 0.082, 0.181)
  r <- qc_score_(d, fake_qc_settings(d))
  box <- r[r$grp == "A" & r$population == "T" & r$metric == "fraction_of_parent", ]
  expect_equal(box$mouse_ID[box$flag], box$mouse_ID[box$value == 0.181])
})

test_that("qc_score_ does not flag populations under min_events", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d))
  rare <- r[r$population == "rare", ]
  expect_gt(nrow(rare), 0)
  expect_false(any(rare$flag))
  # with min_events = 0 the planted outlier's rare fraction is labelled
  r0 <- qc_score_(d, fake_qc_settings(d, min_events = 0))
  rare0 <- r0[r0$population == "rare" & r0$metric == "fraction_of_parent", ]
  expect_true(rare0$flag[rare0$mouse_ID == "m-3"])
})

test_that("qc_score_ handles a box with zero spread", {
  d <- fake_qc_data()
  d$value[d$population == "live" & d$metric == "count"] <- 1e5
  r <- qc_score_(d, fake_qc_settings(d))
  live <- r[r$population == "live" & r$metric == "count", ]
  expect_false(any(live$flag))
  expect_false(anyNA(r$flag))
})

test_that("qc_score_ keeps fraction 0 / 1 and count 0 finite", {
  d <- fake_qc_data()
  d$value[d$mouse_ID == "m-1" & d$population == "live" &
            d$metric == "fraction_of_parent"] <- 1
  d$value[d$mouse_ID == "m-2" & d$population == "T" &
            d$metric %in% c("count", "fraction_of_parent")] <- 0
  r <- qc_score_(d, fake_qc_settings(d))
  expect_false(any(is.infinite(r$t) | is.nan(r$t)))
  expect_false(anyNA(r$flag))
})

test_that("qc_score_ draws the box from all plotted points, excluded ones too", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d), excluded_files = "m-3_spleen.fcs")
  out <- r[r$mouse_ID == "m-3", ]
  expect_true(all(out$excluded))
  expect_false(any(r$excluded[r$mouse_ID != "m-3"]))
  expect_true(any(out$flag))
  r_all <- qc_score_(d, fake_qc_settings(d))
  expect_equal(r$flag, r_all$flag)
})

test_that("qc_score_ does not flag without a group column", {
  d <- fake_qc_data()
  r <- qc_score_(d, fake_qc_settings(d, group_col = NULL))
  expect_gt(nrow(r), 0)
  expect_false(any(r$flag))
  expect_true(all(is.na(r$side)))
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
    population_full_path = paste0("root/", population),
    metric = c(rep("fraction_of_parent", 3), rep("count", 3),
               "fraction_of_parent", "fraction_of_parent", "count"),
    check = c(rep("hierarchy", 3), rep("counts", 3), "hierarchy", "hierarchy",
              "counts"),
    value = 1,
    box_median = 2,
    side = c("low", "high", "low", "low", "high", "low", "high", "low", NA),
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
                             "batch", "population", "metric", "value",
                             "box_median", "side",
                             "excluded"))
  # m-1 before m-2 (same score, mouse_ID order)
  expect_equal(tbl$mouse_ID[1], "m-1")
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

# -- plots ----------------------------------------------------------------------

label_layer_ <- function(p) {
  Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")), p$layers)[[1]]
}

test_that("qc_plot_individuals_ labels only flagged mice", {
  d <- fake_qc_data()
  s <- fake_qc_settings(d)
  p <- qc_plot_individuals_(qc_score_(d, s), s, "hierarchy")
  expect_s3_class(p, "ggplot")
  labs <- ggplot2::ggplot_build(p)$data[[which(vapply(
    p$layers, function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")), logical(1)
  ))[1]]]$label
  expect_equal(unique(labs[labs != ""]), "m-3")
  expect_s3_class(p$facet, "FacetGrid")
})

test_that("qc_plot_individuals_ shows one population's markers for MFI", {
  d <- fake_qc_data()
  s <- fake_qc_settings(d)
  p <- qc_plot_individuals_(qc_score_(d, s), s, "mfi", population = "T")
  expect_s3_class(p, "ggplot")
  expect_equal(unique(as.character(p$data$population)), "T")
})

test_that("qc_plot_batch_ facets by group and uses batch on x", {
  d <- fake_qc_data()
  s <- fake_qc_settings(d)
  p <- qc_plot_batch_(qc_score_(d, s), s, "counts", tissue = "spleen")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetGrid")
  x <- p$data$.batch
  expect_true(is.factor(x))
  expect_equal(nlevels(x), 2L)
  expect_match(levels(x), "^\\d{4}-\\d{2}-\\d{2}$")
  expect_length(Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")),
                       p$layers), 0)
})

test_that("qc plots return NULL without data", {
  d <- fake_qc_data()
  s0 <- fake_qc_settings(d, qc_mfi_populations = character())
  expect_null(qc_plot_individuals_(qc_score_(d, s0), s0, "mfi"))
  s1 <- fake_qc_settings(d, batch_col = NULL)
  expect_null(qc_plot_batch_(qc_score_(d, s1), s1, "counts", tissue = "spleen"))
  s <- fake_qc_settings(d)
  expect_null(qc_plot_batch_(qc_score_(d, s), s, "counts", tissue = "lung"))
})

test_that("qc plots work without ggrepel", {
  local_mocked_bindings(qc_has_ggrepel_ = function() FALSE)
  d <- fake_qc_data()
  s <- fake_qc_settings(d)
  p <- qc_plot_individuals_(qc_score_(d, s), s, "hierarchy")
  expect_true(inherits(label_layer_(p)$geom, "GeomText"))
  expect_no_error(ggplot2::ggplot_build(p))
})

dup_leaf_data_ <- function() {
  d <- fake_qc_data()
  # a second gate also named "T" under another parent, with 3x the events
  other <- d[d$population_full_path == "root/live/T", ]
  other$population_full_path <- "root/other/T"
  other$value[other$metric == "count"] <- 3 * other$value[other$metric == "count"]
  dplyr::bind_rows(d, other)
}

test_that("qc plots keep gates with the same leaf name apart", {
  d <- dup_leaf_data_()
  s <- fake_qc_settings(d)
  p <- qc_plot_individuals_(qc_score_(d, s), s, "counts")
  rows <- unique(as.character(p$data$.row))
  expect_true(all(c("live/T", "other/T") %in% rows))
  expect_true("live" %in% rows)
  b <- qc_plot_batch_(qc_score_(d, s), s, "counts", tissue = "spleen")
  expect_true(all(c("live/T", "other/T") %in% unique(as.character(b$data$.row))))
})

test_that("qc_flag_table_ names the full gate when leaf names repeat", {
  d <- dup_leaf_data_()
  s <- fake_qc_settings(d)
  tbl <- qc_flag_table_(qc_score_(d, s), s)
  expect_true(all(tbl$population[tbl$mouse_ID == "m-3"] %in%
                    c("live/T", "other/T", "live", "rare")))
  expect_true("live/T" %in% tbl$population)
})

test_that("qc_fig_height_ grows with facet rows", {
  expect_equal(qc_fig_height_(1), 3)
  expect_equal(qc_fig_height_(10), 17)
})

test_that("qc_score_ only scores MFI for qc_mfi_populations and positive values", {
  d <- fake_qc_data()
  r0 <- qc_score_(d, fake_qc_settings(d, qc_mfi_populations = character()))
  expect_false("mfi" %in% r0$check)
  d$value[d$mouse_ID == "m-1" & d$metric == "median_CD25"] <- -5
  r <- qc_score_(d, fake_qc_settings(d))
  neg <- r[r$mouse_ID == "m-1" & r$metric == "median_CD25", ]
  expect_true(is.na(neg$t))
  expect_false(neg$flag)
})

# -- empty values / missing stats / keywords ----------------------------------------

# fake_qc_data() with two liver mice and a count_per_g metric that is empty
# for spleen (no organ weights) and for one liver file.
empty_qc_data_ <- function() {
  d <- fake_qc_data(n_per_group = 2)
  d$tissue[d$mouse_ID %in% c("m-1", "m-3")] <- "liver"
  per_g <- dplyr::filter(d, metric == "count") |>
    dplyr::mutate(metric = "count_per_g",
                  value = dplyr::if_else(tissue == "spleen" | mouse_ID == "m-3",
                                         NA_real_, value))
  dplyr::bind_rows(d, per_g)
}

test_that("qc_empty_values_ classifies each metric x tissue", {
  e <- qc_empty_values_(empty_qc_data_())
  st <- function(m, t) e$status[e$metric == m & e$tissue == t]
  expect_equal(st("count_per_g", "spleen"), "all empty")
  expect_equal(st("count_per_g", "liver"), "partly empty")
  expect_equal(st("count", "spleen"), "ok")
  liver <- e[e$metric == "count_per_g" & e$tissue == "liver", ]
  expect_equal(c(liver$n_values, liver$n_empty, liver$n_files_empty), c(6L, 3L, 1L))
})

test_that("qc_empty_values_ marks a metric missing for a tissue as absent", {
  d <- dplyr::filter(empty_qc_data_(), !(metric == "count_per_g" & tissue == "spleen"))
  e <- qc_empty_values_(d)
  absent <- e[e$metric == "count_per_g" & e$tissue == "spleen", ]
  expect_equal(absent$status, "absent")
  expect_equal(absent$n_values, 0L)
})

test_that("qc_empty_grid_ has one row per metric and one column per tissue", {
  g <- qc_empty_grid_(qc_empty_values_(empty_qc_data_()))
  expect_setequal(names(g), c("metric", "liver", "spleen"))
  expect_equal(g$spleen[g$metric == "count_per_g"], "all empty")
  expect_equal(g$liver[g$metric == "count_per_g"], "3/6 empty")
})

test_that("qc_empty_values_ treats empties at 0 events as expected", {
  d <- fake_qc_data(n_per_group = 2)
  rare_count <- d$mouse_ID == "m-2" & d$population == "rare" & d$metric == "count"
  d$value[rare_count] <- 0
  d <- dplyr::bind_rows(
    d,
    dplyr::filter(d, population == "rare", metric == "count") |>
      dplyr::mutate(metric = "median_CD25",
                    value = dplyr::if_else(mouse_ID == "m-2", NA_real_, 500))
  )
  e <- qc_empty_values_(d)
  mfi <- e[e$metric == "median_CD25", ]
  expect_equal(mfi$status, "ok")
  expect_equal(c(mfi$n_empty, mfi$n_no_events, mfi$n_files_empty), c(1L, 1L, 0L))
  expect_equal(qc_empty_grid_(e)$spleen[qc_empty_grid_(e)$metric == "median_CD25"],
               "ok (1 at 0 events)")
})

test_that("qc_missing_stats_ lists gates / stats only some files have", {
  d <- fake_qc_data()
  expect_equal(nrow(qc_missing_stats_(d)), 0)
  d1 <- dplyr::filter(d, !(mouse_ID == "m-4" & metric == "median_CD25"))
  m <- qc_missing_stats_(d1)
  expect_equal(nrow(m), 1)
  expect_equal(c(m$population, m$metrics, m$n_files, m$files),
               c("T", "median_CD25", "15/16", "missing in: m-4"))
  # a stat added to one sample only: metrics collapsed, "only in"
  extra <- dplyr::filter(d, mouse_ID == "m-7", population == "live",
                         metric == "count") |>
    dplyr::slice(c(1, 1)) |>
    dplyr::mutate(metric = c("median_A", "median_B"))
  m2 <- qc_missing_stats_(dplyr::bind_rows(d, extra))
  expect_equal(c(m2$population, m2$metrics, m2$n_files, m2$files),
               c("live", "median_A, median_B", "1/16", "only in: m-7"))
})

test_that("qc_missing_stats_ compares files within a tissue only", {
  d <- empty_qc_data_()
  d <- dplyr::filter(d, !(tissue == "liver" & population == "rare"))
  expect_equal(nrow(qc_missing_stats_(d)), 0)
})

test_that("qc_empty_keywords_ lists files with empty group / batch values", {
  d <- fake_qc_data()
  s <- fake_qc_settings(d)
  expect_equal(nrow(qc_empty_keywords_(d, s)), 0)
  d$grp[d$mouse_ID == "m-2"] <- NA
  d$day[d$mouse_ID == "m-2"] <- NA
  k <- qc_empty_keywords_(d, s)
  expect_equal(k$mouse_ID, c("m-2", "m-2"))
  expect_equal(k$empty, c("day", "grp"))
  expect_equal(nrow(qc_empty_keywords_(d, list())), 0)
})

test_that("qc_dropped_files_ lists files without mouse_ID or tissue", {
  d <- dplyr::bind_rows(
    fake_qc_data(n_per_group = 1),
    tibble::tibble(file_name = c("comp.fcs", "fmo.fcs"),
                   mouse_ID = c(NA, "m-1"), tissue = NA_character_,
                   metric = "count", value = 1)
  )
  expect_equal(qc_dropped_files_(d)$file_name, c("comp.fcs", "fmo.fcs"))
})

test_that("qc_whisker_side_() flags points beyond 1.5 x IQR", {
  expect_equal(qc_whisker_side_(c(1, 2, 2, 3, 2, 50)),
               c(rep(NA_character_, 5), "high"))
  expect_equal(qc_whisker_side_(c(-50, 1, 2, 2, 3, 2))[1], "low")
  expect_equal(qc_whisker_side_(c(NA, 1, 2)), rep(NA_character_, 3))
})

# Shapes of the built point layer, per row of p$data (same order).
point_shapes_ <- function(p) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]
  ggplot2::ggplot_build(p)$data[[i]]$shape
}

test_that("qc plots draw excluded files as crosses", {
  d <- fake_qc_data()
  s <- fake_qc_settings(d)
  sc <- qc_score_(d, s, excluded_files = "m-3_spleen.fcs")
  p <- qc_plot_individuals_(sc, s, "hierarchy")
  sh <- point_shapes_(p)
  # built rows are reordered by panel, so compare counts per shape
  expect_true(all(p$data$.excluded[p$data$mouse_ID == "m-3"] == "TRUE"))
  expect_equal(sum(sh == 4), sum(p$data$mouse_ID == "m-3"))
  expect_equal(sum(sh == 16), sum(p$data$mouse_ID != "m-3"))
  b <- qc_plot_batch_(sc, s, "counts", tissue = "spleen")
  shb <- point_shapes_(b)
  expect_equal(sum(shb == 4), sum(b$data$mouse_ID == "m-3"))
  expect_equal(sum(shb == 16), sum(b$data$mouse_ID != "m-3"))
  p0 <- qc_plot_individuals_(qc_score_(d, s), s, "hierarchy")
  expect_true(all(point_shapes_(p0) == 16))
  # older scored data without the column: all circles
  sc$excluded <- NULL
  expect_true(all(point_shapes_(qc_plot_individuals_(sc, s, "hierarchy")) == 16))
})

# Text layers' fixed colour, named by the non-empty labels each one draws.
label_layers_ <- function(p) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")),
                    logical(1)))
  built <- ggplot2::ggplot_build(p)$data
  stats::setNames(
    vapply(p$layers[i], function(l) l$aes_params$colour, character(1)),
    vapply(i, function(j) paste(sort(unique(setdiff(built[[j]]$label, ""))),
                                collapse = ","), character(1))
  )
}

test_that("outliers get black labels, excluded samples red ones", {
  d <- fake_qc_data()
  s <- fake_qc_settings(d)
  # m-3 is the outlier; exclude m-10 (not an outlier) first
  p <- qc_plot_individuals_(qc_score_(d, s, excluded_files = "m-10_spleen.fcs"),
                            s, "hierarchy")
  expect_equal(label_layers_(p), c("m-3" = "black", "m-10" = "red3"))
  expect_equal(unique(p$data$.label[p$data$.label != ""]), "m-3")
  # every row of the excluded mouse is labelled red, outlier or not
  expect_true(all(p$data$.label_excl[p$data$mouse_ID == "m-10"] == "m-10"))
  expect_true(all(p$data$.label_excl[p$data$mouse_ID != "m-10"] == ""))
  # an excluded outlier gets only the red label
  p2 <- qc_plot_individuals_(qc_score_(d, s, excluded_files = "m-3_spleen.fcs"),
                             s, "hierarchy")
  expect_true(all(p2$data$.label == ""))
  expect_true(all(p2$data$.label_excl[p2$data$mouse_ID == "m-3"] == "m-3"))
  expect_no_error(ggplot2::ggplot_build(p2))
})
