# Cohort overview tables for the analysis documents: meta_summary_table()
# (per group) and meta_mouse_table() (per mouse). The cohort comes from the
# meta workbook's `meta` sheet, which -- unlike the DB -- holds every enrolled
# mouse; DB-style names are accepted too, so the DB can become the source.

#' Cohort summary table per group
#'
#' One row per cohort characteristic, one column per group (plus `not
#' assigned` when an enrolled mouse has no group, and `total`): enrolled mice,
#' sex, age at BMT, strain, cages, BMT dates, dropouts before / at BMT, early
#' ends after BMT, exclusions, end-point day and the number analysed. Meant for
#' the top of an analysis document -- the numbers a methods section quotes.
#'
#' The cohort is the meta workbook's `meta` sheet (`meta_read(path)$meta`),
#' because the DB holds neither mice that dropped out before getting a
#' `mouse_ID` nor (in `db_load_experiment()` output) mice excluded mouse-wide.
#' A mouse is *excluded* when it is excluded mouse-wide (no tissue / panel /
#' condition / file / date restriction) in every import unit of `dat`; it is
#' *analysed* when it has end-point data (a sample on its end date) in at
#' least one import unit and is not excluded -- mice that ended early count.
#'
#' @param workbook The cohort: the workbook `meta` sheet as returned by
#'   [meta_read()] (`mouse_ID`, `id`, `group`, `sex`, `cage`, `line`,
#'   `generation`, `dob`, `bmt_date`, `death_date`, `early`, `reason`), or a
#'   table with DB names (`mouse_id`, `mouse_strain`, `end_date`, ...). `NULL`
#'   uses `dat$meta` (with a message: mice without a DB record are missing).
#' @param dat Optional [db_load_experiment()] list. Needed for the
#'   `Excluded from analysis` and `Analysed` rows (omitted without it) and for
#'   the group order (levels of `dat$meta$group`).
#' @param group_col Group column of `workbook` (default `"group"`).
#' @returns A `gt` table. Its data (`tbl[["_data"]]`) has a `characteristic`
#'   column and one character column per group, `not assigned` and `total`.
#' @examples
#' wb <- data.frame(
#'   mouse_ID = c("m1", "m2", NA), group = c("a", "b", "a"), sex = c("f", "m", "f"),
#'   bmt_date = as.Date("2026-05-12"),
#'   death_date = as.Date(c("2026-06-15", "2026-06-16", "2026-05-20")),
#'   early = c(NA, NA, "yes"), reason = c(NA, NA, "found dead")
#' )
#' meta_summary_table(wb)
#' @export
meta_summary_table <- function(workbook, dat = NULL, group_col = "group") {
  cohort <- cohort_input_(workbook, dat, group_col)
  if (!is.null(dat)) cohort <- analysis_status_(cohort, dat)
  lvls <- group_levels_(cohort, dat)
  cols <- c(lvls, if (anyNA(cohort$group)) "not assigned", "total")
  subsets <- lapply(cols, function(col) {
    if (col == "total") return(cohort)
    if (col == "not assigned") return(cohort[is.na(cohort$group), ])
    cohort[cohort$group %in% col, ]
  })

  cohort <- dplyr::mutate(
    cohort,
    .age_wk_ = as.numeric(.data$bmt_date - .data$dob) / 7,
    .dropped_ = .data$early & (is.na(.data$bmt_date) | .data$end_date <= .data$bmt_date),
    .early_after_ = .data$early & !is.na(.data$bmt_date) & .data$end_date > .data$bmt_date,
    .end_day_ = dplyr::if_else(.data$early, NA_real_,
                               as.numeric(.data$end_date - .data$bmt_date))
  )
  subsets <- lapply(subsets, function(s) cohort[cohort$.row_ %in% s$.row_, ])

  rows <- list(
    "Enrolled, n" = function(s) as.character(nrow(s)),
    "Sex, f / m" = function(s) paste(sum(s$sex %in% "f"), "/", sum(s$sex %in% "m")),
    "Age at BMT, weeks" = function(s) median_range_(s$.age_wk_, digits = 1),
    "Strain / generation" = function(s) {
      v <- unique(stats::na.omit(trimws(paste(
        dplyr::coalesce(s$strain, ""), dplyr::coalesce(s$generation, "")
      ))))
      v <- v[v != ""]
      if (length(v) == 0) "-" else paste(v, collapse = ", ")
    },
    "Cages, n" = function(s) as.character(dplyr::n_distinct(stats::na.omit(s$cage))),
    "BMT date" = function(s) {
      v <- sort(unique(stats::na.omit(s$bmt_date)))
      if (length(v) == 0) "-" else paste(format(v), collapse = ", ")
    },
    "Dropped before / at BMT, n" = function(s) as.character(sum(s$.dropped_)),
    "Ended early after BMT, n" = function(s) as.character(sum(s$.early_after_)),
    "Excluded from analysis, n" = function(s) {
      ids <- sort(s$mouse_id[s$excluded])
      if (length(ids) == 0) "0" else paste0(length(ids), " (", paste(ids, collapse = ", "), ")")
    },
    "End-point, days post-BMT" = function(s) median_range_(s$.end_day_, digits = 0),
    "Analysed, n" = function(s) as.character(sum(s$analysed))
  )
  if (is.null(dat)) rows[c("Excluded from analysis, n", "Analysed, n")] <- NULL
  if (all(is.na(cohort$.age_wk_))) rows[["Age at BMT, weeks"]] <- NULL

  body <- lapply(names(rows), function(label) {
    vals <- vapply(subsets, rows[[label]], character(1))
    tibble::as_tibble(as.list(stats::setNames(c(label, vals), c("characteristic", cols))))
  })
  body <- dplyr::bind_rows(body)

  tbl <- gt::gt(body) |>
    gt::cols_label(characteristic = "") |>
    gt::cols_align("center", columns = !"characteristic")
  for (spec in list(c("Dropped before / at BMT, n", ".dropped_"),
                    c("Ended early after BMT, n", ".early_after_"))) {
    note <- reason_note_(cohort$reason[cohort[[spec[2]]]])
    if (!is.null(note)) {
      tbl <- gt::tab_footnote(
        tbl, footnote = note,
        locations = gt::cells_body(columns = "characteristic",
                                   rows = body$characteristic == spec[1])
      )
    }
  }
  tbl
}

# Internal: the cohort to tabulate, normalised, with a row id.
cohort_input_ <- function(workbook, dat, group_col) {
  if (is.null(workbook)) {
    if (is.null(dat)) stop("Give `workbook` or `dat`.", call. = FALSE)
    message("No workbook given: using dat$meta -- mice without a DB record ",
            "(dropouts without mouse_ID, mice excluded mouse-wide) are not shown.")
    workbook <- dat$meta
  }
  if (group_col != "group") {
    workbook <- dplyr::rename(workbook, group = dplyr::all_of(group_col))
  }
  dplyr::mutate(cohort_normalise_(workbook), .row_ = dplyr::row_number())
}

# Internal: map a workbook `meta` sheet (meta_read() names) or a DB subjects
# table to one shape. Dates become Date; `early` is logical.
cohort_normalise_ <- function(cohort) {
  pick <- function(...) {
    for (nm in c(...)) if (nm %in% names(cohort)) return(cohort[[nm]])
    rep(NA, nrow(cohort))
  }
  early <- pick("early")
  early <- if (is.logical(early)) dplyr::coalesce(early, FALSE) else
    !is.na(early) & !tolower(trimws(as.character(early))) %in% c("", "no", "false")
  tibble::tibble(
    mouse_id = as.character(pick("mouse_ID", "mouse_id")),
    ear_tag = as.character(pick("id", "ID", "mdc_id")),
    group = as.character(pick("group", "treatment_group")),
    sex = as.character(pick("sex")),
    cage = as.character(pick("cage")),
    strain = as.character(pick("line", "mouse_strain", "strain")),
    generation = as.character(pick("generation")),
    dob = as.Date(pick("dob", "DOB")),
    bmt_date = as.Date(pick("bmt_date")),
    end_date = as.Date(pick("death_date", "end_date")),
    early = early,
    reason = as.character(pick("reason"))
  )
}

# Internal: group column order -- dat$meta's factor levels, else the cohort's.
group_levels_ <- function(cohort, dat) {
  present <- unique(stats::na.omit(cohort$group))
  lv <- if (!is.null(dat) && is.factor(dat$meta$group)) levels(dat$meta$group) else sort(present)
  c(intersect(lv, present), setdiff(sort(present), lv))
}

# Internal: the data tibbles of a db_load_experiment() list (not meta,
# exclusions or other helpers), named.
data_units_ <- function(dat) {
  keep <- vapply(seq_along(dat), function(i) {
    nm <- names(dat)[i]
    !nm %in% c("meta", "exclusions", "workbook") && is.data.frame(dat[[i]]) &&
      "mouse_id" %in% names(dat[[i]])
  }, logical(1))
  dat[keep]
}

# Internal: per-mouse analysis status against a db_load_experiment() list.
# Adds in_db, excluded, exclusion_reason, has_endpoint, analysed.
analysis_status_ <- function(cohort, dat) {
  units <- data_units_(dat)
  excl <- dat$exclusions
  if (is.null(excl)) excl <- tibble::tibble(mouse_id = character(), import_unit = character())
  unit_ids <- unique(unlist(lapply(units, function(u) {
    c(u[["import_unit"]], u[["wsp_root"]])
  })))
  unit_ids <- unique(c(unit_ids, excl$import_unit))
  excluded <- if (nrow(excl) > 0) fully_excluded_mice_(excl, unit_ids) else character()
  mouse_wide <- if (nrow(excl) > 0) {
    dplyr::filter(excl, is.na(.data$tissue), is.na(.data$panel), is.na(.data$condition),
                  is.na(.data$file_name), is.na(.data$collected_at))
  } else excl

  ends <- dplyr::select(cohort, "mouse_id", "end_date")
  if (!is.null(dat$meta) && "end_date" %in% names(dat$meta)) {
    ends <- dplyr::bind_rows(
      dplyr::transmute(dat$meta, mouse_id, end_date = as.Date(end_date)),
      ends
    )
  }
  ends <- dplyr::distinct(dplyr::filter(ends, !is.na(.data$mouse_id)), .data$mouse_id, .keep_all = TRUE)

  with_end <- unique(unlist(lapply(units, function(u) {
    if (!"collected_at" %in% names(u)) return(unique(u$mouse_id))
    d <- dplyr::inner_join(dplyr::distinct(u, .data$mouse_id, .data$collected_at),
                           ends, by = "mouse_id")
    unique(d$mouse_id[as.Date(d$collected_at) == d$end_date])
  })))
  in_db <- unique(c(dat$meta$mouse_id, excl$mouse_id,
                    unlist(lapply(units, function(u) u$mouse_id))))

  cohort |>
    dplyr::mutate(
      in_db = !is.na(.data$mouse_id) & .data$mouse_id %in% in_db,
      excluded = !is.na(.data$mouse_id) & .data$mouse_id %in% excluded,
      exclusion_reason = mouse_wide$reason[match(.data$mouse_id, mouse_wide$mouse_id)],
      exclusion_reason = dplyr::if_else(.data$excluded, .data$exclusion_reason, NA_character_),
      has_endpoint = !is.na(.data$mouse_id) & .data$mouse_id %in% with_end,
      analysed = .data$has_endpoint & !.data$excluded
    )
}

# Internal: "median (min-max)" or "-".
median_range_ <- function(x, digits) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  f <- function(v) formatC(v, format = "f", digits = digits)
  paste0(f(stats::median(x)), " (", f(min(x)), "-", f(max(x)), ")")
}

# Internal: footnote text "reason (n), ..." or NULL.
reason_note_ <- function(reasons) {
  reasons <- dplyr::coalesce(reasons, "no reason given")
  if (length(reasons) == 0) return(NULL)
  n <- table(reasons)
  paste(paste0(names(n), " (", as.integer(n), ")"), collapse = ", ")
}
