# Internal histology QC checks on the long table from histo_read_curated().

# Plausibility rules, e.g. "crescents <= gloms". `data` is pivoted wide by
# metric (one row per mouse x tissue x timepoint) and each rule is evaluated
# there; rows where it is FALSE (NA = not a violation) are returned as
# rule, mouse_ID, tissue, timepoint, values ("metric = value" pairs).
histo_qc_rules_ <- function(data, rules) {
  empty <- tibble::tibble(rule = character(), mouse_ID = character(),
                          tissue = character(), timepoint = character(),
                          values = character())
  if (length(rules) == 0L) {
    return(empty)
  }
  data <- dplyr::mutate(data, timepoint = qc_timepoint_(collected_at),
                        mouse_ID = as.character(mouse_ID))
  wide <- data |>
    dplyr::select(mouse_ID, tissue, timepoint, metric, value) |>
    tidyr::pivot_wider(names_from = metric, values_from = value,
                       values_fn = dplyr::first)
  keys <- c("mouse_ID", "tissue", "timepoint")
  metrics <- setdiff(names(wide), keys)
  values <- apply(wide[metrics], 1, function(r) {
    paste(metrics[!is.na(r)], format(r[!is.na(r)], trim = TRUE, digits = 4),
          sep = " = ", collapse = ", ")
  })
  out <- lapply(rules, function(rule) {
    res <- tryCatch(
      eval(parse(text = rule), wide, baseenv()),
      error = function(e) {
        stop("Rule `", rule, "` could not be evaluated (is every metric in ",
             "the data?): ", conditionMessage(e), call. = FALSE)
      }
    )
    if (!is.logical(res) || length(res) != nrow(wide)) {
      stop("Rule `", rule, "` must give one TRUE / FALSE per sample.",
           call. = FALSE)
    }
    bad <- which(!is.na(res) & !res)
    tibble::tibble(rule = rep(rule, length(bad)), wide[bad, keys],
                   values = unname(values[bad]))
  })
  dplyr::bind_rows(empty, out)
}
