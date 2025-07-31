#' Function to create descriptive statistics table with gtsummary
#'
#' @param data dataframe to be analyzed
#' @param group_col contains grouping variable
#' @param tissue_col if the dataframe contains multiple rows from the same animal (such as different tissues) this will subset the data
#' @param cols specify columns to be analyzed. if blank all numeric columns are analyzed
#'
#' @returns
#' @export
#'
#' @examples
create_descriptive_table <- function(data, group_col = "group",
                                     tissue_col = "tissue", cols = NULL) {
  # extract numeric variables to analyze
  if (is.null(cols)) {
    cols <- names(data)[sapply(data, is.numeric)]
    cols <- cols[!cols %in% c("mouse_ID", group_col, tissue_col)]
  }

  # function to create summary table with descriptive statistics
  summary_tbl <- function(input) {
    input |>
      gtsummary::tbl_summary(
        include = where(is.numeric),
        by = "group",
        statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
        digits = all_continuous() ~ 2,
        missing = "no"
      ) |>
      gtsummary::add_p(test = all_continuous() ~ "kruskal.test") %>%
      gtsummary::add_overall() %>%
      gtsummary::modify_header(label ~ paste("**Parameter**")) %>%
      gtsummary::bold_labels() %>%
      gtsummary::italicize_levels() |>
      gtsummary::as_gt() |>
      gt::tab_style(
        style = cell_fill(color = "lightgreen"),
        locations = cells_body(
          columns = p.value,
          rows = p.value < 0.05
        )
      ) |>
      gt::tab_style(
        style = cell_fill(color = "lightcoral"),
        locations = cells_body(
          columns = p.value,
          rows = p.value > 0.05
        )
      )
  }

  ### if tissue column is specified data is subset according to tissue column
  if(!is.null(tissue_col)) {

    # extract specific tissues
    tissues <- as.character(unique(data[[tissue_col]]))

    ## filter function to filter data by tissue
    filter_fun <- function(dat, tis) {
      dat %>%
        filter(.data$tissue == tis)
    }

    ## mapping function to create subsetted list by tissue
    tissue_data <- map(
      .x = tissues,
      .f = ~filter_fun(dat = data, tis = .x)
    ) |>
      set_names(tissues)

    # summarize data by tissue
    ## summarize function to create table for every tissue
    summary_fun <- function(dat, tis) {
      dat[[tis]] |>
        summary_tbl() |>
        tab_caption(caption = md(paste0("**Descriptive Statistics for ", tis, "**")))
    }

    ## mapping function to map summary function over all tissues
    desc_tables <- map(
      .x = tissues,
      .f = ~summary_fun(dat = tissue_data, tis = .x)
    ) |>
      set_names(tissues)
  }
  ### if tissue column is NULL no subsetting is needed and summary table is immediately created###
  else {
    desc_tables <- data |>
      summary_tbl()
  }
  # return descriptive tables
  return(desc_tables)
}



#' Function to create summary table of post-hoc tests. Performs Kruskal-Wallis on all numeric columns and post-hoc Dunn's on all significant columns.
#'
#' @param data dataframe to be analyzed
#' @param group_col contains grouping variable
#' @param tissue_col if the dataframe contains multiple rows from the same animal (such as different tissues) this will subset the data
#'
#' @returns
#' @export
#'
#' @examples
create_posthoc_tables <- function(data, group_col = "group",
                                  tissue_col = "tissue") {

  # summary function to perform kruskal wallis and post-hoc dunns test, creates output in gt format
  posthoc_fun <- function(dat, tis, var, group_col) {
    # create formula for statistical calculations
    form <- reformulate(group_col, response = as.name(var))

    # perform kruskal wallis test to check for significant global differences
    kw <- dat |>
      rstatix::kruskal_test(form)
    kw_p <- kw$p
    kw_stat <- kw$statistic

    # if significant global difference is present perform post-hoc Dunn's test
    if (kw_p < 0.05) {
      dunn <- rstatix::dunn_test(
        data = dat,
        formula = form,
        p.adjust.method = "bonferroni"
      ) |>
        gt::gt() |>
        gt::tab_header(
          title = md(
            if (nzchar(tis)) {
              paste0("**", tis, ":** ", var)
            } else {
              paste0("**", var, "**")
            }
          ),
          subtitle = paste("Kruskal-Wallis p-value:", kw_p)
        ) %>%
        gt::tab_footnote(footnote = "Dunns Post-Hoc Test (Bonferroni-adjusted)") |>
        gt::cols_label(
          group1 = "Group 1",
          group2 = "Group 2",
          p = "Raw p-value",
          p.adj = "Adjusted p-value",
          p.adj.signif = "Significance"
        ) %>%
        gt::fmt_number(columns = vars(p.adj, statistic, p), decimals = 4) %>%
        gt::tab_style(
          style = cell_fill(color = "lightgreen"),
          locations = cells_body(columns = vars(p.adj.signif, p.adj), rows = p.adj < 0.05)
        ) |>
        gt::cols_hide(.y.)
    }
    # return message if KW is not significant
    else {
      dunn <- paste0("Kruskal-Wallis found p-value of ", kw_p, ". No post-hoc testing performed")
    }
    # return results of post-hoc test
    return(dunn)
  }

  ### if tissue column is specified data is subset according to tissue and post hoc testing applied to all numeric columns per tissue
  if(!is.null(tissue_col)) {

    tissues <- as.character(unique(data[[tissue_col]]))

    # function to subset dataframe by tissue and then apply post-hoc function to all populations
    all_pops_fun <- function(dat, tis, group_col) {

      ## subset data according to tissue
      data <- dat |>
        filter(.data$tissue == tis)

      ## extract names of numeric variable columns
      cols <- data |>
        select(where(is.numeric)) |>
        names()

      ## apply post-hoc analysis function to all numeric cols in one tissue
      map(
        .x = cols,
        .f = ~posthoc_fun(dat = data, tis, var = .x, group_col = group_col)
      ) |>
        set_names(cols)
    }

    # appply analysis function for all pops to every tissue
    tab <- map(
      .x = tissues,
      .f = ~all_pops_fun(data, tis = .x, group_col = group_col)
    ) |>
      set_names(tissues)

  }
  ### if tissue column is NULL no subsetting is needed and summary table is immediately created
  else {
    # function to apply post-hoc function to all populations in dataframe
    all_pops_fun <- function(dat, group_col) {

      ## no subsetting of data
      data <- dat

      ## extract names of numeric variable columns
      cols <- data |>
        select(where(is.numeric)) |>
        names()

      ## apply post-hoc function to all numeric cols
      map(
        .x = cols,
        .f = ~posthoc_fun(data, "", var = .x, group_col = group_col)
      ) |>
        set_names(cols)
    }

    tab <- all_pops_fun(data, group_col)

  }
  # return post-hoc tables
  return(tab)
}

