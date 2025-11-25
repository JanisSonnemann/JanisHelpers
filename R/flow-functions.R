# Function to import data from FlowJo workspace
#' Import FlowJo workspace data using fcexpr package
#'
#' @param path path to FlowJo Space
#' @param group group of samples from workspace to include in analysis
#' @param r_stats import stats such as MFIs?
#' @param keywords which keywords to export
#'
#' @return
#' @export
#'
#' @examples
import_workspace <- function(path, group, r_stats, keywords) {

  # Import raw workspace
  ps_raw <- fcexpr::wsx_get_popstats(ws = path, return_stats = r_stats, groups = group)

  # Import keywords
  keys <- fcexpr::wsx_get_keywords(ws = path) |>
    enframe() |>
    rename(
      FileName = name
    ) |>
    unnest(c(FileName, value)) |>
    rename(key = name)

  # Extract channel labels from keywords
  stains <- keys |>
    filter(grepl("\\$P[0-9]+[NS]", key)) |>
    mutate(
      type = case_when(
        grepl("\\$P[0-9]+N", key) ~ "channel",
        grepl("\\$P[0-9]+S", key) ~ "stain"
      ),
      number = str_extract(key, "(?<=P)[0-9]+")
    ) |>
    select(!key) |>
    pivot_wider(
      names_from  = type,
      values_from = value
    ) |>
    select(!number)

  # Merge counts and stats if both are extracted
  if (r_stats == TRUE) {

    ps_clean <- left_join(
      ps_raw[["counts"]],
      ps_raw[["stats"]],
      by = c("FileName", "PopulationFullPath")
    ) |>
      select(FileName, PopulationFullPath, Population, Count, FractionOfParent, statistic, channel, value) |>
      left_join(
        y = stains,
        by = c("FileName", "channel")
      ) |>
      mutate(
        stain = ifelse(stain == "", NA, stain)
      ) |>
      fill(stain, .direction = "downup") |>
      mutate(
        channel = str_remove_all(string = channel, pattern = c("Comp-|-A"))
      ) |>
      pivot_wider(
        names_from = c(statistic, channel, stain),
        names_sep = "_",
        values_from = value
      ) |>
      select(!contains("NA_"))
  }
  # Return only counts if stats are not extracted
  else {
    ps_clean <- ps_raw |>
      select(FileName, PopulationFullPath, Population, Count, FractionOfParent)
  }
  # filter keywords
  keys_clean <- keys |>
    filter(key %in% c(keywords)) %>%
    pivot_wider(
      names_from = key,
      values_from = value
    )

  # Merge workspace data with keywords
  left_join(ps_clean, keys_clean, by = "FileName")
}


#' Function to automate import of FlowJo workspace data automatically. Enables clean import from .wsp with
#' simultaneous writing to excel.
#'
#'
#' @param clean logic, decide if data should be newly imported from .wsp and then saved to excel or
#' simply from previously created Excel
#' @param location name of subdirectory where data is located
#' @param wsp name of workspace file
#' @param group name of group in which samples are organized
#' @param r_stats logic indicating if stats (MFI etc.) should be imported
#' @param keywords character vector of keywords to be imported
#'
#' @return
#' @export
#'
#' @examples
import_fcs_clean <- function(clean, location, wsp, group, r_stats, keywords) {

  if(clean == TRUE) {

    data <- JanisHelpers::import_workspace(
      path = paste0(location, "/", wsp),
      group = group,
      r_stats = r_stats,
      keywords = keywords
    ) |>
      dplyr::tibble()

    file_name <- tools::file_path_sans_ext(wsp)
    file_path <- paste0("/data/", location, "/clean/", file_name, ".xlsx")

    xlsx::write.xlsx(data, file = paste0(here(), file_path))

    print(paste("clean data imported from", wsp, "and written to", file_path))
  }

  if(clean == FALSE) {

    file_name <- tools::file_path_sans_ext(wsp)

    data <- read_excel(
      paste0(here(), "/data/", location, "/clean/", file_name, ".xlsx")
    ) |>
      select(!...1)

    print(paste("imported data from Excel:", file_name))

  }

  data |>
    dplyr::tibble()
}


#' Function to automate import of FlowJo workspace data automatically and work with lists. Enables clean import from .wsp with
#' simultaneous writing to excel.
#'
#' @param clean logic, decide if data should be newly imported from .wsp and then saved to excel or
#' simply from previously created Excel
#' @param path location of wsp
#' @param group name of group in which samples are organized
#' @param r_stats logic indicating if stats (MFI etc.) should be imported
#' @param keywords character vector of keywords to be imported
#'
#' @returns
#' @export
#'
#' @examples
import_fcs <- function(clean, path, group, r_stats, keywords) {

  file_name <- path |>
    basename() |>
    tools::file_path_sans_ext()

  excel_path <- paste0(dirname(path), "/clean/", file_name, ".xlsx")

  if(clean == TRUE) {

    data <- JanisHelpers::import_workspace(
      path = path,
      group = group,
      r_stats = r_stats,
      keywords = keywords
    ) |>
      dplyr::tibble()

    xlsx::write.xlsx(data, file = excel_path)

    print(paste("clean data imported from", file_name, "and written to", excel_path))
  }

  if(clean == FALSE) {

    data <- read_excel(excel_path) |>
      select(!...1)

    print(paste("imported data from Excel:", file_name))

  }

  data |>
    dplyr::tibble()
}


#' Import FACS data from FlowJo workspace in long format
#'
#' @param path path to .wsp file
#' @param group group to be extracted from workspace default = NULL
#' @param r_stats true or false to indicate if statistics such as MFI should be exported
#' @param keywords character vector of keywords to be exported, default = NULL
#'
#' @returns
#' @export
#'
#' @examples
import_workspace_long <- function(path, group = NULL, r_stats = FALSE, keywords = NULL) {

  # Import raw workspace
  ps_raw <- suppressMessages(fcexpr::wsx_get_popstats(ws = path, return_stats = r_stats, groups = group))

  ## extract counts
  counts <- ps_raw[["counts"]] |>
    select(FileName, PopulationFullPath, Population, Count, FractionOfParent) |>
    pivot_longer(
      cols = c(Count, FractionOfParent),
      names_to = "metric",
      values_to = "value"
    )

  # Import keywords
  keys <- fcexpr::wsx_get_keywords(ws = path, return = "data.frame") |>
    enframe() |>
    rename(
      FileName = name
    ) |>
    unnest(c(FileName, value)) |>
    rename(key = name)

  # If stats should be extracted: merge stats and counts
  if (r_stats == TRUE) {

    ## Extract channel labels from keywords
    stains <- keys |>
      filter(grepl("\\$P[0-9]+[NS]", key)) |>
      mutate(
        type = case_when(
          grepl("\\$P[0-9]+N", key) ~ "channel",
          grepl("\\$P[0-9]+S", key) ~ "stain"
        ),
        number = str_extract(key, "(?<=P)[0-9]+")
      ) |>
      select(!key) |>
      pivot_wider(
        names_from  = type,
        values_from = value
      ) |>
      select(!number)

    ## subset stats from list
    stats <- ps_raw[["stats"]] |>
      ## add channel labels
      left_join(
        y = stains,
        by = c("FileName", "channel")
      ) |>
      ## create new column which combined statistic name and channel label
      mutate(
        stain = ifelse(stain == "", NA, stain),
        metric = ifelse(!is.na(statistic), paste0(statistic, "_", stain), NA),
        Population = basename(PopulationFullPath)
      ) |>
      select(FileName, PopulationFullPath, Population, metric, value)

    # combined counts and stats
    df <- bind_rows(
      counts,
      stats
    )

  }
  # Return only counts if stats are not extracted
  else {

    # only take counts from list
    df <- counts

  }
  # extract selected keywords and pivot to wide
  keys_clean <- keys |>
    filter(key %in% c(keywords)) |>
    pivot_wider(
      names_from = key,
      values_from = value
    )

  # Merge workspace data with keywords
  final <- left_join(df, keys_clean, by = "FileName") |>
    arrange(FileName, PopulationFullPath) |>
    relocate(mouse_ID, tissue, PopulationFullPath, Population, metric, value)

  # create summary to make sure every thing was exported correctly
  n_files <- n_distinct(final$FileName)
  n_mice <- n_distinct(final$mouse_ID)
  n_tissues <- n_distinct(final$tissue)

  summary_text <- glue::glue(
    "
Extraction Summary
----------------------------------------------
Extracted groups:          {group}
Number of samples:         {n_files}
Unique mouse_IDs:          {n_mice}
Unique tissues:            {n_tissues}
----------------------------------------------
"
  )
  ## print summary message
  message(summary_text)

  # return final data invisibly
  invisible(final)

}
