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
  # Create path to WS-file in local directory
  path <- here("data", path)

  # Import raw workspace
  ps_raw <- fcexpr::wsx_get_popstats(ws = path, return_stats = r_stats, groups = group)

  # Merge counts and stats if both are extracted
  if (r_stats == TRUE) {
    ps_clean <- left_join(
      ps_raw[["counts"]],
      ps_raw[["stats"]],
      by = c("FileName", "PopulationFullPath")
    ) %>%
      select(FileName, PopulationFullPath, Population, Count, FractionOfParent, statistic, channel, value) %>%
      mutate(
        channel = str_remove_all(string = channel, pattern = c("Comp-|-A"))
      ) %>%
      pivot_wider(
        names_from = c(statistic, channel),
        names_sep = "_",
        values_from = value
      ) %>%
      select(!NA_NA)
    # Return only counts if stats are not extracted
  } else {
    ps_clean <- ps_raw %>%
      select(FileName, PopulationFullPath, Population, Count, FractionOfParent)
  }
  # Extract keywords from workspace
  keys <- fcexpr::wsx_get_keywords(ws = path) %>%
    enframe() %>%
    rename(
      FileName = name
    ) %>%
    unnest(c(FileName, value)) %>%
    rename(key = name) %>%
    filter(key %in% c(keywords)) %>%
    pivot_wider(
      names_from = key,
      values_from = value
    )
  # Merge workspace data with keywords
  left_join(ps_clean, keys, by = "FileName")
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

  excel_path <- paste0(dirname(path), "/clean", file_name, ".xlsx")

  if(clean == TRUE) {

    data <- JanisHelpers::import_workspace(
      path = path,
      group = group,
      r_stats = r_stats,
      keywords = keywords
    ) |>
      dplyr::tibble()

    xlsx::write.xlsx(data, file = excel_path)

    print(paste("clean data imported from", wsp, "and written to", excel_path))
  }

  if(clean == FALSE) {

    data <- read_excel(excel_path) |>
      select(!...1)

    print(paste("imported data from Excel:", file_name))

  }

  data |>
    dplyr::tibble()
}
