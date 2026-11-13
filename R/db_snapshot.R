#' Write a dated snapshot of the database
#'
#' Exports the whole database (`EXPORT DATABASE`, Parquet) to
#' `<root>/<YYYY-MM-DD>/`. A snapshot already written the same day is
#' replaced. The connection stays open. Keep the live `.duckdb` file local
#' and put only snapshots on synced storage (OneDrive sync breaks DuckDB's
#' lock and `.wal` files).
#'
#' Restore into a new file with
#' `DBI::dbExecute(con, "IMPORT DATABASE '<snapshot dir>'")`.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param root Folder holding the dated snapshot folders (created if
#'   missing).
#' @returns Invisibly, the path of the snapshot folder.
#' @examples
#' \dontrun{
#' con <- db_connect("~/databases/lab.duckdb")
#' db_snapshot(con, "data/z_snapshots")
#' }
#' @export
db_snapshot <- function(con, root) {
  dir <- file.path(root, format(Sys.Date(), "%Y-%m-%d"))
  unlink(dir, recursive = TRUE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir, mustWork = FALSE)
  DBI::dbExecute(con, sprintf("EXPORT DATABASE '%s' (FORMAT PARQUET)", dir))
  message("Snapshot written to ", dir)
  invisible(dir)
}
