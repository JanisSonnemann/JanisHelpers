#' update_JanisHelpers_git: update from private git repo
#'
#' @return
#' @export
#'
#' @examples
update_JanisHelpers_git <- function() {
  devtools::install_github("JanisSonnemann/JanisHelpers", auth_token = gh::gh_token())
}

