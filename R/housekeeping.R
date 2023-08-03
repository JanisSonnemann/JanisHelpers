#' update_JanisHelpers: automatically update the package
#'
#' @return
#' @export
#'
#' @examples
update_JanisHelpers <- function(){
  devtools::install("C:/Users/janis/Documents/R/eigene R Pakete/JanisHelpers")
}

update_JanisHelpers_git <- function() {
  devtools::install_github("JanisSonnemann/JanisHelpers", auth_token = gh::gh_token())
}

