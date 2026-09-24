#' Column that holds the two levels of a within-study contrast
#'
#' `modifier` is that column name. The default `"subgroup"` keeps the shipped
#' examples working. `subgroup` is still accepted. Level values are left as
#' stored. A display heading capitalizes only the first letter of the column
#' name, so `dose` is shown as `Dose`.
#'
#' @noRd
NULL

#' Resolve the modifier column name
#'
#' `modifier` wins, then `subgroup`, then `fallback`.
#' @noRd
interactmeta_resolve_modifier_column <- function(modifier = NULL,
                                                subgroup = NULL,
                                                fallback = "subgroup") {
  pick <- function(x) {
    if (is.null(x) || !length(x)) {
      return(NULL)
    }
    x <- as.character(x)[[1L]]
    if (!nzchar(x)) NULL else x
  }
  pick(modifier) %||% pick(subgroup) %||% pick(fallback) %||% "subgroup"
}

#' Display heading for a modifier column
#'
#' Capitalizes the first letter and leaves the rest unchanged.
#' @noRd
interactmeta_modifier_heading <- function(name) {
  name <- as.character(name %||% "subgroup")[[1L]]
  if (!nzchar(name)) {
    name <- "subgroup"
  }
  paste0(toupper(substring(name, 1L, 1L)), substring(name, 2L))
}
