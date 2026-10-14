#' Convert a CAMS effect summary to a tidy data frame
#'
#' @param summary Summary data frame from an effects object.
#' @param id_cols Identifier columns to keep (e.g. `subgroup`, `type`).
#' @noRd
interactmeta_tidy_effects_summary <- function(summary,
                                      id_cols = character()) {
  if (is.null(summary) || !nrow(summary)) {
    return(tibble::tibble())
  }

  out <- as.data.frame(summary)
  if ("prob" %in% names(out)) {
    out$prob <- NULL
  }

  if ("sd" %in% names(out)) {
    names(out)[names(out) == "sd"] <- "std.error"
  }
  if ("lci" %in% names(out)) {
    names(out)[names(out) == "lci"] <- "conf.low"
  }
  if ("uci" %in% names(out)) {
    names(out)[names(out) == "uci"] <- "conf.high"
  }

  stat_cols <- intersect(
    c("estimate", "mean", "std.error", "conf.low", "conf.high"),
    names(out)
  )
  id_keep <- intersect(
    c(id_cols, "component", "parameter", "label", "trial", "subgroup", "interaction", "type"),
    names(out)
  )
  out <- out[, c(id_keep, stat_cols), drop = FALSE]
  tibble::as_tibble(out)
}

#' Tidy a `interactmeta_subgroup_effects` object
#'
#' Returns a [tibble::tibble] with identifying columns (e.g. `subgroup`,
#' `type`) and summary columns (`estimate`, `std.error`, `conf.low`,
#' `conf.high`). Printing uses the ordinary tibble display; for the
#' biostatistical console layout, print the `interactmeta_*` object itself.
#'
#' @param x A `interactmeta_subgroup_effects` object from [subgroups()].
#' @param ... Ignored.
#' @return A [tibble::tibble].
#' @export
#' @importFrom generics tidy
tidy.interactmeta_subgroup_effects <- function(x, ...) {
  id_cols <- if (isTRUE(x$shrinkage)) {
    c("trial", "subgroup", "type")
  } else {
    c("subgroup", "type")
  }
  interactmeta_tidy_effects_summary(x$summary, id_cols = id_cols)
}

#' Tidy an `interactmeta_combined_effect` object
#'
#' @param x An `interactmeta_combined_effect` object from [combined()].
#' @param ... Ignored.
#' @return A [tibble::tibble].
#' @export
#' @importFrom generics tidy
tidy.interactmeta_combined_effect <- function(x, ...) {
  id_cols <- if (isTRUE(x$shrinkage)) {
    c("trial", "type")
  } else {
    "type"
  }
  interactmeta_tidy_effects_summary(x$summary, id_cols = id_cols)
}

#' Tidy a `interactmeta_interaction_effects` object
#'
#' @param x A `interactmeta_interaction_effects` object from [interactions()].
#' @param ... Ignored.
#' @return A [tibble::tibble].
#' @export
#' @importFrom generics tidy
tidy.interactmeta_interaction_effects <- function(x, ...) {
  summ <- x$summary
  if (isTRUE(x$shrinkage) && !is.null(summ)) {
    summ <- comparison_interaction_shrinkage_summary_once(summ)
  }
  id_cols <- if (isTRUE(x$shrinkage)) {
    c("trial", "interaction", "type")
  } else {
    c("interaction", "type")
  }
  interactmeta_tidy_effects_summary(summ, id_cols = id_cols)
}

#' Tidy a `interactmeta_heterogeneity` object
#'
#' @param x A `interactmeta_heterogeneity` object from [heterogeneity()].
#' @param ... Ignored.
#' @return A [tibble::tibble].
#' @export
#' @importFrom generics tidy
tidy.interactmeta_heterogeneity <- function(x, ...) {
  interactmeta_tidy_effects_summary(
    x$summary,
    id_cols = c("parameter", "label", "type")
  )
}
