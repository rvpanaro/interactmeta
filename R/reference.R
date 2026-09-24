#' Resolve and validate the CAMS reference subgroup
#'
#' The reference subgroup is a single source of truth for design columns
#' (`ifrac_<nonref>`, `ifrac2_<nonref> = x - π`), shrinkage
#' (`y_A` = reference, `y_B` = non-reference, `π` = IF of non-reference),
#' and forest / contribution-plot labelling.
#'
#' Resolution order:
#' 1. Explicit `reference` argument (if supplied).
#' 2. Constant `data$reference` / `fit` attribute `interactmeta$reference`.
#'
#' Explicit and stored values must agree when both are present. There is no
#' silent fallback to "first factor level" for fitted objects.
#'
#' @param x A [`interactmeta`], or a data frame with optional `reference` column.
#' @param reference Optional explicit reference label that must agree with `x`.
#' @param subgroup_var Subgroup column name when `x` is a data frame.
#' @return A length-1 character reference label.
#' @keywords internal
#' @noRd
interactmeta_reference <- function(x,
                             reference = NULL,
                             subgroup_var = "subgroup") {
  if (is_comparison_fit(x) || inherits(x, "interactmeta")) {
    st <- attr(x, "interactmeta", exact = TRUE)
    stored <- if (!is.null(st) && !is.null(st$reference) &&
                  length(st$reference) &&
                  !is.na(st$reference[[1L]]) &&
                  nzchar(as.character(st$reference)[[1L]])) {
      as.character(st$reference)[[1L]]
    } else if (!is.null(x$data) && "reference" %in% names(x$data)) {
      interactmeta_reference_from_column(x$data$reference, allow_missing = TRUE)
    } else {
      NULL
    }
    if (is.null(stored) || !nzchar(stored)) {
      stop(
        "Fit has no stored `reference`. Refit with `cam(..., reference = ...)`.",
        call. = FALSE
      )
    }
    if (!is.null(reference)) {
      interactmeta_assert_reference_agrees(stored, reference, where = "fit")
    }
    return(stored)
  }

  if (is.data.frame(x)) {
    col_ref <- if ("reference" %in% names(x)) {
      interactmeta_reference_from_column(x$reference, allow_missing = TRUE)
    } else {
      NULL
    }
    if (!is.null(reference)) {
      ref <- as.character(reference)[[1L]]
      if (!is.null(col_ref)) {
        interactmeta_assert_reference_agrees(col_ref, ref, where = "data")
        return(ref)
      }
      if (subgroup_var %in% names(x)) {
        labs <- unique(as.character(x[[subgroup_var]]))
        labs <- labs[!is.na(labs) & nzchar(labs)]
        labs <- labs[is_ifrac_component_row(labs)]
        if (length(labs) && !(ref %in% labs)) {
          stop(
            "`reference` (", ref, ") must appear in `", subgroup_var, "`.",
            call. = FALSE
          )
        }
      }
      return(ref)
    }
    if (!is.null(col_ref)) {
      return(col_ref)
    }
    stop(
      "Supply `reference` or include a single `reference` column.",
      call. = FALSE
    )
  }

  stop("`x` must be an interactmeta object or a data frame.", call. = FALSE)
}

#' Non-reference subgroup label(s) given reference and design / data
#' @noRd
interactmeta_nonref <- function(x,
                          reference = NULL,
                          subgroup_var = "subgroup") {
  ref <- interactmeta_reference(x, reference = reference, subgroup_var = subgroup_var)
  if (is_comparison_fit(x) || inherits(x, "interactmeta")) {
    st <- attr(x, "interactmeta", exact = TRUE)
    if (!is.null(st$ifrac_vars) && length(st$ifrac_vars)) {
      return(sub("^ifrac_", "", as.character(st$ifrac_vars)))
    }
    if (!is.null(st$ifrac_var) && nzchar(st$ifrac_var)) {
      return(sub("^ifrac_", "", as.character(st$ifrac_var)[[1L]]))
    }
    dat <- st$ifrac_data %||% x$data
  } else {
    dat <- x
  }
  if (is.null(dat) || !subgroup_var %in% names(dat)) {
    stop("Cannot resolve non-reference subgroup label(s).", call. = FALSE)
  }
  labs <- unique(as.character(dat[[subgroup_var]]))
  labs <- labs[!is.na(labs) & nzchar(labs)]
  labs <- labs[!grepl("-", labs, fixed = TRUE)]
  setdiff(labs, ref)
}

#' @noRd
interactmeta_reference_from_column <- function(values, allow_missing = FALSE) {
  vals <- unique(as.character(values))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) {
    if (isTRUE(allow_missing)) {
      return(NULL)
    }
    stop("`reference` column has no non-missing values.", call. = FALSE)
  }
  if (length(vals) > 1L) {
    stop(
      "`reference` must be a single subgroup label; found: ",
      paste(vals, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  vals[[1L]]
}

#' @noRd
interactmeta_assert_reference_agrees <- function(stored, reference, where = "data") {
  ref <- as.character(reference)[[1L]]
  stored <- as.character(stored)[[1L]]
  if (!identical(stored, ref)) {
    stop(
      "`reference = \"", ref, "\"` does not match the stored reference \"",
      stored, "\" on the ", where,
      ". Rebuild with `ifrac(..., reference = \"", ref, "\")` or omit ",
      "`reference` to use the stored value.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Ensure `data$reference` is present and constant
#' @noRd
interactmeta_ensure_reference_column <- function(data, reference) {
  ref <- as.character(reference)[[1L]]
  if ("reference" %in% names(data)) {
    interactmeta_assert_reference_agrees(
      interactmeta_reference_from_column(data$reference),
      ref,
      where = "data"
    )
  }
  data$reference <- ref
  data
}
