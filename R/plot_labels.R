#' Convert a label to a ggplot plotmath object (expression or character)
#'
#' Character strings with plotmath markup (`~`, `[`, `^`, `==`) are parsed;
#' existing expressions and plain text without markup are passed through.
#'
#' @param label Character scalar, expression, or language object.
#' @noRd
interactmeta_as_ggplot_label <- function(label) {
  if (is.null(label)) {
    return(NULL)
  }
  if (inherits(label, "expression")) {
    return(label)
  }
  if (is.language(label)) {
    return(as.expression(label))
  }
  if (!is.character(label) || length(label) != 1L) {
    return(label)
  }
  if (!nzchar(label)) {
    return(label)
  }
  if (grepl("~|\\[|\\]|\\^|==|\\*", label)) {
    parsed <- tryCatch(parse(text = label), error = function(e) NULL)
    if (!is.null(parsed)) {
      return(parsed)
    }
    # Operators such as `>=` in subgroup names are not valid plotmath;
    # never leak the source tildes as literal characters.
    return(gsub("~", " ", label, fixed = TRUE))
  }
  label
}

#' Resolve display subgroup label for overall-IF plots
#' @noRd
overall_ifrac_subgroup_label <- function(ifrac_var, interactmeta = NULL) {
  if (!is.null(interactmeta) && !is.null(interactmeta$ifrac_var) && nzchar(interactmeta$ifrac_var)) {
    lab <- interactmeta_ifrac_subgroup_label(interactmeta$ifrac_var)
    if (!is.null(lab) && nzchar(lab)) {
      return(lab)
    }
  }
  interactmeta_ifrac_subgroup_label(ifrac_var)
}

#' Subgroup label from an ifrac column name (e.g. `ifrac_Older` -> `Older`)
#' @noRd
interactmeta_ifrac_subgroup_label <- function(ifrac_var) {
  if (is.null(ifrac_var) || !nzchar(ifrac_var)) {
    return(NULL)
  }
  sub("^ifrac2?_", "", ifrac_var)
}

#' Default x-axis label for overall IF plots
#'
#' Plain character so subgroup names with operators (`>=65`) are not
#' parsed as plotmath (which would leave `Information~fraction~...` on the
#' axis).
#' @noRd
overall_ifrac_default_xlab <- function(subgroup_label = NULL) {
  if (is.null(subgroup_label) || !nzchar(subgroup_label)) {
    return("Information fraction")
  }
  sprintf("Information fraction (%s)", subgroup_label)
}

#' Default x-axis label for [plot_contribution()] from an `ifrac_*` column name
#' @noRd
bubbleplot_default_xlab <- function(ifrac_var) {
  overall_ifrac_default_xlab(interactmeta_ifrac_subgroup_label(ifrac_var))
}

#' Parse legend labels that may use plotmath strings
#' @noRd
interactmeta_parse_legend_labels <- function(labels, fallback) {
  if (is.null(labels)) {
    return(fallback)
  }
  out <- vapply(
    seq_along(labels),
    function(i) {
      lab <- labels[[i]]
      if (is.na(lab)) {
        return(fallback[[i]])
      }
      parsed <- interactmeta_as_ggplot_label(lab)
      if (inherits(parsed, "expression")) {
        return(parsed[[1]])
      }
      if (is.language(parsed)) {
        return(parsed)
      }
      parsed
    },
    FUN.VALUE = quote(x),
    USE.NAMES = FALSE
  )
  stats::setNames(as.list(out), names(labels))
}
