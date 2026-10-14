#' Display acronym for an analysis family
#' @noRd
interactmeta_analysis_acronym <- function(analysis) {
  switch(
    as.character(analysis)[[1L]],
    cam = "CAMS",
    da = "DA",
    ad = "AD",
    swada = "SWADA",
    conventional = "DA",
    toupper(as.character(analysis)[[1L]])
  )
}

#' Display label for an analysis family
#' @noRd
interactmeta_analysis_display <- function(x) {
  interactmeta_analysis_acronym(interactmeta_fit_analysis(x))
}

#' Inferential qualifier stored on a fit (`engine` / `method`)
#'
#' Bayesian CAMS is `"bayesmeta"`. Frequentist CAMS includes the estimator,
#' e.g. `"metafor (REML)"`.
#' @noRd
interactmeta_inference_display <- function(x) {
  engine <- as.character(x$engine %||% NA_character_)[[1L]]
  method <- as.character(x$method %||% NA_character_)[[1L]]
  inf <- as.character(x$inference %||% NA_character_)[[1L]]
  if (identical(engine, "bayesmeta") || identical(inf, "bayesian") || identical(method, "bayes")) {
    return("bayesmeta")
  }
  if (identical(engine, "metafor") || identical(inf, "frequentist") || method %in% c("REML", "ML")) {
    meth <- if (method %in% c("REML", "ML")) method else "REML"
    return(paste0("metafor (", meth, ")"))
  }
  interactmeta_fit_engine_method(x)
}

#' Approach label for a single fit
#'
#' When every compared fit has a distinct analysis family, the short name
#' (`CAMS`, `DA`, …) is enough. When the same analysis appears more than
#' once, the approach is appended (`CAMS — bayesmeta`).
#' @noRd
interactmeta_fit_display_label <- function(x, unique_analyses = TRUE) {
  base <- interactmeta_analysis_display(x)
  if (isTRUE(unique_analyses)) {
    return(base)
  }
  paste0(base, " \u2014 ", interactmeta_inference_display(x))
}

#' Approach labels for a list of fits
#' @noRd
interactmeta_compare_labels <- function(fits, labels = NULL) {
  if (!is.null(labels)) {
    labels <- as.character(labels)
    if (length(labels) != length(fits)) {
      stop("`labels` must have one entry per fit.", call. = FALSE)
    }
    return(labels)
  }
  analyses <- vapply(fits, interactmeta_fit_analysis, character(1L))
  unique_analyses <- length(unique(analyses)) == length(fits)
  vapply(
    fits,
    interactmeta_fit_display_label,
    character(1L),
    unique_analyses = unique_analyses
  )
}

#' @noRd
interactmeta_compact_fits <- function(...) {
  xs <- list(...)
  xs[!vapply(xs, is.null, logical(1L))]
}

#' CAMS fit (not DA / AD / SWADA / conventional)
#' @noRd
interactmeta_is_cams_fit <- function(x) {
  inherits(x, "interactmeta") && !interactmeta_is_non_cams_analysis(x)
}

#' @noRd
interactmeta_all_cams <- function(fits) {
  length(fits) >= 2L && all(vapply(fits, interactmeta_is_cams_fit, logical(1L)))
}

#' Wrap a comparison table so `plot()` has an S3 method
#' @noRd
interactmeta_compare_result <- function(tab,
                                        class_name,
                                        exponentiate = FALSE,
                                        xlab = NULL) {
  tab <- tibble::as_tibble(tab)
  class(tab) <- unique(c(class_name, "interactmeta_compare", class(tab)))
  attr(tab, "exponentiate") <- isTRUE(exponentiate)
  if (!is.null(xlab)) {
    attr(tab, "xlab") <- as.character(xlab)[[1L]]
  }
  tab
}

#' Harmonize tidy extractor rows onto the comparison table
#' @noRd
interactmeta_compare_tidy_rows <- function(approach, td) {
  td <- as.data.frame(td, stringsAsFactors = FALSE)
  if (!nrow(td)) {
    return(td)
  }
  se <- if ("se" %in% names(td)) {
    td$se
  } else if ("std.error" %in% names(td)) {
    td$std.error
  } else {
    NA_real_
  }
  lo <- if ("ci.lb" %in% names(td)) {
    td$ci.lb
  } else if ("conf.low" %in% names(td)) {
    td$conf.low
  } else {
    NA_real_
  }
  hi <- if ("ci.ub" %in% names(td)) {
    td$ci.ub
  } else if ("conf.high" %in% names(td)) {
    td$conf.high
  } else {
    NA_real_
  }
  keep <- setdiff(
    names(td),
    c("estimate", "se", "std.error", "ci.lb", "ci.ub", "conf.low", "conf.high",
      "mean", "sd")
  )
  out <- data.frame(
    approach = as.character(approach)[[1L]],
    stringsAsFactors = FALSE
  )[rep(1L, nrow(td)), , drop = FALSE]
  if (length(keep)) {
    out <- cbind(out, td[, keep, drop = FALSE])
  }
  out$estimate <- as.numeric(td$estimate)
  out$se <- as.numeric(se)
  out$ci.lb <- as.numeric(lo)
  out$ci.ub <- as.numeric(hi)
  rownames(out) <- NULL
  out
}

#' Stack extractor summaries from two or more fits
#' @noRd
interactmeta_compare_extractor <- function(fits,
                                           builder,
                                           class_name,
                                           caller,
                                           exponentiate = FALSE,
                                           xlab = NULL,
                                           id_cols = character()) {
  if (length(fits) < 2L) {
    stop(caller, "() needs at least two fits.", call. = FALSE)
  }
  labs <- interactmeta_compare_labels(fits)
  interactmeta_warn_incomparable_fits(fits, labs, caller)
  rows <- Map(builder, fits, labs)
  out <- dplyr::bind_rows(rows)
  keep <- intersect(c("approach", id_cols, "estimate", "se", "ci.lb", "ci.ub"), names(out))
  out <- out[, keep, drop = FALSE]
  out$approach <- factor(out$approach, levels = unique(labs))
  interactmeta_compare_result(
    out,
    class_name,
    exponentiate = exponentiate,
    xlab = xlab
  )
}

#' CAMS-only interaction comparison
#' @noRd
interactmeta_compare_interactions_fits <- function(fits,
                                                   exponentiate = FALSE,
                                                   type = "Overall") {
  type <- as.character(type)[[1L]]
  xlab <- if (isTRUE(exponentiate)) "Treatment effect" else "Log treatment effect"
  interactmeta_compare_extractor(
    fits,
    builder = function(fit, lab) {
      ix <- interactions(fit, type = type, exponentiate = exponentiate, summary = TRUE)
      td <- tidy(ix)
      if ("type" %in% names(td)) {
        td <- td[td$type == type, , drop = FALSE]
      }
      if ("interaction" %in% names(td) && !"contrast" %in% names(td)) {
        td$contrast <- td$interaction
      }
      interactmeta_compare_tidy_rows(lab, td)
    },
    class_name = "interactmeta_compare_interactions",
    caller = "compare_interactions",
    exponentiate = exponentiate,
    xlab = xlab,
    id_cols = "contrast"
  )
}

#' CAMS-only subgroup comparison
#' @noRd
interactmeta_compare_subgroups_fits <- function(fits,
                                                exponentiate = TRUE,
                                                type = "Overall") {
  type <- as.character(type)[[1L]]
  xlab <- if (isTRUE(exponentiate)) "Treatment effect" else "Log treatment effect"
  interactmeta_compare_extractor(
    fits,
    builder = function(fit, lab) {
      sg <- subgroups(fit, type = type, exponentiate = exponentiate, summary = TRUE)
      td <- tidy(sg)
      if ("type" %in% names(td)) {
        td <- td[td$type == type, , drop = FALSE]
      }
      interactmeta_compare_tidy_rows(lab, td)
    },
    class_name = "interactmeta_compare_subgroups",
    caller = "compare_subgroups",
    exponentiate = exponentiate,
    xlab = xlab,
    id_cols = "subgroup"
  )
}

#' Side-by-side combined-population effects
#'
#' [combined()] is a CAMS estimand. Pass two or more [cam()] fits (for
#' example Bayesian vs frequentist CAMS). Approach labels keep the analysis
#' family and, when the same family appears more than once, the approach
#' (`CAMS — bayesmeta` vs `CAMS — metafor (REML)`).
#'
#' @param ... Two or more [cam()] fits.
#' @param exponentiate If `TRUE` (default), report on the exponentiated
#'   (ratio) scale.
#' @param type CAMS layer (`"Overall"` default).
#' @return A tibble with class `interactmeta_compare_combined`. Print the
#'   table; [plot()] draws a forest of the compared estimates.
#' @seealso [compare_interactions()], [compare_subgroups()], [combined()]
#' @export
compare_combined <- function(..., exponentiate = TRUE, type = "Overall") {
  fits <- interactmeta_compact_fits(...)
  if (!interactmeta_all_cams(fits)) {
    stop(
      "`compare_combined()` compares two or more CAMS fits. ",
      "combined() is not defined for DA / AD / SWADA.",
      call. = FALSE
    )
  }
  type <- as.character(type)[[1L]]
  xlab <- if (isTRUE(exponentiate)) "Treatment effect" else "Log treatment effect"
  interactmeta_compare_extractor(
    fits,
    builder = function(fit, lab) {
      ov <- combined(fit, type = type, exponentiate = exponentiate, summary = TRUE)
      td <- tidy(ov)
      if ("type" %in% names(td)) {
        td <- td[td$type == type, , drop = FALSE]
      }
      interactmeta_compare_tidy_rows(lab, td)
    },
    class_name = "interactmeta_compare_combined",
    caller = "compare_combined",
    exponentiate = exponentiate,
    xlab = xlab
  )
}

#' Side-by-side heterogeneity summaries
#'
#' Stacks [heterogeneity()] summaries from two or more fits. When two CAMS
#' fits share the analysis family, labels include the inferential approach.
#'
#' @param ... Two or more fitted objects that support [heterogeneity()].
#' @param component Passed to [heterogeneity()] (`"all"` default).
#' @return A tibble with class `interactmeta_compare_heterogeneity`.
#' @seealso [heterogeneity()], [compare_interactions()]
#' @export
compare_heterogeneity <- function(..., component = "all") {
  fits <- interactmeta_compact_fits(...)
  if (length(fits) < 2L) {
    stop("`compare_heterogeneity()` needs at least two fits.", call. = FALSE)
  }
  interactmeta_compare_extractor(
    fits,
    builder = function(fit, lab) {
      het <- heterogeneity(fit, component = component, summary = TRUE)
      td <- tidy(het)
      interactmeta_compare_tidy_rows(lab, td)
    },
    class_name = "interactmeta_compare_heterogeneity",
    caller = "compare_heterogeneity",
    exponentiate = FALSE,
    xlab = "Between-trial SD",
    id_cols = "parameter"
  )
}

#' Side-by-side pooled information fractions
#'
#' Compares the pooled [ifrac()] from two or more CAMS fits. Trial-level
#' \eqn{\pi_j} is a design quantity and is the same across approaches; the
#' compared object is the **pooled** information fraction (integrated over
#' heterogeneity for Bayesian CAMS, plug-in at \eqn{\hat\tau} for frequentist
#' CAMS).
#'
#' @param ... Two or more [cam()] fits.
#' @param ifrac `"interaction"` (default), `"overall"`, or both.
#' @return A tibble with class `interactmeta_compare_ifrac`.
#' @seealso [ifrac()], [compare_weights()]
#' @export
compare_ifrac <- function(..., ifrac = "interaction") {
  dots <- list(...)
  nms <- names(dots)
  if (length(nms)) {
    hit <- nms %in% c("ifrac", "weight")
    if (any(hit)) {
      ifrac <- dots[[which(hit)[[1L]]]]
      dots <- dots[!hit]
    }
  }
  fits <- dots[!vapply(dots, is.null, logical(1L))]
  if (!interactmeta_all_cams(fits)) {
    stop(
      "`compare_ifrac()` compares pooled information fractions from two or more ",
      "CAMS fits. DA / AD / SWADA do not have CAMS information fractions; ",
      "use compare_weights() for those constructions.",
      call. = FALSE
    )
  }
  ifrac <- match.arg(ifrac, c("interaction", "overall"), several.ok = TRUE)
  interactmeta_compare_extractor(
    fits,
    builder = function(fit, lab) {
      rows <- lapply(ifrac, function(w) {
        oi <- interactmeta_eval_ifrac(fit, ifrac = w, summary = FALSE)
        sm <- attr(oi, "summary", exact = TRUE)
        if (is.null(sm)) {
          sm <- overall_ifrac_summary_df(oi, prob = attr(oi, "prob", exact = TRUE) %||% 0.95)
        }
        td <- data.frame(
          ifrac = w,
          estimate = as.numeric(sm$estimate)[[1L]],
          se = NA_real_,
          ci.lb = as.numeric(sm$lci)[[1L]],
          ci.ub = as.numeric(sm$uci)[[1L]],
          stringsAsFactors = FALSE
        )
        if (isTRUE(interactmeta_ifrac_is_plugin(fit))) {
          td$ci.lb <- NA_real_
          td$ci.ub <- NA_real_
        }
        interactmeta_compare_tidy_rows(lab, td)
      })
      dplyr::bind_rows(rows)
    },
    class_name = "interactmeta_compare_ifrac",
    caller = "compare_ifrac",
    exponentiate = FALSE,
    xlab = "Overall information fraction",
    id_cols = "ifrac"
  )
}

#' Plot a comparison table as a forest of estimates
#'
#' @param x A comparison object from [compare_interactions()],
#'   [compare_subgroups()], [compare_combined()], [compare_heterogeneity()],
#'   [compare_ifrac()], or [compare_weights()].
#' @param y Ignored.
#' @param xlim Shared x-axis limits on the displayed scale.
#' @param ... Ignored.
#' @return A **ggplot2** object.
#' @export
#' @method plot interactmeta_compare
plot.interactmeta_compare <- function(x, y, ..., xlim = NULL) {
  if (inherits(x, "interactmeta_compare_weights")) {
    return(interactmeta_plot_compare_weights(x, xlim = xlim))
  }
  interactmeta_plot_compare_forest(x, xlim = xlim)
}

#' @rdname plot.interactmeta_compare
#' @export
#' @method plot interactmeta_compare_interactions
plot.interactmeta_compare_interactions <- function(x, y, ..., xlim = NULL) {
  plot.interactmeta_compare(x, xlim = xlim, ...)
}

#' @rdname plot.interactmeta_compare
#' @export
#' @method plot interactmeta_compare_subgroups
plot.interactmeta_compare_subgroups <- function(x, y, ..., xlim = NULL) {
  plot.interactmeta_compare(x, xlim = xlim, ...)
}

#' @rdname plot.interactmeta_compare
#' @export
#' @method plot interactmeta_compare_combined
plot.interactmeta_compare_combined <- function(x, y, ..., xlim = NULL) {
  plot.interactmeta_compare(x, xlim = xlim, ...)
}

#' @rdname plot.interactmeta_compare
#' @export
#' @method plot interactmeta_compare_heterogeneity
plot.interactmeta_compare_heterogeneity <- function(x, y, ..., xlim = NULL) {
  plot.interactmeta_compare(x, xlim = xlim, ...)
}

#' @rdname plot.interactmeta_compare
#' @export
#' @method plot interactmeta_compare_ifrac
plot.interactmeta_compare_ifrac <- function(x, y, ..., xlim = NULL) {
  plot.interactmeta_compare(x, xlim = xlim, ...)
}

#' @rdname plot.interactmeta_compare
#' @export
#' @method plot interactmeta_compare_weights
plot.interactmeta_compare_weights <- function(x, y, ..., xlim = NULL) {
  interactmeta_plot_compare_weights(x, xlim = xlim)
}

#' Forest of compared estimates / intervals
#' @noRd
interactmeta_plot_compare_forest <- function(x, xlim = NULL) {
  tab <- as.data.frame(x, stringsAsFactors = FALSE)
  exponentiate <- isTRUE(attr(x, "exponentiate"))
  xlab <- attr(x, "xlab", exact = TRUE)
  if (is.null(xlab) || !nzchar(xlab)) {
    xlab <- if (isTRUE(exponentiate)) "Treatment effect" else "Estimate"
  }
  approach_levels <- if (is.factor(x$approach)) {
    levels(x$approach)
  } else {
    unique(as.character(tab$approach))
  }
  tab$approach <- factor(as.character(tab$approach), levels = rev(approach_levels))
  facet_col <- interactmeta_compare_facet_col(tab)
  pal <- interactmeta_approach_colours(approaches = as.character(x$approach))
  const <- interactmeta_plot_constants()
  vline <- if (inherits(x, "interactmeta_compare_ifrac") ||
                inherits(x, "interactmeta_compare_heterogeneity")) {
    NULL
  } else if (isTRUE(exponentiate)) {
    1
  } else {
    0
  }
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$estimate, y = .data$approach, colour = .data$approach))
  if (!is.null(vline)) {
    p <- p + ggplot2::geom_vline(
      xintercept = vline,
      colour = const$null_line,
      linewidth = 0.3
    )
  }
  has_int <- is.finite(tab$ci.lb) & is.finite(tab$ci.ub)
  if (any(has_int)) {
    p <- p + ggplot2::geom_linerange(
      data = tab[has_int, , drop = FALSE],
      ggplot2::aes(xmin = .data$ci.lb, xmax = .data$ci.ub),
      linewidth = const$line_hair
    )
  }
  p <- p +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = pal, breaks = names(pal), guide = "none") +
    ggplot2::labs(x = xlab, y = NULL, colour = NULL) +
    interactmeta_theme_report()
  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste(facet_col, "~ .")),
      scales = "free_y",
      space = "free_y"
    )
  }
  if (isTRUE(exponentiate)) {
    p <- p + ggplot2::scale_x_log10(limits = xlim)
  } else if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  p
}

#' @noRd
interactmeta_compare_facet_col <- function(tab) {
  candidates <- c("subgroup", "parameter", "ifrac", "weight", "contrast")
  for (col in candidates) {
    if (col %in% names(tab) && length(unique(as.character(tab[[col]]))) > 1L) {
      return(col)
    }
  }
  NULL
}

#' Grouped trial-weight comparison
#' @noRd
interactmeta_plot_compare_weights <- function(x, xlim = NULL) {
  tab <- as.data.frame(x, stringsAsFactors = FALSE)
  long <- tidyr::pivot_longer(
    tab,
    cols = -dplyr::all_of("trial"),
    names_to = "approach",
    values_to = "weight"
  )
  long$approach <- factor(long$approach, levels = setdiff(names(tab), "trial"))
  pal <- interactmeta_approach_colours(approaches = as.character(long$approach))
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$weight, y = .data$trial, colour = .data$approach)
  ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = pal, breaks = names(pal)) +
    ggplot2::labs(x = "Weight", y = NULL, colour = NULL) +
    interactmeta_theme_report()
  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  p
}
