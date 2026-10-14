#' CAMS model fit (`interactmeta`)
#'
#' A [cam()] fit with class `c("interactmeta", "cam")`. Fits carry an
#' `"interactmeta"` attribute (plot settings and [ifrac()] data for forest
#' plots).
#'
#' Printing and summarizing follow the style of **bayesmeta**:
#' [print.interactmeta()] is the compact glance (trials, inferential method
#' and computational backend, priors, parameter matrix).
#' [summary.interactmeta()] is the longer report
#' (heterogeneity, overall IF, combined-population / subgroup /
#' interaction). [plot.interactmeta()] draws [plot_contribution()].
#'
#' @param x An `interactmeta` from [cam()].
#' @param y Ignored.
#' @param ... Passed to [plot_contribution()].
#' @name interactmeta-fit
#' @seealso [cam()], [print.interactmeta()], [summary.interactmeta()],
#'   [plot_contribution()]
NULL

#' @export
#' @rdname interactmeta-fit
plot.interactmeta <- function(x, y, ...) {
  plot_contribution(x, ...)
}

#' Footnote name for the interval actually used.
#'
#' Bayesian CAMS names `"shortest"` or `"central"`. Other Bayesian engines
#' keep equal-tailed draw quantiles.
#' @noRd
interactmeta_credible_interval_label <- function(x) {
  eng <- if (inherits(x, "interactmeta") && !is.null(x$engine)) {
    as.character(x$engine)[[1L]]
  } else {
    interactmeta_origin_engine(x)
  }
  if (!identical(eng, "bayesmeta")) {
    return("central, equal-tailed credible intervals")
  }
  it <- if (inherits(x, "interactmeta")) x$interval.type else NULL
  if (is.null(it) && is.list(x) && inherits(x$fit, "interactmeta")) {
    it <- x$fit$interval.type
  }
  it <- interactmeta_match_interval_type(it)
  if (identical(it, "central")) {
    "central credible intervals"
  } else {
    "shortest credible intervals"
  }
}

#' Print a CAMS fit (bayesmeta-style)
#'
#' Compact overview of a [cam()] fit: trials / subgroups, inferential method
#' and computational backend, priors (when Bayesian), and the marginal
#' parameter matrix. Derived CAMS quantities (heterogeneity, overall IF,
#' combined-population / subgroup / interaction effects) are in
#' [summary.interactmeta()].
#'
#' @param x An `interactmeta` from [cam()].
#' @param digits Decimal places for the summary matrix (default `2`).
#' @param prob Coverage for the fit's interval type (default `0.95`).
#'   Bayesian CAMS uses `"shortest"` unless `interval.type = "central"`
#'   was passed to the fit.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @seealso [summary.interactmeta()], [cam()]
#' @export
print.interactmeta <- function(x, digits = 2, prob = 0.95, ...) {
  check_comparison_fit(x)
  meta <- comparison_fit_print_meta(x)
  cat(" 'interactmeta' object.\n")
  comparison_fit_print_header(meta)
  comparison_fit_print_priors(x, meta)
  mat <- comparison_fit_parameters_matrix(x, meta = meta, prob = prob)
  comparison_fit_print_parameters(mat, meta = meta, digits = digits)
  if (!isTRUE(meta$mle) && !identical(meta$engine, "metafor")) {
    cat("\n(quoted intervals are ", interactmeta_credible_interval_label(x), ".)\n", sep = "")
  } else if (identical(meta$engine, "metafor")) {
    cat("\n(quoted intervals are Wald confidence intervals from the stored rma covariance.)\n")
  }
  cat(
    "Use summary() for heterogeneity, overall IF, and derived effects.\n"
  )
  invisible(x)
}

#' Summarize a CAMS fit (bayesmeta-style)
#'
#' Longer companion to [print.interactmeta()]. Returns a `summary.interactmeta`
#' list and, when printed, adds heterogeneity, overall information fraction,
#' and marginalized combined-population / subgroup / interaction effects.
#'
#' @param object An `interactmeta` from [cam()].
#' @param digits Decimal places used by [print.summary.interactmeta()].
#' @param prob Coverage for intervals (default `0.95`).
#' @param exponentiate If `TRUE`, report marginalized effects on the
#'   multiplicative scale (e.g. OR / IRR).
#' @param ... Ignored.
#' @return An object of class `summary.interactmeta`.
#' @seealso [print.interactmeta()], [cam()]
#' @export
summary.interactmeta <- function(object,
                               digits = 2,
                               prob = 0.95,
                               exponentiate = FALSE,
                               ...) {
  check_comparison_fit(object)
  meta <- comparison_fit_print_meta(object)
  param <- comparison_fit_parameters_matrix(object, meta = meta, prob = prob)
  het <- tryCatch(
    heterogeneity(object, component = "all", prob = prob, summary = TRUE),
    error = function(e) NULL
  )
  ov <- tryCatch(
    combined(object, summary = TRUE, prob = prob, exponentiate = exponentiate),
    error = function(e) NULL
  )
  sg <- tryCatch(
    subgroups(object, summary = TRUE, prob = prob, exponentiate = exponentiate),
    error = function(e) NULL
  )
  int <- tryCatch(
    interactions(object, summary = TRUE, prob = prob, exponentiate = exponentiate),
    error = function(e) NULL
  )
  oi <- tryCatch(
    ifrac(object, summary = FALSE, prob = prob),
    error = function(e) NULL
  )
  out <- list(
    fit = object,
    meta = meta,
    data = comparison_fit_summary_data(object),
    parameters = param,
    heterogeneity = het,
    overall = ov,
    combined = ov,
    subgroups = sg,
    interaction = int,
    overall_ifrac = oi,
    digits = as.integer(digits)[[1L]],
    prob = prob,
    exponentiate = isTRUE(exponentiate),
    call = match.call(expand.dots = FALSE)
  )
  class(out) <- "summary.interactmeta"
  out
}

#' @rdname summary.interactmeta
#' @param x A `summary.interactmeta` object.
#' @export
print.summary.interactmeta <- function(x, digits = NULL, ...) {
  if (is.null(digits)) {
    digits <- x$digits %||% 2L
  }
  digits <- as.integer(digits)[[1L]]
  prob <- x$prob %||% 0.95
  meta <- x$meta
  cat("Summary of 'interactmeta' object.\n")
  comparison_fit_print_header(meta)

  comparison_fit_print_priors(x$fit, meta)
  comparison_fit_print_parameters(x$parameters, meta = meta, digits = digits)
  if (identical(meta$engine, "bayesmeta") && !isTRUE(meta$mle)) {
    cat(
      "\n(bayesmeta approach: overall via bmr, interaction via bayesmeta; ",
      "derived Ecological = theta - gamma from bayesmeta::convolve(). ",
      "quoted intervals are ",
      interactmeta_credible_interval_label(x$fit),
      ".)\n",
      sep = ""
    )
  } else if (identical(meta$engine, "metafor")) {
    cat(
      "\n(quoted intervals are Wald confidence intervals from the stored rma covariance.)\n"
    )
  } else if (!isTRUE(meta$mle)) {
    cat(
      "\n(quoted intervals are central, equal-tailed credible intervals.)\n"
    )
  }

  if (!is.null(x$heterogeneity) && !is.null(x$heterogeneity$summary)) {
    cat("\n")
    comparison_print_effects_table(
      x$heterogeneity$summary,
      digits = digits,
      prob = prob,
      id_cols = c("parameter", "type"),
      title = "heterogeneity:",
      footnote = FALSE,
      style = interactmeta_interval_style(x$fit)
    )
  }

  if (!is.null(x$overall_ifrac)) {
    cat("\n")
    if (inherits(x$overall_ifrac, "interactmeta_overall_ifrac")) {
      print(x$overall_ifrac, digits = digits)
    } else {
      cat("overall information fraction:\n")
      print(x$overall_ifrac, digits = digits)
    }
    if (isTRUE(meta$mle) || identical(meta$engine, "metafor")) {
      cat("(plug-in at fitted heterogeneity)\n")
    }
  }

  style <- interactmeta_interval_style(x$fit)
  scale_note <- if (isTRUE(x$exponentiate)) {
    " (exponentiated)"
  } else {
    ""
  }
  if (!is.null(x$overall) && !is.null(x$overall$summary)) {
    cat("\n")
    comparison_print_effects_table(
      x$overall$summary,
      digits = digits,
      prob = prob,
      id_cols = c("trial", "type"),
      title = paste0("combined-population effect", scale_note, ":"),
      footnote = FALSE,
      style = style
    )
  }
  if (!is.null(x$subgroups) && !is.null(x$subgroups$summary)) {
    cat("\n")
    comparison_print_effects_table(
      x$subgroups$summary,
      digits = digits,
      prob = prob,
      id_cols = c("trial", "subgroup", "type"),
      title = paste0(
        interactmeta_modifier_heading(
          attr(x$fit, "interactmeta", exact = TRUE)$modifier %||%
            attr(x$fit, "interactmeta", exact = TRUE)$subgroup_var
        ),
        " effects",
        scale_note,
        ":"
      ),
      footnote = FALSE,
      style = style
    )
  }
  if (!is.null(x$interaction) && !is.null(x$interaction$summary)) {
    cat("\n")
    int_sum <- x$interaction$summary
    if (isTRUE(x$interaction$shrinkage)) {
      int_sum <- comparison_interaction_shrinkage_summary_once(int_sum)
    }
    comparison_print_effects_table(
      int_sum,
      digits = digits,
      prob = prob,
      id_cols = c("trial", "interaction", "type"),
      title = paste0("interaction effect", scale_note, ":"),
      footnote = FALSE,
      style = style
    )
  }

  if (identical(style, "frequentist")) {
    cat("\n(CI = confidence interval; Pred = prediction interval where shown)\n")
  } else {
    cat("\n(CrI = credible interval; Pred = prediction interval where shown)\n")
  }
  invisible(x)
}

# ---- biostatistical print tables --------------------------------------------

#' Format a numeric CrI as "(lci, uci)"
#' @noRd
comparison_format_cri <- function(lci, uci, digits = 2) {
  paste0(
    "(",
    format(round(as.numeric(lci), digits), nsmall = digits, trim = TRUE),
    ", ",
    format(round(as.numeric(uci), digits), nsmall = digits, trim = TRUE),
    ")"
  )
}

#' Fold Prediction rows into a `95% Pred` column on Overall rows
#'
#' Keeps Overall median / mean / sd / CrI; attaches Prediction interval as an
#' extra column so Prediction is not printed as its own line.
#' @noRd
comparison_collapse_pooled_prediction <- function(summary, prob = 0.95) {
  if (is.null(summary) || !nrow(summary) || !"type" %in% names(summary)) {
    return(list(summary = summary, has_prediction = FALSE))
  }
  df <- as.data.frame(summary, stringsAsFactors = FALSE)
  types <- unique(as.character(df$type))
  if (!all(c("Overall", "Prediction") %in% types)) {
    return(list(summary = summary, has_prediction = FALSE))
  }

  # Normalize interval column names
  if ("conf.low" %in% names(df) && !"lci" %in% names(df)) {
    names(df)[names(df) == "conf.low"] <- "lci"
  }
  if ("conf.high" %in% names(df) && !"uci" %in% names(df)) {
    names(df)[names(df) == "conf.high"] <- "uci"
  }

  id_cols <- setdiff(
    intersect(
      c("trial", "subgroup", "interaction", "parameter", "label"),
      names(df)
    ),
    "type"
  )
  pooled <- df[as.character(df$type) == "Overall", , drop = FALSE]
  pred <- df[as.character(df$type) == "Prediction", , drop = FALSE]
  if (!nrow(pooled) || !nrow(pred)) {
    return(list(summary = summary, has_prediction = FALSE))
  }

  if (length(id_cols)) {
    key_p <- do.call(paste, c(pooled[id_cols], sep = "\r"))
    key_r <- do.call(paste, c(pred[id_cols], sep = "\r"))
    ord <- match(key_p, key_r)
    if (anyNA(ord)) {
      return(list(summary = summary, has_prediction = FALSE))
    }
    pred <- pred[ord, , drop = FALSE]
  } else if (nrow(pooled) != 1L || nrow(pred) != 1L) {
    return(list(summary = summary, has_prediction = FALSE))
  }

  out <- pooled
  out$type <- NULL
  if (!all(c("lci", "uci") %in% names(pred))) {
    return(list(summary = summary, has_prediction = FALSE))
  }
  out$pred_lci <- as.numeric(pred$lci)
  out$pred_uci <- as.numeric(pred$uci)
  list(summary = out, has_prediction = TRUE)
}

#' Build a plain data.frame for biostatistical console printing
#'
#' Avoids tibble metadata. Bayesian columns: median, mean, sd, and a
#' formatted credible-interval string. Frequentist columns: estimate, se,
#' and a formatted confidence-interval string. When Overall and Prediction are
#' both present, Prediction is folded into a `95% Pred` column.
#' @noRd
comparison_effects_print_df <- function(summary,
                                  digits = 2,
                                  prob = 0.95,
                                  id_cols = NULL,
                                  collapse_prediction = TRUE,
                                  style = c("bayesian", "frequentist")) {
  style <- match.arg(style)
  if (is.null(summary) || !nrow(summary)) {
    return(NULL)
  }
  has_prediction <- FALSE
  if (isTRUE(collapse_prediction)) {
    collapsed <- comparison_collapse_pooled_prediction(summary, prob = prob)
    summary <- collapsed$summary
    has_prediction <- isTRUE(collapsed$has_prediction)
  }
  df <- as.data.frame(summary, stringsAsFactors = FALSE)

  if (identical(style, "bayesian")) {
    if ("estimate" %in% names(df) && !"median" %in% names(df)) {
      names(df)[names(df) == "estimate"] <- "median"
    }
    if ("std.error" %in% names(df) && !"sd" %in% names(df)) {
      names(df)[names(df) == "std.error"] <- "sd"
    }
  } else {
    if ("median" %in% names(df) && !"estimate" %in% names(df)) {
      names(df)[names(df) == "median"] <- "estimate"
    }
    if ("sd" %in% names(df) && !"se" %in% names(df)) {
      names(df)[names(df) == "sd"] <- "se"
    }
    if ("std.error" %in% names(df) && !"se" %in% names(df)) {
      names(df)[names(df) == "std.error"] <- "se"
    }
  }
  if ("conf.low" %in% names(df) && !"lci" %in% names(df)) {
    names(df)[names(df) == "conf.low"] <- "lci"
  }
  if ("conf.high" %in% names(df) && !"uci" %in% names(df)) {
    names(df)[names(df) == "conf.high"] <- "uci"
  }

  default_ids <- c(
    "parameter", "label", "trial", "subgroup", "interaction", "type"
  )
  if (is.null(id_cols)) {
    id_cols <- intersect(default_ids, names(df))
  } else {
    id_cols <- intersect(as.character(id_cols), names(df))
  }
  # type is dropped after Overall/Prediction collapse
  if (isTRUE(has_prediction)) {
    id_cols <- setdiff(id_cols, "type")
  }

  # Prefer compact parameter over long prose label when both requested
  if (all(c("label", "parameter") %in% id_cols)) {
    id_cols <- setdiff(id_cols, "label")
  }

  if (length(id_cols)) {
    out_ids <- df[id_cols]
    keep_id <- vapply(out_ids, function(col) {
      any(nzchar(as.character(col)) & !is.na(col))
    }, logical(1))
    out_ids <- out_ids[keep_id]
  } else {
    out_ids <- data.frame(row.names = seq_len(nrow(df)))
  }

  interval_lab <- if (identical(style, "frequentist")) {
    if (isTRUE(all.equal(as.numeric(prob), 0.95))) "95% CI" else sprintf("%g%% CI", 100 * as.numeric(prob))
  } else {
    if (isTRUE(all.equal(as.numeric(prob), 0.95))) "95% CrI" else sprintf("%g%% CrI", 100 * as.numeric(prob))
  }
  pred_lab <- if (isTRUE(all.equal(as.numeric(prob), 0.95))) {
    "95% Pred"
  } else {
    sprintf("%g%% Pred", 100 * as.numeric(prob))
  }

  out <- out_ids
  fmt_num <- function(v) {
    v <- as.numeric(v)
    out <- format(round(v, digits), nsmall = digits, trim = TRUE)
    out[!is.finite(v)] <- ""
    out
  }
  if (identical(style, "frequentist")) {
    if ("estimate" %in% names(df)) {
      out[["estimate"]] <- fmt_num(df$estimate)
    }
    if ("se" %in% names(df) && any(is.finite(as.numeric(df$se)))) {
      out[["se"]] <- fmt_num(df$se)
    }
  } else {
    if ("median" %in% names(df)) {
      out[["median"]] <- fmt_num(df$median)
    }
    if ("mean" %in% names(df)) {
      out[["mean"]] <- fmt_num(df$mean)
    }
    if ("sd" %in% names(df) && any(is.finite(as.numeric(df$sd)))) {
      out[["sd"]] <- fmt_num(df$sd)
    }
  }
  has_interval <- all(c("lci", "uci") %in% names(df)) &&
    any(is.finite(as.numeric(df$lci)) & is.finite(as.numeric(df$uci)))
  if (has_interval) {
    cri <- comparison_format_cri(df$lci, df$uci, digits = digits)
    cri[!is.finite(as.numeric(df$lci)) | !is.finite(as.numeric(df$uci))] <- ""
    out[[interval_lab]] <- cri
  }
  if (isTRUE(has_prediction) && all(c("pred_lci", "pred_uci") %in% names(df))) {
    out[[pred_lab]] <- comparison_format_cri(df$pred_lci, df$pred_uci, digits = digits)
  }
  if (!ncol(out)) {
    return(NULL)
  }
  attr(out, "has_prediction") <- has_prediction
  attr(out, "has_interval") <- has_interval
  attr(out, "style") <- style
  out
}

#' Print a CAMS effects / heterogeneity summary without tibbles
#' @noRd
comparison_print_effects_table <- function(summary,
                                     digits = 2,
                                     prob = 0.95,
                                     id_cols = NULL,
                                     title = NULL,
                                     footnote = TRUE,
                                     collapse_prediction = TRUE,
                                     style = c("bayesian", "frequentist")) {
  style <- match.arg(style)
  tab <- comparison_effects_print_df(
    summary,
    digits = digits,
    prob = prob,
    id_cols = id_cols,
    collapse_prediction = collapse_prediction,
    style = style
  )
  if (!is.null(title) && nzchar(title)) {
    cat(title, "\n", sep = "")
  }
  if (is.null(tab) || !nrow(tab)) {
    cat("(no summary)\n")
    return(invisible(NULL))
  }
  comparison_print_aligned_df(tab)
  if (isTRUE(footnote)) {
    has_pred <- isTRUE(attr(tab, "has_prediction"))
    has_int <- isTRUE(attr(tab, "has_interval"))
    if (identical(style, "frequentist")) {
      if (has_pred) {
        cat("(CI = confidence interval; Pred = prediction interval)\n")
      } else if (has_int) {
        cat("(Wald confidence intervals)\n")
      }
    } else if (has_pred) {
      cat("(CrI = credible interval; Pred = prediction interval)\n")
    } else if (has_int) {
      it <- attr(summary, "interval.type", exact = TRUE)
      if (identical(it, "shortest")) {
        cat("(shortest credible intervals)\n")
      } else if (identical(it, "central")) {
        cat("(central credible intervals)\n")
      } else {
        cat("(central equal-tailed credible intervals)\n")
      }
    }
  }
  invisible(NULL)
}

#' Console-aligned print for a character/numeric data.frame (no tibble chrome)
#' @noRd
comparison_print_aligned_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  nms <- names(df)
  cols <- lapply(seq_along(df), function(j) {
    as.character(df[[j]])
  })
  widths <- vapply(seq_along(cols), function(j) {
    max(nchar(c(nms[[j]], cols[[j]]), type = "width"), na.rm = TRUE)
  }, integer(1))
  # Left-align identifiers; right-align numeric-ish columns
  id_like <- nms %in% c(
    "parameter", "label", "trial", "subgroup", "interaction", "type"
  )
  pad <- function(x, w, right = TRUE) {
    if (isTRUE(right)) {
      format(x, width = w, justify = "right")
    } else {
      format(x, width = w, justify = "left")
    }
  }
  header <- vapply(seq_along(nms), function(j) {
    pad(nms[[j]], widths[[j]], right = !id_like[[j]])
  }, character(1))
  cat(paste(header, collapse = "  "), "\n", sep = "")
  for (i in seq_len(nrow(df))) {
    row <- vapply(seq_along(cols), function(j) {
      pad(cols[[j]][[i]], widths[[j]], right = !id_like[[j]])
    }, character(1))
    cat(paste(row, collapse = "  "), "\n", sep = "")
  }
  invisible(NULL)
}

#' Resolve digits / prob from an effects object or `...`
#' @noRd
comparison_print_digits_prob <- function(x, digits = NULL, prob = NULL, ...) {
  dots <- list(...)
  if (is.null(digits)) {
    digits <- dots$digits %||% 2L
  }
  if (is.null(prob)) {
    prob <- dots$prob %||% attr(x, "prob", exact = TRUE) %||% x$prob %||% 0.95
  }
  list(digits = as.integer(digits)[[1L]], prob = as.numeric(prob)[[1L]])
}

#' @noRd
comparison_fit_print_meta <- function(x) {
  st <- attr(x, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% "trial"
  subgroup_var <- st$subgroup_var %||% "subgroup"
  dat <- x$data
  if (is.null(dat) || !nrow(dat)) {
    dat <- st$data_fit
  }
  trials <- if (!is.null(dat) && trial_var %in% names(dat)) {
    unique(as.character(dat[[trial_var]]))
  } else {
    character()
  }
  subgroups <- as.character(x$design$subgroup_levels %||% character())
  if (!length(subgroups) && !is.null(dat) && subgroup_var %in% names(dat)) {
    subgroups <- unique(as.character(dat[[subgroup_var]]))
  }
  list(
    engine = x$engine %||% st$engine %||% NA_character_,
    method = x$method %||% st$method %||% NA_character_,
    algorithm = x$algorithm %||% NA_character_,
    mle = comparison_is_mle(x),
    ecological = st$ecological_slope %||% st$ecological %||% TRUE,
    ecological_slope = st$ecological_slope %||% st$ecological %||% TRUE,
    reference = st$reference,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    modifier = st$modifier %||% subgroup_var,
    trials = trials,
    subgroups = subgroups,
    n_trials = length(trials),
    n_subgroups = length(subgroups),
    n_draws = if (!is.null(x$draws)) nrow(x$draws) else NA_integer_,
    formula = x$formula
  )
}

#' @noRd
comparison_fit_print_header <- function(meta) {
  heading <- interactmeta_modifier_heading(meta$modifier %||% meta$subgroup_var)
  is_subgroup <- identical(tolower(heading), "subgroup")
  approach <- switch(
    tolower(as.character(meta$engine)[[1L]]),
    bayesmeta = "Bayesian (bayesmeta)",
    metafor = "Frequentist (metafor)",
    as.character(meta$engine)[[1L]]
  )
  count_heading <- if (is_subgroup) {
    if (identical(meta$n_subgroups, 1L)) "subgroup" else "subgroups"
  } else {
    paste0(heading, if (identical(meta$n_subgroups, 1L)) " level" else " levels")
  }
  levels_heading <- if (is_subgroup) "Subgroups" else paste0(heading, " levels")
  cat(paste0(
    "\n",
    meta$n_trials,
    " trials, ",
    meta$n_subgroups,
    " ",
    count_heading
  ))
  if (!is.null(meta$reference) && nzchar(as.character(meta$reference)[[1L]])) {
    cat(paste0(" (reference: ", as.character(meta$reference)[[1L]], ")"))
  }
  cat(":\n")
  labs <- meta$trials
  if (length(labs) > 10L) {
    cat(paste(c(labs[seq_len(10L)], "..."), collapse = ", "))
  } else if (length(labs)) {
    cat(paste(labs, collapse = ", "))
  } else {
    cat("(none)")
  }
  cat("\n")
  if (length(meta$subgroups)) {
    cat(paste0(levels_heading, ": ", paste(meta$subgroups, collapse = ", "), "\n"))
  }
  cat(paste0(
    "\napproach: ",
    approach,
    "\nmethod: ",
    meta$method
  ))
  if (!is.na(meta$algorithm) && nzchar(meta$algorithm)) {
    cat(paste0(" (", meta$algorithm, ")"))
  }
  if (isTRUE(meta$mle)) {
    cat("\n")
  } else if (identical(meta$engine, "metafor")) {
    cat("\n")
  } else if (is.finite(meta$n_draws)) {
    cat(paste0("\n", comparison_fit_draws_label(meta), ": ", meta$n_draws, "\n"))
  } else {
    cat("\n")
  }
}

#' Print label for `$draws` rows (quadrature nodes vs MCMC draws)
#' @noRd
comparison_fit_draws_label <- function(meta) {
  alg <- as.character(meta$algorithm %||% "")[[1L]]
  eng <- as.character(meta$engine %||% "")[[1L]]
  if (identical(alg, "quadrature") || identical(eng, "bayesmeta")) {
    "quadrature nodes"
  } else {
    "posterior draws"
  }
}

#' Print-facing name for an `ifrac2_*` design / RE term
#'
#' Maps `ifrac2_Older` → `Interaction_Older` for formula and heterogeneity
#' print labels. Internal draw / design columns stay `ifrac2_*`.
#' @noRd
comparison_interaction_print_name <- function(ifrac2_name) {
  sub("^ifrac2_", "Interaction_", as.character(ifrac2_name))
}

#' Print-facing name for a fixed-effect `b_ifrac_*` coefficient
#'
#' Maps `b_ifrac_Older` / `ifrac_Older` → `Ecological_Older` in the marginal
#' posterior summary (composition / IF slope). Internal draw columns stay
#' `b_ifrac_*`.
#' @noRd
comparison_ifrac_fe_print_name <- function(ifrac_coef) {
  nm <- as.character(ifrac_coef)[[1L]]
  nm <- sub("^b_", "", nm)
  sub("^ifrac_", "Ecological_", nm)
}

#' Print-facing name for a subgroup fixed-effect coefficient
#'
#' Maps `b_subgroupOlder` / `subgroupOlder` → `Interaction_Older`.
#' Internal draw columns stay `b_subgroup*`.
#' @noRd
comparison_subgroup_fe_print_name <- function(subgroup_coef) {
  nm <- as.character(subgroup_coef)[[1L]]
  sg <- sub("^b_subgroup", "", nm)
  sg <- sub("^subgroup", "", sg)
  paste0("Interaction_", sg)
}

#' Rewrite `ifrac2_*` tokens to `Interaction_*` in a formula for print
#' @noRd
comparison_formula_print_text <- function(formula) {
  if (is.null(formula)) {
    return(NULL)
  }
  txt <- if (inherits(formula, "formula") || is.language(formula)) {
    paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  } else {
    as.character(formula)[[1L]]
  }
  gsub("ifrac2_", "Interaction_", txt, fixed = TRUE)
}

#' Parameter matrix for print / summary (bayesmeta native table when available)
#' @noRd
comparison_fit_parameters_matrix <- function(x, meta = NULL, prob = 0.95) {
  if (is.null(meta)) {
    meta <- comparison_fit_print_meta(x)
  }
  if (identical(meta$engine, "metafor")) {
    interactmeta_metafor_parameters_df(x, prob = prob)
  } else if (identical(meta$engine, "bayesmeta")) {
    interactmeta_bayesmeta_param_summary(x, prob = prob) %||%
      comparison_fit_param_summary_matrix(x, prob = prob)
  } else {
    comparison_fit_param_summary_matrix(x, prob = prob)
  }
}

#' Frequentist CAMS parameter table from stored rma objects
#' @noRd
interactmeta_metafor_parameters_df <- function(fit, prob = 0.95) {
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  if (is.null(ov) || is.null(int)) {
    return(NULL)
  }
  z <- stats::qnorm(1 - (1 - as.numeric(prob)[[1L]]) / 2)
  row_est <- function(parameter, estimate, se, lci = NA_real_, uci = NA_real_) {
    estimate <- as.numeric(estimate)[[1L]]
    se <- as.numeric(se)[[1L]]
    if (!is.finite(lci) && is.finite(estimate) && is.finite(se)) {
      lci <- estimate - z * se
      uci <- estimate + z * se
    }
    data.frame(
      parameter = parameter,
      estimate = estimate,
      se = se,
      lci = as.numeric(lci)[[1L]],
      uci = as.numeric(uci)[[1L]],
      stringsAsFactors = FALSE
    )
  }
  tau_interval <- function(model) {
    ci <- tryCatch(
      suppressWarnings(stats::confint(model, level = prob)),
      error = function(e) NULL
    )
    random <- ci$random %||% NULL
    if (is.null(random) || is.null(dim(random)) ||
        !"tau" %in% rownames(random) ||
        !all(c("ci.lb", "ci.ub") %in% colnames(random))) {
      return(c(lci = NA_real_, uci = NA_real_))
    }
    c(
      lci = as.numeric(random["tau", "ci.lb"]),
      uci = as.numeric(random["tau", "ci.ub"])
    )
  }
  tau_row <- function(parameter, tau2, se_tau2, interval) {
    tau <- sqrt(max(as.numeric(tau2)[[1L]], 0))
    se <- NA_real_
    se_tau2 <- as.numeric(se_tau2)[[1L]]
    if (is.finite(se_tau2) && is.finite(tau) && tau > 0) {
      se <- se_tau2 / (2 * tau)
    }
    data.frame(
      parameter = parameter,
      estimate = tau,
      se = se,
      lci = unname(interval[["lci"]]),
      uci = unname(interval[["uci"]]),
      stringsAsFactors = FALSE
    )
  }
  alpha <- as.numeric(ov$beta[[1L]])
  theta <- as.numeric(ov$beta[[2L]])
  gamma <- as.numeric(int$beta[[1L]])
  eco <- theta - gamma
  se_eco <- sqrt(max(as.numeric(ov$se[[2L]])^2 + as.numeric(int$se[[1L]])^2, 0))
  sg <- grep("^b_subgroup", names(as.data.frame(fit$draws)), value = TRUE)
  int_lab <- if (length(sg)) {
    comparison_subgroup_fe_print_name(sg[[1L]])
  } else {
    "Interaction"
  }
  ifrac_b <- grep("^b_ifrac_", names(as.data.frame(fit$draws)), value = TRUE)
  eco_lab <- if (length(ifrac_b)) {
    comparison_ifrac_fe_print_name(ifrac_b[[1L]])
  } else {
    "Ecological"
  }
  rbind(
    tau_row("tau", ov$tau2, ov$se.tau2, tau_interval(ov)),
    tau_row("tau_gamma", int$tau2, int$se.tau2, tau_interval(int)),
    row_est("Intercept", alpha, ov$se[[1L]], ov$ci.lb[[1L]], ov$ci.ub[[1L]]),
    row_est("theta", theta, ov$se[[2L]], ov$ci.lb[[2L]], ov$ci.ub[[2L]]),
    row_est(eco_lab, eco, se_eco),
    row_est(int_lab, gamma, int$se[[1L]], int$ci.lb[[1L]], int$ci.ub[[1L]])
  )
}

#' Print the fitted-parameter matrix (shared by print and summary)
#' @noRd
comparison_fit_print_parameters <- function(mat, meta, digits) {
  if (identical(meta$engine, "metafor") && is.data.frame(mat)) {
    meth <- as.character(meta$method %||% "REML")[[1L]]
    cat("\n", meth, " estimates:\n", sep = "")
    comparison_print_effects_table(
      mat,
      digits = digits,
      prob = 0.95,
      id_cols = "parameter",
      title = NULL,
      footnote = FALSE,
      collapse_prediction = FALSE,
      style = "frequentist"
    )
    return(invisible(NULL))
  }
  if (is.null(mat) || !length(mat)) {
    return(invisible(NULL))
  }
  if (isTRUE(meta$mle)) {
    cat("\nML estimates:\n")
    print(round(mat["median", , drop = FALSE], digits = digits))
    cat("\n(reporting information fractions use the tau MLE as a plug-in.)\n")
  } else if (identical(meta$engine, "bayesmeta")) {
    cat("\nmarginal posterior summary (overall | interaction partition):\n")
    show <- round(mat, digits = digits)
    attr(show, "interval.type") <- NULL
    print(show)
  } else {
    cat("\nmarginal posterior summary:\n")
    print(round(mat, digits = digits))
  }
  invisible(NULL)
}

#' @noRd
comparison_fit_print_priors <- function(fit, meta) {
  if (identical(meta$engine, "metafor")) {
    cat("\npriors: none (frequentist; Bayesian prior arguments are ignored)\n")
    return(invisible(NULL))
  }
  if (isTRUE(meta$mle)) {
    cat("\npriors: soft priors off (likelihood + RE hierarchy; Stan optimize)\n")
    return(invisible(NULL))
  }
  if (inherits(fit$priors, "cam_priors")) {
    cat("\n")
    print(fit$priors)
    if (identical(meta$engine, "bayesmeta") &&
          identical(fit$priors$fixed$family, "normal")) {
      cat(
        "bayesmeta: independent Normal prior on (alpha, theta) with beta.prior.cov = diag(sd^2); ",
        "mu.prior for gamma separately. ",
        "Independent normals on (alpha, theta) are not independent normals on (alpha, beta, gamma), ",
        "because theta = beta + gamma.\n"
      )
    }
    return(invisible(NULL))
  }
  pr <- fit$prior
  if (is.null(pr) || !nrow(pr)) {
    pr <- tryCatch(default_comparison_prior("ifrac2", group = "trial"), error = function(e) NULL)
  }
  cat("\npriors:\n")
  if (isFALSE(meta$ecological_slope %||% meta$ecological)) {
    cat("  Intercept, b_subgroup ~ normal(0, 10)\n")
    cat("  b_ifrac omitted (ecological_slope = FALSE; identically 0)\n")
  } else {
    cat("  (alpha, theta) and gamma ~ independent normal(0, 10); ")
    cat("stored b_ifrac = theta - gamma\n")
  }
  if (!is.null(pr) && nrow(pr)) {
    for (i in seq_len(nrow(pr))) {
      cat(paste0("  ", pr$class[[i]], " ~ ", pr$prior[[i]], "\n"))
    }
  } else {
    cat("  sd_trial[*] ~ improper uniform (bayesmeta default)\n")
  }
  invisible(NULL)
}

#' Data matrix like bayesmeta summary (y, sigma by trial-subgroup)
#' @noRd
comparison_fit_summary_data <- function(x) {
  st <- attr(x, "interactmeta", exact = TRUE)
  dat <- x$data
  if (is.null(dat) || !nrow(dat)) {
    dat <- st$data_fit
  }
  if (is.null(dat) || !nrow(dat)) {
    return(NULL)
  }
  trial_var <- st$trial_var %||% "trial"
  subgroup_var <- st$subgroup_var %||% "subgroup"
  need <- c("yi", "vi")
  if (!all(need %in% names(dat))) {
    return(NULL)
  }
  lab <- if (all(c(trial_var, subgroup_var) %in% names(dat))) {
    paste(as.character(dat[[trial_var]]), as.character(dat[[subgroup_var]]), sep = " / ")
  } else {
    as.character(seq_len(nrow(dat)))
  }
  out <- cbind(
    y = as.numeric(dat$yi),
    sigma = sqrt(as.numeric(dat$vi))
  )
  rownames(out) <- make.unique(lab)
  out
}

#' Approximate posterior mode from a sample (density peak)
#' @noRd
comparison_sample_mode <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  if (length(unique(x)) <= 2L) {
    return(stats::median(x))
  }
  d <- tryCatch(stats::density(x, n = 512L), error = function(e) NULL)
  if (is.null(d)) {
    return(stats::median(x))
  }
  d$x[[which.max(d$y)]]
}

#' Bayesmeta-style summary matrix for one numeric draw vector
#' @noRd
comparison_draw_summary_vec <- function(x, prob = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  lower_lab <- if (isTRUE(all.equal(prob, 0.95))) {
    "95% lower"
  } else {
    sprintf("%g%% lower", 100 * prob)
  }
  upper_lab <- if (isTRUE(all.equal(prob, 0.95))) {
    "95% upper"
  } else {
    sprintf("%g%% upper", 100 * prob)
  }
  row_names <- c("mode", "median", "mean", "sd", lower_lab, upper_lab)
  empty <- stats::setNames(rep(NA_real_, 6L), row_names)
  if (!length(x)) {
    return(empty)
  }
  mom <- comparison_draw_moments(x, prob = prob)
  mode_v <- if (length(x) == 1L) {
    x[[1L]]
  } else {
    comparison_sample_mode(x)
  }
  stats::setNames(
    c(mode_v, mom$median, mom$mean, mom$sd, mom$lci, mom$uci),
    row_names
  )
}

#' Select key draw columns and build bayesmeta-style summary matrix
#' @noRd
comparison_fit_param_summary_matrix <- function(x, prob = 0.95) {
  post <- comparison_posterior_df(x)
  cols <- character()
  labels <- character()

  if ("sd_trial__Intercept" %in% names(post)) {
    cols <- c(cols, "sd_trial__Intercept")
    labels <- c(labels, "tau")
  }
  gamma_cols <- grep("^sd_trial__ifrac2_", names(post), value = TRUE)
  if (length(gamma_cols)) {
    for (gc in gamma_cols) {
      lab <- sub("^sd_trial__ifrac2_", "tau_gamma[", gc)
      lab <- paste0(lab, "]")
      if (length(gamma_cols) == 1L) {
        lab <- "tau_gamma"
      }
      cols <- c(cols, gc)
      labels <- c(labels, lab)
    }
  }
  if ("b_Intercept" %in% names(post)) {
    cols <- c(cols, "b_Intercept")
    labels <- c(labels, "Intercept")
  }
  ifrac_b <- grep("^b_ifrac_", names(post), value = TRUE)
  for (bc in ifrac_b) {
    cols <- c(cols, bc)
    labels <- c(labels, comparison_ifrac_fe_print_name(bc))
  }
  sg_b <- grep("^b_subgroup", names(post), value = TRUE)
  for (bc in sg_b) {
    cols <- c(cols, bc)
    labels <- c(labels, comparison_subgroup_fe_print_name(bc))
  }

  if (!length(cols)) {
    return(matrix(NA_real_, nrow = 6L, ncol = 0L))
  }

  mat <- vapply(
    cols,
    function(nm) comparison_draw_summary_vec(post[[nm]], prob = prob),
    numeric(6L)
  )
  colnames(mat) <- make.unique(labels)
  mat
}
