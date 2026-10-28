# ---- print_summary.R ----
#' CAMS model fit (`interactmeta`)
#'
#' A [cams()] fit with class `"interactmeta"`. Fits carry an
#' `"interactmeta"` attribute (plot settings and [ifrac()] data for forest
#' plots).
#'
#' Printing and summarizing follow the style of **bayesmeta**:
#' [print.interactmeta()] is the compact glance (trials, `method` as in
#' `cams()`, parameter matrix).
#' [summary.interactmeta()] is the longer report
#' (heterogeneity, overall IF, combined-population / subgroup /
#' interaction). [plot.interactmeta()] draws contribution, heterogeneity,
#' or forest plots via `type`.
#'
#' @param x An `interactmeta` from [cams()].
#' @param y Ignored.
#' @param type `"contribution"` (default), `"heterogeneity"`, or `"forest"`.
#' @param ... Passed to the selected plot method.
#' @name interactmeta-fit
#' @seealso [cams()], [print.interactmeta()], [summary.interactmeta()],
#'   [forestplot()]
NULL

#' @export
#' @rdname interactmeta-fit
plot.interactmeta <- function(x,
                              y,
                              type = c("contribution", "heterogeneity", "forest"),
                              ...) {
  type <- match.arg(type)
  switch(
    type,
    contribution = plot_contribution(x, ...),
    heterogeneity = plot_heterogeneity(x, ...),
    forest = forestplot(x, ...)
  )
}

#' Footnote name for the interval actually used.
#'
#' Bayesian CAMS names `"shortest"` or `"central"`. Other Bayesian engines
#' use central draw quantiles.
#' @noRd
interactmeta_credible_interval_label <- function(x) {
  rule <- interactmeta_credible_interval_rule(x)
  paste0(rule, " credible intervals")
}

#' One-line Bayesian interval footnote for [print.interactmeta()] / [summary.interactmeta()]
#' @noRd
interactmeta_bayesian_interval_footnote <- function(x) {
  fit <- if (inherits(x, "summary.interactmeta")) x$fit else x
  paste0(
    "(Quoted intervals are ",
    interactmeta_credible_interval_label(fit),
    "; overall via bmr, interaction via bayesmeta.)"
  )
}

#' One-word interval rule for footnotes (`central` or `shortest`).
#' @noRd
interactmeta_credible_interval_rule <- function(x) {
  if (!interactmeta_is_bayesian_fit(x)) {
    return("central")
  }
  it <- if (inherits(x, "interactmeta")) {
    x$interval.type
  } else if (is.list(x) && inherits(x$fit, "interactmeta")) {
    x$fit$interval.type
  } else {
    NULL
  }
  it <- interactmeta_match_interval_type(it)
  if (identical(it, "shortest")) "shortest" else "central"
}

#' Print a CAMS fit (bayesmeta-style)
#'
#' Compact overview of a [cams()] fit: trials / subgroups, the same `method`
#' label as in `cams(..., method = )`, and the marginal parameter matrix. Derived CAMS quantities (heterogeneity, overall IF,
#' combined-population / subgroup / interaction effects) are in
#' [summary.interactmeta()].
#'
#' @param x An `interactmeta` from [cams()].
#' @param digits Decimal places for the summary matrix (default `2`).
#' @param prob Coverage for the fit's interval type (default `0.95`).
#'   Bayesian CAMS uses `"central"` by default; pass `interval.type = "shortest"`
#'   to [cams()] for the alternative rule.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @seealso [summary.interactmeta()], [cams()]
#' @export
print.interactmeta <- function(x, digits = 2, prob = 0.95, ...) {
  check_comparison_fit(x)
  interactmeta_ensure_reporting(x, prob = prob, warm_marginals = FALSE)
  interactmeta_sync_reporting_from_sidecar(x)
  meta <- comparison_fit_print_meta(x)
  cat(" 'interactmeta' object.\n")
  comparison_fit_print_header(meta)
  rep <- interactmeta_reporting_sidecar(x) %||% x$reporting
  mat <- if (is.list(rep) && is.matrix(rep$parameters)) {
    rep$parameters
  } else {
    comparison_fit_parameters_matrix(x, meta = meta, prob = prob)
  }
  comparison_fit_print_parameters(mat, meta = meta, digits = digits)
  if (!isTRUE(meta$mle) && interactmeta_is_bayesian_fit(x)) {
    cat("\n", interactmeta_bayesian_interval_footnote(x), "\n", sep = "")
  } else if (!interactmeta_is_bayesian_fit(x)) {
    cat("\n(quoted intervals are Wald confidence intervals from the stored rma covariance.)\n")
  }
  cat(
    "Use summary() for heterogeneity, IF, and derived effects.\n"
  )
  invisible(x)
}

#' Summarize a CAMS fit (bayesmeta-style)
#'
#' Longer companion to [print.interactmeta()]. Returns a `summary.interactmeta`
#' list and, when printed, adds heterogeneity, overall information fraction,
#' and marginalized combined-population / subgroup / interaction effects.
#'
#' @param object An `interactmeta` from [cams()].
#' @param digits Decimal places used by [print.summary.interactmeta()].
#' @param prob Coverage for intervals (default `0.95`).
#' @param exponentiate If `TRUE`, report marginalized effects on the
#'   multiplicative scale (e.g. OR / IRR).
#' @param ... Ignored.
#' @return An object of class `summary.interactmeta`.
#' @seealso [print.interactmeta()], [cams()]
#' @export
summary.interactmeta <- function(object,
                               digits = 2,
                               prob = 0.95,
                               exponentiate = FALSE,
                               ...) {
  check_comparison_fit(object)
  interactmeta_ensure_reporting(object, prob = prob, warm_marginals = TRUE)
  interactmeta_sync_reporting_from_sidecar(object)
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
    modifiers(object, summary = TRUE, prob = prob, exponentiate = exponentiate),
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
  fit <- x$fit
  cat("Summary of 'interactmeta' object.\n")
  comparison_fit_print_header(meta)
  comparison_fit_print_parameters(x$parameters, meta = meta, digits = digits)
  bayes <- interactmeta_is_bayesian_fit(fit)
  if (bayes && !isTRUE(meta$mle)) {
    cat("\n", interactmeta_bayesian_interval_footnote(fit), "\n", sep = "")
  } else if (!bayes) {
    cat(
      "\n(quoted intervals are Wald confidence intervals from the stored rma covariance.)\n"
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
      style = interactmeta_interval_style(x$fit)
    )
  }

  if (!is.null(x$overall_ifrac)) {
    cat("\n")
    print(x$overall_ifrac, digits = digits)
  }

  style <- interactmeta_interval_style(x$fit)
  if (!is.null(x$overall) && !is.null(x$overall$summary)) {
    cat("\n")
    comparison_print_effects_table(
      x$overall$summary,
      digits = digits,
      prob = prob,
      id_cols = if (isTRUE(x$overall$shrinkage)) {
        c("trial", "type")
      } else {
        "combined"
      },
      title = NULL,
      style = style
    )
  }
  if (!is.null(x$subgroups) && !is.null(x$subgroups$summary)) {
    cat("\n")
    comparison_print_effects_table(
      x$subgroups$summary,
      digits = digits,
      prob = prob,
      id_cols = c("trial", interactmeta_level_column(), "type"),
      title = NULL,
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
      title = NULL,
      style = style
    )
  }

  invisible(x)
}

# ---- biostatistical print tables --------------------------------------------

#' @noRd
comparison_interval_prob_tag <- function(prob) {
  prob <- as.numeric(prob)[[1L]]
  if (isTRUE(all.equal(prob, 0.95))) {
    "95%"
  } else {
    sprintf("%g%%", 100 * prob)
  }
}

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

#' Fold Prediction rows into a `Prediction CrI` column on Overall rows
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
#' both present, Prediction is folded into a `Prediction CrI` column.
#' @noRd
comparison_effects_print_df <- function(summary,
                                  digits = 2,
                                  prob = 0.95,
                                  id_cols = NULL,
                                  collapse_prediction = TRUE,
                                  interval_scope = c("marginal", "single"),
                                  style = c("bayesian", "frequentist")) {
  style <- match.arg(style)
  interval_scope <- match.arg(interval_scope)
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
    "parameter", "label", "IF", "combined", "trial", "subgroup", "interaction", "type"
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

  tag <- comparison_interval_prob_tag(prob)
  interval_lab <- if (identical(interval_scope, "single")) {
    if (identical(style, "frequentist")) {
      paste0(tag, " CI")
    } else {
      paste0(tag, " CrI")
    }
  } else if (identical(style, "frequentist")) {
    paste0("Overall ", tag, " CI")
  } else {
    paste0("Overall ", tag, " CrI")
  }
  pred_lab <- if (identical(style, "frequentist")) {
    paste0("Prediction ", tag, " CI")
  } else {
    paste0("Prediction ", tag, " CrI")
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
                                     collapse_prediction = TRUE,
                                     interval_scope = c("marginal", "single"),
                                     style = c("bayesian", "frequentist")) {
  style <- match.arg(style)
  tab <- comparison_effects_print_df(
    summary,
    digits = digits,
    prob = prob,
    id_cols = id_cols,
    collapse_prediction = collapse_prediction,
    interval_scope = interval_scope,
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
  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
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
    backend = interactmeta_fit_backend(x),
    method = x$method %||% st$method %||% NA_character_,
    inference = x$inference %||% NA_character_,
    algorithm = x$algorithm %||% NA_character_,
    ecological_slope = st$ecological_slope %||% TRUE,
    reference = st$reference,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    modifier = st$modifier %||% subgroup_var,
    trials = trials,
    subgroups = subgroups,
    n_trials = length(trials),
    n_subgroups = length(subgroups),
    n_support = if (!is.null(x$support)) nrow(x$support) else NA_integer_,
    formula = x$formula
  )
}

#' @noRd
comparison_fit_print_header <- function(meta) {
  heading <- interactmeta_modifier_heading(meta$modifier %||% meta$subgroup_var)
  is_subgroup <- identical(tolower(heading), "subgroup")
  method_label <- interactmeta_method_display(
    list(
      method = meta$method,
      inference = meta$inference,
      bayesmeta = if (identical(meta$backend, "bayesmeta")) list() else NULL,
      metafor = if (identical(meta$backend, "metafor")) list() else NULL
    )
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
  cat(paste0("\nmethod: ", method_label, "\n"))
  alg <- as.character(meta$algorithm %||% "")[[1L]]
  inference <- as.character(meta$inference %||% "")[[1L]]
  if (is.finite(meta$n_support) &&
        meta$n_support > 1L &&
        identical(inference, "bayesian") &&
        !identical(alg, "quadrature")) {
    cat(paste0("support rows: ", meta$n_support, "\n"))
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

#' Parameter matrix for print / summary (bayesmeta native table when available)
#' @noRd
comparison_fit_parameters_matrix <- function(x, meta = NULL, prob = 0.95) {
  if (is.null(meta)) {
    meta <- comparison_fit_print_meta(x)
  }
  if (!interactmeta_is_bayesian_fit(x)) {
    interactmeta_metafor_param_summary(x, prob = prob)
  } else if (interactmeta_is_bayesian_fit(x)) {
    interactmeta_bayesmeta_param_summary(x, prob = prob) %||%
      comparison_fit_param_summary_matrix(x, prob = prob)
  } else {
    comparison_fit_param_summary_matrix(x, prob = prob)
  }
}

#' Row names for frequentist fit parameter matrices (metafor Model Results)
#' @noRd
interactmeta_metafor_param_stat_rows <- function() {
  c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
}

#' Tau interval from a stored `rma` object
#' @noRd
interactmeta_metafor_tau_confint <- function(model, prob = 0.95) {
  ci <- tryCatch(
    suppressWarnings(stats::confint(model, level = prob)),
    error = function(e) NULL
  )
  random <- ci$random %||% NULL
  if (is.null(random) || is.null(dim(random)) ||
      !"tau" %in% rownames(random) ||
      !all(c("ci.lb", "ci.ub") %in% colnames(random))) {
    return(c(ci.lb = NA_real_, ci.ub = NA_real_))
  }
  c(
    ci.lb = as.numeric(random["tau", "ci.lb"]),
    ci.ub = as.numeric(random["tau", "ci.ub"])
  )
}

#' Heterogeneity rows (`tau^2`, `tau`) in metafor summary wording
#' @noRd
interactmeta_metafor_heterogeneity_df <- function(model, model_label, prob = 0.95) {
  tau2 <- as.numeric(model$tau2)[[1L]]
  se_tau2 <- as.numeric(model$se.tau2)[[1L]]
  tau <- sqrt(max(tau2, 0))
  tau_ci <- interactmeta_metafor_tau_confint(model, prob = prob)
  dplyr::bind_rows(
    data.frame(
      model = model_label,
      term = "tau^2",
      estimate = tau2,
      se = se_tau2,
      zval = NA_real_,
      pval = NA_real_,
      ci.lb = NA_real_,
      ci.ub = NA_real_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      model = model_label,
      term = "tau",
      estimate = tau,
      se = NA_real_,
      zval = NA_real_,
      pval = NA_real_,
      ci.lb = tau_ci[["ci.lb"]],
      ci.ub = tau_ci[["ci.ub"]],
      stringsAsFactors = FALSE
    )
  )
}

#' Fixed-effect rows from `coef(summary(rma))` with native term names
#' @noRd
interactmeta_metafor_coef_df <- function(model, model_label) {
  cs <- tryCatch(coef(summary(model)), error = function(e) NULL)
  if (is.null(cs) || !nrow(cs)) {
    return(NULL)
  }
  req <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
  if (!all(req %in% colnames(cs))) {
    return(NULL)
  }
  data.frame(
    model = model_label,
    term = rownames(cs),
    cs[, req, drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Long parameter table from stored `rma` objects (metafor aliases)
#' @noRd
interactmeta_metafor_parameters_df <- function(fit, prob = 0.95) {
  prob <- as.numeric(prob)[[1L]]
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  if (is.null(ov) || is.null(int)) {
    return(NULL)
  }
  pieces <- list(
    interactmeta_metafor_heterogeneity_df(ov, "overall", prob = prob),
    interactmeta_metafor_coef_df(ov, "overall"),
    interactmeta_metafor_heterogeneity_df(int, "interaction", prob = prob),
    interactmeta_metafor_coef_df(int, "interaction")
  )
  pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
  if (!length(pieces)) {
    return(NULL)
  }
  dplyr::bind_rows(pieces)
}

#' Public CAMS column labels for the fit parameter matrix (Bayes / metafor)
#' @noRd
interactmeta_cams_param_column_labels <- function(fit) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  nonref <- st$coef_map$ifrac_orig[[1L]] %||%
    (st$ifrac_vars[[1L]] %||% "Older")
  nonref <- sub("^ifrac_", "", as.character(nonref)[[1L]])
  list(
    ecological = paste0("Ecological_", nonref),
    interaction = paste0("Interaction_", nonref)
  )
}

#' One metafor Model Results row as a six-stat vector
#' @noRd
interactmeta_metafor_df_stat_vec <- function(df, model, term) {
  stat_rows <- interactmeta_metafor_param_stat_rows()
  i <- which(df$model == model & df$term == term)
  if (!length(i)) {
    return(stats::setNames(rep(NA_real_, length(stat_rows)), stat_rows))
  }
  stats::setNames(as.numeric(df[i[1L], stat_rows, drop = TRUE]), stat_rows)
}

#' Ecological (beta = theta - gamma) Wald stats from stored `rma` objects
#' @noRd
interactmeta_metafor_ecological_stat_vec <- function(fit, prob = 0.95) {
  stat_rows <- interactmeta_metafor_param_stat_rows()
  lin <- interactmeta_metafor_lincomb(fit, L_ov = c(0, 1), L_int = -1, prob = prob)
  if (is.null(lin)) {
    return(stats::setNames(rep(NA_real_, length(stat_rows)), stat_rows))
  }
  se <- as.numeric(lin$sd)[[1L]]
  est <- as.numeric(lin$estimate)[[1L]]
  zval <- if (is.finite(se) && se > 0) {
    est / se
  } else {
    NA_real_
  }
  pval <- if (is.finite(zval)) {
    2 * stats::pnorm(-abs(zval))
  } else {
    NA_real_
  }
  stats::setNames(
    c(est, se, zval, pval, lin$lci, lin$uci),
    stat_rows
  )
}

#' Frequentist CAMS parameter matrix (metafor Model Results layout)
#' @noRd
interactmeta_metafor_param_summary <- function(fit, prob = 0.95) {
  prob <- as.numeric(prob)[[1L]]
  rep <- fit$reporting
  if (is.list(rep) && !is.null(rep$parameters) &&
        isTRUE(all.equal(as.numeric(rep$prob), prob)) &&
        is.matrix(rep$parameters)) {
    return(rep$parameters)
  }
  df <- interactmeta_metafor_parameters_df(fit, prob = prob)
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }
  stat_rows <- interactmeta_metafor_param_stat_rows()
  labs <- interactmeta_cams_param_column_labels(fit)
  ifrac_term <- df$term[df$model == "overall" & grepl("^ifrac", df$term)][1L]
  if (is.na(ifrac_term)) {
    ifrac_term <- "ifrac"
  }
  mat <- cbind(
    tau = interactmeta_metafor_df_stat_vec(df, "overall", "tau"),
    tau_gamma = interactmeta_metafor_df_stat_vec(df, "interaction", "tau"),
    Intercept = interactmeta_metafor_df_stat_vec(df, "overall", "intrcpt"),
    Ecological = interactmeta_metafor_ecological_stat_vec(fit, prob = prob),
    Interaction = interactmeta_metafor_df_stat_vec(df, "interaction", "intrcpt")
  )
  colnames(mat)[4:5] <- c(labs$ecological, labs$interaction)
  rownames(mat) <- stat_rows
  attr(mat, "interval.type") <- "central"
  mat
}

#' Print the fitted-parameter matrix (shared by print and summary)
#' @noRd
comparison_fit_print_parameters <- function(mat, meta, digits) {
  if (is.null(mat) || !length(mat)) {
    return(invisible(NULL))
  }
  is_bayes <- identical(meta$inference, "bayesian")
  if (is_bayes) {
    cat("\nmarginal posterior summary:\n")
  } else {
    cat("\n", interactmeta_method_display(list(
      method = meta$method,
      inference = meta$inference
    )), " model results:\n", sep = "")
  }
  show <- round(mat, digits = digits)
  attr(show, "interval.type") <- NULL
  print(show)
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
  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
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
  post <- comparison_support_df(x)
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

# ---- tidy.R ----
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
#' @param x A `interactmeta_subgroup_effects` object from [modifiers()].
#' @param ... Ignored.
#' @return A [tibble::tibble].
#' @export
#' @importFrom generics tidy
tidy.interactmeta_subgroup_effects <- function(x, ...) {
  lvl_col <- interactmeta_modifier_summary_col(x$summary)
  id_cols <- if (isTRUE(x$shrinkage)) {
    c("trial", lvl_col, "type")
  } else {
    c(lvl_col, "type")
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
    "combined"
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
