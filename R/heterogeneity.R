#' Between-trial heterogeneity from a CAMS fit
#'
#' Summarizes between-trial standard deviations in the CAMS random-effects
#' structure. For Bayesian CAMS this uses the DIRECT support for
#' \eqn{\tau_m} and \eqn{\tau_g} (not IID MCMC draws). For frequentist CAMS
#' it uses the stored `rma` heterogeneity estimates.
#'
#' @param fit A [`interactmeta`] from [cam()].
#' @param component Which SD parameters to return: `"intercept"`
#'   (trial intercept, treatment-effect heterogeneity at the
#'   reference information fraction), `"slope"` (trial
#'   SDs of `ifrac2_*` slopes), or `"all"` (default).
#' @param trial_var Name of the trial grouping variable; default from the fit's
#'   `"interactmeta"` attribute.
#' @param prob Coverage probability for `lci` / `uci` when `summary = TRUE`.
#' @param summary If `TRUE` (default), include a summary data frame.
#'
#' @return A `interactmeta_heterogeneity` object (list) with elements `draws` (long format
#'   with columns `parameter`, `type`, `value`) and, when requested, `summary`
#'   (`parameter`, `type`, `label`, `estimate`, `mean`, `sd`, `lci`, `uci`,
#'   `prob`). Attributes: `prob`, `trial_var`, `reference`, `component`.
#'
#' @examples
#' \dontrun{
#' library(interactmeta)
#' fit <- cam(es, reference = "Younger")
#' heterogeneity(fit)
#' heterogeneity(fit, component = "intercept")
#' }
#'
#' @seealso [ifrac()], [cam()], [plot_heterogeneity()]
#' @export
heterogeneity <- function(fit,
                          component = c("all", "intercept", "slope"),
                          trial_var = NULL,
                          prob = 0.95,
                          summary = TRUE) {
  if (interactmeta_is_non_cams_analysis(fit)) {
    return(interactmeta_extract_heterogeneity(fit, component = component, prob = prob))
  }
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cam()`.", call. = FALSE)
  }

  component <- match.arg(component)
  specs <- comparison_heterogeneity_specs(fit, trial_var = trial_var, component = component)
  if (!nrow(specs)) {
    stop("No heterogeneity parameters matched `component`.", call. = FALSE)
  }

  draws <- comparison_heterogeneity_draws(fit, specs)
  st <- attr(fit, "interactmeta", exact = TRUE)
  out <- list(draws = draws)
  if (isTRUE(summary)) {
    out$summary <- comparison_heterogeneity_summary(draws, specs = specs, prob = prob)
    out$summary <- interactmeta_metafor_fill_heterogeneity(fit, out$summary)
    out$summary <- interactmeta_bayesmeta_fill_heterogeneity(fit, out$summary, prob = prob)
  }
  attr(out, "prob") <- prob
  attr(out, "trial_var") <- specs$trial_var[[1L]]
  attr(out, "reference") <- st$reference
  attr(out, "component") <- component
  class(out) <- c("interactmeta_heterogeneity", "list")
  out <- interactmeta_attach_origin(out, fit, quantity = "heterogeneity", fitted_or_derived = "fitted")
  attr(out, "interactmeta") <- fit
  out
}

#' Tau bounds from the stored fit's `$post.interval` (not node quantiles).
#' @noRd
interactmeta_bayesmeta_fill_heterogeneity <- function(fit, summary, prob = 0.95) {
  if (!identical(fit$engine, "bayesmeta") || is.null(fit$bayesmeta$overall) ||
        is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  it <- interactmeta_match_interval_type(
    fit$interval.type %||% fit$bayesmeta$overall$interval.type
  )
  prob <- as.numeric(prob)[[1L]]
  for (i in seq_len(nrow(summary))) {
    ty <- as.character(summary$type[[i]])
    obj <- if (identical(ty, "slope")) {
      fit$bayesmeta$interaction
    } else {
      fit$bayesmeta$overall
    }
    loc <- interactmeta_native_location(
      obj,
      prob = prob,
      interval.type = it,
      margin = "tau"
    )
    sm <- obj$summary
    summary$estimate[[i]] <- loc$median
    summary$mean[[i]] <- as.numeric(sm["mean", "tau"])[[1L]]
    summary$sd[[i]] <- as.numeric(sm["sd", "tau"])[[1L]]
    summary$lci[[i]] <- loc$lower
    summary$uci[[i]] <- loc$upper
    summary$prob[[i]] <- prob
  }
  attr(summary, "interval.type") <- it
  summary
}

#' Fill frequentist heterogeneity SE from stored rma tau2 (no profile CI)
#' @noRd
interactmeta_metafor_fill_heterogeneity <- function(fit, summary) {
  if (!identical(fit$engine, "metafor") || is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  if (is.null(ov) || is.null(int)) {
    return(summary)
  }
  for (i in seq_len(nrow(summary))) {
    ty <- as.character(summary$type[[i]])
    tau2 <- if (identical(ty, "slope")) int$tau2 else ov$tau2
    se_tau2 <- if (identical(ty, "slope")) int$se.tau2 else ov$se.tau2
    tau <- sqrt(max(as.numeric(tau2)[[1L]], 0))
    summary$estimate[[i]] <- tau
    summary$mean[[i]] <- tau
    se <- NA_real_
    if (is.finite(as.numeric(se_tau2)[[1L]]) && tau > 0) {
      se <- as.numeric(se_tau2)[[1L]] / (2 * tau)
    }
    summary$sd[[i]] <- se
  }
  summary
}

#' Format heterogeneity summaries as a one-line character string
#'
#' @param x A `interactmeta_heterogeneity` object from [heterogeneity()], or its
#'   `summary` data frame.
#' @param prob Coverage probability associated with `lci` / `uci` (kept for
#'   callers; not printed in the footnote text). When `NULL`, taken from `x`
#'   (or `0.95`).
#' @param digits Decimal places for displayed SDs.
#' @param include_ci If `TRUE` (default), append each parameter's credibility
#'   interval as `(lci, uci)`.
#'
#' @return A single character string (components separated by `" · "`).
#'
#' @seealso [heterogeneity()]
#' @keywords internal
#' @export
heterogeneity_footnote <- function(x, prob = NULL, digits = 2, include_ci = TRUE) {
  if (inherits(x, "interactmeta_heterogeneity")) {
    summ <- x$summary
    if (is.null(summ)) {
      stop("`x` has no `summary`; call `heterogeneity(..., summary = TRUE)`.", call. = FALSE)
    }
    if (is.null(prob)) {
      prob <- attr(x, "prob", exact = TRUE) %||% 0.95
    }
  } else if (is.data.frame(x)) {
    summ <- x
    if (is.null(prob)) {
      prob <- if ("prob" %in% names(summ)) summ$prob[[1L]] else 0.95
    }
  } else {
    stop("`x` must be a `interactmeta_heterogeneity` object or summary data frame.", call. = FALSE)
  }
  force(prob)

  if (!nrow(summ)) {
    return("")
  }

  # Overall (intercept) before interaction (slope) when both are present.
  if ("type" %in% names(summ)) {
    ord <- order(match(as.character(summ$type), c("intercept", "slope"), nomatch = 99L))
    summ <- summ[ord, , drop = FALSE]
  }

  fmt <- function(v) format(round(as.numeric(v), digits), nsmall = digits, trim = TRUE)
  parts <- vapply(seq_len(nrow(summ)), function(i) {
    out <- paste0(summ$label[[i]], ": ", fmt(summ$estimate[[i]]))
    if (isTRUE(include_ci) &&
          "lci" %in% names(summ) &&
          "uci" %in% names(summ) &&
          is.finite(summ$lci[[i]]) &&
          is.finite(summ$uci[[i]])) {
      out <- paste0(
        out,
        " (",
        fmt(summ$lci[[i]]),
        ", ",
        fmt(summ$uci[[i]]),
        ")"
      )
    }
    out
  }, character(1))
  paste(parts, collapse = " \u00b7 ")
}

#' @noRd
comparison_heterogeneity_interaction_label <- function(ifrac2_name, reference) {
  # Footnotes name the component only; do not suffix with (nonref - ref).
  heterogeneity_tau_label("interaction")
}

#' Footnote labels for CAMS heterogeneity SDs
#' @noRd
heterogeneity_tau_label <- function(which = c("overall", "interaction")) {
  which <- match.arg(which)
  if (identical(which, "interaction")) {
    "Interaction heterogeneity (\u03c4_\u03b3)"
  } else {
    "Heterogeneity (\u03c4)"
  }
}

#' @export
print.interactmeta_heterogeneity <- function(x, digits = 2, ...) {
  opts <- comparison_print_digits_prob(x, digits = digits, ...)
  comparison_print_effects_table(
    x$summary,
    digits = opts$digits,
    prob = opts$prob,
    id_cols = c("parameter", "type"),
    title = "Between-trial heterogeneity (trial SD parameters)",
    style = interactmeta_interval_style(x)
  )
  if (is.null(x$summary)) {
    cat("Draws: ", nrow(x$draws), " rows\n", sep = "")
  }
  invisible(x)
}

#' @noRd
comparison_heterogeneity_specs <- function(fit, trial_var = NULL, component = "all") {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st)) {
    stop(
      "Fit has no `\"interactmeta\"` attribute; refit with `cam()`.",
      call. = FALSE
    )
  }

  trial_var <- trial_var %||% st$trial_var
  if (is.null(trial_var) || !nzchar(trial_var)) {
    trial_var <- "trial"
  }

  ifrac2_std <- st$coef_map$ifrac2_std
  ifrac2_orig <- st$coef_map$ifrac2_orig
  if (is.null(ifrac2_std)) {
    ifrac2_std <- st$ifrac2_vars
    ifrac2_orig <- ifrac2_std
  }
  if (is.null(ifrac2_std) || !length(ifrac2_std)) {
    n <- st$n_nonref %||% 0L
    ifrac2_std <- comparison_std_ifrac2_names(n)
    ifrac2_orig <- ifrac2_std
  }

  specs <- data.frame(
    parameter = "Overall",
    hypothesis = sprintf("%s__Intercept = 0", trial_var),
    label = heterogeneity_tau_label("overall"),
    type = "intercept",
    trial_var = trial_var,
    stringsAsFactors = FALSE
  )

  ref <- st$reference
  if (is.null(ref) || !nzchar(ref)) {
    ref <- ""
  }

  if (length(ifrac2_orig)) {
    re_coef <- vapply(
      ifrac2_orig,
      comparison_heterogeneity_re_coef,
      character(1),
      fit = fit,
      ifrac2_std = ifrac2_std
    )
    slope_specs <- data.frame(
      parameter = comparison_interaction_print_name(ifrac2_orig),
      hypothesis = sprintf("%s__%s = 0", trial_var, re_coef),
      label = vapply(
        ifrac2_orig,
        comparison_heterogeneity_interaction_label,
        character(1),
        reference = ref
      ),
      type = "slope",
      trial_var = trial_var,
      stringsAsFactors = FALSE
    )
    specs <- rbind(specs, slope_specs)
  }

  if (identical(component, "intercept")) {
    specs <- specs[specs$type == "intercept", , drop = FALSE]
  } else if (identical(component, "slope")) {
    specs <- specs[specs$type == "slope", , drop = FALSE]
  }

  specs
}

#' Random-effect coefficient name for trial SD columns
#' @noRd
comparison_heterogeneity_re_coef <- function(ifrac2_orig, fit, ifrac2_std = NULL) {
  if (inherits(fit, "interactmeta")) {
    post <- names(fit$draws)
    if (paste0("sd_trial__", ifrac2_orig) %in% post) {
      return(ifrac2_orig)
    }
    if (!is.null(ifrac2_std)) {
      hit <- match(ifrac2_orig, ifrac2_std)
      if (!is.na(hit) && paste0("sd_trial__", ifrac2_std[[hit]]) %in% post) {
        return(ifrac2_std[[hit]])
      }
    }
    return(ifrac2_orig)
  }
  ifrac2_orig
}

#' @noRd
comparison_sd_samples <- function(fit, hypothesis) {
  # hypothesis like "trial__Intercept = 0" or "trial__ifrac2_Older = 0"
  param <- sub("\\s*=\\s*0\\s*$", "", hypothesis)
  col <- paste0("sd_", param)
  post <- comparison_posterior_df(fit)
  if (col %in% names(post)) {
    return(as.numeric(post[[col]]))
  }
  stop(sprintf("Could not find SD draws column `%s`.", col), call. = FALSE)
}

#' @noRd
comparison_heterogeneity_draws <- function(fit, specs) {
  post <- comparison_posterior_df(fit)
  pieces <- lapply(seq_len(nrow(specs)), function(i) {
    vals <- comparison_sd_samples(fit, specs$hypothesis[[i]])
    n <- length(vals)
    df <- data.frame(
      .chain = 1L,
      .iteration = seq_len(n),
      .draw = seq_len(n),
      parameter = specs$parameter[[i]],
      type = specs$type[[i]],
      value = vals,
      stringsAsFactors = FALSE
    )
    if (".weight" %in% names(post) && length(post$.weight) == n) {
      df$.weight <- post$.weight
    }
    df
  })
  dplyr::bind_rows(pieces)
}

#' @noRd
comparison_heterogeneity_summary <- function(draws, specs, prob = 0.95) {
  summ <- comparison_median_qi_grouped(
    draws = draws,
    value_col = "value",
    group_cols = c("parameter", "type"),
    prob = prob
  )
  labels <- specs[, c("parameter", "label"), drop = FALSE]
  dplyr::left_join(summ, labels, by = "parameter")
}
