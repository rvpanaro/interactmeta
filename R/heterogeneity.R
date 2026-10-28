# ---- heterogeneity.R ----
#' Between-trial heterogeneity from a CAMS fit
#'
#' Summarizes between-trial standard deviations in the CAMS random-effects
#' structure. For Bayesian CAMS this uses the DIRECT support for
#' \eqn{\tau} and \eqn{\tau_\gamma}. For frequentist CAMS
#' it uses the stored `rma` heterogeneity estimates (delta-method SE on
#' `\eqn{\hat\tau}`; `lci` / `uci` from [stats::confint()] profile intervals
#' on each orthogonal `rma` block).
#'
#' @param fit A [`interactmeta`] from [cams()].
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
#' library(interact)
#' fit <- cams(es, reference = "Lower")
#' heterogeneity(fit)
#' heterogeneity(fit, component = "intercept")
#' }
#'
#' @seealso [ifrac()], [cams()], [plot_heterogeneity()]
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
    stop("`fit` must be an interactmeta object from `cams()`.", call. = FALSE)
  }

  component <- match.arg(component)
  specs <- comparison_heterogeneity_specs(fit, trial_var = trial_var, component = component)
  if (!nrow(specs)) {
    stop("No heterogeneity parameters matched `component`.", call. = FALSE)
  }

  support <- comparison_heterogeneity_draws(fit, specs)
  st <- attr(fit, "interactmeta", exact = TRUE)
  out <- list(support = support)
  if (isTRUE(summary)) {
    out$summary <- comparison_heterogeneity_summary(support, specs = specs, prob = prob)
    out$summary <- interactmeta_metafor_fill_heterogeneity(
      fit,
      out$summary,
      prob = prob
    )
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
  if (!interactmeta_is_bayesian_fit(fit) || is.null(fit$bayesmeta$overall) ||
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

#' Profile / likelihood-based CI bounds for tau from a stored rma fit
#' @noRd
interactmeta_metafor_rma_tau_ci <- function(model, prob = 0.95) {
  ci <- tryCatch(
    suppressWarnings(stats::confint(model, level = as.numeric(prob)[[1L]])),
    error = function(e) NULL
  )
  random <- ci$random %||% NULL
  if (is.null(random) ||
        is.null(dim(random)) ||
        !"tau" %in% rownames(random) ||
        !all(c("ci.lb", "ci.ub") %in% colnames(random))) {
    return(c(lci = NA_real_, uci = NA_real_))
  }
  c(
    lci = as.numeric(random["tau", "ci.lb"])[[1L]],
    uci = as.numeric(random["tau", "ci.ub"])[[1L]]
  )
}

#' Fill frequentist heterogeneity SE and tau CIs from stored rma fits
#' @noRd
interactmeta_metafor_fill_heterogeneity <- function(fit, summary, prob = 0.95) {
  if (!!interactmeta_is_bayesian_fit(fit) || is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  if (is.null(ov) || is.null(int)) {
    return(summary)
  }
  prob <- as.numeric(prob)[[1L]]
  for (i in seq_len(nrow(summary))) {
    ty <- as.character(summary$type[[i]])
    obj <- if (identical(ty, "slope")) int else ov
    tau2 <- obj$tau2
    se_tau2 <- obj$se.tau2
    tau <- sqrt(max(as.numeric(tau2)[[1L]], 0))
    summary$estimate[[i]] <- tau
    summary$mean[[i]] <- tau
    se <- NA_real_
    if (is.finite(as.numeric(se_tau2)[[1L]]) && tau > 0) {
      se <- as.numeric(se_tau2)[[1L]] / (2 * tau)
    }
    summary$sd[[i]] <- se
    interval <- interactmeta_metafor_rma_tau_ci(obj, prob = prob)
    summary$lci[[i]] <- interval[["lci"]]
    summary$uci[[i]] <- interval[["uci"]]
    summary$prob[[i]] <- prob
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

  # Combined (intercept) before interaction (slope) when both are present.
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
    paste0("Interaction heterogeneity (", "\u03c4", "_", "\u03b3", ")")
  } else {
    paste0("Heterogeneity (", "\u03c4", ")")
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
    cat("Draws: ", nrow(x$support), " rows\n", sep = "")
  }
  invisible(x)
}

#' @noRd
comparison_heterogeneity_specs <- function(fit, trial_var = NULL, component = "all") {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st)) {
    stop(
      "Fit has no `\"interactmeta\"` attribute; refit with `cams()`.",
      call. = FALSE
    )
  }

  trial_var <- trial_var %||% st$trial_var
  if (is.null(trial_var) || !nzchar(trial_var)) {
    trial_var <- interactmeta_study_column()
  }
  re_group <- interactmeta_sd_re_group()

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
    parameter = "Combined",
    hypothesis = sprintf("%s__Intercept = 0", re_group),
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
      hypothesis = sprintf("%s__%s = 0", re_group, re_coef),
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
    post <- names(fit$support)
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
  post <- comparison_support_df(fit)
  if (col %in% names(post)) {
    return(as.numeric(post[[col]]))
  }
  stop(sprintf("Could not find SD draws column `%s`.", col), call. = FALSE)
}

#' @noRd
comparison_heterogeneity_draws <- function(fit, specs) {
  post <- comparison_support_df(fit)
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

# ---- plot_heterogeneity.R ----
#' Heterogeneity plots
#'
#' Posterior distributions of CAMS between-trial standard deviations.
#' Bayesian CAMS uses base-R tau-marginals ([bayesmeta::plot.bayesmeta]
#' `which = 4` style) from the stored `bmr` / `bayesmeta` objects (`tau`,
#' `tau_gamma`). `"both"` is those two `which` plots side by side (`par(mfrow)`),
#' with a shared y-axis. Frequentist CAMS uses the same two panels as
#' profile (restricted) log-likelihoods of `τ` from the stored `rma`
#' objects ([metafor::profile.rma]), also side by side with shared
#' x- and y-axis limits (0 to the larger upper `τ` bound). The region
#' between the profile and the profile cutoff is shaded; a vertical line
#' marks the REML/ML `\hat\tau`; a dotted horizontal line marks the cutoff at
#' `max(ℓ) − qchisq(prob, 1)/2` for the heterogeneity summary `prob`
#' (default 95%), with a small in-panel label.
#'
#' @param x A `interactmeta_heterogeneity` object from [heterogeneity()], or a
#'   [`interactmeta`] from [cams()] (in which case [heterogeneity()] is called).
#' @param which Which components to draw: `"both"` (default; overall /
#'   combined and interaction **side by side** with a shared y-axis),
#'   `"overall"` / `"combined"` (`tau`), or `"interaction"` (`tau_gamma`).
#' @param probs Quantiles to annotate on ggplot densities (default
#'   `c(0.5, 0.75, 0.9)`). Ignored by the bayesmeta and metafor paths.
#' @param bins Number of equal-width histogram bins for continuous samples (default
#'   `40`). Ignored by the bayesmeta and metafor paths and when tau is a
#'   discrete support.
#' @param digits Decimal places for ggplot quantile labels (default `2`).
#' @param colours Length-2 character vector `c(overall, interaction)` for the
#'   ggplot path. Defaults to density slate tokens (not subgroup red/blue).
#' @param base_size Base font size passed to [interactmeta_theme_report()].
#' @param title Optional ggplot title (default `NULL`). Ignored by the
#'   bayesmeta and metafor paths.
#' @param xlab Optional shared x-axis title. The default (`NULL`) labels
#'   each panel `"heterogeneity τ"` (`tau`) or `"heterogeneity τ_γ"`
#'   (`tau_gamma`). A supplied value is used on
#'   ggplot panels only.
#' @param printplot If `TRUE` (default), print the ggplot. Ignored by the
#'   bayesmeta and metafor paths (always drawn).
#' @param prior Overlay the tau prior density (dashed line). On the bayesmeta
#'   path this matches [bayesmeta::plot.bayesmeta]'s `prior`. On the ggplot
#'   fallback, uses `fit$priors` when available (default `TRUE`).
#' @param ... Passed to [heterogeneity()] when `x` is a fit (`prob`,
#'   `component`, …). Ignored for `interactmeta_heterogeneity` objects except that
#'   unused names are dropped quietly.
#'
#' @return For Bayesian or frequentist CAMS with stored native objects,
#'   `x` invisibly after drawing base-R plots. Otherwise a **ggplot2**
#'   object (invisibly when `printplot = TRUE`).
#'
#' @examples
#' \dontrun{
#' fit <- cams(es, reference = "Lower")
#' plot_heterogeneity(fit)
#' plot_heterogeneity(fit, which = "combined")
#' }
#'
#' @seealso [heterogeneity()], [bayesmeta::plot.bayesmeta],
#'   [metafor::profile.rma], `vignette("cams")`
#' @keywords internal
plot_heterogeneity <- function(x,
                               which = c("both", "overall", "interaction", "combined"),
                               probs = c(0.5, 0.75, 0.9),
                               bins = 40L,
                               digits = 2,
                               colours = NULL,
                               base_size = interactmeta_base_size(default = 11L),
                               title = NULL,
                               xlab = NULL,
                               printplot = TRUE,
                               prior = TRUE,
                               ...) {
  which <- plot_heterogeneity_which(which)
  het <- plot_heterogeneity_resolve(x, ...)
  objs <- heterogeneity_native_tau_objects(het, which = which)
  if (!is.null(objs)) {
    heterogeneity_plot_bayesmeta(objs, prior = isTRUE(prior))
    return(invisible(het))
  }
  rma_objs <- heterogeneity_native_rma_objects(het, which = which)
  if (!is.null(rma_objs)) {
    heterogeneity_plot_metafor_profile(
      rma_objs,
      profile_ci = attr(het, "prob", exact = TRUE)
    )
    return(invisible(het))
  }
  p <- plot_heterogeneity_ggplot(
    het,
    which = which,
    probs = probs,
    bins = bins,
    digits = digits,
    colours = colours,
    base_size = base_size,
    title = title,
    xlab = xlab,
    prior = isTRUE(prior)
  )
  if (isTRUE(printplot)) {
    print(p)
    return(invisible(p))
  }
  p
}

#' @rdname plot_heterogeneity
#' @method plot interactmeta_heterogeneity
#' @export
plot.interactmeta_heterogeneity <- function(x,
                                            which = c("both", "overall", "interaction", "combined"),
                                            prior = TRUE,
                                            ...) {
  plot_heterogeneity(x, which = which, prior = prior, ...)
}

#' @noRd
plot_heterogeneity_which <- function(which) {
  which <- as.character(which)[[1L]]
  if (which %in% c("combined", "combined_population")) {
    which <- "overall"
  }
  match.arg(which, c("both", "overall", "interaction"))
}

#' Native bayesmeta / bmr objects for tau plots (no refit)
#' @noRd
heterogeneity_native_tau_objects <- function(het, which = "both") {
  fit <- attr(het, "interactmeta", exact = TRUE)
  if (is.null(fit) || !interactmeta_is_bayesian_fit(fit)) {
    return(NULL)
  }
  bm <- fit$bayesmeta
  if (is.null(bm)) {
    return(NULL)
  }
  objs <- list()
  if (which %in% c("both", "overall") && !is.null(bm$overall)) {
    objs$tau <- bm$overall
  }
  if (which %in% c("both", "interaction") && !is.null(bm$interaction)) {
    objs$tau_gamma <- bm$interaction
  }
  if (!length(objs)) {
    NULL
  } else {
    objs
  }
}

#' Native metafor rma objects for tau profile plots (no refit)
#' @noRd
heterogeneity_native_rma_objects <- function(het, which = "both") {
  fit <- attr(het, "interactmeta", exact = TRUE)
  if (is.null(fit) || interactmeta_is_bayesian_fit(fit)) {
    return(NULL)
  }
  mf <- fit$metafor
  if (is.null(mf)) {
    return(NULL)
  }
  objs <- list()
  if (which %in% c("both", "overall") && !is.null(mf$overall)) {
    objs$tau <- mf$overall
  }
  if (which %in% c("both", "interaction") && !is.null(mf$interaction)) {
    objs$tau_gamma <- mf$interaction
  }
  if (!length(objs)) {
    NULL
  } else {
    objs
  }
}

#' Side-by-side bayesmeta tau-marginal plots (shared ylim)
#' @noRd
heterogeneity_plot_bayesmeta <- function(objs, prior = FALSE, ...) {
  check_bayesmeta()
  force(list(...))
  panels <- lapply(names(objs), function(nm) {
    d <- heterogeneity_bayesmeta_tau_prep(objs[[nm]], prior = prior)
    d$xlab <- heterogeneity_tau_xlab(nm)
    d
  })
  ylim <- heterogeneity_plot_bayesmeta_ylim(panels, prior = isTRUE(prior))
  n <- length(panels)
  if (n > 1L) {
    old_par <- graphics::par(mfrow = n2mfrow_traceplot(n), oma = c(0, 0, 0, 0))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  for (d in panels) {
    heterogeneity_plot_bayesmeta_tau_draw(d, ylim = ylim)
  }
  invisible(objs)
}

#' @noRd
heterogeneity_plot_bayesmeta_ylim <- function(panels, prior = FALSE) {
  ymax <- max(vapply(panels, function(d) {
    y <- d$maxdens
    if (isTRUE(prior) && !is.null(d$pdens)) {
      y <- c(y, d$pdens[is.finite(d$pdens)])
    }
    max(y, na.rm = TRUE)
  }, numeric(1)))
  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- 1
  }
  c(0, ymax)
}

#' Panel x-axis label: `tau` vs `tau_gamma`
#' @noRd
heterogeneity_tau_xlab <- function(nm) {
  nm <- as.character(nm)[[1L]]
  if (nm %in% c("tau", "Combined", "overall", "combined")) {
    return(expression("heterogeneity " * tau))
  }
  if (nm %in% c("tau_gamma", "Interaction", "interaction")) {
    return(expression("heterogeneity " * tau[gamma]))
  }
  expression("heterogeneity " * tau)
}

#' ggplot facet labeller for panel-specific tau x-axis labels
#' @noRd
heterogeneity_tau_facet_labeller <- function(labels) {
  list(panel = lapply(as.character(labels$panel), heterogeneity_tau_xlab))
}

#' @noRd
heterogeneity_bayesmeta_tau_prep <- function(x, prior = FALSE) {
  tau_hi <- tryCatch(
    as.numeric(x$qposterior(tau = 0.995))[[1L]],
    error = function(e) NA_real_
  )
  if (!is.finite(tau_hi) || tau_hi <= 0) {
    tau_hi <- 1
  }
  taurange <- c(0, tau_hi * 1.1)
  tau <- seq(
    max(c(0, taurange[[1L]] - 0.1 * diff(taurange))),
    taurange[[2L]] + 0.1 * diff(taurange),
    length.out = 200L
  )
  dens <- as.numeric(x$dposterior(tau = tau))
  maxdens <- max(dens[is.finite(dens)], na.rm = TRUE)
  if (!is.finite(maxdens) || maxdens <= 0) {
    maxdens <- 1
  }
  dens[!is.finite(dens)] <- 10 * maxdens
  lo <- as.numeric(x$summary["95% lower", "tau"])[[1L]]
  med <- as.numeric(x$summary["median", "tau"])[[1L]]
  hi <- as.numeric(x$summary["95% upper", "tau"])[[1L]]
  d_lo <- NA_real_
  d_hi <- NA_real_
  indi <- NULL
  if (is.finite(lo) && is.finite(hi) && hi > lo) {
    indi <- (tau >= lo) & (tau <= hi)
    d_lo <- as.numeric(x$dposterior(tau = lo))[[1L]]
    d_hi <- as.numeric(x$dposterior(tau = hi))[[1L]]
    if (!is.finite(d_lo)) {
      d_lo <- 10 * maxdens
    }
    if (!is.finite(d_hi)) {
      d_hi <- 10 * maxdens
    }
  }
  d_med <- NA_real_
  if (is.finite(med)) {
    d_med <- as.numeric(x$dposterior(tau = med))[[1L]]
  }
  pdens <- NULL
  if (isTRUE(prior)) {
    prior_ok <- isTRUE(x$tau.prior.proper) || is.function(x$dprior)
    if (isTRUE(prior_ok)) {
      pdens <- tryCatch(
        as.numeric(x$dprior(tau = tau)),
        error = function(e) NULL
      )
      if (!is.null(pdens) && length(pdens) != length(tau)) {
        pdens <- NULL
      }
    }
  }
  list(
    x = x,
    tau = tau,
    dens = dens,
    taurange = taurange,
    maxdens = maxdens,
    lo = lo,
    med = med,
    hi = hi,
    d_lo = d_lo,
    d_hi = d_hi,
    d_med = d_med,
    indi = indi,
    pdens = pdens,
    prior = isTRUE(prior),
    xlab = heterogeneity_tau_xlab("tau")
  )
}

#' @noRd
heterogeneity_plot_bayesmeta_tau_draw <- function(d, ylim = NULL) {
  if (is.null(ylim)) {
    ylim <- c(0, d$maxdens)
  }
  graphics::plot(
    c(d$taurange[[1L]], d$taurange[[2L]]),
    ylim,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  graphics::polygon(
    c(0, d$tau, max(d$tau)),
    c(0, d$dens, 0),
    border = NA,
    col = grDevices::grey(0.9)
  )
  if (!is.null(d$indi) && any(d$indi)) {
    graphics::polygon(
      c(rep(d$lo, 2L), d$tau[d$indi], rep(d$hi, 2L)),
      c(0, d$d_lo, d$dens[d$indi], d$d_hi, 0),
      border = NA,
      col = grDevices::grey(0.8)
    )
  }
  if (is.finite(d$med) && is.finite(d$d_med)) {
    graphics::lines(rep(d$med, 2L), c(0, d$d_med), col = grDevices::grey(0.6))
  }
  graphics::lines(d$tau, d$dens, col = "black")
  graphics::abline(h = 0, v = 0, col = grDevices::grey(0.4))
  if (!is.null(d$pdens)) {
    graphics::lines(d$tau, d$pdens, col = "black", lty = "dashed")
  }
  xlab <- d$xlab
  if (is.null(xlab)) {
    xlab <- heterogeneity_tau_xlab("tau")
  }
  graphics::mtext(
    side = 1,
    line = graphics::par("mgp")[[1L]],
    xlab
  )
  graphics::mtext(
    side = 2,
    line = graphics::par("mgp")[[2L]],
    expression("marginal posterior density")
  )
  graphics::axis(1)
  graphics::box()
  invisible(d$x)
}

#' Profile-likelihood confidence levels for tau (chi-squared, df = 1)
#' @noRd
heterogeneity_metafor_profile_levels <- function(het_prob = NULL) {
  x <- as.numeric(het_prob)
  x <- x[is.finite(x) & x > 0 & x < 1]
  if (length(x) > 1L) {
    return(sort(unique(x)))
  }
  prob <- if (length(x)) x[[1L]] else 0.95
  prob
}

#' @noRd
heterogeneity_metafor_profile_cuts <- function(maxll, levels) {
  if (!is.finite(maxll) || !length(levels)) {
    return(numeric())
  }
  stats::setNames(
    maxll - 0.5 * stats::qchisq(as.numeric(levels), df = 1L),
    paste0(round(100 * levels), "%")
  )
}

#' Side-by-side metafor tau profile-likelihood plots (shared xlim and ylim)
#' @noRd
heterogeneity_plot_metafor_profile <- function(objs,
                                               steps = 20L,
                                               profile_ci = NULL) {
  profile_ci <- heterogeneity_metafor_profile_levels(het_prob = profile_ci)
  panels <- lapply(names(objs), function(nm) {
    d <- heterogeneity_metafor_tau_prep(
      objs[[nm]],
      steps = steps,
      profile_ci = profile_ci
    )
    if (is.null(d)) {
      return(NULL)
    }
    d$xlab <- heterogeneity_tau_xlab(nm)
    d
  })
  panels <- Filter(Negate(is.null), panels)
  if (!length(panels)) {
    stop(
      "Cannot profile tau for this frequentist CAMS fit.",
      call. = FALSE
    )
  }
  ylim <- heterogeneity_plot_metafor_ylim(panels)
  xlim <- heterogeneity_plot_metafor_xlim(panels)
  n <- length(panels)
  if (n > 1L) {
    old_par <- graphics::par(mfrow = n2mfrow_traceplot(n), oma = c(0, 0, 0, 0))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  for (d in panels) {
    heterogeneity_plot_metafor_tau_draw(d, ylim = ylim, xlim = xlim)
  }
  invisible(objs)
}

#' @noRd
heterogeneity_plot_metafor_xlim <- function(panels) {
  hi <- max(vapply(panels, function(d) d$xlim[[2L]], numeric(1)))
  if (!is.finite(hi) || hi <= 0) {
    hi <- 1
  }
  c(0, hi)
}

#' @noRd
heterogeneity_plot_metafor_ylim <- function(panels) {
  lo <- min(vapply(panels, function(d) d$ylim[[1L]], numeric(1)))
  hi <- max(vapply(panels, function(d) d$ylim[[2L]], numeric(1)))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    c(-1, 0)
  } else {
    c(lo, hi)
  }
}

#' Profile (restricted) log-likelihood of tau from a stored rma fit
#' @noRd
heterogeneity_metafor_tau_prep <- function(x,
                                           steps = 20L,
                                           profile_ci = NULL) {
  pr <- heterogeneity_metafor_profile_tau(x, steps = steps)
  if (is.null(pr)) {
    return(NULL)
  }
  tau <- if (!is.null(pr$tau)) {
    as.numeric(pr$tau)
  } else {
    sqrt(pmax(as.numeric(pr$tau2), 0))
  }
  ll <- as.numeric(pr$ll)
  ok <- is.finite(tau) & is.finite(ll)
  tau <- tau[ok]
  ll <- ll[ok]
  if (!length(tau)) {
    return(NULL)
  }
  vc <- sqrt(max(as.numeric(x$tau2)[[1L]], 0))
  if (!is.finite(vc)) {
    vc <- NA_real_
  }
  maxll <- as.numeric(pr$maxll)[[1L]]
  if (!is.finite(maxll)) {
    maxll <- max(ll)
  }
  if (is.finite(vc) && min(tau) <= vc && max(tau) >= vc) {
    pos <- which(tau >= vc)[[1L]]
    tau <- c(tau[seq_len(pos - 1L)], vc, tau[pos:length(tau)])
    ll <- c(ll[seq_len(pos - 1L)], maxll, ll[pos:length(ll)])
  }
  xlim <- c(0, max(c(as.numeric(pr$xlim), tau, vc), na.rm = TRUE))
  if (!all(is.finite(xlim)) || xlim[[2L]] <= xlim[[1L]]) {
    xlim <- c(0, max(tau, 1))
  }
  ylim <- as.numeric(pr$ylim)
  if (length(ylim) != 2L || !all(is.finite(ylim))) {
    ylim <- range(c(maxll, ll), na.rm = TRUE) + c(-0.1, 0.1)
  }
  if (is.null(profile_ci)) {
    profile_ci <- heterogeneity_metafor_profile_levels(NULL)
  } else {
    profile_ci <- heterogeneity_metafor_profile_levels(het_prob = profile_ci)
  }
  profile_cuts <- heterogeneity_metafor_profile_cuts(maxll, profile_ci)
  ylim[[1L]] <- min(
    c(ylim[[1L]], ll, unname(profile_cuts)),
    na.rm = TRUE
  )
  ylim[[2L]] <- max(
    c(ylim[[2L]], ll, maxll, unname(profile_cuts)),
    na.rm = TRUE
  )
  list(
    x = x,
    tau = tau,
    ll = ll,
    vc = vc,
    maxll = maxll,
    profile_ci = profile_ci,
    profile_cuts = profile_cuts,
    xlim = xlim,
    ylim = ylim,
    ylab = heterogeneity_metafor_profile_ylab(pr$method %||% x$method),
    xlab = heterogeneity_tau_xlab("tau")
  )
}

#' @noRd
heterogeneity_metafor_profile_ylab <- function(method) {
  if (identical(as.character(method %||% "")[[1L]], "REML")) {
    "restricted log-likelihood"
  } else {
    "log-likelihood"
  }
}

#' @noRd
heterogeneity_metafor_profile_tau <- function(x, steps = 20L) {
  if (!inherits(x, "rma")) {
    return(NULL)
  }
  steps <- max(5L, as.integer(steps)[[1L]])
  withCallingHandlers(
    tryCatch(
      stats::profile(
        x,
        plot = FALSE,
        progbar = FALSE,
        sqrt = TRUE,
        steps = steps
      ),
      error = function(e) NULL
    ),
    warning = function(w) {
      if (grepl(
        "larger than the log-likelihood",
        conditionMessage(w),
        ignore.case = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Polygon(s) for profile-likelihood set {tau : ell(tau) >= y_cut}
#' @noRd
heterogeneity_metafor_cutoff_polygons <- function(tau, ll, y_cut) {
  if (!is.finite(y_cut) || length(tau) < 2L) {
    return(list())
  }
  o <- order(tau)
  tau <- tau[o]
  ll <- ll[o]
  n <- length(tau)
  cross_tau <- function(i) {
    t1 <- tau[i]
    t2 <- tau[i + 1L]
    y1 <- ll[i]
    y2 <- ll[i + 1L]
    if (!is.finite(y1) || !is.finite(y2) || y1 == y2) {
      return(t1)
    }
    t1 + (y_cut - y1) * (t2 - t1) / (y2 - y1)
  }
  in_set <- ll >= y_cut
  if (!any(in_set)) {
    return(list())
  }
  polys <- list()
  i <- 1L
  while (i <= n) {
    if (!in_set[i]) {
      i <- i + 1L
      next
    }
    t0 <- if (i > 1L && !in_set[i - 1L]) {
      cross_tau(i - 1L)
    } else {
      tau[i]
    }
    j <- i
    while (j < n && in_set[j + 1L]) {
      j <- j + 1L
    }
    t1 <- if (j < n && !in_set[j + 1L]) {
      cross_tau(j)
    } else {
      tau[j]
    }
    if (t1 < t0) {
      tmp <- t0
      t0 <- t1
      t1 <- tmp
    }
    tx <- unique(c(t0, tau[tau >= t0 & tau <= t1], t1))
    ty_top <- stats::approx(tau, ll, xout = tx, rule = 2)$y
    ty_top <- pmax(ty_top, y_cut)
    polys[[length(polys) + 1L]] <- list(
      x = c(tx, rev(tx)),
      y = c(ty_top, rep(y_cut, length(tx)))
    )
    i <- j + 1L
  }
  polys
}

#' @noRd
heterogeneity_plot_metafor_tau_draw <- function(d, ylim = NULL, xlim = NULL) {
  if (is.null(ylim)) {
    ylim <- d$ylim
  }
  if (is.null(xlim)) {
    xlim <- d$xlim
  }
  tau <- d$tau
  ll <- d$ll
  graphics::plot(
    xlim,
    ylim,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  ok <- is.finite(tau) & is.finite(ll)
  tau <- tau[ok]
  ll <- ll[ok]
  cuts <- d$profile_cuts
  if (length(cuts)) {
    cuts <- cuts[is.finite(cuts)]
  } else {
    cuts <- numeric()
  }
  if (length(cuts)) {
    for (i in seq_along(cuts)) {
      y_cut <- unname(cuts[[i]])
      lab <- names(cuts)[[i]]
      if (!nzchar(lab)) {
        lab <- "profile CI"
      }
      shade <- heterogeneity_metafor_cutoff_polygons(tau, ll, y_cut)
      for (pg in shade) {
        graphics::polygon(
          pg$x,
          pg$y,
          border = NA,
          col = grDevices::grey(0.88)
        )
      }
      graphics::abline(
        h = y_cut,
        lty = 3,
        col = grDevices::grey(0.72),
        lwd = 0.75
      )
      graphics::text(
        x = xlim[[2L]],
        y = y_cut,
        labels = paste0(lab, " cutoff"),
        adj = c(1, 0.5),
        cex = graphics::par("cex") * 0.72,
        col = grDevices::grey(0.45),
        xpd = NA
      )
    }
  }
  if (length(tau)) {
    graphics::lines(tau, ll, col = "black", lwd = 1.4)
  }
  if (is.finite(d$vc)) {
    graphics::abline(v = d$vc, lty = 1, col = grDevices::grey(0.55), lwd = 0.9)
    if (is.finite(d$maxll)) {
      graphics::points(d$vc, d$maxll, pch = 19, cex = 0.65, col = "black")
    }
  }
  graphics::abline(v = xlim[[1L]], col = grDevices::grey(0.4), lwd = 0.5)
  xlab <- d$xlab
  if (is.null(xlab)) {
    xlab <- heterogeneity_tau_xlab("tau")
  }
  graphics::mtext(
    side = 1,
    line = graphics::par("mgp")[[1L]],
    xlab
  )
  graphics::mtext(
    side = 2,
    line = graphics::par("mgp")[[2L]],
    d$ylab
  )
  graphics::axis(1)
  graphics::box()
  invisible(d$x)
}

#' Side-by-side heterogeneity ggplot used by [plot_heterogeneity()]
#' @noRd
plot_heterogeneity_ggplot <- function(het,
                                      which = c("both", "overall", "interaction", "combined"),
                                      probs = c(0.5, 0.75, 0.9),
                                      bins = 40L,
                                      digits = 2,
                                      colours = NULL,
                                      base_size = interactmeta_base_size(default = 11L),
                                      title = NULL,
                                      xlab = NULL,
                                      prior = TRUE) {
  if (!inherits(het, "interactmeta_heterogeneity")) {
    stop("`het` must be a `interactmeta_heterogeneity` object from `heterogeneity()`.", call. = FALSE)
  }
  which <- plot_heterogeneity_which(which)
  probs <- sort(unique(as.numeric(probs)))
  probs <- probs[is.finite(probs) & probs > 0 & probs < 1]
  if (!length(probs)) {
    stop("`probs` must contain values in (0, 1).", call. = FALSE)
  }
  bins <- max(5L, as.integer(bins)[[1L]])
  digits <- max(0L, as.integer(digits)[[1L]])
  colours <- plot_heterogeneity_colours(colours)
  const <- interactmeta_plot_constants()

  pieces <- plot_heterogeneity_apply_colours(
    plot_heterogeneity_pieces(het, which = which),
    colours
  )
  hist_df <- dplyr::bind_rows(lapply(pieces, function(pc) {
    plot_heterogeneity_hist(
      pc$values,
      panel = pc$panel,
      bins = bins,
      fill = pc$fill,
      weight = pc$weight
    )
  }))
  dens_df <- dplyr::bind_rows(lapply(pieces, function(pc) {
    plot_heterogeneity_density(
      pc$values,
      panel = pc$panel,
      colour = pc$line,
      weight = pc$weight
    )
  }))
  ann_df <- dplyr::bind_rows(lapply(pieces, function(pc) {
    plot_heterogeneity_annotations(
      pc$values,
      panel = pc$panel,
      probs = probs,
      digits = digits,
      colour = pc$line,
      weight = pc$weight
    )
  }))
  fit <- attr(het, "interactmeta", exact = TRUE)
  prior_df <- plot_heterogeneity_prior_df(
    pieces = pieces,
    fit = fit,
    prior = prior,
    x_hi = max(c(
      hist_df$xmax,
      dens_df$x,
      ann_df$x,
      if (length(pieces)) {
        vapply(pieces, function(pc) max(pc$values, na.rm = TRUE), numeric(1))
      } else {
        numeric()
      }
    ), na.rm = TRUE)
  )
  panel_levels <- vapply(pieces, function(pc) pc$panel, character(1))
  if (nrow(hist_df)) {
    hist_df$panel <- factor(hist_df$panel, levels = panel_levels)
  } else {
    hist_df$panel <- factor(character(), levels = panel_levels)
  }
  if (nrow(dens_df)) {
    dens_df$panel <- factor(dens_df$panel, levels = panel_levels)
  } else {
    dens_df$panel <- factor(character(), levels = panel_levels)
  }
  if (nrow(ann_df)) {
    ann_df$panel <- factor(ann_df$panel, levels = panel_levels)
  } else {
    ann_df$panel <- factor(character(), levels = panel_levels)
  }
  if (nrow(prior_df)) {
    prior_df$panel <- factor(prior_df$panel, levels = panel_levels)
  } else {
    prior_df$panel <- factor(character(), levels = panel_levels)
  }

  xmax <- max(c(hist_df$xmax, dens_df$x, ann_df$x, prior_df$x), na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) {
    xmax <- 1
  }
  ymax <- max(c(hist_df$ymax, dens_df$y, ann_df$y, prior_df$y), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- 1
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = const$axis_line,
      linewidth = const$line_hair
    ) +
    ggplot2::geom_rect(
      data = hist_df,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$fill
      ),
      colour = NA,
      alpha = 0.45
    ) +
    ggplot2::geom_line(
      data = dens_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$panel,
        colour = .data$colour
      ),
      linewidth = const$line_mid
    )
  if (nrow(prior_df)) {
    p <- p + ggplot2::geom_line(
      data = prior_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$panel
      ),
      colour = "black",
      linetype = "dashed",
      linewidth = const$line_mid
    )
  }
  p <- p +
    ggplot2::geom_segment(
      data = ann_df[ann_df$vline, , drop = FALSE],
      ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = .data$y),
      colour = const$null_line,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = ann_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y_lab,
        label = .data$label,
        colour = .data$colour
      ),
      vjust = -0.35,
      size = base_size / 3.6,
      fontface = "plain"
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.06))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::coord_cartesian(
      xlim = c(0, xmax * 1.05),
      ylim = c(0, ymax * 1.2),
      clip = "off"
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    interactmeta_theme_report(base_size = base_size, variant = "density") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 10, 8, 8, unit = "pt")
    )
  if (!is.null(xlab)) {
    p <- p + ggplot2::labs(x = xlab)
    if (length(panel_levels) > 1L) {
      p <- p + ggplot2::facet_wrap(~panel, ncol = 2L)
    }
  } else if (length(panel_levels) > 1L) {
    p <- p +
      ggplot2::facet_wrap(
        ~panel,
        ncol = 2L,
        strip.position = "bottom",
        labeller = heterogeneity_tau_facet_labeller
      ) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.placement = "outside"
      )
  } else if (length(panel_levels) == 1L) {
    p <- p + ggplot2::labs(x = heterogeneity_tau_xlab(panel_levels[[1L]]))
  }

  p
}

# --- internals --------------------------------------------------------------

#' @noRd
plot_heterogeneity_resolve <- function(x, ...) {
  if (inherits(x, "interactmeta_construction_heterogeneity")) {
    stop(
      "plot_heterogeneity() requires the CAMS information-fraction construction. ",
      "DA, AD, and SWADA heterogeneity is not the CAMS intercept/slope pair, ",
      "and a Bayesian posterior does not add CAMS information-fraction functionality.",
      call. = FALSE
    )
  }
  interactmeta_abort_non_cams(x, "plot_heterogeneity()")
  if (inherits(x, "interactmeta_heterogeneity")) {
    return(x)
  }
  if (inherits(x, "interactmeta")) {
    dots <- list(...)
    # Drop plot-only args that may arrive via plot.interactmeta_heterogeneity → ...
    drop <- c(
      "which", "probs", "bins", "digits", "colours", "base_size",
      "title", "xlab", "printplot", "prior"
    )
    dots <- dots[setdiff(names(dots), drop)]
    return(do.call(heterogeneity, c(list(fit = x), dots)))
  }
  stop(
    "`x` must be a `interactmeta_heterogeneity` from `heterogeneity()` or an `interactmeta` from `cams()`.",
    call. = FALSE
  )
}

#' @noRd
plot_heterogeneity_colours <- function(colours = NULL) {
  const <- interactmeta_plot_constants()
  default <- c(
    overall = const$density_fill,
    interaction = const$density_point
  )
  if (is.null(colours)) {
    return(default)
  }
  colours <- as.character(colours)
  if (length(colours) == 1L) {
    colours <- rep(colours, 2L)
  }
  if (length(colours) < 2L) {
    stop("`colours` must be length 2 (overall, interaction).", call. = FALSE)
  }
  stats::setNames(colours[1:2], c("overall", "interaction"))
}

#' @noRd
plot_heterogeneity_prior_label <- function(type) {
  if (identical(type, "intercept")) {
    "tau"
  } else {
    "tau_gamma"
  }
}

#' Tau prior curves for the ggplot heterogeneity path
#' @noRd
plot_heterogeneity_prior_df <- function(pieces,
                                        fit,
                                        prior = TRUE,
                                        x_hi = NULL,
                                        n = 200L) {
  empty <- data.frame(
    x = numeric(),
    y = numeric(),
    panel = character(),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(prior) || is.null(fit) || !length(pieces)) {
    return(empty)
  }
  priors <- fit$priors
  if (is.null(priors)) {
    return(empty)
  }
  tau_ov <- priors$tau
  tau_int <- priors$tau_interaction
  if (is.null(tau_ov) && inherits(priors, "cam_priors")) {
    tau_ov <- priors$tau
    tau_int <- priors$tau_interaction
  }
  if (is.null(tau_ov) && is.list(priors$common)) {
    tau_ov <- priors$common$tau
    tau_int <- priors$common$tau_interaction
  }
  n <- as.integer(n)[[1L]]
  rows <- lapply(pieces, function(pc) {
    pr <- if (identical(pc$type, "intercept")) tau_ov else tau_int
    if (is.null(pr)) {
      return(NULL)
    }
    hi <- x_hi
    if (!is.finite(hi) || hi <= 0) {
      hi <- max(pc$values, na.rm = TRUE)
    }
    if (!is.finite(hi) || hi <= 0) {
      hi <- 1
    }
    tau <- seq(0, hi * 1.1, length.out = n)
    y <- tryCatch(
      cam_prior_tau_density(pr, tau),
      error = function(e) NULL
    )
    if (is.null(y) || !length(y)) {
      return(NULL)
    }
    y[!is.finite(y)] <- NA_real_
    data.frame(
      x = tau,
      y = as.numeric(y),
      panel = pc$panel,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(empty)
  }
  dplyr::bind_rows(rows)
}

#' @noRd
plot_heterogeneity_pieces <- function(het, which = "both") {
  draws <- het$support
  if (is.null(draws) || !nrow(draws)) {
    stop("`het` has no posterior draws.", call. = FALSE)
  }
  types <- unique(as.character(draws$type))
  want_int <- which %in% c("both", "overall")
  want_slp <- which %in% c("both", "interaction")
  has_int <- want_int && "intercept" %in% types
  has_slp <- want_slp && "slope" %in% types

  pieces <- list()
  if (has_int) {
    ii <- draws$type == "intercept"
    vals <- draws$value[ii]
    pieces[[length(pieces) + 1L]] <- list(
      panel = "tau",
      type = "intercept",
      values = vals,
      weight = plot_heterogeneity_piece_weight(draws, ii, length(vals)),
      fill = NA_character_,
      line = NA_character_,
      title = plot_heterogeneity_prior_label("intercept")
    )
  }
  if (has_slp) {
    sl <- draws[draws$type == "slope", , drop = FALSE]
    # One interaction-block parameter (first if several).
    param <- unique(as.character(sl$parameter))[[1L]]
    keep <- sl$parameter == param
    vals <- sl$value[keep]
    pieces[[length(pieces) + 1L]] <- list(
      panel = "tau_gamma",
      type = "slope",
      values = vals,
      weight = plot_heterogeneity_piece_weight(sl, keep, length(vals)),
      fill = NA_character_,
      line = NA_character_,
      title = plot_heterogeneity_prior_label("slope")
    )
  }
  if (!length(pieces)) {
    stop(
      "No matching heterogeneity draws for `which = \"", which, "\"`.",
      call. = FALSE
    )
  }
  pieces
}

#' @noRd
plot_heterogeneity_apply_colours <- function(pieces, colours) {
  colours <- plot_heterogeneity_colours(colours)
  lapply(pieces, function(pc) {
    key <- if (identical(pc$type, "intercept")) "overall" else "interaction"
    pc$line <- colours[[key]]
    pc$fill <- colours[[key]]
    pc
  })
}

#' @noRd
plot_heterogeneity_piece_weight <- function(df, i, n) {
  if (!is.data.frame(df) || !".weight" %in% names(df)) {
    return(NULL)
  }
  w <- as.numeric(df$.weight[i])
  if (length(w) != n || !any(is.finite(w) & w > 0)) {
    return(NULL)
  }
  w
}

#' @noRd
plot_heterogeneity_normalize_weight <- function(values, weight) {
  n <- length(values)
  if (is.null(weight) || length(weight) != n) {
    return(rep(1 / max(n, 1L), n))
  }
  w <- pmax(as.numeric(weight), 0)
  w[!is.finite(w)] <- 0
  if (!sum(w)) {
    return(rep(1 / max(n, 1L), n))
  }
  w / sum(w)
}

#' Discrete tau support (quadrature nodes) versus continuous samples
#' @noRd
plot_heterogeneity_is_discrete <- function(values, bins = 40L) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  n <- length(values)
  if (n < 2L) {
    return(FALSE)
  }
  n_u <- length(unique(values))
  n_u <= max(2L * as.integer(bins)[[1L]], 80L) && n_u < 0.5 * n
}

#' @noRd
plot_heterogeneity_hist <- function(values,
                                    panel = "tau",
                                    bins = 40L,
                                    fill = "grey",
                                    weight = NULL) {
  empty <- data.frame(
    xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric(),
    panel = character(), fill = character(), stringsAsFactors = FALSE
  )
  values <- as.numeric(values)
  ok <- is.finite(values) & values >= 0
  values <- values[ok]
  if (!length(values)) {
    return(empty)
  }
  weight <- plot_heterogeneity_normalize_weight(values, if (is.null(weight)) NULL else weight[ok])
  if (length(unique(values)) == 1L) {
    return(empty)
  }
  if (plot_heterogeneity_is_discrete(values, bins = bins)) {
    agg <- stats::aggregate(list(w = weight), by = list(x = values), FUN = sum)
    o <- order(agg$x)
    u <- agg$x[o]
    w <- agg$w[o] / sum(agg$w)
    edges <- plot_heterogeneity_node_edges(u)
    width <- pmax(edges$right - edges$left, .Machine$double.eps)
    y <- w / width
    return(data.frame(
      xmin = edges$left,
      xmax = edges$right,
      ymin = 0,
      ymax = y,
      panel = panel,
      fill = fill,
      stringsAsFactors = FALSE
    ))
  }
  br <- graphics::hist(values, breaks = bins, plot = FALSE)
  data.frame(
    xmin = br$breaks[-length(br$breaks)],
    xmax = br$breaks[-1L],
    ymin = 0,
    ymax = br$density,
    panel = panel,
    fill = fill,
    stringsAsFactors = FALSE
  )
}

#' Voronoi edges on [0, Inf) for unique support nodes
#' @noRd
plot_heterogeneity_node_edges <- function(u) {
  u <- sort(unique(as.numeric(u)))
  if (length(u) == 1L) {
    half <- max(u[[1L]] / 2, 0.02)
    return(list(left = max(0, u[[1L]] - half), right = u[[1L]] + half))
  }
  mid <- (u[-length(u)] + u[-1L]) / 2
  left0 <- max(0, u[[1L]] - (mid[[1L]] - u[[1L]]))
  rightk <- u[[length(u)]] + (u[[length(u)]] - mid[[length(mid)]])
  list(left = c(left0, mid), right = c(mid, rightk))
}

#' @noRd
plot_heterogeneity_safe_density <- function(values, n = 512L, weight = NULL) {
  values <- as.numeric(values)
  ok <- is.finite(values) & values >= 0
  values <- values[ok]
  if (length(values) < 2L || isTRUE(stats::sd(values) <= 0)) {
    return(NULL)
  }
  w <- NULL
  if (!is.null(weight)) {
    w <- plot_heterogeneity_normalize_weight(values, weight[ok])
    if (length(unique(w)) <= 1L) {
      w <- NULL
    }
  }
  tryCatch(
    if (is.null(w)) {
      stats::density(values, from = 0, n = as.integer(n)[[1L]])
    } else {
      stats::density(values, from = 0, n = as.integer(n)[[1L]], weights = w)
    },
    error = function(e) NULL
  )
}

#' @noRd
plot_heterogeneity_density <- function(values,
                                       panel = "tau",
                                       colour = "black",
                                       n = 512L,
                                       weight = NULL) {
  d <- plot_heterogeneity_safe_density(values, n = n, weight = weight)
  if (is.null(d)) {
    return(data.frame(
      x = numeric(), y = numeric(), panel = character(), colour = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    x = d$x,
    y = d$y,
    panel = panel,
    colour = colour,
    stringsAsFactors = FALSE
  )
}

#' @noRd
plot_heterogeneity_annotations <- function(values,
                                           panel = "tau",
                                           probs = c(0.5, 0.75, 0.9),
                                           digits = 2,
                                           colour = "black",
                                           weight = NULL) {
  values <- as.numeric(values)
  ok <- is.finite(values) & values >= 0
  values <- values[ok]
  empty <- data.frame(
    x = numeric(), y = numeric(), y_lab = numeric(), label = character(),
    panel = character(), colour = character(), vline = logical(),
    stringsAsFactors = FALSE
  )
  if (!length(values)) {
    return(empty)
  }
  if (!is.null(weight)) {
    weight <- weight[ok]
  }
  qs <- if (!is.null(weight) && length(weight) == length(values)) {
    unname(interactmeta_weighted_quantile(values, weight, probs))
  } else {
    stats::quantile(values, probs = probs, names = FALSE, type = 7)
  }
  keep <- !duplicated(round(qs, digits + 4L))
  qs <- qs[keep]
  if (!length(qs)) {
    return(empty)
  }
  d <- plot_heterogeneity_safe_density(values, n = 512L, weight = weight)
  y_at <- if (is.null(d)) {
    rep(1, length(qs))
  } else {
    vapply(qs, function(q) {
      stats::approx(d$x, d$y, xout = q, rule = 2)$y
    }, numeric(1))
  }
  data.frame(
    x = as.numeric(qs),
    y = y_at,
    y_lab = y_at,
    label = format(round(qs, digits), nsmall = digits, trim = TRUE),
    panel = panel,
    colour = colour,
    vline = c(TRUE, rep(FALSE, max(0L, length(qs) - 1L))),
    stringsAsFactors = FALSE
  )
}
