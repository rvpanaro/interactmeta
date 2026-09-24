#' Single-study MAP predictive reproducibility check
#'
#' Directional posterior predictive assessment of whether a **target** study's
#' overall effect and treatment-by-subgroup interaction are compatible with a
#' replication of a **source-only** CAMS fit.
#'
#' The primary procedure uses the **same-information-fraction** replication
#' assumption: the intended replication has
#' \eqn{\pi_{\mathrm{target}}=\pi_{\mathrm{source}}}, so
#' \eqn{\beta(\pi_{\mathrm{target}}-\pi_{\mathrm{source}})=0}. Neither
#' \eqn{\beta} nor study-specific information fractions enter the primary test.
#' If the observed target has a different information fraction, any resulting
#' composition difference is counted as part of marginal non-reproducibility.
#'
#' Parameters \eqn{\tau_m} and \eqn{\tau_g} must come from a fit that does
#' **not** include the target study. A joint fit that includes the target is
#' rejected with an informative error.
#'
#' @section Overall effect (primary):
#' \deqn{D_m=m_{\mathrm{target}}-m_{\mathrm{source}},}
#' \deqn{Q_m=\frac{D_m^2}{s_{m,\mathrm{source}}^2+s_{m,\mathrm{target}}^2+2\tau_m^2}.}
#'
#' @section Interaction:
#' \deqn{D_g=g_{\mathrm{target}}-g_{\mathrm{source}},}
#' \deqn{Q_g=\frac{D_g^2}{s_{g,\mathrm{source}}^2+s_{g,\mathrm{target}}^2+2\tau_g^2}.}
#'
#' Draw-level \eqn{Q} statistics are compared to a \eqn{\chi^2_1} reference;
#' reported p-values average the chi-square tails over the posterior
#' (\eqn{p_m}, \eqn{p_g}). Reject when \eqn{p < } `level`.
#'
#' @param fit_source A [`interactmeta`] from [cam()] fitted on the **source
#'   study only**.
#' @param target_data A data frame (or [`metafor::escalc`] / [ifrac()] result)
#'   containing the **target** study's subgroup-level `yi` / `vi`. If `ifrac_*`
#'   design columns are missing, they are built automatically via [ifrac()]
#'   when needed (e.g. for `method = "composition_decomposition"`).
#' @param source Character ID of the source study. Defaults to the single trial
#'   present in `fit_source`.
#' @param target Character ID of the target study. Defaults to the single trial
#'   present in `target_data` when unambiguous.
#' @param level Significance level for the reject flags (default `0.05`).
#' @param method `"same_information_fraction"` (default; primary test) or
#'   `"composition_decomposition"` (optional diagnostic that adjusts the overall
#'   difference by \eqn{\beta(\pi_{\mathrm{target}}-\pi_{\mathrm{source}})}).
#' @param beta_draws Optional numeric vector of composition
#'   slopes \eqn{\beta}, used only when
#'   `method = "composition_decomposition"`. When `NULL` in that mode, uses
#'   rows from `fit_source$draws` (typically prior-driven in a one-study fit)
#'   and issues a warning.
#' @param ifrac Optional named numeric vector with elements for `source` and/or
#'   `target`, used only when `method = "composition_decomposition"`.
#'
#' @return A `interactmeta_reproducibility` list with elements `source`, `target`,
#'   `assumption`, `method`, `level`, `overall`, and `interaction`.
#'   Parameter-dependent `q` objects are numeric vectors aligned to
#'   `$draws` rows.
#'
#' @seealso [cam()], [observed_effects()]
#' @examples
#' \dontrun{
#' library(interactmeta)
#' es <- escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
#' fit_source <- bayescam(
#'   subset(es, trial == "DEFINE"),
#'   reference = "Younger"
#' )
#' reproducibility(
#'   fit_source,
#'   target_data = subset(es, trial == "FREEDOMS"),
#'   source = "DEFINE",
#'   target = "FREEDOMS",
#'   level = 0.05
#' )
#' }
#' @export
reproducibility <- function(fit_source,
                            target_data,
                            source = NULL,
                            target = NULL,
                            level = 0.05,
                            method = c(
                              "same_information_fraction",
                              "composition_decomposition"
                            ),
                            beta_draws = NULL,
                            ifrac = NULL) {
  check_comparison_fit(fit_source)

  method <- match.arg(method)
  level <- comparison_reproducibility_level(level)
  st <- attr(fit_source, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% "trial"
  subgroup_var <- st$subgroup_var %||% "subgroup"
  ref <- as.character(st$reference %||% "")[[1L]]
  if (!nzchar(ref)) {
    stop("Fit metadata is missing the reference subgroup.", call. = FALSE)
  }

  source <- comparison_reproducibility_resolve_source(
    fit_source = fit_source,
    source = source,
    trial_var = trial_var
  )
  comparison_reproducibility_check_source_only(
    fit_source = fit_source,
    source = source,
    target = target,
    trial_var = trial_var
  )

  target_dat <- comparison_reproducibility_prepare_target_data(
    target_data = target_data,
    fit_source = fit_source,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    need_ifrac = identical(method, "composition_decomposition")
  )
  target <- comparison_reproducibility_resolve_target(
    target_data = target_dat,
    target = target,
    source = source,
    trial_var = trial_var
  )

  parts <- comparison_reproducibility_extract_source(
    fit_source = fit_source,
    source = source,
    method = method,
    beta_draws = beta_draws
  )

  target_obs <- comparison_reproducibility_observed_from_data(
    data = target_dat,
    trial = target,
    reference = parts$reference,
    nonref = parts$nonref,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    ifrac_var = if (identical(method, "composition_decomposition")) {
      parts$ifrac_var
    } else {
      NULL
    }
  )

  if (identical(method, "composition_decomposition")) {
    pi <- comparison_reproducibility_pi(
      source = source,
      target = target,
      pi_source = parts$pi_source,
      pi_target = target_obs$pi,
      ifrac = ifrac
    )
    out <- comparison_reproducibility_compute(
      m_source = parts$m_source,
      se_m_source = parts$se_m_source,
      m_target = target_obs$m,
      se_m_target = target_obs$se_m,
      g_source = parts$g_source,
      se_g_source = parts$se_g_source,
      g_target = target_obs$g,
      se_g_target = target_obs$se_g,
      tau_m = parts$tau_m,
      tau_g = parts$tau_g,
      source = source,
      target = target,
      level = level,
      method = method,
      pi_source = pi[[source]],
      pi_target = pi[[target]],
      beta = parts$beta
    )
    out$pi <- pi
  } else {
    out <- comparison_reproducibility_compute(
      m_source = parts$m_source,
      se_m_source = parts$se_m_source,
      m_target = target_obs$m,
      se_m_target = target_obs$se_m,
      g_source = parts$g_source,
      se_g_source = parts$se_g_source,
      g_target = target_obs$g,
      se_g_target = target_obs$se_g,
      tau_m = parts$tau_m,
      tau_g = parts$tau_g,
      source = source,
      target = target,
      level = level,
      method = method
    )
  }

  out$reference <- parts$reference
  out$nonref <- parts$nonref
  out$ifrac_var <- parts$ifrac_var
  class(out) <- c("interactmeta_reproducibility", "list")
  out
}

#' @rdname reproducibility
#' @param x A `interactmeta_reproducibility` object.
#' @param digits Decimal places for printed p-values (default `2`).
#' @param ... Ignored.
#' @export
print.interactmeta_reproducibility <- function(x, digits = 2, ...) {
  cat(comparison_reproducibility_text(x, digits = digits), sep = "")
  invisible(x)
}

#' @rdname reproducibility
#' @param object A `interactmeta_reproducibility` object.
#' @export
summary.interactmeta_reproducibility <- function(object, digits = 2, ...) {
  overall <- list(
    observed_difference = object$overall$observed_difference,
    q_median = stats::median(as.numeric(object$overall$q)),
    p_value = object$overall$p_value,
    reject = object$overall$reject
  )
  if (!is.null(object$overall$composition_adjustment)) {
    overall$composition_adjustment_median <- stats::median(
      as.numeric(object$overall$composition_adjustment)
    )
    overall$adjusted_difference_median <- stats::median(
      as.numeric(object$overall$adjusted_difference)
    )
  }
  out <- list(
    source = object$source,
    target = object$target,
    assumption = object$assumption,
    method = object$method,
    level = object$level,
    overall = overall,
    interaction = list(
      observed_difference = object$interaction$observed_difference,
      q_median = stats::median(as.numeric(object$interaction$q)),
      p_value = object$interaction$p_value,
      reject = object$interaction$reject
    ),
    digits = as.integer(digits)[[1L]],
    text = comparison_reproducibility_text(object, digits = digits)
  )
  class(out) <- c("summary.interactmeta_reproducibility", "list")
  out
}

#' @rdname reproducibility
#' @export
print.summary.interactmeta_reproducibility <- function(x, ...) {
  cat(x$text, sep = "")
  invisible(x)
}

# ---- internals --------------------------------------------------------------

#' @noRd
comparison_reproducibility_level <- function(level) {
  level <- as.numeric(level)[[1L]]
  if (!is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a probability in (0, 1).", call. = FALSE)
  }
  level
}

#' @noRd
comparison_reproducibility_study_id <- function(x, what) {
  x <- as.character(x)
  if (length(x) != 1L || !nzchar(x) || is.na(x)) {
    stop("`", what, "` must be a single non-missing study name.", call. = FALSE)
  }
  x
}

#' @noRd
comparison_reproducibility_fit_trials <- function(fit_source, trial_var) {
  dat <- fit_source$data
  if (is.null(dat) || !nrow(dat)) {
    stop("Fit has no model data; refit with `cam()`.", call. = FALSE)
  }
  if (!trial_var %in% names(dat)) {
    stop("Fit data is missing trial column `", trial_var, "`.", call. = FALSE)
  }
  unique(as.character(dat[[trial_var]]))
}

#' @noRd
comparison_reproducibility_resolve_source <- function(fit_source, source, trial_var) {
  avail <- comparison_reproducibility_fit_trials(fit_source, trial_var)
  if (is.null(source)) {
    if (length(avail) != 1L) {
      stop(
        "`source` must be supplied when `fit_source` contains more than one study.",
        call. = FALSE
      )
    }
    return(avail[[1L]])
  }
  source <- comparison_reproducibility_study_id(source, "source")
  if (!source %in% avail) {
    stop(
      "Source study `",
      source,
      "` not found in `fit_source`. Available: ",
      paste(avail, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  source
}

#' Ensure fit_source is source-only (target must not inform the predictive distribution)
#' @noRd
comparison_reproducibility_check_source_only <- function(fit_source,
                                                  source,
                                                  target,
                                                  trial_var) {
  avail <- comparison_reproducibility_fit_trials(fit_source, trial_var)
  target_chr <- if (!is.null(target)) {
    comparison_reproducibility_study_id(target, "target")
  } else {
    NA_character_
  }

  if (!is.na(target_chr) && target_chr %in% avail) {
    stop(
      "The reproducibility check requires a CAMS fit based on the source study only.\n",
      "The target study must not contribute to the posterior predictive distribution.",
      call. = FALSE
    )
  }

  if (length(avail) != 1L || !identical(avail[[1L]], source)) {
    stop(
      "The reproducibility check requires a CAMS fit based on the source study only.\n",
      "The target study must not contribute to the posterior predictive distribution.\n",
      "Found studies in `fit_source`: ",
      paste(avail, collapse = ", "),
      "; expected only `",
      source,
      "`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @noRd
comparison_reproducibility_prepare_target_data <- function(target_data,
                                                    fit_source,
                                                    trial_var,
                                                    subgroup_var,
                                                    need_ifrac = FALSE) {
  if (inherits(target_data, "interactmeta")) {
    stop(
      "`target_data` must be raw study data, not a CAMS fit. ",
      "Pass the target study's `yi` / `vi` data frame (or escalc / ifrac result).",
      call. = FALSE
    )
  }
  st <- attr(fit_source, "interactmeta", exact = TRUE)
  ref <- as.character(st$reference)[[1L]]
  ifrac_var <- st$ifrac_var %||% st$ifrac_vars[[1L]]
  ifrac_var <- if (!is.null(ifrac_var)) as.character(ifrac_var)[[1L]] else NA_character_

  has_core <- all(c("yi", "vi", trial_var, subgroup_var) %in% names(target_data)) ||
    all(c("yi", "vi", "trial", "subgroup") %in% names(target_data))
  needs_ifrac_cols <- isTRUE(need_ifrac) &&
    (is.na(ifrac_var) || !ifrac_var %in% names(target_data) || !is_ifrac_output(target_data))

  if (!has_core || needs_ifrac_cols) {
    target_data <- prepare_ifrac_input(
      target_data,
      reference = ref,
      trial = if (trial_var %in% names(target_data)) trial_var else "trial",
      subgroup = if (subgroup_var %in% names(target_data)) {
        subgroup_var
      } else {
        "subgroup"
      }
    )
  } else {
    target_data <- as.data.frame(target_data)
  }

  if (!trial_var %in% names(target_data) && "trial" %in% names(target_data)) {
    trial_var <- "trial"
  }
  if (!subgroup_var %in% names(target_data) && "subgroup" %in% names(target_data)) {
    subgroup_var <- "subgroup"
  }
  need <- c("yi", "vi", trial_var, subgroup_var)
  miss <- setdiff(need, names(target_data))
  if (length(miss)) {
    stop(
      "`target_data` must include columns: ",
      paste(need, collapse = ", "),
      ". Missing: ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (subgroup_var %in% names(target_data)) {
    keep <- is_ifrac_component_row(target_data[[subgroup_var]])
    if (any(!keep)) {
      target_data <- target_data[keep, , drop = FALSE]
    }
  }
  target_data
}

#' @noRd
comparison_reproducibility_resolve_target <- function(target_data,
                                               target,
                                               source,
                                               trial_var) {
  avail <- unique(as.character(target_data[[trial_var]]))
  avail <- avail[nzchar(avail) & !is.na(avail)]
  if (!length(avail)) {
    stop("`target_data` contains no study IDs.", call. = FALSE)
  }

  if (is.null(target)) {
    cand <- setdiff(avail, source)
    if (length(avail) == 1L) {
      target <- avail[[1L]]
    } else if (length(cand) == 1L) {
      target <- cand[[1L]]
    } else {
      stop(
        "`target` must be supplied when `target_data` contains multiple studies.",
        call. = FALSE
      )
    }
  } else {
    target <- comparison_reproducibility_study_id(target, "target")
  }

  if (identical(target, source)) {
    stop("`source` and `target` must name distinct studies.", call. = FALSE)
  }
  if (!target %in% avail) {
    stop(
      "Target study `",
      target,
      "` not found in `target_data`. Available: ",
      paste(avail, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  target
}

#' Source observations + predictive parameters from the source-only fit
#' @noRd
comparison_reproducibility_extract_source <- function(fit_source,
                                               source,
                                               method = "same_information_fraction",
                                               beta_draws = NULL) {
  st <- attr(fit_source, "interactmeta", exact = TRUE)
  ref <- as.character(st$reference)[[1L]]
  ifrac_var <- st$ifrac_var %||% st$ifrac_vars[[1L]]
  if (!is.null(ifrac_var) && nzchar(as.character(ifrac_var)[[1L]])) {
    ifrac_var <- as.character(ifrac_var)[[1L]]
  } else {
    ifrac_var <- NULL
  }
  trial_var <- st$trial_var %||% "trial"
  post <- as.data.frame(fit_source$draws)

  ov <- observed_effects(fit_source, type = "overall")
  int <- observed_effects(fit_source, type = "interaction")
  if (!source %in% names(ov$estimate) || !source %in% names(int$estimate)) {
    stop(
      "Observed effects for source study `",
      source,
      "` are missing from `fit_source`.",
      call. = FALSE
    )
  }
  m_source <- unname(ov$estimate[[source]])
  se_m_source <- sqrt(pmax(unname(ov$variance[[source]]), 0))
  g_source <- unname(int$estimate[[source]])
  se_g_source <- sqrt(pmax(unname(int$variance[[source]]), 0))
  if (!is.finite(m_source) || !is.finite(se_m_source)) {
    stop("Observed overall effect / SE is missing for the source study.", call. = FALSE)
  }
  if (!is.finite(g_source) || !is.finite(se_g_source)) {
    stop("Observed interaction effect / SE is missing for the source study.", call. = FALSE)
  }

  if (!"sd_trial__Intercept" %in% names(post)) {
    stop("Fit `$draws` are missing overall-block heterogeneity (`sd_trial__Intercept`).", call. = FALSE)
  }
  slope <- grep("^sd_trial__ifrac2_", names(post), value = TRUE)
  if (!length(slope)) {
    stop("Fit `$draws` are missing interaction-block heterogeneity (`sd_trial__ifrac2_*`).", call. = FALSE)
  }
  tau_m <- as.numeric(post$sd_trial__Intercept)
  tau_g <- as.numeric(post[[slope[[1L]]]])
  nonref <- comparison_reproducibility_nonref(fit_source = fit_source)

  out <- list(
    m_source = m_source,
    se_m_source = se_m_source,
    g_source = g_source,
    se_g_source = se_g_source,
    tau_m = tau_m,
    tau_g = tau_g,
    reference = ref,
    nonref = nonref,
    ifrac_var = ifrac_var,
    pi_source = NA_real_,
    beta = NULL
  )

  if (identical(method, "composition_decomposition")) {
    if (is.null(ifrac_var)) {
      stop("Fit metadata is missing the information-fraction column name.", call. = FALSE)
    }
    out$pi_source <- comparison_trial_information_fractions(
      data = fit_source$data,
      trial_var = trial_var,
      ifrac_var = ifrac_var,
      trials = source
    )[[source]]
    out$beta <- comparison_reproducibility_beta(
      beta_draws = beta_draws,
      post = post,
      nonref = nonref,
      n_tau = length(tau_m)
    )
  }
  out
}

#' @noRd
comparison_reproducibility_nonref <- function(fit_source) {
  st <- attr(fit_source, "interactmeta", exact = TRUE)
  ref <- as.character(st$reference)[[1L]]
  sg <- unique(as.character(fit_source$data[[st$subgroup_var %||% "subgroup"]]))
  nonrefs <- setdiff(sg, ref)
  if (!length(nonrefs)) {
    stop("No non-reference subgroup found in `fit_source`.", call. = FALSE)
  }
  if (length(nonrefs) > 1L) {
    stop(
      "reproducibility() currently supports a single non-reference subgroup; ",
      "found: ",
      paste(nonrefs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  nonrefs[[1L]]
}

#' Resolve beta draws (external or from source-only fit)
#' @noRd
comparison_reproducibility_beta <- function(beta_draws, post, nonref, n_tau) {
  n_tau <- as.integer(n_tau)[[1L]]
  if (is.null(beta_draws)) {
    col <- paste0("b_ifrac_", nonref)
    if (!col %in% names(post)) {
      alt <- grep("^b_ifrac_", names(post), value = TRUE)
      if (length(alt) == 1L) {
        col <- alt[[1L]]
      } else {
        stop(
          "Fit `$draws` are missing the composition slope (`b_ifrac_*`). ",
          "Supply `beta_draws` explicitly.",
          call. = FALSE
        )
      }
    }
    warning(
      "With a source-only fit, the composition slope (`beta`) is generally not ",
      "identified from between-study composition variation. Values of `beta` ",
      "from `fit_source` are typically prior-driven. Supply `beta_draws` from ",
      "an externally justified prior or estimate if available.",
      call. = FALSE
    )
    beta <- as.numeric(post[[col]])
  } else if (is.numeric(beta_draws)) {
    beta <- as.numeric(beta_draws)
  } else {
    stop("`beta_draws` must be NULL or a numeric vector.", call. = FALSE)
  }

  n_beta <- length(beta)
  if (n_beta == 1L && n_tau > 1L) {
    beta <- rep(beta, n_tau)
  } else if (n_beta != n_tau) {
    stop(
      "`beta_draws` has ",
      n_beta,
      " value(s) but `fit_source` has ",
      n_tau,
      " row(s) in `$draws`. Lengths must match (or supply a single value).",
      call. = FALSE
    )
  }
  beta
}

#' Observed overall / interaction / IF for one trial from raw data
#' @noRd
comparison_reproducibility_observed_from_data <- function(data,
                                                   trial,
                                                   reference,
                                                   nonref,
                                                   trial_var,
                                                   subgroup_var,
                                                   ifrac_var = NULL) {
  d <- data[as.character(data[[trial_var]]) == trial, , drop = FALSE]
  if (!nrow(d)) {
    stop("No rows for target study `", trial, "` in `target_data`.", call. = FALSE)
  }
  blocks <- comparison_trial_effect_blocks(
    data = d,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = reference,
    nonref = nonref,
    ifrac_var = ifrac_var
  )
  m <- unname(blocks$overall$yi[[1L]])
  se_m <- sqrt(pmax(unname(blocks$overall$vi[[1L]]), 0))
  g <- unname(blocks$interaction$yi[[1L]])
  se_g <- sqrt(pmax(unname(blocks$interaction$vi[[1L]]), 0))
  if (!is.finite(m) || !is.finite(se_m)) {
    stop("Observed overall effect / SE is missing for the target study.", call. = FALSE)
  }
  if (!is.finite(g) || !is.finite(se_g)) {
    stop("Observed interaction effect / SE is missing for the target study.", call. = FALSE)
  }

  pi <- NA_real_
  if (!is.null(ifrac_var) && ifrac_var %in% names(d)) {
    pi <- unname(comparison_trial_information_fractions(
      data = d,
      trial_var = trial_var,
      ifrac_var = ifrac_var,
      trials = trial
    )[[trial]])
  } else if (!is.null(blocks$pi) && length(blocks$pi) && is.finite(blocks$pi[[trial]])) {
    pi <- unname(blocks$pi[[trial]])
  }

  list(m = m, se_m = se_m, g = g, se_g = se_g, pi = pi)
}

#' Resolve source/target information fractions (composition diagnostic only)
#' @noRd
comparison_reproducibility_pi <- function(source,
                                    target,
                                    pi_source,
                                    pi_target,
                                    ifrac = NULL) {
  pi <- stats::setNames(c(as.numeric(pi_source), as.numeric(pi_target)), c(source, target))

  if (!is.null(ifrac)) {
    nm <- names(ifrac)
    ifrac <- as.numeric(ifrac)
    if (!length(ifrac) || any(!is.finite(ifrac))) {
      stop("`ifrac` must contain finite information fractions.", call. = FALSE)
    }
    if (is.null(nm) || !all(nzchar(nm))) {
      stop(
        "`ifrac` must be a named numeric vector with names matching `source` ",
        "and/or `target`.",
        call. = FALSE
      )
    }
    names(ifrac) <- nm
    extra <- setdiff(names(ifrac), c(source, target))
    if (length(extra)) {
      stop(
        "`ifrac` names must match `source` / `target`. Unknown: ",
        paste(extra, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    pi[names(ifrac)] <- ifrac
  }

  if (any(!is.finite(pi))) {
    miss <- names(pi)[!is.finite(pi)]
    stop(
      "Could not retrieve finite information fractions for ",
      paste(miss, collapse = ", "),
      ". Ensure `target_data` includes sample sizes (`n`) or `ifrac_*` ",
      "design columns, or supply `ifrac` explicitly ",
      "(e.g. `ifrac = c(",
      source,
      " = 0.3, ",
      target,
      " = 0.5)`).",
      call. = FALSE
    )
  }
  pi
}

#' Core MAP predictive reproducibility calculations (testable without a fit)
#' @noRd
comparison_reproducibility_compute <- function(m_source,
                                         se_m_source,
                                         m_target,
                                         se_m_target,
                                         g_source,
                                         se_g_source,
                                         g_target,
                                         se_g_target,
                                         tau_m,
                                         tau_g,
                                         source,
                                         target,
                                         level = 0.05,
                                         method = c(
                                           "same_information_fraction",
                                           "composition_decomposition"
                                         ),
                                         pi_source = NULL,
                                         pi_target = NULL,
                                         beta = NULL) {
  method <- match.arg(method)
  source <- comparison_reproducibility_study_id(source, "source")
  target <- comparison_reproducibility_study_id(target, "target")
  if (identical(source, target)) {
    stop("`source` and `target` must name distinct studies.", call. = FALSE)
  }
  level <- comparison_reproducibility_level(level)

  m_source <- as.numeric(m_source)[[1L]]
  m_target <- as.numeric(m_target)[[1L]]
  se_m_source <- as.numeric(se_m_source)[[1L]]
  se_m_target <- as.numeric(se_m_target)[[1L]]
  g_source <- as.numeric(g_source)[[1L]]
  g_target <- as.numeric(g_target)[[1L]]
  se_g_source <- as.numeric(se_g_source)[[1L]]
  se_g_target <- as.numeric(se_g_target)[[1L]]

  delta_m <- m_target - m_source
  tau_m <- as.numeric(tau_m)
  tau_g <- as.numeric(tau_g)
  if (!is.null(beta)) {
    beta <- as.numeric(beta)
  }
  var_m <- se_m_source^2 + se_m_target^2 + 2 * tau_m^2

  if (identical(method, "same_information_fraction")) {
    q_m <- (delta_m^2) / var_m
    composition_adjustment <- NULL
    adjusted_difference <- NULL
    assumption <- "same information fraction"
  } else {
    if (is.null(beta) || is.null(pi_source) || is.null(pi_target)) {
      stop(
        "`composition_decomposition` requires `beta`, `pi_source`, and `pi_target`.",
        call. = FALSE
      )
    }
    pi_source <- as.numeric(pi_source)[[1L]]
    pi_target <- as.numeric(pi_target)[[1L]]
    composition_adjustment <- beta * (pi_target - pi_source)
    adjusted_difference <- delta_m - composition_adjustment
    q_m <- adjusted_difference^2 / var_m
    assumption <- "composition decomposition"
  }

  p_m_draw <- as.numeric(
    stats::pchisq(
      as.numeric(q_m),
      df = 1,
      lower.tail = FALSE
    )
  )
  p_m <- mean(p_m_draw)

  delta_g <- g_target - g_source
  var_g <- se_g_source^2 + se_g_target^2 + 2 * tau_g^2
  q_g <- (delta_g^2) / var_g
  p_g_draw <- as.numeric(
    stats::pchisq(
      as.numeric(q_g),
      df = 1,
      lower.tail = FALSE
    )
  )
  p_g <- mean(p_g_draw)

  overall <- list(
    observed_difference = delta_m,
    q = q_m,
    p_value_draws = p_m_draw,
    p_value = p_m,
    reject = isTRUE(p_m < level),
    sampling_var_sum = se_m_source^2 + se_m_target^2,
    se = stats::setNames(c(se_m_source, se_m_target), c(source, target)),
    estimate = stats::setNames(c(m_source, m_target), c(source, target))
  )
  if (!is.null(composition_adjustment)) {
    overall$composition_adjustment <- composition_adjustment
    overall$adjusted_difference <- adjusted_difference
  }

  list(
    source = source,
    target = target,
    assumption = assumption,
    method = method,
    level = level,
    overall = overall,
    interaction = list(
      observed_difference = delta_g,
      q = q_g,
      p_value_draws = p_g_draw,
      p_value = p_g,
      reject = isTRUE(p_g < level),
      sampling_var_sum = se_g_source^2 + se_g_target^2,
      se = stats::setNames(c(se_g_source, se_g_target), c(source, target)),
      estimate = stats::setNames(c(g_source, g_target), c(source, target))
    )
  )
}

#' Narrative print text for reproducibility results
#' @noRd
comparison_reproducibility_text <- function(x, digits = 2) {
  fmt_p <- function(p) {
    format(round(as.numeric(p), digits), nsmall = digits, trim = TRUE)
  }
  msg_for <- function(reject) {
    if (isTRUE(reject)) {
      "Evidence of excess non-reproducibility.\n"
    } else {
      "No evidence of excess non-reproducibility.\n"
    }
  }

  overall_label <- if (identical(x$method, "composition_decomposition")) {
    "Overall-effect reproducibility (composition-adjusted MAP PPP)"
  } else {
    "Overall-effect reproducibility (same-IF MAP PPP)"
  }

  paste0(
    "CAMS single-study MAP predictive reproducibility (",
    x$source,
    " -> ",
    x$target,
    ")\n",
    "Assumption: ",
    x$assumption %||% "same information fraction",
    "\n",
    "\n",
    overall_label,
    "\n",
    "p_m = ",
    fmt_p(x$overall$p_value),
    "\n",
    msg_for(x$overall$reject),
    "\n",
    "Interaction reproducibility (MAP PPP)\n",
    "p_g = ",
    fmt_p(x$interaction$p_value),
    "\n",
    msg_for(x$interaction$reject)
  )
}
