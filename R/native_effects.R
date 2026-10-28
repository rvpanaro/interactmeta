#' Engine-native summaries for default bayesmeta marginal effects
#'
#' When the estimand is a single `bmr` / `bayesmeta` margin (combined
#' \eqn{M(\pi)} at one \eqn{\pi}, interaction \eqn{\gamma}), use stored
#' `$pred.interval` / `$post.interval` instead of re-summarising the full
#' DIRECT product grid.
#' @noRd
interactmeta_native_marginal_eligible <- function(fit,
                                                  quantity,
                                                  at,
                                                  ifrac,
                                                  type,
                                                  prob,
                                                  max_support,
                                                  summary,
                                                  exponentiate,
                                                  shrinkage) {
  if (!isTRUE(summary)) {
    return(FALSE)
  }
  if (isTRUE(shrinkage) || isTRUE(exponentiate)) {
    return(FALSE)
  }
  if (!is.null(max_support)) {
    return(FALSE)
  }
  if (interactmeta_is_bayesian_fit(fit)) {
    bm <- fit$bayesmeta
    if (is.null(bm$overall) || is.null(bm$interaction)) {
      return(FALSE)
    }
  } else {
    mf <- fit$metafor
    if (is.null(mf$overall) || is.null(mf$interaction)) {
      return(FALSE)
    }
  }
  type <- as.character(type)
  if (length(type) != 1L || !identical(type[[1L]], "Overall")) {
    return(FALSE)
  }
  prob <- as.numeric(prob)[[1L]]
  if (!is.finite(prob) || prob <= 0 || prob >= 1) {
    return(FALSE)
  }
  if (is.numeric(at) && length(at) > 1L) {
    return(FALSE)
  }
  if (identical(quantity, "interactions")) {
    return(TRUE)
  }
  if (identical(quantity, "combined")) {
    omega <- interactmeta_native_marginal_omega(fit, at, ifrac)
    return(is.finite(omega))
  }
  FALSE
}

#' @noRd
interactmeta_native_marginal_omega <- function(fit, at, ifrac) {
  ifrac <- match.arg(ifrac, c("interaction", "overall"))
  if (is.null(at)) {
    return(as.numeric(interactmeta_eval_ifrac(fit, ifrac = ifrac))[[1L]])
  }
  if (inherits(at, "interactmeta_overall_ifrac")) {
    summ <- attr(at, "summary", exact = TRUE)
    if (is.data.frame(summ) && nrow(summ)) {
      return(as.numeric(summ$estimate[[1L]])[[1L]])
    }
  }
  as.numeric(at)[[1L]]
}

#' @noRd
interactmeta_native_marginal_interval_type <- function(fit) {
  interactmeta_match_interval_type(
    fit$interval.type %||% fit$bayesmeta$overall$interval.type
  )
}

#' @noRd
interactmeta_try_native_combined <- function(fit,
                                             at,
                                             ifrac,
                                             prob,
                                             contrast_sep) {
  if (!interactmeta_native_marginal_eligible(
    fit,
    quantity = "combined",
    at = at,
    ifrac = ifrac,
    type = "Overall",
    prob = prob,
    max_support = NULL,
    summary = TRUE,
    exponentiate = FALSE,
    shrinkage = FALSE
  )) {
    return(NULL)
  }
  if (!interactmeta_is_bayesian_fit(fit)) {
    return(interactmeta_try_native_combined_metafor(
      fit, at = at, ifrac = ifrac, prob = prob, contrast_sep = contrast_sep
    ))
  }
  ifrac <- match.arg(ifrac, c("interaction", "overall"))
  source <- if (is.null(at)) {
    "fit"
  } else {
    "user"
  }
  oi <- if (is.null(at)) {
    interactmeta_eval_ifrac(fit, ifrac = ifrac)
  } else {
    interactmeta_native_marginal_omega(fit, at, ifrac)
  }
  omega <- as.numeric(oi)[[1L]]
  it <- interactmeta_native_marginal_interval_type(fit)
  loc <- interactmeta_bayesmeta_population_location(
    fit,
    pi = omega,
    gamma_scale = 0,
    re_scale = 0,
    predictive = FALSE,
    prob = prob,
    interval.type = it
  )
  gm <- interactmeta_bmr_population_linpred_gm(fit$bayesmeta$overall, omega)
  design <- comparison_marginal_design(fit, contrast_sep = contrast_sep)
  summary <- data.frame(
    combined = comparison_combined_population_label(),
    type = "Overall",
    estimate = loc$median,
    mean = interactmeta_gm_mean(gm),
    sd = sqrt(max(interactmeta_gm_var(gm), 0)),
    lci = loc$lower,
    uci = loc$upper,
    prob = prob,
    stringsAsFactors = FALSE
  )
  attr(summary, "interval.type") <- it
  out <- list(
    support = data.frame(
      type = character(),
      effect = numeric(),
      stringsAsFactors = FALSE
    ),
    overall_ifrac = oi,
    overall_ifrac_source = source,
    ifrac = ifrac,
    reference = design$reference,
    prob = prob,
    exponentiate = FALSE,
    shrinkage = FALSE,
    summary = summary
  )
  class(out) <- c("interactmeta_combined_effect", class(out))
  interactmeta_attach_origin(out, fit, quantity = "combined")
}

#' @noRd
interactmeta_try_native_interactions <- function(fit,
                                                 at,
                                                 ifrac,
                                                 prob,
                                                 contrast_sep) {
  if (!interactmeta_native_marginal_eligible(
    fit,
    quantity = "interactions",
    at = at,
    ifrac = ifrac,
    type = "Overall",
    prob = prob,
    max_support = NULL,
    summary = TRUE,
    exponentiate = FALSE,
    shrinkage = FALSE
  )) {
    return(NULL)
  }
  if (!interactmeta_is_bayesian_fit(fit)) {
    return(interactmeta_try_native_interactions_metafor(
      fit, at = at, ifrac = ifrac, prob = prob, contrast_sep = contrast_sep
    ))
  }
  ifrac <- match.arg(ifrac, c("interaction", "overall"))
  it <- interactmeta_native_marginal_interval_type(fit)
  int <- fit$bayesmeta$interaction
  loc <- interactmeta_native_location(
    int,
    prob = prob,
    interval.type = it,
    margin = "effect"
  )
  gm <- interactmeta_gm_from_bayesmeta(int)
  int_design <- comparison_interaction_design(fit, contrast_sep = contrast_sep)
  if (!length(int_design$interaction_levels)) {
    return(NULL)
  }
  lbl <- int_design$interaction_levels[[1L]]
  source <- if (is.null(at)) {
    "fit"
  } else {
    "user"
  }
  oi <- if (is.null(at)) {
    interactmeta_eval_ifrac(fit, ifrac = ifrac)
  } else {
    interactmeta_native_marginal_omega(fit, at, ifrac)
  }
  summary <- data.frame(
    interaction = lbl,
    type = "Overall",
    estimate = loc$median,
    mean = interactmeta_gm_mean(gm),
    sd = sqrt(max(interactmeta_gm_var(gm), 0)),
    lci = loc$lower,
    uci = loc$upper,
    prob = prob,
    stringsAsFactors = FALSE
  )
  attr(summary, "interval.type") <- it
  out <- list(
    support = data.frame(
      interaction = character(),
      type = character(),
      effect = numeric(),
      stringsAsFactors = FALSE
    ),
    overall_ifrac = oi,
    overall_ifrac_source = source,
    ifrac = ifrac,
    reference = int_design$reference,
    prob = prob,
    exponentiate = FALSE,
    shrinkage = FALSE,
    summary = summary
  )
  class(out) <- c("interactmeta_interaction_effects", class(out))
  interactmeta_attach_origin(out, fit, quantity = "interaction")
}

#' @noRd
interactmeta_metafor_wald_summary <- function(est, se, prob) {
  z <- stats::qnorm(1 - (1 - prob) / 2)
  data.frame(
    estimate = est,
    mean = est,
    sd = se,
    lci = est - z * se,
    uci = est + z * se,
    prob = prob,
    stringsAsFactors = FALSE
  )
}

#' @noRd
interactmeta_try_native_combined_metafor <- function(fit, at, ifrac, prob, contrast_sep) {
  ifrac <- match.arg(ifrac, c("interaction", "overall"))
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  alpha <- as.numeric(ov$beta[[1L]])
  theta <- as.numeric(ov$beta[[2L]])
  gamma <- as.numeric(int$beta[[1L]])
  beta <- theta - gamma
  omega <- interactmeta_native_marginal_omega(fit, at, ifrac)
  est <- alpha + omega * (beta + gamma)
  se <- sqrt(max(as.numeric(ov$se[[1L]])^2 + (omega^2) * (
    as.numeric(ov$se[[2L]])^2 + as.numeric(int$se[[1L]])^2
  ), 0))
  design <- comparison_marginal_design(fit, contrast_sep = contrast_sep)
  summary <- interactmeta_metafor_wald_summary(est, se, prob)
  summary <- cbind(
    data.frame(combined = comparison_combined_population_label(), type = "Overall"),
    summary
  )
  oi <- if (is.null(at)) interactmeta_eval_ifrac(fit, ifrac = ifrac) else at
  out <- list(
    support = data.frame(type = character(), effect = numeric(), stringsAsFactors = FALSE),
    overall_ifrac = oi,
    overall_ifrac_source = if (is.null(at)) "fit" else "user",
    ifrac = ifrac,
    reference = design$reference,
    prob = prob,
    exponentiate = FALSE,
    shrinkage = FALSE,
    summary = summary
  )
  class(out) <- c("interactmeta_combined_effect", class(out))
  interactmeta_attach_origin(out, fit, quantity = "combined")
}

#' @noRd
interactmeta_try_native_interactions_metafor <- function(fit, at, ifrac, prob, contrast_sep) {
  ifrac <- match.arg(ifrac, c("interaction", "overall"))
  int <- fit$metafor$interaction
  est <- as.numeric(int$beta[[1L]])
  se <- as.numeric(int$se[[1L]])
  int_design <- comparison_interaction_design(fit, contrast_sep = contrast_sep)
  lbl <- int_design$interaction_levels[[1L]]
  summary <- interactmeta_metafor_wald_summary(est, se, prob)
  summary <- cbind(data.frame(interaction = lbl, type = "Overall"), summary)
  oi <- if (is.null(at)) interactmeta_eval_ifrac(fit, ifrac = ifrac) else at
  out <- list(
    support = data.frame(
      interaction = character(),
      type = character(),
      effect = numeric(),
      stringsAsFactors = FALSE
    ),
    overall_ifrac = oi,
    overall_ifrac_source = if (is.null(at)) "fit" else "user",
    ifrac = ifrac,
    reference = int_design$reference,
    prob = prob,
    exponentiate = FALSE,
    shrinkage = FALSE,
    summary = summary
  )
  class(out) <- c("interactmeta_interaction_effects", class(out))
  interactmeta_attach_origin(out, fit, quantity = "interaction")
}
