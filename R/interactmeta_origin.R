#' Engine tagging and derived-effect helpers for one CAMS fit
#'
#' CAMS downstream quantities are transformations of one selected engine's
#' fit. SWADA / DA / AD are not CAMS fits and must not borrow IF quantities.
#' @noRd
NULL

#' Classes that are not CAMS fits
#' @noRd
interactmeta_non_cams_classes <- function() {
  c(
    "interactmeta_swada",
    "interactmeta_da",
    "interactmeta_ad",
    "interactmeta_conventional"
  )
}

#' Informative error when a CAMS/IF method is applied to SWADA / DA / AD
#' @noRd
interactmeta_abort_non_cams <- function(object, feature) {
  if (!interactmeta_is_non_cams_analysis(object)) {
    return(invisible(NULL))
  }
  interactmeta_require_capability(object, "ifrac", feature = feature)
}

#' Native backend object stored on a fit (no refit)
#'
#' For Bayesian CAMS (`engine = "bayesmeta"`), returns the list of stored
#' univariate fits: `$overall` ([bayesmeta::bmr()] of \eqn{m_j} on \eqn{\pi_j}),
#' `$interaction` ([bayesmeta::bayesmeta()] of \eqn{g_j}), and `$trial_blocks`.
#' Those objects are the `bayesmeta` / DIRECT backend, not CAMS estimands.
#' Plot them with the native `bayesmeta` methods
#' (`vignette("further")`). For frequentist CAMS the
#' analogue is `$metafor`.
#'
#' @param fit A [cam()] / [da()] / [ad()] / [swada()] / [conventional()] object.
#' @return For Bayesian CAMS, the list with `$overall`, `$interaction`, and
#'   `$trial_blocks`. For frequentist CAMS, `fit$metafor`. For DA / AD /
#'   SWADA / conventional, the stored `$fits` list. `NULL` if none is stored.
#' @seealso [traceplot.interactmeta()], `vignette("further")`
#' @keywords internal
#' @export
native_fit <- function(fit) {
  analysis <- interactmeta_fit_analysis(fit)
  if (!is.na(analysis) && analysis %in% c("da", "ad", "swada", "conventional")) {
    return(fit$fits)
  }
  check_comparison_fit(fit)
  eng <- as.character(fit$engine %||% "")[[1L]]
  if (identical(eng, "bayesmeta")) {
    return(fit$bayesmeta)
  }
  if (identical(eng, "metafor")) {
    return(fit$metafor)
  }
  if (identical(eng, "brms")) {
    return(fit$brmsfit)
  }
  NULL
}

#' Engine / quantity tag on a derived CAMS object
#' @noRd
interactmeta_origin <- function(x) {
  attr(x, "interactmeta_origin", exact = TRUE)
}

#' @noRd
interactmeta_attach_origin <- function(x,
                                      fit,
                                      quantity,
                                      fitted_or_derived = "derived",
                                      uncertainty = NULL) {
  eng <- as.character(fit$engine %||% NA_character_)[[1L]]
  meth <- as.character(fit$method %||% NA_character_)[[1L]]
  algo <- as.character(fit$algorithm %||% NA_character_)[[1L]]
  if (is.null(uncertainty)) {
    uncertainty <- interactmeta_default_uncertainty(eng, meth)
  }
  attr(x, "interactmeta_origin") <- list(
    fit_engine = eng,
    fit_method = meth,
    fit_algorithm = algo,
    quantity = quantity,
    fitted_or_derived = fitted_or_derived,
    uncertainty = uncertainty,
    call = fit$call
  )
  x
}

#' Inferential interval style for print / summary tables
#'
#' `"frequentist"` is `engine = "metafor"` only. Bayesian engines keep
#' `"bayesian"` even when `$draws` happens to have one row.
#' @noRd
interactmeta_interval_style <- function(x) {
  eng <- interactmeta_origin_engine(x)
  if (identical(eng, "metafor")) {
    "frequentist"
  } else {
    "bayesian"
  }
}

#' Engine stored on a fit or on a derived object's origin
#' @noRd
interactmeta_origin_engine <- function(x) {
  if (inherits(x, "interactmeta") && !is.null(x$engine)) {
    return(as.character(x$engine)[[1L]])
  }
  org <- attr(x, "interactmeta_origin", exact = TRUE)
  if (!is.null(org) && !is.null(org$fit_engine)) {
    return(as.character(org$fit_engine)[[1L]])
  }
  if (is.list(x) && !is.null(x$fit) && inherits(x$fit, "interactmeta")) {
    return(as.character(x$fit$engine)[[1L]])
  }
  if (is.list(x) && !is.null(x$engine)) {
    return(as.character(x$engine)[[1L]])
  }
  NA_character_
}

#' Whether IF pooling is a plug-in at fitted heterogeneity (not a posterior)
#' @noRd
interactmeta_ifrac_is_plugin <- function(x) {
  eng <- interactmeta_origin_engine(x)
  if (identical(eng, "metafor")) {
    return(TRUE)
  }
  comparison_is_mle(x)
}

#' @noRd
interactmeta_default_uncertainty <- function(engine, method) {
  if (identical(engine, "bayesmeta")) {
    return(paste(
      "equal-tail interval from the DIRECT-node Gaussian mixture",
      "(.weight and within-node .cond_*); native post.interval / convolve",
      "remain on $latent"
    ))
  }
  if (identical(engine, "metafor")) {
    return("Delta method from the selected metafor::rma coefficient covariance")
  }
  if (identical(method, "mle")) {
    return("Point estimate from a removed cmdstanr optimize() fit; interval unavailable")
  }
  paste(
    "Equal-tailed quantile of posterior draws from the selected",
    engine,
    "fit"
  )
}

#' Mixture / weighted moments for a scalar derived from one fit
#'
#' `cond_sd` is the within-node standard deviation (0 for MCMC draws).
#' @noRd
interactmeta_mixture_moments <- function(mean,
                                         weight = NULL,
                                         cond_sd = NULL,
                                         prob = 0.95,
                                         exponentiate = FALSE) {
  mean <- as.numeric(mean)
  n <- length(mean)
  empty <- list(
    median = NA_real_,
    mean = NA_real_,
    sd = NA_real_,
    lci = NA_real_,
    uci = NA_real_,
    estimate = NA_real_,
    prob = prob
  )
  if (!n) {
    return(empty)
  }
  if (is.null(weight) || length(weight) != n) {
    weight <- rep(1 / n, n)
  }
  w <- as.numeric(weight)
  w[!is.finite(w) | w < 0] <- 0
  if (!sum(w)) {
    w <- rep(1 / n, n)
  }
  w <- w / sum(w)
  if (is.null(cond_sd) || length(cond_sd) != n) {
    cond_sd <- rep(0, n)
  }
  s <- pmax(as.numeric(cond_sd), 0)
  s[!is.finite(s)] <- 0

  ok <- is.finite(mean)
  if (!any(ok)) {
    return(empty)
  }
  mean <- mean[ok]
  w <- w[ok]
  s <- s[ok]
  w <- w / sum(w)

  if (length(mean) == 1L && s[[1L]] <= 0) {
    v <- mean[[1L]]
    if (isTRUE(exponentiate)) {
      v <- exp(v)
    }
    return(list(
      median = v,
      mean = v,
      sd = NA_real_,
      lci = NA_real_,
      uci = NA_real_,
      estimate = v,
      prob = prob
    ))
  }

  if (all(s <= 0) && length(unique(w)) == 1L) {
    if (isTRUE(exponentiate)) {
      mean <- exp(mean)
    }
    return(comparison_draw_moments(mean, prob = prob))
  }

  gm <- interactmeta_gm(mean = mean, sd = s, weight = w)
  sm <- interactmeta_gm_summary(gm, prob = prob)
  if (isTRUE(exponentiate)) {
    ln <- interactmeta_gm_lognormal_moments(gm)
    return(list(
      median = exp(sm$median),
      mean = ln$mean,
      sd = ln$sd,
      lci = exp(sm$lower),
      uci = exp(sm$upper),
      estimate = exp(sm$median),
      prob = prob
    ))
  }
  list(
    median = sm$median,
    mean = sm$mean,
    sd = sm$sd,
    lci = sm$lower,
    uci = sm$upper,
    estimate = sm$median,
    prob = prob
  )
}

#' Conditional SD of a linear combination of CAMS primitives on one draws table
#'
#' Uses the stored overall-block joint moments `Var(alpha)`, `Var(theta)`,
#' and `Cov(alpha, theta)`, plus independent `Var(gamma)`. This is how
#' population `M* = alpha + pi* theta` keeps the covariance; do not replace
#' it with an independent convolution of the `alpha` and `theta` mixture
#' marginals.
#' Coefficients \eqn{(a_\alpha,a_\theta,a_\gamma)} apply to
#' \eqn{(\alpha,\theta,\gamma)}, **not** \eqn{(\alpha, b_{\mathrm{ifrac}}, \gamma)}.
#' Stored `.cond_var_theta` is \eqn{\mathrm{Var}(\theta)} from the overall
#' `bmr` block. Use \eqn{(1, \pi, 0)} for \eqn{M=\alpha+\pi\theta},
#' \eqn{(1, \pi, -\pi)} for \eqn{\mu_A}, and \eqn{(1, \pi, 1-\pi)} for \eqn{\mu_B}.
#' @noRd
interactmeta_lincomb_cond_sd <- function(post,
                                         a_alpha,
                                         a_theta,
                                         a_gamma) {
  n <- nrow(post)
  rec <- function(a) {
    a <- as.numeric(a)
    if (length(a) == 1L) {
      a <- rep(a, n)
    }
    if (length(a) != n) {
      stop("Linear-combination coefficients must recycle to the draws table.", call. = FALSE)
    }
    a
  }
  a_alpha <- rec(a_alpha)
  a_theta <- rec(a_theta)
  a_gamma <- rec(a_gamma)
  va <- if (".cond_var_alpha" %in% names(post)) {
    pmax(as.numeric(post$.cond_var_alpha), 0)
  } else {
    rep(0, n)
  }
  vt <- if (".cond_var_theta" %in% names(post)) {
    pmax(as.numeric(post$.cond_var_theta), 0)
  } else {
    rep(0, n)
  }
  cat <- if (".cond_cov_alpha_theta" %in% names(post)) {
    as.numeric(post$.cond_cov_alpha_theta)
  } else {
    rep(0, n)
  }
  vg <- if (".cond_sd_gamma" %in% names(post)) {
    pmax(as.numeric(post$.cond_sd_gamma), 0)^2
  } else {
    rep(0, n)
  }
  va[!is.finite(va)] <- 0
  vt[!is.finite(vt)] <- 0
  cat[!is.finite(cat)] <- 0
  vg[!is.finite(vg)] <- 0
  sqrt(pmax(
    a_alpha^2 * va + a_theta^2 * vt + 2 * a_alpha * a_theta * cat + a_gamma^2 * vg,
    0
  ))
}

#' Extra predictive SD from new-trial random effects (same fit's tau)
#' @noRd
interactmeta_prediction_re_sd <- function(post, omega, x = NULL) {
  n <- nrow(post)
  omega <- as.numeric(omega)
  if (length(omega) == 1L) {
    omega <- rep(omega, n)
  }
  tau_m <- if ("sd_trial__Intercept" %in% names(post)) {
    pmax(as.numeric(post$sd_trial__Intercept), 0)
  } else {
    rep(0, n)
  }
  slope <- grep("^sd_trial__ifrac2_", names(post), value = TRUE)
  tau_g <- if (length(slope)) {
    pmax(as.numeric(post[[slope[[1L]]]]), 0)
  } else {
    rep(0, n)
  }
  if (is.null(x)) {
    sqrt(pmax(tau_m^2 + (omega * tau_g)^2, 0))
  } else {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, n)
    }
    sqrt(pmax(tau_m^2 + ((x - omega) * tau_g)^2, 0))
  }
}

#' Whether this posterior table carries quadrature conditional moments
#' @noRd
interactmeta_has_cond_moments <- function(post) {
  is.data.frame(post) &&
    ".cond_var_alpha" %in% names(post) &&
    ".cond_sd_gamma" %in% names(post)
}

#' Fill metafor Wald intervals from the stored rma objects (no refit)
#' @noRd
interactmeta_metafor_lincomb <- function(fit,
                                         L_ov,
                                         L_int,
                                         pred_var = 0,
                                         exponentiate = FALSE,
                                         prob = 0.95) {
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  if (is.null(ov) || is.null(int)) {
    return(NULL)
  }
  b_ov <- as.numeric(ov$beta)
  b_int <- as.numeric(int$beta)[[1L]]
  L_ov <- as.numeric(L_ov)
  L_int <- as.numeric(L_int)[[1L]]
  if (length(b_ov) != length(L_ov)) {
    return(NULL)
  }

  if (isFALSE(fit$ecological_slope %||% fit$metafor$ecological_slope)) {
    # theta = gamma: fold the interaction loading onto the shared gamma.
    L_alpha <- L_ov[[1L]]
    L_gamma <- L_ov[[2L]] + L_int
    est <- L_alpha * b_ov[[1L]] + L_gamma * b_ov[[2L]]
    vb <- as.matrix(ov$vb)
    var <- L_alpha^2 * vb[1L, 1L] +
      L_gamma^2 * vb[2L, 2L] +
      2 * L_alpha * L_gamma * vb[1L, 2L] +
      as.numeric(pred_var)[[1L]]
  } else {
    est <- sum(L_ov * b_ov) + L_int * b_int
    vb_ov <- as.matrix(ov$vb)
    vb_int <- as.matrix(int$vb)
    var <- as.numeric(t(L_ov) %*% vb_ov %*% L_ov) +
      (L_int^2) * as.numeric(vb_int[[1L]]) +
      as.numeric(pred_var)[[1L]]
  }
  se <- sqrt(max(var, 0))
  z <- stats::qnorm(1 - (1 - prob) / 2)
  lo <- est - z * se
  hi <- est + z * se
  if (isTRUE(exponentiate)) {
    est_y <- exp(est)
    return(list(
      estimate = est_y,
      mean = est_y,
      sd = est_y * se,
      lci = exp(lo),
      uci = exp(hi),
      prob = prob
    ))
  }
  list(
    estimate = est,
    mean = est,
    sd = se,
    lci = lo,
    uci = hi,
    prob = prob
  )
}

#' @noRd
interactmeta_metafor_omega <- function(fit, overall_ifrac = NULL) {
  if (inherits(overall_ifrac, "interactmeta_overall_ifrac")) {
    return(stats::median(as.numeric(overall_ifrac), na.rm = TRUE))
  }
  if (is.numeric(overall_ifrac) && length(overall_ifrac)) {
    return(stats::median(as.numeric(overall_ifrac), na.rm = TRUE))
  }
  oi <- tryCatch(ifrac(fit, summary = FALSE), error = function(e) NULL)
  if (is.null(oi)) {
    return(NA_real_)
  }
  stats::median(as.numeric(oi), na.rm = TRUE)
}

#' Replace NA intervals on a metafor effect summary using the stored rma vcov
#' @noRd
interactmeta_metafor_fill_effects <- function(fit,
                                              summary,
                                              quantity,
                                              omega,
                                              exponentiate = FALSE,
                                              prob = 0.95,
                                              delta_ifrac2 = 1) {
  if (!identical(fit$engine, "metafor") || is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  if (all(is.finite(summary$lci)) && all(is.finite(summary$uci))) {
    return(summary)
  }
  omega <- as.numeric(omega)[[1L]]
  if (!is.finite(omega)) {
    return(summary)
  }
  tau_m2 <- max(as.numeric(fit$metafor$overall$tau2), 0)
  tau_g2 <- max(as.numeric(fit$metafor$interaction$tau2), 0)
  for (i in seq_len(nrow(summary))) {
    ty <- if ("type" %in% names(summary)) {
      as.character(summary$type[[i]])
    } else {
      "Pooled"
    }
    pred <- identical(ty, "Prediction")
    lin <- NULL
    if (identical(quantity, "interaction")) {
      L_ov <- c(0, 0)
      L_int <- 1
      pvar <- if (pred) tau_g2 * (as.numeric(delta_ifrac2)[[1L]]^2) else 0
      lin <- interactmeta_metafor_lincomb(
        fit, L_ov, L_int, pred_var = pvar,
        exponentiate = exponentiate, prob = prob
      )
    } else if (identical(quantity, "combined")) {
      L_ov <- c(1, omega)
      L_int <- 0
      pvar <- if (pred) tau_m2 + (omega^2) * tau_g2 else 0
      lin <- interactmeta_metafor_lincomb(
        fit, L_ov, L_int, pred_var = pvar,
        exponentiate = exponentiate, prob = prob
      )
    } else if (identical(quantity, "subgroups") && "subgroup" %in% names(summary)) {
      st <- attr(fit, "interactmeta", exact = TRUE)
      ref <- as.character(st$reference %||% "")[[1L]]
      sg <- as.character(summary$subgroup[[i]])
      if (identical(sg, ref)) {
        L_ov <- c(1, omega)
        L_int <- -omega
        pvar <- if (pred) tau_m2 + (omega^2) * tau_g2 else 0
      } else {
        L_ov <- c(1, omega)
        L_int <- 1 - omega
        pvar <- if (pred) tau_m2 + ((1 - omega)^2) * tau_g2 else 0
      }
      lin <- interactmeta_metafor_lincomb(
        fit, L_ov, L_int, pred_var = pvar,
        exponentiate = exponentiate, prob = prob
      )
    }
    if (is.null(lin)) {
      next
    }
    summary$estimate[[i]] <- lin$estimate
    summary$mean[[i]] <- lin$mean
    summary$sd[[i]] <- lin$sd
    summary$lci[[i]] <- lin$lci
    summary$uci[[i]] <- lin$uci
    summary$prob[[i]] <- lin$prob
  }
  summary
}

#' Attach quadrature weights / within-tau SDs onto effect-draw rows
#' @noRd
interactmeta_decorate_effect_draws <- function(df, post, cond_sd = NULL) {
  if (is.data.frame(post) && ".weight" %in% names(post)) {
    if (nrow(df) == nrow(post)) {
      df$.weight <- post$.weight
    } else if (nrow(df) > 0L && nrow(post) > 0L && (nrow(df) %% nrow(post) == 0L)) {
      df$.weight <- rep(post$.weight, length.out = nrow(df))
    }
  }
  if (!is.null(cond_sd)) {
    df$.cond_sd <- as.numeric(cond_sd)
  }
  df
}

#' Summarize derived estimands from the same fit (no refit)
#' @noRd
interactmeta_summarize_effects <- function(draws,
                                           group_cols,
                                           prob,
                                           exponentiate,
                                           fit,
                                           quantity,
                                           omega,
                                           delta_ifrac2 = 1) {
  mix <- is.data.frame(draws) && ".cond_sd" %in% names(draws)
  if (mix) {
    # Public Bayesian intervals: equal-tail mixture on DIRECT nodes (fast).
    # Native post.interval / convolve stay on $latent only.
    summ <- comparison_effects_summary(
      draws,
      group_cols = group_cols,
      prob = prob,
      exponentiate = exponentiate
    )
    draws <- comparison_transform_effect_draws(draws, exponentiate)
  } else {
    draws <- comparison_transform_effect_draws(draws, exponentiate)
    summ <- comparison_effects_summary(draws, group_cols = group_cols, prob = prob)
    summ <- interactmeta_metafor_fill_effects(
      fit = fit,
      summary = summ,
      quantity = quantity,
      omega = omega,
      exponentiate = exponentiate,
      prob = prob,
      delta_ifrac2 = delta_ifrac2
    )
  }
  list(draws = draws, summary = summ)
}
