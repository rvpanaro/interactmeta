#' Bayesian DA and AD via independent bayesmeta random-effects models
#'
#' DA: two independent univariate Bayesian MAs, then μ_B − μ_A by
#' bayesmeta::convolve(). The Gaussian mixture is a check. AD: one Bayesian
#' MA of within-trial g_j, with intervals from $post.interval. Neither uses CAMS IF.
#' @noRd
NULL

#' @noRd
interactmeta_construction_tau_prior <- function(analysis, tau.prior = NULL) {
  if (!is.null(tau.prior)) {
    return(normalize_cam_prior(tau.prior, role = "tau"))
  }
  if (identical(analysis, "ad")) {
    cam_default_tau_interaction_prior()
  } else {
    cam_default_tau_prior()
  }
}

#' @noRd
interactmeta_construction_mu_prior <- function(mu.prior = NULL) {
  if (is.null(mu.prior)) {
    return(cam_default_fixed_prior())
  }
  normalize_cam_prior(mu.prior, role = "fixed")
}

#' @noRd
interactmeta_bayesmeta_univariate <- function(y,
                                              vi,
                                              labels,
                                              tau.prior,
                                              mu.prior,
                                              engine_args = list()) {
  check_bayesmeta()
  tau_bm <- prior_to_bayesmeta(tau.prior)
  fe <- prior_to_bayesmeta_fixed(mu.prior)
  args <- if (is.null(engine_args)) list() else engine_args
  if ("convolve_derived" %in% names(args)) {
    stop(
      "`convolve_derived` cannot change reported intervals. Independent sums use bayesmeta::convolve().",
      call. = FALSE
    )
  }
  interval.type <- interactmeta_match_interval_type(args$interval.type)
  args$interval.type <- NULL
  args[c("method", "engine", "chains", "iter", "warmup")] <- NULL
  if (!anyNA(fe$mu.prior)) {
    args$mu.prior <- unname(fe$mu.prior)
  }
  do.call(
    bayesmeta::bayesmeta,
    c(
      list(
        y = as.numeric(y),
        sigma = interactmeta_sigma_from_vi(vi),
        labels = as.character(labels),
        tau.prior = tau_bm,
        interval.type = interval.type
      ),
      args
    )
  )
}

#' @noRd
interactmeta_bayesmeta_tau_summary <- function(fit, prob = 0.95) {
  tau <- as.numeric(fit$support[, "tau"])
  w <- as.numeric(fit$support[, "weight"])
  interactmeta_discrete_summary(interactmeta_discrete(tau, w), prob = prob)
}

#' @noRd
interactmeta_bayesmeta_tau2_mean <- function(fit) {
  tau <- as.numeric(fit$support[, "tau"])
  w <- as.numeric(fit$support[, "weight"])
  w <- w / sum(w)
  sum(w * tau^2)
}

#' @noRd
interactmeta_bayesmeta_re_weights <- function(fit) {
  tau2 <- interactmeta_bayesmeta_tau2_mean(fit)
  sig2 <- as.numeric(fit$sigma)^2
  w <- 1 / (sig2 + tau2)
  w / sum(w)
}

#' @noRd
interactmeta_bayesmeta_mu_row <- function(fit, label, method = "bayes", prob = 0.95) {
  it <- interactmeta_match_interval_type(fit$interval.type)
  loc <- interactmeta_native_location(fit, prob = prob, interval.type = it, margin = "effect")
  tau_loc <- interactmeta_native_location(fit, prob = prob, interval.type = it, margin = "tau")
  data.frame(
    label = label,
    estimate = loc$median,
    se = as.numeric(fit$summary["sd", "mu"])[[1L]],
    ci.lb = loc$lower,
    ci.ub = loc$upper,
    tau2 = interactmeta_bayesmeta_tau2_mean(fit),
    tau = tau_loc$median,
    method = method,
    stringsAsFactors = FALSE
  )
}

#' Bayesian DA: independent RE-MA per subgroup, then μ_B − μ_A
#' @noRd
interactmeta_da_from_pairs_bayesmeta <- function(pairs,
                                                 tau.prior,
                                                 mu.prior,
                                                 engine_args,
                                                 call) {
  ref <- attr(pairs, "reference")
  nonref <- attr(pairs, "nonref")
  fit_ref <- interactmeta_bayesmeta_univariate(
    y = pairs$y_ref,
    vi = pairs$v_ref,
    labels = pairs$trial,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    engine_args = engine_args
  )
  fit_non <- interactmeta_bayesmeta_univariate(
    y = pairs$y_nonref,
    vi = pairs$v_nonref,
    labels = pairs$trial,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    engine_args = engine_args
  )
  gm_ref <- interactmeta_gm_from_bayesmeta(fit_ref)
  gm_non <- interactmeta_gm_from_bayesmeta(fit_non)
  # Exact mixture is a check. The published difference is convolve().
  gm_da <- interactmeta_gm_linear_sum(gm_non, 1, gm_ref, -1)
  grid_delta <- fit_non$delta %||% fit_ref$delta
  grid_eps <- fit_non$epsilon %||% fit_ref$epsilon
  conv <- interactmeta_convolve_sum(
    dens1 = interactmeta_bayesmeta_mu_density(fit_non),
    dens2 = interactmeta_negate_density(interactmeta_bayesmeta_mu_density(fit_ref)),
    cdf1 = interactmeta_bayesmeta_mu_cdf(fit_non),
    cdf2 = interactmeta_negate_cdf(interactmeta_bayesmeta_mu_cdf(fit_ref)),
    delta = grid_delta,
    epsilon = grid_eps
  )
  it <- interactmeta_match_interval_type(fit_non$interval.type)
  tol <- fit_non$tol.uniroot %||% fit_ref$tol.uniroot
  bounds <- interactmeta_interval_from_distribution(
    conv$quantile, conv$cdf, prob = 0.95, interval.type = it, tol = tol
  )
  sm_da <- list(
    median = as.numeric(conv$quantile(0.5))[[1L]],
    mean = interactmeta_gm_mean(gm_da),
    sd = sqrt(max(interactmeta_gm_var(gm_da), 0)),
    lower = bounds[[1L]],
    upper = bounds[[2L]]
  )
  sm_ref <- interactmeta_bayesmeta_mu_row(fit_ref, ref)
  sm_non <- interactmeta_bayesmeta_mu_row(fit_non, nonref)
  list(
    reference = ref,
    nonref = nonref,
    pairs = pairs,
    subgroups = rbind(sm_ref, sm_non),
    interaction = data.frame(
      label = paste0(nonref, " - ", ref, " (DA)"),
      estimate = sm_da$median,
      se = sm_da$sd,
      ci.lb = sm_da$lower,
      ci.ub = sm_da$upper,
      tau2 = NA_real_,
      tau = NA_real_,
      method = "bayes",
      stringsAsFactors = FALSE
    ),
    weights = data.frame(
      trial = pairs$trial,
      weight_ref = interactmeta_bayesmeta_re_weights(fit_ref),
      weight_nonref = interactmeta_bayesmeta_re_weights(fit_non),
      stringsAsFactors = FALSE
    ),
    fits = list(reference = fit_ref, nonref = fit_non),
    posterior = list(
      mu_reference = interactmeta_posterior_1d(
        gm = gm_ref,
        name = "mu_reference",
        published = interactmeta_native_published(fit_ref, "effect"),
        interval.type = it,
        tol.uniroot = tol
      ),
      mu_nonref = interactmeta_posterior_1d(
        gm = gm_non,
        name = "mu_nonref",
        published = interactmeta_native_published(fit_non, "effect"),
        interval.type = it,
        tol.uniroot = tol
      ),
      interaction = interactmeta_posterior_1d(
        gm = gm_da,
        convolve = conv,
        name = "gamma_DA",
        note = "Independent Bayesian DA: published posterior is convolve(mu_B, -mu_A); gm is the mixture check.",
        interval.type = it,
        tol.uniroot = tol
      )
    ),
    priors = list(tau = tau.prior, fixed = mu.prior),
    method = "bayes",
    call = call
  )
}

#' Bayesian AD: RE-MA of g_j = y_B − y_A
#' @noRd
interactmeta_ad_from_pairs_bayesmeta <- function(pairs,
                                                 tau.prior,
                                                 mu.prior,
                                                 engine_args,
                                                 call) {
  ref <- attr(pairs, "reference")
  nonref <- attr(pairs, "nonref")
  fit_g <- interactmeta_bayesmeta_univariate(
    y = pairs$g,
    vi = pairs$vg,
    labels = pairs$trial,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    engine_args = engine_args
  )
  sm <- interactmeta_bayesmeta_mu_row(fit_g, paste0(nonref, " - ", ref))
  gm <- interactmeta_gm_from_bayesmeta(fit_g)
  list(
    reference = ref,
    nonref = nonref,
    pairs = pairs,
    interaction = sm,
    weights = data.frame(
      trial = pairs$trial,
      weight = interactmeta_bayesmeta_re_weights(fit_g),
      stringsAsFactors = FALSE
    ),
    fits = list(interaction = fit_g),
    posterior = list(
      interaction = interactmeta_posterior_1d(
        gm = gm,
        name = "gamma_AD",
        published = interactmeta_native_published(fit_g, "effect"),
        interval.type = fit_g$interval.type
      )
    ),
    priors = list(tau = tau.prior, fixed = mu.prior),
    method = "bayes",
    call = call
  )
}
