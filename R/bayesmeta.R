# ---- interactmeta_bayesmeta.R ----
#' Fit CAMS via orthogonal bayesmeta / bmr decomposition
#'
#' After the information-fraction rotation, binary CAMS separates into two
#' independent models:
#' 1. overall meta-regression of trial-level IV-pooled effects on IF (`bmr`)
#' 2. interaction RE meta-analysis of within-trial contrasts (`bayesmeta`)
#'
#' Fixed effects are recovered as
#' \eqn{\alpha} = overall intercept, \eqn{\gamma} = interaction mean,
#' \eqn{\beta = \theta - \gamma} where \eqn{\theta} is the IF slope of the
#' overall meta-regression (\eqn{\alpha + \pi(\beta+\gamma)}).
#'
#' Published **engine-native** intervals on `$latent` use one rule
#' (`interval.type`, default `"central"`): margins from `$post.interval`,
#' independent-block sums from `bayesmeta::convolve()`. Public extractors and
#' [predict_cams_by_if()] default to central intervals on the DIRECT-node
#' Gaussian mixture; [plot_contribution()] ribbon limits default to central
#' (equal-tailed) mixture intervals. `$support` is the Cartesian product
#' of DIRECT tau-nodes with `.weight = w_i * v_k` (conditional means;
#' within-node variance is in `.cond_*` and `$latent`). `$mu` remains the
#' **fixed-effect** linear predictor. No `rposterior()` calls.
#'
#' @noRd
.interactmeta_fit_bayesmeta <- function(data_fit,
                                        design,
                                        trial_var = "study",
                                        subgroup_var = "level",
                                        priors = NULL,
                                        fit_args = list(),
                                        seed = 42L,
                                        silent = 2,
                                        lazy_latent = TRUE,
                                        ...) {
  force(seed)
  check_bayesmeta()
  if (design$n_nonref != 1L) {
    stop(
      "`engine = \"bayesmeta\"` currently supports binary subgroups only ",
      "(one non-reference `ifrac_*` column).",
      call. = FALSE
    )
  }
  if (is.null(priors)) {
    priors <- cam_priors_resolve()
  }

  blocks <- interactmeta_bayesmeta_trial_blocks(
    data_fit = data_fit,
    design = design,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )

  fit_args <- interactmeta_merge_engine_args(fit_args, list(...))
  if ("convolve_derived" %in% names(fit_args)) {
    stop(
      "`convolve_derived` cannot change reported intervals. Derived CAMS posteriors always come from bayesmeta::convolve() of independent blocks.",
      call. = FALSE
    )
  }
  interval.type <- interactmeta_match_interval_type(fit_args$interval.type)
  fit_args$interval.type <- NULL
  fit_args[c(
    "chains", "iter", "warmup", "thin", "refresh",
    "iter_warmup", "iter_sampling", "parallel_chains",
    "method", "cores", "control",
    "tau.prior", "tau.prior.overall", "tau.prior.interaction",
    "tau.interaction.prior", "fixed.prior", "beta.prior.mean",
    "beta.prior.sd", "beta.prior.cov", "mu.prior",
    "lazy_latent", "precompute"
  )] <- NULL
  n_support <- fit_args$n_support
  if (!is.null(n_support)) {
    n_support <- as.integer(n_support)[[1L]]
    fit_args$n_support <- NULL
  }

  tau_overall <- prior_to_bayesmeta(priors$tau)
  tau_int <- prior_to_bayesmeta(priors$tau_interaction)
  fe <- prior_to_bayesmeta_fixed(priors$fixed)

  bmr_args <- fit_args
  bm_args <- fit_args
  if (!is.null(fe$beta.prior.mean)) {
    bmr_args$beta.prior.mean <- unname(fe$beta.prior.mean)
    bmr_args$beta.prior.sd <- unname(fe$beta.prior.sd)
    bmr_args$beta.prior.cov <- fe$beta.prior.cov
  }
  if (!anyNA(fe$mu.prior)) {
    bm_args$mu.prior <- fe$mu.prior
  }

  bmr_call <- list(
    y = blocks$overall$yi,
    sigma = interactmeta_sigma_from_vi(blocks$overall$vi),
    labels = as.character(blocks$overall$trial),
    X = cbind(
      intercept = 1,
      ifrac = blocks$overall$ifrac
    ),
    tau.prior = tau_overall,
    interval.type = interval.type
  )
  bm_call <- list(
    y = blocks$interaction$yi,
    sigma = interactmeta_sigma_from_vi(blocks$interaction$vi),
    labels = as.character(blocks$interaction$trial),
    tau.prior = tau_int,
    interval.type = interval.type
  )
  if (!is.null(n_support)) {
    bmr_call$n.support <- n_support
    bm_call$n.support <- n_support
  }
  overall_fit <- do.call(bayesmeta::bmr, c(bmr_call, bmr_args))
  interaction_fit <- do.call(bayesmeta::bayesmeta, c(bm_call, bm_args))

  support <- interactmeta_bayesmeta_draws_from_support(
    overall_fit = overall_fit,
    interaction_fit = interaction_fit,
    design = design
  )
  mu <- interactmeta_bayesmeta_mu_matrix(support, data_fit, design)
  latent <- interactmeta_bayesmeta_latent(
    overall_fit = overall_fit,
    interaction_fit = interaction_fit,
    blocks = blocks,
    interval.type = interval.type,
    trial_marginals = if (isTRUE(lazy_latent)) "none" else "full"
  )

  structure(
    list(
      support = support,
      mu = mu,
      latent = latent,
      data = data_fit,
      formula = stats::as.formula(design$formula),
      bayesmeta = list(
        overall = overall_fit,
        interaction = interaction_fit,
        trial_blocks = blocks,
        ecological_slope = TRUE
      ),
      interval.type = interval.type,
      prior = design$prior,
      priors = priors,
      algorithm = "quadrature",
      ecological_slope = TRUE
    ),
    class = c("interactmeta", "list")
  )
}

#' Trial-level overall / interaction blocks for bayesmeta CAMS
#' @noRd
interactmeta_bayesmeta_trial_blocks <- function(data_fit,
                                                design,
                                                trial_var,
                                                subgroup_var) {
  ref_fit <- design$coef_map$subgroup_fit[[1L]]
  nonref_fit <- design$coef_map$subgroup_fit[[2L]]
  ifrac_std <- design$coef_map$ifrac_std[[1L]]

  blocks <- comparison_trial_effect_blocks(
    data = data_fit,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = ref_fit,
    nonref = nonref_fit,
    ifrac_var = ifrac_std,
    min_vi = 0
  )
  blocks$reference <- design$reference
  blocks$nonref <- design$non_ref_labels[[1L]]
  blocks$ifrac_var <- design$ifrac_vars[[1L]]
  blocks
}

#' Deterministic CAMS draws from independent DIRECT tau-node product
#'
#' Cartesian product of overall and interaction support nodes with
#' `.weight = w_i v_k`. Rows store **conditional means** given each tau node
#' (within-tau coefficient variance lives in `$latent`). No CDF pairing and
#' no `rposterior()`.
#' @noRd
interactmeta_bayesmeta_draws_from_support <- function(overall_fit,
                                                      interaction_fit,
                                                      design) {
  ov <- interactmeta_bmr_support_table(overall_fit)
  int <- interactmeta_bayesmeta_support_table(interaction_fit)

  need_ov <- c("weight", "tau", "intercept", "ifrac")
  missing_ov <- setdiff(need_ov, names(ov))
  if (length(missing_ov)) {
    stop(
      "Unexpected `bmr` support columns; missing: ",
      paste(missing_ov, collapse = ", "),
      call. = FALSE
    )
  }
  need_int <- c("weight", "tau", "mu")
  missing_int <- setdiff(need_int, names(int))
  if (length(missing_int)) {
    stop(
      "Unexpected `bayesmeta` support columns; missing: ",
      paste(missing_int, collapse = ", "),
      call. = FALSE
    )
  }

  n_ov <- nrow(ov)
  n_int <- nrow(int)
  i_ov <- rep(seq_len(n_ov), each = n_int)
  i_int <- rep(seq_len(n_int), times = n_ov)
  alpha <- ov$intercept[i_ov]
  theta <- ov$ifrac[i_ov]
  gamma <- int$mu[i_int]
  beta <- theta - gamma
  tau_m <- pmax(as.numeric(ov$tau[i_ov]), 0)
  tau_g <- pmax(as.numeric(int$tau[i_int]), 0)
  w <- as.numeric(ov$weight[i_ov]) * as.numeric(int$weight[i_int])
  w <- w / sum(w)

  cond <- interactmeta_bayesmeta_cond_moments(
    overall_fit = overall_fit,
    interaction_fit = interaction_fit,
    ov = ov,
    int = int
  )
  var_alpha <- cond$var_alpha[i_ov]
  var_theta <- cond$var_theta[i_ov]
  cov_at <- cond$cov_alpha_theta[i_ov]
  sd_gamma <- cond$sd_gamma[i_int]

  sg_nonref <- design$coef_map$subgroup_fit[[2L]]
  ifrac2_std <- design$coef_map$ifrac2_std[[1L]]
  n <- length(alpha)

  out <- data.frame(
    .chain = 1L,
    .iteration = seq_len(n),
    .draw = seq_len(n),
    .weight = w,
    .cond_var_alpha = var_alpha,
    .cond_var_theta = var_theta,
    .cond_cov_alpha_theta = cov_at,
    .cond_sd_gamma = sd_gamma,
    b_Intercept = as.numeric(alpha),
    b_ifrac_1 = as.numeric(beta),
    sd_trial__Intercept = tau_m,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[[paste0("b_subgroup", sg_nonref)]] <- as.numeric(gamma)
  out[[paste0("sd_trial__", ifrac2_std)]] <- tau_g
  out
}

#' Within-tau moments from the native bmr / bayesmeta objects (no refit)
#' @noRd
interactmeta_bayesmeta_cond_moments <- function(overall_fit, interaction_fit, ov, int) {
  n_ov <- nrow(ov)
  n_int <- nrow(int)
  var_alpha <- rep(0, n_ov)
  var_theta <- rep(0, n_ov)
  cov_at <- rep(0, n_ov)
  cov_arr <- overall_fit$support$covariance
  if (!is.null(cov_arr) && length(dim(cov_arr)) == 3L) {
    nms <- colnames(overall_fit$support$mean)
    i_a <- match("intercept", nms)
    if (is.na(i_a)) {
      i_a <- 1L
    }
    i_t <- match("ifrac", nms)
    if (is.na(i_t)) {
      i_t <- min(2L, length(nms))
    }
    var_alpha <- pmax(as.numeric(cov_arr[, i_a, i_a]), 0)
    var_theta <- pmax(as.numeric(cov_arr[, i_t, i_t]), 0)
    cov_at <- as.numeric(cov_arr[, i_a, i_t])
    if (length(var_alpha) != n_ov) {
      var_alpha <- rep(var_alpha[[1L]] %||% 0, n_ov)
      var_theta <- rep(var_theta[[1L]] %||% 0, n_ov)
      cov_at <- rep(cov_at[[1L]] %||% 0, n_ov)
    }
  }

  sd_gamma <- if ("sd" %in% names(int)) {
    pmax(as.numeric(int$sd), 0)
  } else {
    rep(0, n_int)
  }
  if (!length(sd_gamma) || all(!is.finite(sd_gamma))) {
    gm <- tryCatch(
      interactmeta_gm_from_bayesmeta(interaction_fit),
      error = function(e) NULL
    )
    if (!is.null(gm) && length(gm$sd) == n_int) {
      sd_gamma <- pmax(as.numeric(gm$sd), 0)
    } else {
      sd_gamma <- rep(0, n_int)
    }
  }
  sd_gamma[!is.finite(sd_gamma)] <- 0
  list(
    var_alpha = var_alpha,
    var_theta = var_theta,
    cov_alpha_theta = cov_at,
    sd_gamma = sd_gamma
  )
}

#' Grid settings stored on the fitted objects (not a second pair of literals).
#' @noRd
interactmeta_bayesmeta_grid <- function(overall_fit, interaction_fit) {
  delta <- overall_fit$delta %||% interaction_fit$delta
  epsilon <- overall_fit$epsilon %||% interaction_fit$epsilon
  if (is.null(delta) || is.null(epsilon)) {
    stop("Fitted bayesmeta objects are missing `$delta` or `$epsilon`.", call. = FALSE)
  }
  if (!is.null(interaction_fit$delta) && !isTRUE(all.equal(delta, interaction_fit$delta))) {
    stop("Combined and interaction fits stored different `$delta` values.", call. = FALSE)
  }
  if (!is.null(interaction_fit$epsilon) && !isTRUE(all.equal(epsilon, interaction_fit$epsilon))) {
    stop("Combined and interaction fits stored different `$epsilon` values.", call. = FALSE)
  }
  tol <- overall_fit$tol.uniroot %||% interaction_fit$tol.uniroot %||%
    (2^16 * .Machine$double.eps)
  list(delta = delta, epsilon = epsilon, tol.uniroot = tol)
}

#' @noRd
interactmeta_native_published <- function(object,
                                         margin = c("effect", "tau", "shrink"),
                                         which = NULL,
                                         predict = FALSE,
                                         individual = FALSE) {
  list(
    kind = "native",
    object = object,
    margin = match.arg(margin),
    which = which,
    predict = isTRUE(predict),
    individual = individual
  )
}

#' Marginal 1-D posteriors for derived CAMS quantities
#'
#' Cross-block independence is used. Within-block trial effects that share
#' hyperparameters are **not** given a joint posterior here.
#'
#' `$mu_A[[j]]` / `$mu_B[[j]]` are **trial-specific** latent subgroup
#' posteriors
#' \eqn{\mu_{Aj}=m_j-\pi_j g_j}, \eqn{\mu_{Bj}=m_j+(1-\pi_j)g_j}
#' at that trial's design information fraction \eqn{\pi_j}. They are
#' **not** population subgroup effects at a reporting IF \eqn{\pi^\star}.
#' Those are \eqn{\mu_A(\pi^\star)=\alpha+\pi^\star(\theta-\gamma)} and
#' \eqn{\mu_B(\pi^\star)=\alpha+\pi^\star(\theta-\gamma)+\gamma}, returned
#' by [modifiers()]. Do not form \eqn{M^\star=\alpha+\pi^\star\theta} by
#' independently combining the `alpha` and `theta` mixture marginals.
#' @noRd
interactmeta_bayesmeta_latent <- function(overall_fit,
                                          interaction_fit,
                                          blocks,
                                          interval.type = NULL,
                                          trial_marginals = c("full", "none")) {
  trial_marginals <- match.arg(trial_marginals)
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% overall_fit$interval.type
  )
  grid <- interactmeta_bayesmeta_grid(overall_fit, interaction_fit)
  trials <- as.character(blocks$overall$trial)
  pi <- as.numeric(blocks$overall$ifrac)
  names(pi) <- trials

  intercept_nm <- interactmeta_bmr_coef_name(overall_fit, "intercept")
  slope_nm <- interactmeta_bmr_coef_name(overall_fit, "slope")
  gm_alpha <- interactmeta_gm_from_bmr(overall_fit, intercept_nm)
  gm_theta <- interactmeta_gm_from_bmr(overall_fit, slope_nm)
  gm_gamma <- interactmeta_gm_tag_block(
    interactmeta_gm_from_bayesmeta(interaction_fit),
    "interaction"
  )
  gm_beta <- interactmeta_gm_linear_sum(gm_theta, 1, gm_gamma, -1)

  conv_beta <- interactmeta_convolve_sum(
    dens1 = interactmeta_bmr_coef_density(overall_fit, slope_nm),
    dens2 = interactmeta_negate_density(
      interactmeta_bayesmeta_mu_density(interaction_fit)
    ),
    cdf1 = interactmeta_bmr_coef_cdf(overall_fit, slope_nm),
    cdf2 = interactmeta_negate_cdf(
      interactmeta_bayesmeta_mu_cdf(interaction_fit)
    ),
    delta = grid$delta,
    epsilon = grid$epsilon
  )

  one <- function(gm, name, note, published = NULL, convolve = NULL) {
    interactmeta_posterior_1d(
      gm = gm,
      convolve = convolve,
      name = name,
      note = note,
      published = published,
      interval.type = interval.type,
      tol.uniroot = grid$tol.uniroot
    )
  }

  m <- vector("list", length(trials))
  g <- vector("list", length(trials))
  mu_A <- vector("list", length(trials))
  mu_B <- vector("list", length(trials))
  names(m) <- names(g) <- names(mu_A) <- names(mu_B) <- trials
  lazy_trials <- identical(trial_marginals, "none")

  if (!lazy_trials) {
  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    pij <- pi[[tr]]
    gm_m <- interactmeta_gm_from_bmr_shrink(overall_fit, tr)
    gm_g <- interactmeta_gm_from_bayesmeta_individual(interaction_fit, tr)
    m[[tr]] <- one(
      gm = gm_m,
      name = paste0("m[", tr, "]"),
      note = "marginal over overall-block hyperparameters; interval from bmr$shrink.interval",
      published = interactmeta_native_published(overall_fit, "shrink", which = tr)
    )
    g[[tr]] <- one(
      gm = gm_g,
      name = paste0("g[", tr, "]"),
      note = "marginal over interaction-block hyperparameters; interval from bayesmeta$post.interval(individual=)",
      published = interactmeta_native_published(
        interaction_fit, "effect", individual = tr
      )
    )

    gm_A <- interactmeta_gm_linear_sum(gm_m, 1, gm_g, -pij)
    gm_B <- interactmeta_gm_linear_sum(gm_m, 1, gm_g, 1 - pij)
    dens_m <- interactmeta_bmr_shrink_density(overall_fit, tr)
    cdf_m <- interactmeta_bmr_shrink_cdf(overall_fit, tr)
    dens_g <- interactmeta_bayesmeta_individual_density(interaction_fit, tr)
    cdf_g <- interactmeta_bayesmeta_individual_cdf(interaction_fit, tr)
    conv_A <- NULL
    pub_A <- NULL
    if (abs(pij) >= 1e-15) {
      conv_A <- interactmeta_convolve_sum(
        dens1 = dens_m,
        dens2 = interactmeta_scale_density(dens_g, -pij),
        cdf1 = cdf_m,
        cdf2 = interactmeta_scale_cdf(cdf_g, -pij),
        delta = grid$delta,
        epsilon = grid$epsilon
      )
    } else {
      pub_A <- interactmeta_native_published(overall_fit, "shrink", which = tr)
    }
    conv_B <- NULL
    pub_B <- NULL
    cB <- 1 - pij
    if (abs(cB) >= 1e-15) {
      conv_B <- interactmeta_convolve_sum(
        dens1 = dens_m,
        dens2 = interactmeta_scale_density(dens_g, cB),
        cdf1 = cdf_m,
        cdf2 = interactmeta_scale_cdf(cdf_g, cB),
        delta = grid$delta,
        epsilon = grid$epsilon
      )
    } else {
      pub_B <- interactmeta_native_published(overall_fit, "shrink", which = tr)
    }
    mu_A[[tr]] <- one(
      gm = gm_A,
      convolve = conv_A,
      published = pub_A,
      name = paste0("mu_A[", tr, "]"),
      note = "trial-specific latent; mu_Aj = m_j - pi_j * g_j via convolve() (not modifiers() at pi*)"
    )
    mu_B[[tr]] <- one(
      gm = gm_B,
      convolve = conv_B,
      published = pub_B,
      name = paste0("mu_B[", tr, "]"),
      note = "trial-specific latent; mu_Bj = m_j + (1 - pi_j) * g_j via convolve() (not modifiers() at pi*)"
    )
  }
  }

  list(
    independence = "cross_block",
    interval.type = interval.type,
    delta = grid$delta,
    epsilon = grid$epsilon,
    tol.uniroot = grid$tol.uniroot,
    note = paste(
      "Marginal posteriors only. Combined and interaction blocks are independent.",
      "Trials within a block share hyperparameters and are not mutually independent.",
      "Published intervals on one stored fit use that object's $post.interval",
      "(or bmr$shrink.interval). Published intervals for independent-block sums",
      "use bayesmeta::convolve(). The Gaussian mixture is a check, not a second interval.",
      "alpha and theta are not convolved: they share the overall-block support.",
      "latent$mu_A[[j]] / $mu_B[[j]] are trial-specific (m_j, g_j, pi_j),",
      "not population mu_A(pi*) / mu_B(pi*); those come from modifiers()."
    ),
    alpha = one(
      gm = gm_alpha,
      name = "alpha",
      note = "interval from bmr$post.interval(which.beta = intercept)",
      published = interactmeta_native_published(overall_fit, "effect", which = intercept_nm)
    ),
    theta = one(
      gm = gm_theta,
      name = "theta",
      note = "interval from bmr$post.interval(which.beta = slope). Not convolved with alpha.",
      published = interactmeta_native_published(overall_fit, "effect", which = slope_nm)
    ),
    gamma = one(
      gm = gm_gamma,
      name = "gamma",
      note = "interval from bayesmeta$post.interval",
      published = interactmeta_native_published(interaction_fit, "effect")
    ),
    beta = one(
      gm = gm_beta,
      convolve = conv_beta,
      name = "beta",
      note = "published posterior is convolve(theta, -gamma); gm is the exact mixture check"
    ),
    tau_m = interactmeta_discrete(overall_fit$support$tau, overall_fit$support$weight),
    tau_g = interactmeta_discrete(
      as.numeric(interaction_fit$support[, "tau"]),
      as.numeric(interaction_fit$support[, "weight"])
    ),
    m = m,
    g = g,
    mu_A = mu_A,
    mu_B = mu_B,
    pi = pi,
    trials = trials,
    cache = new.env(parent = emptyenv()),
    .lazy_trials = lazy_trials
  )
}

#' Whether per-trial `$latent$m` / `$g` Gaussian mixtures are materialized
#' @noRd
interactmeta_bayesmeta_trial_latent_ready <- function(fit) {
  if (!interactmeta_is_bayesian_fit(fit) || !is.list(fit$latent)) {
    return(FALSE)
  }
  lat <- fit$latent
  if (isTRUE(lat$.lazy_trials)) {
    return(FALSE)
  }
  trials <- lat$trials
  if (is.null(trials) || !length(trials)) {
    trials <- names(lat$m)
  }
  if (is.null(trials) || !length(trials)) {
    return(FALSE)
  }
  tr <- trials[[1L]]
  gm_m <- lat$m[[tr]]$gm
  gm_g <- lat$g[[tr]]$gm
  is.list(gm_m) &&
    length(gm_m$mean) &&
    is.list(gm_g) &&
    length(gm_g$mean)
}

#' Build or extend `$latent` after a lazy [cams()] fit
#' @param full_trials If `TRUE`, ensure per-trial convolve marginals exist.
#' @noRd
interactmeta_ensure_bayesmeta_latent <- function(fit, full_trials = FALSE) {
  if (!inherits(fit, "interactmeta") || !interactmeta_is_bayesian_fit(fit)) {
    return(invisible(fit))
  }
  lat <- fit$latent
  bm <- fit$bayesmeta
  if (is.null(bm$overall) || is.null(bm$interaction) || is.null(bm$trial_blocks)) {
    return(invisible(fit))
  }
  if (isTRUE(full_trials)) {
    if (is.list(lat) && !isTRUE(lat$.lazy_trials)) {
      return(invisible(fit))
    }
    mode <- "full"
  } else {
    if (is.list(lat) && !is.null(lat$beta)) {
      return(invisible(fit))
    }
    mode <- "none"
  }
  interval.type <- interactmeta_match_interval_type(
    fit$interval.type %||% if (is.list(lat)) lat$interval.type else NULL
  )
  fit$latent <- interactmeta_bayesmeta_latent(
    overall_fit = bm$overall,
    interaction_fit = bm$interaction,
    blocks = bm$trial_blocks,
    interval.type = interval.type,
    trial_marginals = mode
  )
  fit
}

#' Gaussian mixture for \eqn{M(\pi)=\alpha+\pi\theta} on the bmr DIRECT support
#' @noRd
interactmeta_bmr_population_linpred_gm <- function(overall_fit, pi) {
  pi <- as.numeric(pi)[[1L]]
  if (!is.finite(pi)) {
    stop("`pi` must be finite.", call. = FALSE)
  }
  tab <- interactmeta_bmr_support_table(overall_fit)
  effect <- as.numeric(tab$intercept) + pi * as.numeric(tab$ifrac)
  w <- as.numeric(tab$weight)
  w <- w / sum(w)
  nms <- colnames(overall_fit$support$mean)
  i_a <- match("intercept", nms)
  if (is.na(i_a)) {
    i_a <- 1L
  }
  i_t <- match("ifrac", nms)
  if (is.na(i_t)) {
    i_t <- min(2L, length(nms))
  }
  cond_sd <- rep(0, length(effect))
  cov_arr <- overall_fit$support$covariance
  if (!is.null(cov_arr) && length(dim(cov_arr)) == 3L) {
    va <- pmax(as.numeric(cov_arr[, i_a, i_a]), 0)
    vt <- pmax(as.numeric(cov_arr[, i_t, i_t]), 0)
    cov_at <- as.numeric(cov_arr[, i_a, i_t])
    cond_sd <- sqrt(pmax(va + (pi^2) * vt + 2 * pi * cov_at, 0))
  }
  interactmeta_gm(mean = effect, sd = cond_sd, weight = w)
}

#' Normalize bmr$support (list) to a data frame
#' @noRd
interactmeta_bmr_support_table <- function(overall_fit) {
  sup <- overall_fit$support
  if (is.null(sup)) {
    stop("`bmr` fit has no `$support` grid.", call. = FALSE)
  }
  if (is.list(sup) && !is.data.frame(sup) && !is.matrix(sup)) {
    mean_mat <- as.matrix(sup$mean)
    out <- data.frame(
      weight = as.numeric(sup$weight),
      tau = as.numeric(sup$tau),
      stringsAsFactors = FALSE
    )
    for (nm in colnames(mean_mat)) {
      out[[nm]] <- as.numeric(mean_mat[, nm])
    }
    return(out)
  }
  as.data.frame(sup)
}

#' Normalize bayesmeta$support (matrix) to a data frame
#' @noRd
interactmeta_bayesmeta_support_table <- function(interaction_fit) {
  sup <- interaction_fit$support
  if (is.null(sup)) {
    stop("`bayesmeta` fit has no `$support` grid.", call. = FALSE)
  }
  out <- as.data.frame(sup)
  if (!"mu" %in% names(out) && "mean" %in% names(out)) {
    out$mu <- out$mean
  }
  out
}

#' Pair overall support rows to interaction rows by cumulative weight
#'
#' Comonotonic coupling of independent blocks. Kept only as a diagnostic of
#' the previous reconstruction; the fit path uses a product of weights.
#' @noRd
interactmeta_support_pair_indices <- function(w_ov, w_int) {
  w_ov <- as.numeric(w_ov)
  w_int <- as.numeric(w_int)
  w_ov <- w_ov / sum(w_ov)
  w_int <- w_int / sum(w_int)
  cdf_ov <- cumsum(w_ov) - w_ov / 2
  cdf_int <- cumsum(w_int)
  vapply(cdf_ov, function(p) {
    which.min(abs(cdf_int - p))[[1L]]
  }, integer(1))
}

#' Fixed-effect linear predictor matrix (nodes × observations)
#' @noRd
interactmeta_bayesmeta_mu_matrix <- function(draws, data_fit, design) {
  nonref_fit <- design$coef_map$subgroup_fit[[2L]]
  ifrac_std <- design$coef_map$ifrac_std[[1L]]
  b_sg <- paste0("b_subgroup", nonref_fit)
  alpha <- draws$b_Intercept
  beta <- draws$b_ifrac_1
  gamma <- draws[[b_sg]]
  pi <- as.numeric(data_fit[[ifrac_std]])
  sg <- as.character(data_fit[[design$subgroup_var]])
  is_nonref <- sg == nonref_fit

  n_obs <- nrow(data_fit)
  n_support <- length(alpha)
  mu <- matrix(NA_real_, nrow = n_support, ncol = n_obs)
  for (j in seq_len(n_obs)) {
    mu[, j] <- alpha + beta * pi[[j]] + gamma * as.numeric(is_nonref[[j]])
  }
  mu
}

#' Parameter matrix from each stored fit's `$post.interval` and from convolve().
#'
#' Mean, sd, and mode of a one-fit margin come from that object's `$summary`
#' (those rows are not the probability interval). Bounds are `$post.interval`
#' at `prob`, not the copied 95% rows. The ecological column is the shortest
#' or central interval on `convolve(theta, -gamma)`.
#' @noRd
interactmeta_bayesmeta_param_summary <- function(fit, prob = 0.95) {
  interactmeta_ensure_bayesmeta_latent(fit, full_trials = FALSE)
  prob <- as.numeric(prob)[[1L]]
  rep <- interactmeta_reporting_sidecar(fit) %||% fit$reporting
  if (is.list(rep) && !is.null(rep$parameters) &&
        isTRUE(all.equal(as.numeric(rep$prob), prob))) {
    return(rep$parameters)
  }
  bm <- fit$bayesmeta
  if (is.null(bm$overall) || is.null(bm$interaction)) {
    return(NULL)
  }
  ov_sum <- bm$overall$summary
  int_sum <- bm$interaction$summary
  if (is.null(ov_sum) || is.null(int_sum)) {
    return(NULL)
  }
  if (is.null(fit$latent$beta) || is.null(fit$latent$beta$convolve)) {
    stop(
      "Bayesian ecological interval requires latent$beta from bayesmeta::convolve(theta, -gamma).",
      call. = FALSE
    )
  }
  interval.type <- interactmeta_match_interval_type(
    fit$interval.type %||% fit$latent$interval.type %||% bm$overall$interval.type
  )
  intercept_nm <- interactmeta_bmr_coef_name(bm$overall, "intercept")
  native_col <- function(summary, object, column, margin, which = NULL) {
    bounds <- interactmeta_post_interval(
      object,
      prob = prob,
      interval.type = interval.type,
      margin = margin,
      which = which
    )
    med <- interactmeta_native_quantile(
      object,
      0.5,
      margin = margin,
      which = which
    )[[1L]]
    c(
      mode = as.numeric(summary["mode", column])[[1L]],
      median = med,
      mean = as.numeric(summary["mean", column])[[1L]],
      sd = as.numeric(summary["sd", column])[[1L]],
      lower = bounds[[1L]],
      upper = bounds[[2L]]
    )
  }

  beta_sm <- interactmeta_posterior_1d_summary(
    fit$latent$beta,
    prob = prob,
    interval.type = interval.type
  )
  beta_row <- c(
    mode = NA_real_,
    median = beta_sm$median,
    mean = beta_sm$mean,
    sd = beta_sm$sd,
    lower = beta_sm$lower,
    upper = beta_sm$upper
  )

  alpha_row <- native_col(ov_sum, bm$overall, intercept_nm, "effect", intercept_nm)
  gamma_row <- native_col(int_sum, bm$interaction, "mu", "effect")
  tau_m <- native_col(ov_sum, bm$overall, "tau", "tau")
  tau_g <- native_col(int_sum, bm$interaction, "tau", "tau")

  st <- attr(fit, "interactmeta", exact = TRUE)
  nonref <- st$coef_map$ifrac_orig[[1L]] %||%
    (st$ifrac_vars[[1L]] %||% "Older")
  nonref <- sub("^ifrac_", "", as.character(nonref)[[1L]])
  eco_nm <- paste0("Ecological_", nonref)
  int_nm <- paste0("Interaction_", nonref)

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
  mat <- cbind(
    tau = tau_m,
    tau_gamma = tau_g,
    Intercept = alpha_row,
    Ecological = beta_row,
    Interaction = gamma_row
  )
  colnames(mat)[4:5] <- c(eco_nm, int_nm)
  rownames(mat) <- c("mode", "median", "mean", "sd", lower_lab, upper_lab)
  attr(mat, "interval.type") <- interval.type
  if (isTRUE(all.equal(prob, 0.95))) {
    if (is.null(fit$reporting) || !is.list(fit$reporting)) {
      fit$reporting <- list()
    }
    fit$reporting$prob <- prob
    fit$reporting$parameters <- mat
    fit$reporting$interval.type <- interval.type
  }
  mat
}

#' Weighted sample quantiles (deterministic)
#' @noRd
interactmeta_weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) {
    return(stats::setNames(rep(NA_real_, length(probs)), probs))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  cw <- cumsum(w)
  vapply(probs, function(p) {
    x[[which.max(cw >= p)]]
  }, numeric(1))
}

#' @noRd
interactmeta_native_location <- function(object,
                                        prob,
                                        interval.type,
                                        margin = "effect",
                                        which = NULL,
                                        predict = FALSE,
                                        individual = FALSE) {
  bounds <- interactmeta_post_interval(
    object,
    prob = prob,
    interval.type = interval.type,
    margin = margin,
    which = which,
    predict = predict,
    individual = individual
  )
  list(
    median = interactmeta_native_quantile(
      object, 0.5, margin = margin, which = which,
      predict = predict, individual = individual
    )[[1L]],
    lower = bounds[[1L]],
    upper = bounds[[2L]]
  )
}

#' Linear predictor `alpha + pi * theta` on the joint bmr support.
#' @noRd
interactmeta_bmr_linpred_fns <- function(overall_fit, pi, mean = TRUE) {
  x <- c(1, as.numeric(pi)[[1L]])
  list(
    density = function(z, log = FALSE) {
      as.numeric(overall_fit$dpredict(theta = z, x = x, mean = mean, log = log))
    },
    cdf = function(z) {
      as.numeric(overall_fit$ppredict(theta = z, x = x, mean = mean))
    },
    quantile = function(p) {
      as.numeric(overall_fit$qpredict(p = p, x = x, mean = mean))
    }
  )
}

#' Marginal law of `scale * u`, `u | tau ~ N(0, tau^2)`, from his tau support.
#' @noRd
interactmeta_zero_re_fns <- function(interaction_fit, scale) {
  sup <- as.data.frame(interaction_fit$support)
  tau <- pmax(as.numeric(sup[, "tau"]), 0)
  w <- as.numeric(sup[, "weight"])
  w <- w / sum(w)
  sd <- interactmeta_gm_sd_eval(abs(as.numeric(scale)[[1L]]) * tau)
  logw <- log(w)
  list(
    density = function(x, log = FALSE) {
      x <- as.numeric(x)
      out <- vapply(x, function(z) {
        interactmeta_log_sum_exp(logw + stats::dnorm(z, mean = 0, sd = sd, log = TRUE))
      }, numeric(1))
      if (isTRUE(log)) out else exp(out)
    },
    cdf = function(x) {
      x <- as.numeric(x)
      mat <- matrix(NA_real_, nrow = length(x), ncol = length(w))
      for (j in seq_along(w)) {
        mat[, j] <- stats::pnorm(x, mean = 0, sd = sd[[j]])
      }
      as.numeric(mat %*% w)
    }
  )
}

#' @noRd
interactmeta_bayesmeta_population_location <- function(fit,
                                                      pi,
                                                      gamma_scale = 0,
                                                      re_scale = 0,
                                                      predictive = FALSE,
                                                      prob = 0.95,
                                                      interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% fit$interval.type
  )
  pi <- as.numeric(pi)[[1L]]
  gamma_scale <- as.numeric(gamma_scale)[[1L]]
  re_scale <- as.numeric(re_scale)[[1L]]
  predictive <- isTRUE(predictive)
  prob <- as.numeric(prob)[[1L]]

  ov <- fit$bayesmeta$overall
  int <- fit$bayesmeta$interaction
  grid <- interactmeta_bayesmeta_grid(ov, int)
  pure <- abs(gamma_scale) < 1e-15 && abs(re_scale) < 1e-15
  if (pure) {
    x <- c(1, pi)
    bounds <- ov$pred.interval(
      level = prob,
      x = x,
      mean = !predictive,
      method = interval.type
    )
    med <- as.numeric(ov$qpredict(p = 0.5, x = x, mean = !predictive))[[1L]]
    return(list(
      median = med,
      lower = unname(as.numeric(bounds)[[1L]]),
      upper = unname(as.numeric(bounds)[[2L]])
    ))
  }
  key <- sprintf(
    "pop|%s|%d|%.8g|%.8g|%.8g|%.8g",
    interval.type, as.integer(predictive), pi, gamma_scale, re_scale, prob
  )
  cache <- fit$latent$cache
  if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  loc <- interactmeta_bmr_linpred_fns(ov, pi, mean = !predictive)
  if (abs(gamma_scale) >= 1e-15) {
    if (predictive) {
      g_dens <- function(z, log = FALSE) {
        as.numeric(int$dposterior(mu = z, predict = TRUE, log = log))
      }
      g_cdf <- function(z) as.numeric(int$pposterior(mu = z, predict = TRUE))
    } else {
      g_dens <- interactmeta_bayesmeta_mu_density(int)
      g_cdf <- interactmeta_bayesmeta_mu_cdf(int)
    }
    dens2 <- interactmeta_scale_density(g_dens, gamma_scale)
    cdf2 <- interactmeta_scale_cdf(g_cdf, gamma_scale)
  } else {
    re <- interactmeta_zero_re_fns(int, re_scale)
    dens2 <- re$density
    cdf2 <- re$cdf
  }
  conv <- interactmeta_convolve_sum(
    dens1 = loc$density,
    dens2 = dens2,
    cdf1 = loc$cdf,
    cdf2 = cdf2,
    delta = grid$delta,
    epsilon = grid$epsilon
  )
  bounds <- interactmeta_interval_from_distribution(
    conv$quantile, conv$cdf, prob, interval.type, grid$tol.uniroot
  )
  out <- list(
    median = as.numeric(conv$quantile(0.5))[[1L]],
    lower = bounds[[1L]],
    upper = bounds[[2L]]
  )
  if (is.environment(cache)) {
    assign(key, out, envir = cache)
  }
  out
}

# ---- interactmeta_quadrature.R ----
#' Independent-block quadrature helpers for binary bayesmeta CAMS
#'
#' DIRECT support from `bmr` / `bayesmeta` is a mixture of normals (plus a
#' discrete distribution on tau). Independent blocks combine by a product of
#' mixture weights. That Gaussian-mixture product is an exact check on sums of
#' independent blocks. The published 1-D law of those sums is
#' `bayesmeta::convolve()` of his density and cdf, and the published interval
#' is his shortest or central rule on that convolution (or `$post.interval`
#' when the margin lives on one stored fit).
#'
#' These objects are **marginal** posteriors. Trials that share a block also
#' share hyperparameters and are not mutually independent.
#'
#' @noRd
NULL

interactmeta_log_sum_exp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(-Inf)
  }
  m <- max(x)
  m + log(sum(exp(x - m)))
}

interactmeta_gm <- function(mean, sd, weight) {
  mean <- as.numeric(mean)
  sd <- pmax(as.numeric(sd), 0)
  weight <- as.numeric(weight)
  n <- length(mean)
  if (length(sd) != n || length(weight) != n) {
    stop("Gaussian mixture `mean`, `sd`, and `weight` must have equal length.", call. = FALSE)
  }
  ok <- is.finite(mean) & is.finite(sd) & is.finite(weight) & weight > 0
  if (!any(ok)) {
    stop("Gaussian mixture has no finite positive-weight components.", call. = FALSE)
  }
  weight <- weight[ok]
  weight <- weight / sum(weight)
  structure(
    list(mean = mean[ok], sd = sd[ok], weight = weight),
    class = "interactmeta_gm"
  )
}

interactmeta_gm_sd_eval <- function(sd, min_sd = 1e-12) {
  pmax(as.numeric(sd), min_sd)
}

interactmeta_gm_density <- function(gm, x, log = FALSE) {
  x <- as.numeric(x)
  sd <- interactmeta_gm_sd_eval(gm$sd)
  logw <- log(gm$weight)
  out <- vapply(x, function(z) {
    interactmeta_log_sum_exp(logw + stats::dnorm(z, mean = gm$mean, sd = sd, log = TRUE))
  }, numeric(1))
  if (isTRUE(log)) {
    out
  } else {
    exp(out)
  }
}

interactmeta_gm_cdf <- function(gm, x) {
  x <- as.numeric(x)
  sd <- interactmeta_gm_sd_eval(gm$sd)
  n <- length(x)
  k <- length(gm$weight)
  cdf <- matrix(NA_real_, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    cdf[, j] <- stats::pnorm(x, mean = gm$mean[[j]], sd = sd[[j]])
  }
  as.numeric(cdf %*% gm$weight)
}

interactmeta_gm_quantile <- function(gm, probs) {
  probs <- as.numeric(probs)
  mu <- sum(gm$weight * gm$mean)
  v <- sum(gm$weight * (gm$sd^2 + (gm$mean - mu)^2))
  s <- sqrt(max(v, 0))
  span <- max(s, 1e-8)
  vapply(probs, function(p) {
    if (!is.finite(p)) {
      return(NA_real_)
    }
    if (p <= 0) {
      return(-Inf)
    }
    if (p >= 1) {
      return(Inf)
    }
    lo <- mu - 8 * span
    hi <- mu + 8 * span
    f_lo <- interactmeta_gm_cdf(gm, lo) - p
    f_hi <- interactmeta_gm_cdf(gm, hi) - p
    tries <- 0L
    while (f_lo > 0 && tries < 20L) {
      lo <- lo - (hi - lo + span)
      f_lo <- interactmeta_gm_cdf(gm, lo) - p
      tries <- tries + 1L
    }
    tries <- 0L
    while (f_hi < 0 && tries < 20L) {
      hi <- hi + (hi - lo + span)
      f_hi <- interactmeta_gm_cdf(gm, hi) - p
      tries <- tries + 1L
    }
    stats::uniroot(
      function(z) interactmeta_gm_cdf(gm, z) - p,
      lower = lo,
      upper = hi,
      tol = 1e-8
    )$root
  }, numeric(1))
}

interactmeta_gm_mean <- function(gm) {
  sum(gm$weight * gm$mean)
}

interactmeta_gm_var <- function(gm) {
  mu <- interactmeta_gm_mean(gm)
  sum(gm$weight * (gm$sd^2 + (gm$mean - mu)^2))
}

interactmeta_gm_summary <- function(gm, prob = 0.95,
                                    interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(interval.type)
  mu <- interactmeta_gm_mean(gm)
  s <- sqrt(max(interactmeta_gm_var(gm), 0))
  median <- interactmeta_gm_quantile(gm, 0.5)
  bounds <- interactmeta_interval_from_distribution(
    quantile = function(p) interactmeta_gm_quantile(gm, p),
    cdf = function(x) interactmeta_gm_cdf(gm, x),
    prob = prob,
    interval.type = interval.type
  )
  qs <- c(bounds[[1L]], median[[1L]], bounds[[2L]])
  list(
    mean = mu,
    sd = s,
    median = unname(median[[1L]]),
    lower = unname(qs[[1L]]),
    upper = unname(qs[[3L]]),
    q = stats::setNames(qs, c("lower", "median", "upper")),
    interval.type = interval.type
  )
}

#' Mean and SD of exp(X) when X is a Gaussian mixture (mixture of lognormals)
#' @noRd
interactmeta_gm_lognormal_moments <- function(gm) {
  mean_y <- sum(gm$weight * exp(gm$mean + 0.5 * gm$sd^2))
  mean_y2 <- sum(gm$weight * exp(2 * gm$mean + 2 * gm$sd^2))
  list(
    mean = mean_y,
    sd = sqrt(max(mean_y2 - mean_y^2, 0))
  )
}

#' Linear combination aX + bY of **independent** Gaussian mixtures
#'
#' The product of mixture weights assumes `gm_x` and `gm_y` come from
#' independent blocks. Valid for cross-block CAMS quantities such as
#' `theta - gamma` and trial-level `m_j ± c g_j`.
#'
#' Do **not** use this to form `M_star = alpha + pi_star * theta`. `alpha`
#' and `theta` share the overall-block DIRECT support and can be
#' correlated; that linear combination must keep the paired support
#' (and `Cov(alpha, theta)`), not a product of marginal mixtures.
#' @noRd
interactmeta_gm_linear_sum <- function(gm_x, a, gm_y, b) {
  interactmeta_gm_assert_independent(gm_x, gm_y)
  a <- as.numeric(a)[[1L]]
  b <- as.numeric(b)[[1L]]
  if (!is.finite(a) || !is.finite(b)) {
    stop("Linear combination coefficients must be finite.", call. = FALSE)
  }
  ix <- seq_along(gm_x$weight)
  iy <- seq_along(gm_y$weight)
  i <- rep(ix, each = length(iy))
  j <- rep(iy, times = length(ix))
  mean <- a * gm_x$mean[i] + b * gm_y$mean[j]
  sd <- sqrt((a * gm_x$sd[i])^2 + (b * gm_y$sd[j])^2)
  weight <- gm_x$weight[i] * gm_y$weight[j]
  interactmeta_gm(mean = mean, sd = sd, weight = weight)
}

#' Tag a Gaussian mixture with its CAMS block (`"overall"` / `"interaction"`)
#' @noRd
interactmeta_gm_tag_block <- function(gm, block) {
  attr(gm, "block") <- block
  gm
}

#' Stop if both mixtures are tagged as the same CAMS block
#' @noRd
interactmeta_gm_assert_independent <- function(gm_x, gm_y) {
  bx <- attr(gm_x, "block", exact = TRUE)
  by <- attr(gm_y, "block", exact = TRUE)
  if (is.null(bx) || is.null(by) || !identical(bx, by)) {
    return(invisible(NULL))
  }
  stop(
    "`interactmeta_gm_linear_sum()` assumes independent mixtures. ",
    "Do not combine two components from the same CAMS block ",
    "(e.g. `alpha` and `theta` from the overall `bmr` support). ",
    "Form `M* = alpha + pi* theta` from the joint support, not by ",
    "convolving the two marginals.",
    call. = FALSE
  )
}

#' Discrete (point-mass) distribution on a support, e.g. tau nodes
#' @noRd
interactmeta_discrete <- function(x, weight) {
  x <- as.numeric(x)
  weight <- as.numeric(weight)
  ok <- is.finite(x) & is.finite(weight) & weight > 0
  weight <- weight[ok]
  weight <- weight / sum(weight)
  list(x = x[ok], weight = weight)
}

interactmeta_discrete_summary <- function(dist, prob = 0.95) {
  mu <- sum(dist$weight * dist$x)
  s <- sqrt(max(sum(dist$weight * (dist$x - mu)^2), 0))
  qs <- interactmeta_weighted_quantile(
    dist$x,
    dist$weight,
    probs = c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2)
  )
  list(
    mean = mu,
    sd = s,
    median = unname(qs[[2L]]),
    lower = unname(qs[[1L]]),
    upper = unname(qs[[3L]]),
    q = qs
  )
}

#' Density / CDF of Z = c * X, with Jacobian. c = 0 is degenerate at 0.
#' @noRd
interactmeta_call_density <- function(dens, x, log = FALSE) {
  if ("log" %in% names(formals(dens))) {
    dens(x, log = log)
  } else if (isTRUE(log)) {
    log(dens(x))
  } else {
    dens(x)
  }
}

interactmeta_scale_density <- function(dens, c) {
  c <- as.numeric(c)[[1L]]
  if (!is.finite(c)) {
    stop("Scale coefficient must be finite.", call. = FALSE)
  }
  if (abs(c) < 1e-15) {
    stop("Scale coefficient is zero; use a degenerate point mass at 0.", call. = FALSE)
  }
  function(x, log = FALSE) {
    z <- as.numeric(x) / c
    if (isTRUE(log)) {
      as.numeric(interactmeta_call_density(dens, z, log = TRUE) - log(abs(c)))
    } else {
      as.numeric(interactmeta_call_density(dens, z, log = FALSE) / abs(c))
    }
  }
}

interactmeta_scale_cdf <- function(cdf, c) {
  c <- as.numeric(c)[[1L]]
  if (!is.finite(c)) {
    stop("Scale coefficient must be finite.", call. = FALSE)
  }
  if (abs(c) < 1e-15) {
    stop("Scale coefficient is zero; use a degenerate point mass at 0.", call. = FALSE)
  }
  if (c > 0) {
    function(x) as.numeric(cdf(as.numeric(x) / c))
  } else {
    function(x) as.numeric(1 - cdf(as.numeric(x) / c))
  }
}

interactmeta_negate_density <- function(dens) {
  interactmeta_scale_density(dens, -1)
}

interactmeta_negate_cdf <- function(cdf) {
  interactmeta_scale_cdf(cdf, -1)
}

#' Normal-normal shrinkage conditional on hyperparameters
#'
#' tau = 0 is handled without 1/tau^2: the latent equals mu0 with variance 0
#' given the hyperparameters.
#' @noRd
interactmeta_normal_normal_shrink <- function(y, sigma2, mu0, tau) {
  y <- as.numeric(y)
  sigma2 <- as.numeric(sigma2)
  mu0 <- as.numeric(mu0)
  tau <- as.numeric(tau)
  n <- max(length(y), length(sigma2), length(mu0), length(tau))
  y <- rep_len(y, n)
  sigma2 <- rep_len(sigma2, n)
  mu0 <- rep_len(mu0, n)
  tau <- rep_len(tau, n)
  mean <- mu0
  sd <- rep(0, n)
  pos <- is.finite(tau) & tau > 0 & is.finite(sigma2) & sigma2 > 0
  if (any(pos)) {
    prec_y <- 1 / sigma2[pos]
    prec_0 <- 1 / (tau[pos]^2)
    post_prec <- prec_y + prec_0
    sd[pos] <- sqrt(1 / post_prec)
    mean[pos] <- (prec_y * y[pos] + prec_0 * mu0[pos]) / post_prec
  }
  zero_sig <- is.finite(sigma2) & sigma2 <= 0 & is.finite(tau) & tau > 0
  if (any(zero_sig)) {
    mean[zero_sig] <- y[zero_sig]
    sd[zero_sig] <- 0
  }
  cbind(mean = mean, sd = sd)
}

#' Extract mixture nodes from a bayesmeta() support matrix (mu | tau)
#'
#' Untagged on purpose: independent univariate `bayesmeta` fits (DA
#' `mu_B - mu_A`) may be combined with `interactmeta_gm_linear_sum()`.
#' The CAMS interaction mean is tagged `"interaction"` where it is
#' constructed in `interactmeta_bayesmeta_latent()`.
#' @noRd
interactmeta_gm_from_bayesmeta <- function(fit) {
  sup <- as.data.frame(fit$support)
  if (!"mu" %in% names(sup) && "mean" %in% names(sup)) {
    sup$mu <- sup$mean
  }
  interactmeta_gm(
    mean = sup$mu,
    sd = sup$sd,
    weight = sup$weight
  )
}

#' Extract mixture nodes for one bmr() coefficient
#' @noRd
interactmeta_bmr_coef_name <- function(fit, role = c("intercept", "slope")) {
  role <- match.arg(role)
  nms <- colnames(fit$support$mean)
  if (is.null(nms) || !length(nms)) {
    stop("`bmr` support$mean has no column names.", call. = FALSE)
  }
  if (identical(role, "intercept")) {
    hit <- intersect(c("intercept", "Intercept", "(Intercept)"), nms)
    if (length(hit)) {
      return(hit[[1L]])
    }
    return(nms[[1L]])
  }
  hit <- intersect(c("ifrac", "slope"), nms)
  if (length(hit)) {
    return(hit[[1L]])
  }
  if (length(nms) < 2L) {
    stop("No slope column in `bmr` support$mean.", call. = FALSE)
  }
  nms[[2L]]
}

interactmeta_gm_from_bmr <- function(fit, which) {
  mu <- as.numeric(fit$support$mean[, which])
  var <- as.numeric(fit$support$covariance[, which, which])
  interactmeta_gm_tag_block(
    interactmeta_gm(
      mean = mu,
      sd = sqrt(pmax(var, 0)),
      weight = fit$support$weight
    ),
    "overall"
  )
}

#' Mixture for latent m_j from bmr$shrink.moments
#' @noRd
interactmeta_gm_from_bmr_shrink <- function(fit, which) {
  mom <- fit$shrink.moments(which = which)
  interactmeta_gm_tag_block(
    interactmeta_gm(
      mean = mom[, "mean"],
      sd = mom[, "sd"],
      weight = fit$support$weight
    ),
    "overall"
  )
}

#' Mixture for latent g_j from bayesmeta cond.moment(individual=)
#' @noRd
interactmeta_gm_from_bayesmeta_individual <- function(fit, which) {
  tau <- as.numeric(fit$support[, "tau"])
  w <- as.numeric(fit$support[, "weight"])
  cm <- fit$cond.moment(tau = tau, individual = which)
  if (is.null(dim(cm))) {
    stop("Unexpected `cond.moment()` result for an individual study.", call. = FALSE)
  }
  if (length(dim(cm)) == 3L) {
    cm <- cm[, , 1L]
  }
  interactmeta_gm_tag_block(
    interactmeta_gm(
      mean = cm[, "mean"],
      sd = cm[, "sd"],
      weight = w
    ),
    "interaction"
  )
}

#' Unary density/CDF wrappers for convolve()
#' @noRd
interactmeta_bmr_coef_density <- function(fit, which) {
  force(fit)
  force(which)
  function(x, log = FALSE) {
    as.numeric(fit$dposterior(beta = x, which.beta = which, log = log))
  }
}

interactmeta_bmr_coef_cdf <- function(fit, which) {
  force(fit)
  force(which)
  function(x) {
    as.numeric(fit$pposterior(beta = x, which.beta = which))
  }
}

interactmeta_bayesmeta_mu_density <- function(fit) {
  force(fit)
  function(x, log = FALSE) {
    as.numeric(fit$dposterior(mu = x, log = log))
  }
}

interactmeta_bayesmeta_mu_cdf <- function(fit) {
  force(fit)
  function(x) {
    as.numeric(fit$pposterior(mu = x))
  }
}

interactmeta_bmr_shrink_density <- function(fit, which) {
  force(fit)
  force(which)
  function(x, log = FALSE) {
    as.numeric(fit$dshrink(theta = x, which = which, log = log))
  }
}

interactmeta_bmr_shrink_cdf <- function(fit, which) {
  force(fit)
  force(which)
  function(x) {
    as.numeric(fit$pshrink(theta = x, which = which))
  }
}

interactmeta_bayesmeta_individual_density <- function(fit, which) {
  force(fit)
  force(which)
  function(x, log = FALSE) {
    as.numeric(fit$dposterior(mu = x, individual = which, log = log))
  }
}

interactmeta_bayesmeta_individual_cdf <- function(fit, which) {
  force(fit)
  force(which)
  function(x) {
    as.numeric(fit$pposterior(mu = x, individual = which))
  }
}

#' Convolution of two independent continuous distributions via bayesmeta::convolve
#'
#' Inputs are density/CDF of independent random variables X and Y.
#' Output is a DIRECT approximation of the law of X+Y. Settings `delta`
#' (grid spacing) and `epsilon` (tail mass) are numerical accuracy controls,
#' not model parameters.
#' @noRd
interactmeta_convolve_sum <- function(dens1,
                                      dens2,
                                      cdf1,
                                      cdf2,
                                      delta = 0.01,
                                      epsilon = 1e-4) {
  bayesmeta::convolve(
    dens1 = dens1,
    dens2 = dens2,
    cdf1 = cdf1,
    cdf2 = cdf2,
    delta = delta,
    epsilon = epsilon
  )
}

#' `"shortest"` or `"central"`. `"evidentiary"` is a single-model bayesmeta rule.
#' @noRd
interactmeta_match_interval_type <- function(interval.type = NULL) {
  if (is.null(interval.type) || !length(interval.type) ||
        !nzchar(as.character(interval.type)[[1L]])) {
    return("central")
  }
  interval.type <- as.character(interval.type)[[1L]]
  if (identical(interval.type, "evidentiary")) {
    stop(
      "`interval.type = \"evidentiary\"` is defined for a single bayesmeta model, not for a derived CAMS sum.",
      call. = FALSE
    )
  }
  if (!interval.type %in% c("shortest", "central")) {
    stop(
      "`interval.type` must be \"shortest\" or \"central\".",
      call. = FALSE
    )
  }
  interval.type
}

#' Stop before ifrac() replaces a non-positive variance with a placeholder.
#'
#' A supplied `vi <= 0` is not a missing subgroup. Metafor still uses the
#' existing drop-and-placeholder path.
#' @noRd
interactmeta_reject_nonpositive_sampling_vi <- function(x) {
  if (!is.data.frame(x) && !inherits(x, "escalc")) {
    return(invisible(NULL))
  }
  if (!"vi" %in% names(x)) {
    return(invisible(NULL))
  }
  interactmeta_sigma_from_vi(x$vi)
  invisible(NULL)
}

#' Sampling standard deviations passed to bmr() / bayesmeta().
#'
#' Variances must be finite and strictly positive. No `1e-12` floor.
#' @noRd
interactmeta_sigma_from_vi <- function(vi) {
  vi <- as.numeric(vi)
  if (!length(vi) || any(!is.finite(vi)) || any(vi <= 0)) {
    stop(
      "Sampling variances passed to bayesmeta must be finite and strictly positive.",
      call. = FALSE
    )
  }
  sqrt(vi)
}

#' bmr `$post.interval` / `$qposterior` index a coefficient by position.
#'
#' A character name reaches `sum(which.beta)` inside his argument check
#' (`&` does not short-circuit), so the name is resolved here first.
#' @noRd
interactmeta_bmr_which_index <- function(object, which) {
  if (is.null(which) || is.numeric(which)) {
    return(which)
  }
  nms <- colnames(object$support$mean)
  if (is.null(nms) || !length(nms)) {
    nms <- setdiff(colnames(object$summary), "tau")
  }
  hit <- match(as.character(which)[[1L]], nms)
  if (is.na(hit)) {
    stop(
      "No bmr coefficient named ",
      deparse(as.character(which)[[1L]]),
      ".",
      call. = FALSE
    )
  }
  hit
}

#' Call `$post.interval` with that object's own formals.
#'
#' `bmr` uses `beta.level` / `which.beta` (and `shrink.interval` for a
#' study margin). `bayesmeta` uses `mu.level` / `individual`. Tau uses
#' `tau.level` on either object. Effect intervals are not pinned at 0;
#' his tau method may be, and that stays inside his function.
#' @noRd
interactmeta_post_interval <- function(object,
                                      prob,
                                      interval.type = NULL,
                                      margin = c("effect", "tau", "shrink"),
                                      which = NULL,
                                      predict = FALSE,
                                      individual = FALSE) {
  margin <- match.arg(margin)
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% object$interval.type
  )
  prob <- as.numeric(prob)[[1L]]
  if (!is.finite(prob) || prob <= 0 || prob >= 1) {
    stop("`prob` must lie in (0, 1).", call. = FALSE)
  }
  if (identical(margin, "shrink") && is.function(object$shrink.interval)) {
    fn <- object$shrink.interval
    fml <- names(formals(fn))
    args <- list()
    if ("level" %in% fml) {
      args$level <- prob
    }
    if ("which" %in% fml) {
      args$which <- which
    }
    if ("method" %in% fml) {
      args$method <- interval.type
    }
    res <- do.call(fn, args)
    return(unname(as.numeric(res)[1:2]))
  }
  fn <- object$post.interval
  if (!is.function(fn)) {
    stop("Stored fit has no `$post.interval`.", call. = FALSE)
  }
  fml <- names(formals(fn))
  args <- list()
  if ("method" %in% fml) {
    args$method <- interval.type
  }
  if (identical(margin, "tau")) {
    args$tau.level <- prob
  } else if ("beta.level" %in% fml) {
    args$beta.level <- prob
    args$which.beta <- interactmeta_bmr_which_index(object, which)
  } else {
    args$mu.level <- prob
    if (isTRUE(predict) && "predict" %in% fml) {
      args$predict <- TRUE
    }
    if (!identical(individual, FALSE) && "individual" %in% fml) {
      args$individual <- individual
    }
  }
  res <- do.call(fn, args)
  unname(as.numeric(res)[1:2])
}

#' His quantile at one or more probabilities, using that object's formals.
#' @noRd
interactmeta_native_quantile <- function(object,
                                        p,
                                        margin = c("effect", "tau", "shrink"),
                                        which = NULL,
                                        predict = FALSE,
                                        individual = FALSE) {
  margin <- match.arg(margin)
  p <- as.numeric(p)
  if (identical(margin, "shrink") && is.function(object$qshrink)) {
    return(as.numeric(vapply(p, function(pp) {
      object$qshrink(p = pp, which = which)
    }, numeric(1))))
  }
  fn <- object$qposterior
  fml <- names(formals(fn))
  vapply(p, function(pp) {
    args <- list()
    if (identical(margin, "tau")) {
      args$tau.p <- pp
    } else if ("beta.p" %in% fml) {
      args$beta.p <- pp
      args$which.beta <- interactmeta_bmr_which_index(object, which)
    } else {
      args$mu.p <- pp
      if (isTRUE(predict) && "predict" %in% fml) {
        args$predict <- TRUE
      }
      if (!identical(individual, FALSE) && "individual" %in% fml) {
        args$individual <- individual
      }
    }
    as.numeric(do.call(fn, args))[[1L]]
  }, numeric(1))
}

#' Central or shortest interval from a posterior quantile function.
#'
#' Central is `quantile(c((1-p)/2, 1-(1-p)/2))`. Shortest minimizes
#' `q(a+p)-q(a)` over the lower-tail probability `a`. Effect parameters are
#' not forced to a boundary at 0 (his tau `$post.interval` may be; call that
#' instead of this helper).
#' @noRd
interactmeta_interval_from_distribution <- function(quantile,
                                                   cdf,
                                                   prob,
                                                   interval.type = NULL,
                                                   tol = NULL) {
  interval.type <- interactmeta_match_interval_type(interval.type)
  prob <- as.numeric(prob)[[1L]]
  if (is.null(tol) || !is.finite(tol) || tol <= 0) {
    tol <- 2^16 * .Machine$double.eps
  }
  if (identical(interval.type, "central")) {
    return(unname(as.numeric(quantile(c((1 - prob) / 2, 1 - (1 - prob) / 2)))))
  }
  force(cdf)
  tail_max <- 1 - prob
  tail_min <- tail_max / 50
  width <- function(a) {
    bounds <- as.numeric(quantile(c(a, a + prob)))
    bounds[[2L]] - bounds[[1L]]
  }
  tail <- stats::optimize(
    width,
    lower = tail_min,
    upper = tail_max,
    tol = tol
  )$minimum
  unname(as.numeric(quantile(c(tail, tail + prob))))
}

#' One-dimensional posterior (mixture check and/or published convolve / native)
#'
#' `gm` is the exact Gaussian-mixture check. It does not supply the published
#' interval when `convolve` or `published` is set. `published` is a list
#' `kind` (`"native"` or `"convolve"`), plus the stored fit and margin
#' selectors for `"native"`.
#' @noRd
interactmeta_posterior_1d <- function(gm = NULL,
                                      convolve = NULL,
                                      density = NULL,
                                      cdf = NULL,
                                      quantile = NULL,
                                      name = NULL,
                                      note = NULL,
                                      published = NULL,
                                      interval.type = NULL,
                                      tol.uniroot = NULL) {
  if (is.null(density)) {
    if (!is.null(convolve)) {
      density <- convolve$density
    } else if (!is.null(gm)) {
      density <- function(x, log = FALSE) interactmeta_gm_density(gm, x, log = log)
    }
  }
  if (is.null(cdf)) {
    if (!is.null(convolve)) {
      cdf <- convolve$cdf
    } else if (!is.null(gm)) {
      cdf <- function(x) interactmeta_gm_cdf(gm, x)
    }
  }
  if (is.null(quantile)) {
    if (!is.null(convolve)) {
      quantile <- convolve$quantile
    } else if (!is.null(published) && identical(published$kind, "native")) {
      quantile <- function(p) {
        interactmeta_native_quantile(
          published$object,
          p,
          margin = published$margin %||% "effect",
          which = published$which,
          predict = isTRUE(published$predict),
          individual = published$individual %||% FALSE
        )
      }
    } else if (!is.null(gm)) {
      quantile <- function(p) interactmeta_gm_quantile(gm, p)
    }
  }
  structure(
    list(
      name = name,
      gm = gm,
      convolve = convolve,
      density = density,
      cdf = cdf,
      quantile = quantile,
      note = note,
      published = published,
      interval.type = if (is.null(interval.type)) {
        NULL
      } else {
        interactmeta_match_interval_type(interval.type)
      },
      tol.uniroot = tol.uniroot
    ),
    class = "interactmeta_posterior_1d"
  )
}

#' Published summary of a 1-D posterior.
#'
#' Bounds come from `$post.interval` or from the shortest/central rule on
#' `bayesmeta::convolve()`. The Gaussian mixture, when stored, supplies the
#' mean and sd only.
#' @noRd
interactmeta_posterior_1d_summary <- function(x, prob = 0.95, interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% x$interval.type
  )
  tol <- x$tol.uniroot %||% (2^16 * .Machine$double.eps)
  use_published <- !is.null(x$published) || !is.null(x$convolve)
  if (!use_published) {
    if (!is.null(x$gm)) {
      return(interactmeta_gm_summary(
        x$gm,
        prob = prob,
        interval.type = interval.type
      ))
    }
    if (is.null(x$quantile)) {
      stop("1-D posterior has neither a Gaussian mixture nor a quantile function.", call. = FALSE)
    }
  }
  med <- NA_real_
  bounds <- c(NA_real_, NA_real_)
  if (!is.null(x$published) && identical(x$published$kind, "native")) {
    pub <- x$published
    med <- interactmeta_native_quantile(
      pub$object,
      0.5,
      margin = pub$margin %||% "effect",
      which = pub$which,
      predict = isTRUE(pub$predict),
      individual = pub$individual %||% FALSE
    )[[1L]]
    bounds <- interactmeta_post_interval(
      pub$object,
      prob = prob,
      interval.type = interval.type,
      margin = pub$margin %||% "effect",
      which = pub$which,
      predict = isTRUE(pub$predict),
      individual = pub$individual %||% FALSE
    )
  } else if (!is.null(x$convolve)) {
    med <- as.numeric(x$convolve$quantile(0.5))[[1L]]
    bounds <- interactmeta_interval_from_distribution(
      quantile = x$convolve$quantile,
      cdf = x$convolve$cdf,
      prob = prob,
      interval.type = interval.type,
      tol = tol
    )
  } else if (!is.null(x$quantile)) {
    qs <- x$quantile(c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2))
    med <- unname(qs[[2L]])
    bounds <- unname(qs[c(1L, 3L)])
  }
  mean_v <- NA_real_
  sd_v <- NA_real_
  if (!is.null(x$gm)) {
    mean_v <- interactmeta_gm_mean(x$gm)
    sd_v <- sqrt(max(interactmeta_gm_var(x$gm), 0))
  }
  list(
    mean = mean_v,
    sd = sd_v,
    median = med,
    lower = bounds[[1L]],
    upper = bounds[[2L]],
    q = stats::setNames(
      c(bounds[[1L]], med, bounds[[2L]]),
      c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2)
    )
  )
}

#' Weighted Pearson correlation of two series with the same weights
#' @noRd
interactmeta_weighted_cor <- function(x, y, w) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  w <- as.numeric(w)
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  cxy <- sum(w * (x - mx) * (y - my))
  sx <- sqrt(max(sum(w * (x - mx)^2), 0))
  sy <- sqrt(max(sum(w * (y - my)^2), 0))
  if (sx <= 0 || sy <= 0) {
    return(NA_real_)
  }
  cxy / (sx * sy)
}

#' Inverse CAMS rotation (trial-level)
#'
#' `mu_Aj = m_j - pi_j * g_j`
#' `mu_Bj = m_j + (1 - pi_j) * g_j`
#' These are trial-specific latent modifier-level effects at that trial's
#' `pi_j`, not population `mu_A(pi*)` / `mu_B(pi*)` from `modifiers()`.
#' @noRd
interactmeta_muA_from_mg <- function(m, g, pi) {
  as.numeric(m) - as.numeric(pi) * as.numeric(g)
}

interactmeta_muB_from_mg <- function(m, g, pi) {
  as.numeric(m) + (1 - as.numeric(pi)) * as.numeric(g)
}
