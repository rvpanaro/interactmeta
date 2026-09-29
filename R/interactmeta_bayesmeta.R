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
#' (`interval.type`, default `"shortest"`): margins from `$post.interval`,
#' independent-block sums from `bayesmeta::convolve()`. Public extractors and
#' [predict_cams_by_if()] / [plot_contribution()] instead use **equal-tail**
#' quantiles of the DIRECT-node mixture (fast). The mixture is not a second
#' native interval on `$latent`. `$draws` is the Cartesian product of DIRECT
#' tau-nodes with `.weight = w_i * v_k` (conditional means; within-node
#' variance is in `.cond_*` and `$latent`). `$mu` remains the **fixed-effect**
#' linear predictor. No `rposterior()` calls.
#'
#' @noRd
.interactmeta_fit_bayesmeta <- function(data_fit,
                                        design,
                                        trial_var = "trial",
                                        subgroup_var = "subgroup",
                                        priors = NULL,
                                        engine_args = list(),
                                        seed = 42L,
                                        silent = 2,
                                        ecological = TRUE,
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

  engine_args <- interactmeta_merge_engine_args(engine_args, list(...))
  if ("convolve_derived" %in% names(engine_args)) {
    stop(
      "`convolve_derived` cannot change reported intervals. Derived CAMS posteriors always come from bayesmeta::convolve() of independent blocks.",
      call. = FALSE
    )
  }
  interval.type <- interactmeta_match_interval_type(engine_args$interval.type)
  engine_args$interval.type <- NULL
  engine_args[c(
    "chains", "iter", "warmup", "thin", "refresh",
    "iter_warmup", "iter_sampling", "parallel_chains", "n_draws",
    "method", "cores", "control",
    "tau.prior", "tau.prior.overall", "tau.prior.interaction",
    "tau.interaction.prior", "fixed.prior", "beta.prior.mean",
    "beta.prior.sd", "beta.prior.cov", "mu.prior"
  )] <- NULL

  tau_overall <- prior_to_bayesmeta(priors$tau)
  tau_int <- prior_to_bayesmeta(priors$tau_interaction)
  fe <- prior_to_bayesmeta_fixed(priors$fixed)

  if (!isTRUE(ecological)) {
    return(interactmeta_fit_bayesmeta_constrained(
      data_fit = data_fit,
      design = design,
      blocks = blocks,
      priors = priors,
      fe = fe,
      tau_overall = tau_overall,
      tau_int = tau_int,
      engine_args = engine_args,
      interval.type = interval.type,
      silent = silent
    ))
  }

  bmr_args <- engine_args
  bm_args <- engine_args
  if (!is.null(fe$beta.prior.mean)) {
    bmr_args$beta.prior.mean <- unname(fe$beta.prior.mean)
    bmr_args$beta.prior.sd <- unname(fe$beta.prior.sd)
    bmr_args$beta.prior.cov <- fe$beta.prior.cov
  }
  if (!anyNA(fe$mu.prior)) {
    bm_args$mu.prior <- fe$mu.prior
  }

  overall_fit <- do.call(
    bayesmeta::bmr,
    c(
      list(
        y = blocks$overall$yi,
        sigma = interactmeta_sigma_from_vi(blocks$overall$vi),
        labels = as.character(blocks$overall$trial),
        X = cbind(
          intercept = 1,
          ifrac = blocks$overall$ifrac
        ),
        tau.prior = tau_overall,
        interval.type = interval.type
      ),
      bmr_args
    )
  )
  interaction_fit <- do.call(
    bayesmeta::bayesmeta,
    c(
      list(
        y = blocks$interaction$yi,
        sigma = interactmeta_sigma_from_vi(blocks$interaction$vi),
        labels = as.character(blocks$interaction$trial),
        tau.prior = tau_int,
        interval.type = interval.type
      ),
      bm_args
    )
  )

  draws <- interactmeta_bayesmeta_draws_from_support(
    overall_fit = overall_fit,
    interaction_fit = interaction_fit,
    design = design
  )
  mu <- interactmeta_bayesmeta_mu_matrix(draws, data_fit, design)
  latent <- interactmeta_bayesmeta_latent(
    overall_fit = overall_fit,
    interaction_fit = interaction_fit,
    blocks = blocks,
    interval.type = interval.type
  )

  if (!isTRUE(as.integer(silent)[[1L]] >= 2L)) {
    message(
      "interactmeta bayesmeta approach: orthogonal overall (bmr) + interaction ",
      "(bayesmeta) fits via independent product quadrature (",
      nrow(draws),
      " nodes)."
    )
  }

  structure(
    list(
      draws = draws,
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
      engine = "bayesmeta",
      interval.type = interval.type,
      prior = design$prior,
      priors = priors,
      algorithm = "quadrature",
      method = "bayes",
      ecological_slope = TRUE
    ),
    class = c("interactmeta", "cam", "cams", "list")
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
    stop("Overall and interaction fits stored different `$delta` values.", call. = FALSE)
  }
  if (!is.null(interaction_fit$epsilon) && !isTRUE(all.equal(epsilon, interaction_fit$epsilon))) {
    stop("Overall and interaction fits stored different `$epsilon` values.", call. = FALSE)
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
#' by [subgroups()]. Do not form \eqn{M^\star=\alpha+\pi^\star\theta} by
#' independently combining the `alpha` and `theta` mixture marginals.
#' @noRd
interactmeta_bayesmeta_latent <- function(overall_fit,
                                          interaction_fit,
                                          blocks,
                                          interval.type = "shortest") {
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
      note = "trial-specific latent; mu_Aj = m_j - pi_j * g_j via convolve() (not subgroups() at pi*)"
    )
    mu_B[[tr]] <- one(
      gm = gm_B,
      convolve = conv_B,
      published = pub_B,
      name = paste0("mu_B[", tr, "]"),
      note = "trial-specific latent; mu_Bj = m_j + (1 - pi_j) * g_j via convolve() (not subgroups() at pi*)"
    )
  }

  list(
    independence = "cross_block",
    interval.type = interval.type,
    delta = grid$delta,
    epsilon = grid$epsilon,
    tol.uniroot = grid$tol.uniroot,
    note = paste(
      "Marginal posteriors only. Overall and interaction blocks are independent.",
      "Trials within a block share hyperparameters and are not mutually independent.",
      "Published intervals on one stored fit use that object's $post.interval",
      "(or bmr$shrink.interval). Published intervals for independent-block sums",
      "use bayesmeta::convolve(). The Gaussian mixture is a check, not a second interval.",
      "alpha and theta are not convolved: they share the overall-block support.",
      "latent$mu_A[[j]] / $mu_B[[j]] are trial-specific (m_j, g_j, pi_j),",
      "not population mu_A(pi*) / mu_B(pi*); those come from subgroups()."
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
    cache = new.env(parent = emptyenv())
  )
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

#' Weighted posterior expectation over a support-derived draws table
#' @noRd
interactmeta_support_expectation <- function(draws, metric) {
  w <- draws$.weight
  if (is.null(w)) {
    w <- rep(1 / nrow(draws), nrow(draws))
  }
  w <- as.numeric(w)
  w <- w / sum(w)
  val <- metric(draws)
  sum(as.numeric(val) * w)
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
  n_draws <- length(alpha)
  mu <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
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

#' Mean and sd of a node-wise normal mixture, including conditional variance.
#'
#' Does not compute an interval. Published bounds come from `$post.interval`
#' or `bayesmeta::convolve()`.
#' @noRd
interactmeta_node_mean_sd <- function(mean, weight = NULL, cond_sd = NULL, exponentiate = FALSE) {
  mean <- as.numeric(mean)
  n <- length(mean)
  if (!n) {
    return(list(mean = NA_real_, sd = NA_real_))
  }
  if (is.null(weight) || length(weight) != n) {
    weight <- rep(1 / n, n)
  }
  if (is.null(cond_sd) || length(cond_sd) != n) {
    cond_sd <- rep(0, n)
  }
  gm <- interactmeta_gm(mean = mean, sd = pmax(as.numeric(cond_sd), 0), weight = weight)
  if (isTRUE(exponentiate)) {
    ln <- interactmeta_gm_lognormal_moments(gm)
    return(list(mean = ln$mean, sd = ln$sd))
  }
  list(
    mean = interactmeta_gm_mean(gm),
    sd = sqrt(max(interactmeta_gm_var(gm), 0))
  )
}

#' Grouped mean / sd for quadrature nodes. Intervals are filled by the caller.
#' @noRd
interactmeta_node_moment_grouped <- function(draws, group_cols, prob = 0.95, exponentiate = FALSE) {
  draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(d, ...) {
      mom <- interactmeta_node_mean_sd(
        d$effect,
        weight = if (".weight" %in% names(d)) d$.weight else NULL,
        cond_sd = if (".cond_sd" %in% names(d)) d$.cond_sd else NULL,
        exponentiate = exponentiate
      )
      tibble::tibble(
        estimate = NA_real_,
        mean = mom$mean,
        sd = mom$sd,
        lci = NA_real_,
        uci = NA_real_,
        prob = prob
      )
    }) |>
    dplyr::ungroup()
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

  if (isFALSE(fit$ecological_slope %||% fit$bayesmeta$ecological_slope)) {
    return(interactmeta_constrained_population_location(
      fit = fit,
      pi = pi,
      gamma_scale = gamma_scale,
      re_scale = re_scale,
      predictive = predictive,
      prob = prob,
      interval.type = interval.type
    ))
  }

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

#' Replace effect bounds with native bayesmeta / convolve intervals.
#'
#' Engine-native path kept for `$latent` / diagnostics. Public extractors use
#' equal-tail mixture summaries instead (`interactmeta_summarize_effects`).
#' Interaction uses `$post.interval`. Combined `M*` uses `bmr$pred.interval`
#' on the joint `(alpha, theta)` support. Subgroup effects convolve that
#' linear predictor with the independent gamma block. Alpha and theta are
#' not convolved with each other.
#' @noRd
interactmeta_bayesmeta_fill_effect_summary <- function(fit,
                                                      summary,
                                                      quantity,
                                                      prob,
                                                      exponentiate = FALSE,
                                                      omega = NULL) {
  if (!identical(fit$engine, "bayesmeta") || is.null(fit$bayesmeta$overall) ||
        is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  it <- interactmeta_match_interval_type(
    fit$interval.type %||% fit$bayesmeta$overall$interval.type
  )
  prob <- as.numeric(prob)[[1L]]
  pi_star <- if (is.null(omega)) NA_real_ else stats::median(as.numeric(omega), na.rm = TRUE)
  st <- attr(fit, "interactmeta", exact = TRUE)
  ref <- as.character(st$reference %||% "")[[1L]]
  for (i in seq_len(nrow(summary))) {
    ty <- if ("type" %in% names(summary)) as.character(summary$type[[i]]) else "Pooled"
    predictive <- identical(ty, "Prediction")
    loc <- NULL
    if (identical(quantity, "interaction")) {
      loc <- interactmeta_native_location(
        fit$bayesmeta$interaction,
        prob = prob,
        interval.type = it,
        predict = predictive
      )
    } else if (identical(quantity, "combined")) {
      if (!is.finite(pi_star)) {
        next
      }
      loc <- interactmeta_bayesmeta_population_location(
        fit,
        pi = pi_star,
        gamma_scale = 0,
        re_scale = if (predictive) pi_star else 0,
        predictive = predictive,
        prob = prob,
        interval.type = it
      )
    } else if (identical(quantity, "subgroups") && "subgroup" %in% names(summary)) {
      if (!is.finite(pi_star)) {
        next
      }
      sg <- as.character(summary$subgroup[[i]])
      loc <- interactmeta_bayesmeta_population_location(
        fit,
        pi = pi_star,
        gamma_scale = if (identical(sg, ref)) -pi_star else (1 - pi_star),
        re_scale = 0,
        predictive = predictive,
        prob = prob,
        interval.type = it
      )
    }
    if (is.null(loc)) {
      next
    }
    med <- loc$median
    lo <- loc$lower
    hi <- loc$upper
    if (isTRUE(exponentiate)) {
      med <- exp(med)
      lo <- exp(lo)
      hi <- exp(hi)
    }
    summary$estimate[[i]] <- med
    summary$lci[[i]] <- lo
    summary$uci[[i]] <- hi
  }
  attr(summary, "interval.type") <- it
  summary
}

#' Trial-level shrinkage bounds from his shrinkage interval or from convolve().
#' @noRd
interactmeta_bayesmeta_fill_shrinkage_summary <- function(fit,
                                                         summary,
                                                         type,
                                                         prob = 0.95,
                                                         exponentiate = FALSE) {
  if (!identical(fit$engine, "bayesmeta") || is.null(fit$latent) ||
        is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  it <- interactmeta_match_interval_type(fit$interval.type %||% fit$latent$interval.type)
  ref <- as.character(attr(fit, "interactmeta", exact = TRUE)$reference %||% "")[[1L]]
  for (i in seq_len(nrow(summary))) {
    tr <- as.character(summary$trial[[i]])
    post <- if (identical(type, "subgroups")) {
      sg <- as.character(summary$subgroup[[i]])
      if (identical(sg, ref)) fit$latent$mu_A[[tr]] else fit$latent$mu_B[[tr]]
    } else if (identical(type, "interaction")) {
      fit$latent$g[[tr]]
    } else if (identical(type, "overall")) {
      fit$latent$m[[tr]]
    } else {
      NULL
    }
    if (is.null(post)) {
      next
    }
    sm <- interactmeta_posterior_1d_summary(post, prob = prob, interval.type = it)
    med <- sm$median
    lo <- sm$lower
    hi <- sm$upper
    if (isTRUE(exponentiate)) {
      med <- exp(med)
      lo <- exp(lo)
      hi <- exp(hi)
    }
    summary$estimate[[i]] <- med
    summary$lci[[i]] <- lo
    summary$uci[[i]] <- hi
  }
  attr(summary, "interval.type") <- it
  summary
}
