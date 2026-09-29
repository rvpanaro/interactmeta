#' Joint (alpha, gamma) moments under ecological_slope = FALSE (beta = 0)
#'
#' Constrains \eqn{\theta = \gamma} so
#' \eqn{m_j = \alpha + \pi_j\gamma + \varepsilon_{mj}} and
#' \eqn{g_j = \gamma + \varepsilon_{gj}} with independent RE SDs
#' \eqn{\tau_m} and \eqn{\tau_g}. Across-trial heterogeneity on both
#' blocks remains free.
#' @noRd
interactmeta_joint_alpha_gamma_moments <- function(m,
                                                   vm,
                                                   pi,
                                                   g,
                                                   vg,
                                                   tau_m,
                                                   tau_g,
                                                   prior_mean = c(alpha = 0, gamma = 0),
                                                   prior_sd = c(alpha = Inf, gamma = Inf)) {
  m <- as.numeric(m)
  g <- as.numeric(g)
  vm <- as.numeric(vm)
  vg <- as.numeric(vg)
  pi <- as.numeric(pi)
  k <- length(m)
  if (length(g) != k || length(vm) != k || length(vg) != k || length(pi) != k) {
    stop("Constrained CAMS blocks must have equal trial length.", call. = FALSE)
  }
  tau_m <- max(as.numeric(tau_m)[[1L]], 0)
  tau_g <- max(as.numeric(tau_g)[[1L]], 0)
  vm_t <- vm + tau_m^2
  vg_t <- vg + tau_g^2
  if (any(!is.finite(vm_t) | vm_t <= 0) || any(!is.finite(vg_t) | vg_t <= 0)) {
    stop("Non-finite or non-positive block variances in constrained CAMS.", call. = FALSE)
  }

  # X' V^{-1} X for beta = (alpha, gamma)
  w_m <- 1 / vm_t
  w_g <- 1 / vg_t
  prec11 <- sum(w_m)
  prec12 <- sum(w_m * pi)
  prec22 <- sum(w_m * pi^2) + sum(w_g)
  xy1 <- sum(w_m * m)
  xy2 <- sum(w_m * pi * m) + sum(w_g * g)

  prior_mean <- as.numeric(prior_mean)[c(1L, 2L)]
  prior_sd <- as.numeric(prior_sd)[c(1L, 2L)]
  if (length(prior_mean) != 2L || length(prior_sd) != 2L) {
    stop("Constrained fixed prior needs length-2 mean and sd.", call. = FALSE)
  }
  pprec <- ifelse(is.finite(prior_sd) & prior_sd > 0, 1 / prior_sd^2, 0)
  prec11 <- prec11 + pprec[[1L]]
  prec22 <- prec22 + pprec[[2L]]
  xy1 <- xy1 + pprec[[1L]] * prior_mean[[1L]]
  xy2 <- xy2 + pprec[[2L]] * prior_mean[[2L]]

  det <- prec11 * prec22 - prec12^2
  if (!is.finite(det) || det <= 0) {
    stop("Constrained (alpha, gamma) precision is singular.", call. = FALSE)
  }
  cov <- matrix(
    c(prec22, -prec12, -prec12, prec11) / det,
    nrow = 2L,
    dimnames = list(c("alpha", "gamma"), c("alpha", "gamma"))
  )
  mean <- as.numeric(cov %*% c(xy1, xy2))
  names(mean) <- c("alpha", "gamma")
  list(mean = mean, cov = cov)
}

#' Build stacked trial data for constrained metafor CAMS (beta = 0)
#' @noRd
interactmeta_metafor_constrained_data <- function(blocks) {
  ov <- blocks$overall
  int <- blocks$interaction
  k <- nrow(ov)
  tr <- as.character(ov$trial)
  data.frame(
    yi = c(ov$yi, int$yi),
    vi = c(ov$vi, int$vi),
    trial = factor(c(tr, tr), levels = unique(tr)),
    block = factor(
      c(rep("overall", k), rep("interaction", k)),
      levels = c("overall", "interaction")
    ),
    X_alpha = c(rep(1, k), rep(0, k)),
    X_gamma = c(ov$ifrac, rep(1, k)),
    stringsAsFactors = FALSE
  )
}

#' Frequentist constrained CAMS (theta = gamma, beta = 0)
#'
#' Prefers a joint \code{rma.mv} DIAG model; falls back to a two-step REML/ML
#' fit (interaction first, then overall with IF offset) when the joint
#' optimizer does not converge (common for very small \eqn{k}).
#' @noRd
interactmeta_metafor_fit_constrained <- function(blocks,
                                                 method,
                                                 engine_args = list()) {
  dat <- interactmeta_metafor_constrained_data(blocks)
  args <- c(
    list(
      yi = dat$yi,
      V = dat$vi,
      mods = ~ X_alpha + X_gamma - 1,
      random = ~ block | trial,
      struct = "DIAG",
      method = method,
      data = dat,
      slab = paste(dat$trial, dat$block, sep = ":")
    ),
    engine_args
  )
  fit <- tryCatch(do.call(metafor::rma.mv, args), error = function(e) e)
  if (!inherits(fit, "error")) {
    cf <- stats::coef(fit)
    alpha <- as.numeric(cf[["X_alpha"]])
    gamma <- as.numeric(cf[["X_gamma"]])
    vb <- as.matrix(stats::vcov(fit))
    se_alpha <- sqrt(max(as.numeric(vb["X_alpha", "X_alpha"]), 0))
    se_gamma <- sqrt(max(as.numeric(vb["X_gamma", "X_gamma"]), 0))
    tau2 <- as.numeric(fit$tau2)
    if (length(tau2) < 2L) {
      tau2 <- c(tau2[[1L]], 0)
    }
    nms <- names(fit$tau2)
    if (!is.null(nms) && length(nms) >= 2L) {
      i_ov <- match("overall", nms)
      i_int <- match("interaction", nms)
      if (is.na(i_ov)) i_ov <- 1L
      if (is.na(i_int)) i_int <- 2L
      tau2_m <- max(as.numeric(tau2[[i_ov]]), 0)
      tau2_g <- max(as.numeric(tau2[[i_int]]), 0)
    } else {
      tau2_m <- max(as.numeric(tau2[[1L]]), 0)
      tau2_g <- max(as.numeric(tau2[[2L]]), 0)
    }
    return(list(
      joint = fit,
      alpha = alpha,
      gamma = gamma,
      theta = gamma,
      beta = 0,
      se_alpha = se_alpha,
      se_gamma = se_gamma,
      tau2_m = tau2_m,
      tau2_g = tau2_g,
      data = dat,
      method_path = "rma.mv"
    ))
  }

  # Two-step fallback: interaction RE-MA, then overall intercept with offset.
  ov <- blocks$overall
  int <- blocks$interaction
  interaction_fit <- interactmeta_rma(
    yi = int$yi,
    vi = int$vi,
    method = method,
    slab = as.character(int$trial),
    engine_args = engine_args
  )
  gamma <- as.numeric(interaction_fit$beta[[1L]])
  se_gamma <- as.numeric(interaction_fit$se[[1L]])
  y_adj <- as.numeric(ov$yi) - as.numeric(ov$ifrac) * gamma
  overall_fit <- interactmeta_rma(
    yi = y_adj,
    vi = ov$vi,
    method = method,
    slab = as.character(ov$trial),
    engine_args = engine_args
  )
  alpha <- as.numeric(overall_fit$beta[[1L]])
  se_alpha <- as.numeric(overall_fit$se[[1L]])
  list(
    joint = NULL,
    overall_step = overall_fit,
    interaction_step = interaction_fit,
    alpha = alpha,
    gamma = gamma,
    theta = gamma,
    beta = 0,
    se_alpha = se_alpha,
    se_gamma = se_gamma,
    tau2_m = max(as.numeric(overall_fit$tau2), 0),
    tau2_g = max(as.numeric(interaction_fit$tau2), 0),
    data = dat,
    method_path = "two_step"
  )
}

#' Bayesian constrained CAMS (beta = 0, theta = gamma)
#'
#' Interaction block is still a univariate bayesmeta fit (gives \eqn{\tau_g}
#' DIRECT support). Overall \eqn{\tau_m} support comes from bayesmeta on
#' offset-adjusted overall effects \eqn{m_j - \pi_j\hat\gamma}. At each
#' product node, \eqn{(\alpha,\gamma)} is the joint Gaussian posterior under
#' the stacked likelihood with free \eqn{\tau_m} and \eqn{\tau_g}. Setting
#' \eqn{\beta=0} does **not** remove across-trial heterogeneity
#' (\eqn{\tau_m}, \eqn{\tau_g} remain free).
#' @noRd
interactmeta_fit_bayesmeta_constrained <- function(data_fit,
                                                   design,
                                                   blocks,
                                                   priors,
                                                   fe,
                                                   tau_overall,
                                                   tau_int,
                                                   engine_args,
                                                   interval.type,
                                                   silent = 2) {
  bm_args <- engine_args
  if (!anyNA(fe$mu.prior)) {
    bm_args$mu.prior <- fe$mu.prior
  }
  ov <- blocks$overall
  int <- blocks$interaction
  sm <- interactmeta_sigma_from_vi(ov$vi)
  sg <- interactmeta_sigma_from_vi(int$vi)

  interaction_fit <- do.call(
    bayesmeta::bayesmeta,
    c(
      list(
        y = int$yi,
        sigma = sg,
        labels = as.character(int$trial),
        tau.prior = tau_int,
        interval.type = interval.type
      ),
      bm_args
    )
  )
  gamma_hat <- as.numeric(interaction_fit$summary["median", "mu"])[[1L]]
  if (!is.finite(gamma_hat)) {
    gamma_hat <- as.numeric(interaction_fit$summary["mean", "mu"])[[1L]]
  }

  y_adj <- as.numeric(ov$yi) - as.numeric(ov$ifrac) * gamma_hat
  overall_fit <- do.call(
    bayesmeta::bayesmeta,
    c(
      list(
        y = y_adj,
        sigma = sm,
        labels = as.character(ov$trial),
        tau.prior = tau_overall,
        interval.type = interval.type
      ),
      bm_args
    )
  )

  prior_mean <- c(alpha = 0, gamma = 0)
  prior_sd <- c(alpha = Inf, gamma = Inf)
  if (!anyNA(fe$mu.prior) && is.finite(fe$mu.prior[["sd"]])) {
    prior_mean[] <- as.numeric(fe$mu.prior[["mean"]])
    prior_sd[] <- as.numeric(fe$mu.prior[["sd"]])
  }

  draws <- interactmeta_bayesmeta_draws_constrained(
    overall_fit = overall_fit,
    interaction_fit = interaction_fit,
    blocks = blocks,
    design = design,
    prior_mean = prior_mean,
    prior_sd = prior_sd
  )
  mu <- interactmeta_bayesmeta_mu_matrix(draws, data_fit, design)
  latent <- interactmeta_bayesmeta_latent_constrained(
    interaction_fit = interaction_fit,
    interval.type = interval.type
  )

  if (!isTRUE(as.integer(silent)[[1L]] >= 2L)) {
    message(
      "interactmeta bayesmeta approach: constrained CAMS ",
      "(ecological_slope = FALSE; beta = 0, theta = gamma) with joint ",
      "(alpha, gamma) at product quadrature nodes (",
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
        ecological_slope = FALSE,
        gamma_offset = gamma_hat
      ),
      engine = "bayesmeta",
      interval.type = interval.type,
      prior = design$prior,
      priors = priors,
      algorithm = "quadrature",
      method = "bayes",
      ecological_slope = FALSE
    ),
    class = c("interactmeta", "cam", "cams", "list")
  )
}

#' Product DIRECT draws with joint (alpha, gamma) and beta = 0
#' @noRd
interactmeta_bayesmeta_draws_constrained <- function(overall_fit,
                                                     interaction_fit,
                                                     blocks,
                                                     design,
                                                     prior_mean,
                                                     prior_sd) {
  ov <- interactmeta_bayesmeta_support_table(overall_fit)
  int <- interactmeta_bayesmeta_support_table(interaction_fit)
  need <- c("weight", "tau")
  if (length(setdiff(need, names(ov))) || length(setdiff(need, names(int)))) {
    stop("Constrained bayesmeta support is missing weight/tau.", call. = FALSE)
  }
  n_ov <- nrow(ov)
  n_int <- nrow(int)
  i_ov <- rep(seq_len(n_ov), each = n_int)
  i_int <- rep(seq_len(n_int), times = n_ov)
  tau_m <- pmax(as.numeric(ov$tau[i_ov]), 0)
  tau_g <- pmax(as.numeric(int$tau[i_int]), 0)
  w <- as.numeric(ov$weight[i_ov]) * as.numeric(int$weight[i_int])
  w <- w / sum(w)

  m <- as.numeric(blocks$overall$yi)
  vm <- as.numeric(blocks$overall$vi)
  pi <- as.numeric(blocks$overall$ifrac)
  g <- as.numeric(blocks$interaction$yi)
  vg <- as.numeric(blocks$interaction$vi)

  n <- length(w)
  alpha <- rep(NA_real_, n)
  gamma <- rep(NA_real_, n)
  var_alpha <- rep(NA_real_, n)
  var_gamma <- rep(NA_real_, n)
  cov_ag <- rep(NA_real_, n)
  for (j in seq_len(n)) {
    mom <- interactmeta_joint_alpha_gamma_moments(
      m = m,
      vm = vm,
      pi = pi,
      g = g,
      vg = vg,
      tau_m = tau_m[[j]],
      tau_g = tau_g[[j]],
      prior_mean = prior_mean,
      prior_sd = prior_sd
    )
    alpha[[j]] <- mom$mean[["alpha"]]
    gamma[[j]] <- mom$mean[["gamma"]]
    var_alpha[[j]] <- mom$cov["alpha", "alpha"]
    var_gamma[[j]] <- mom$cov["gamma", "gamma"]
    cov_ag[[j]] <- mom$cov["alpha", "gamma"]
  }

  sg_nonref <- design$coef_map$subgroup_fit[[2L]]
  ifrac2_std <- design$coef_map$ifrac2_std[[1L]]
  out <- data.frame(
    .chain = 1L,
    .iteration = seq_len(n),
    .draw = seq_len(n),
    .weight = w,
    .cond_var_alpha = var_alpha,
    .cond_var_theta = var_gamma,
    .cond_cov_alpha_theta = cov_ag,
    .cond_sd_gamma = sqrt(pmax(var_gamma, 0)),
    b_Intercept = as.numeric(alpha),
    b_ifrac_1 = 0,
    sd_trial__Intercept = tau_m,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[[paste0("b_subgroup", sg_nonref)]] <- as.numeric(gamma)
  out[[paste0("sd_trial__", ifrac2_std)]] <- tau_g
  out
}

#' Minimal latent bundle when beta is fixed at 0
#' @noRd
interactmeta_bayesmeta_latent_constrained <- function(interaction_fit,
                                                      interval.type = "shortest") {
  interval.type <- interactmeta_match_interval_type(interval.type)
  gamma_iv <- as.numeric(interaction_fit$post.interval(
    mu.level = 0.95,
    method = interval.type
  ))
  list(
    ecological_slope = FALSE,
    beta = list(
      median = 0,
      lower = 0,
      upper = 0,
      note = paste(
        "ecological_slope = FALSE constrains beta = 0;",
        "across-trial tau_m and tau_g remain free"
      )
    ),
    gamma = list(
      median = as.numeric(interaction_fit$summary["median", "mu"]),
      lower = gamma_iv[[1L]],
      upper = gamma_iv[[2L]],
      note = "interaction-block bayesmeta post.interval"
    ),
    theta = list(
      note = "theta = gamma under ecological_slope = FALSE"
    ),
    alpha = list(
      note = "joint with gamma at each DIRECT product node; see draws"
    ),
    interval.type = interval.type
  )
}

#' Population location under ecological_slope = FALSE from DIRECT draws
#'
#' With \eqn{\beta=0} and \eqn{\theta=\gamma}, estimands reduce to linear
#' forms in the joint \eqn{(\alpha,\gamma)} node moments. Across-trial
#' \eqn{\tau_m} and \eqn{\tau_g} remain in the mixture.
#' @noRd
interactmeta_constrained_population_location <- function(fit,
                                                         pi,
                                                         gamma_scale = 0,
                                                         re_scale = 0,
                                                         predictive = FALSE,
                                                         prob = 0.95,
                                                         interval.type = "shortest") {
  force(interval.type)
  post <- as.data.frame(fit$draws)
  w <- if (".weight" %in% names(post)) {
    as.numeric(post$.weight)
  } else {
    rep(1, nrow(post))
  }
  w <- w / sum(w)
  alpha <- as.numeric(post$b_Intercept)
  sg <- grep("^b_subgroup", names(post), value = TRUE)
  if (!length(sg)) {
    stop("Constrained CAMS draws lack a subgroup interaction coefficient.", call. = FALSE)
  }
  gamma <- as.numeric(post[[sg[[1L]]]])
  coef_g <- as.numeric(pi)[[1L]] + as.numeric(gamma_scale)[[1L]]
  mean <- alpha + coef_g * gamma
  va <- if (".cond_var_alpha" %in% names(post)) as.numeric(post$.cond_var_alpha) else 0
  vg <- if (".cond_var_theta" %in% names(post)) as.numeric(post$.cond_var_theta) else 0
  cag <- if (".cond_cov_alpha_theta" %in% names(post)) {
    as.numeric(post$.cond_cov_alpha_theta)
  } else {
    0
  }
  var <- pmax(va + coef_g^2 * vg + 2 * coef_g * cag, 0)
  if (isTRUE(predictive) && abs(as.numeric(re_scale)[[1L]]) > 0) {
    tau_m <- as.numeric(post$sd_trial__Intercept)
    var <- var + (as.numeric(re_scale)[[1L]]^2) * pmax(tau_m, 0)^2
  }
  sd <- sqrt(var)
  qs <- interactmeta_gm_quantile(
    interactmeta_gm(mean, sd, w),
    probs = c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2)
  )
  list(
    median = as.numeric(qs[[2L]]),
    lower = as.numeric(qs[[1L]]),
    upper = as.numeric(qs[[3L]])
  )
}
