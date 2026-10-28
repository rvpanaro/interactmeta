cams_latent_example_data <- function() {
  data.frame(
    study = rep(paste0("T", 1:4), each = 2),
    level = factor(rep(c("Younger", "Older"), 4), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.28, 0.08, 0.22, 0.15, 0.30),
    vi = c(0.04, 0.02, 0.05, 0.03, 0.03, 0.04, 0.02, 0.025)
  )
}

test_that("trial-effect blocks match the CAMS IV rotation and are orthogonal", {
  yA <- 0.10
  yB <- 0.30
  vA <- 0.04
  vB <- 0.01
  pi <- (1 / vB) / (1 / vA + 1 / vB)
  dat <- data.frame(
    trial = "T",
    level = c("Younger", "Older"),
    yi = c(yA, yB),
    vi = c(vA, vB),
    ifrac_Older = pi
  )
  blocks <- interactmeta:::comparison_trial_effect_blocks(
    data = dat,
    reference = "Younger",
    nonref = "Older",
    ifrac_var = "ifrac_Older"
  )
  expect_equal(blocks$overall$yi, (1 - pi) * yA + pi * yB, tolerance = 1e-12)
  expect_equal(blocks$overall$vi, 1 / (1 / vA + 1 / vB), tolerance = 1e-12)
  expect_equal(blocks$overall$vi, (1 - pi)^2 * vA + pi^2 * vB, tolerance = 1e-12)
  expect_equal(blocks$interaction$yi, yB - yA, tolerance = 1e-12)
  expect_equal(blocks$interaction$vi, vA + vB, tolerance = 1e-12)
  expect_equal(blocks$overall$ifrac, pi, tolerance = 1e-12)
  expect_equal(-(1 - pi) * vA + pi * vB, 0, tolerance = 1e-12)
})

test_that("inverse CAMS transformation recovers subgroup effects", {
  m <- 1.2
  g <- -0.4
  pi <- 0.35
  mu_A <- interactmeta:::interactmeta_muA_from_mg(m, g, pi)
  mu_B <- interactmeta:::interactmeta_muB_from_mg(m, g, pi)
  expect_equal((1 - pi) * mu_A + pi * mu_B, m, tolerance = 1e-12)
  expect_equal(mu_B - mu_A, g, tolerance = 1e-12)
  expect_equal(mu_A, m - pi * g, tolerance = 1e-12)
  expect_equal(mu_B, m + (1 - pi) * g, tolerance = 1e-12)
})

test_that("normal-normal latent conditional handles tau = 0, small, and large", {
  y <- 1
  sigma2 <- 0.25
  mu0 <- 0.2

  z <- interactmeta:::interactmeta_normal_normal_shrink(y, sigma2, mu0, tau = 0)
  expect_equal(unname(z[, "mean"]), mu0)
  expect_equal(unname(z[, "sd"]), 0)

  tiny <- interactmeta:::interactmeta_normal_normal_shrink(y, sigma2, mu0, tau = 1e-8)
  expect_equal(unname(tiny[, "mean"]), mu0, tolerance = 1e-6)
  expect_lt(unname(tiny[, "sd"]), 1e-7)

  huge <- interactmeta:::interactmeta_normal_normal_shrink(y, sigma2, mu0, tau = 1e6)
  expect_equal(unname(huge[, "mean"]), y, tolerance = 1e-8)
  expect_equal(unname(huge[, "sd"]), sqrt(sigma2), tolerance = 1e-8)

  mid <- interactmeta:::interactmeta_normal_normal_shrink(y, sigma2, mu0, tau = 0.5)
  prec_y <- 1 / sigma2
  prec_0 <- 1 / 0.5^2
  expect_equal(unname(mid[, "mean"]), (prec_y * y + prec_0 * mu0) / (prec_y + prec_0))
  expect_equal(unname(mid[, "sd"]), sqrt(1 / (prec_y + prec_0)))
})

test_that("CDF pairing induces comonotonic dependence unlike the product", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
    reference = "Younger",
        silent = 2
  ))
  ov <- interactmeta:::interactmeta_bmr_support_table(fit$bayesmeta$overall)
  int <- interactmeta:::interactmeta_bayesmeta_support_table(fit$bayesmeta$interaction)
  idx <- interactmeta:::interactmeta_support_pair_indices(ov$weight, int$weight)
  w <- as.numeric(ov$weight)
  w <- w / sum(w)
  corr_pair <- interactmeta:::interactmeta_weighted_cor(ov$ifrac, int$mu[idx], w)

  i <- rep(seq_len(nrow(ov)), each = nrow(int))
  k <- rep(seq_len(nrow(int)), times = nrow(ov))
  wp <- as.numeric(ov$weight[i]) * as.numeric(int$weight[k])
  corr_prod <- interactmeta:::interactmeta_weighted_cor(ov$ifrac[i], int$mu[k], wp)

  expect_lt(abs(corr_prod), 1e-10)
  expect_gt(corr_pair, 0.2)

  beta_pair <- ov$ifrac - int$mu[idx]
  sd_pair <- sqrt(sum(w * (beta_pair - sum(w * beta_pair))^2))
  sd_mix <- sqrt(interactmeta:::interactmeta_gm_var(fit$latent$beta$gm))
  expect_lt(sd_pair, sd_mix)
})

test_that("bayesmeta latent posteriors exist and beta matches convolve", {
  skip_if_not_installed("bayesmeta")

  dat <- cams_latent_example_data()
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
        silent = 2
  ))
  elapsed <- proc.time()[["elapsed"]] - t0

  ov_n <- length(fit$bayesmeta$overall$support$weight)
  int_n <- nrow(as.data.frame(fit$bayesmeta$interaction$support))
  expect_equal(nrow(fit$support), ov_n * int_n)
  expect_true(abs(sum(fit$support$.weight) - 1) < 1e-8)
  expect_equal(fit$latent$independence, "cross_block")
  expect_true(grepl("Marginal posteriors only", fit$latent$note, fixed = TRUE))

  sm_mix <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$beta)
  expect_true(!is.null(fit$latent$beta$convolve))
  bounds <- interactmeta:::interactmeta_interval_from_distribution(
    fit$latent$beta$convolve$quantile,
    fit$latent$beta$convolve$cdf,
    prob = 0.95,
    interval.type = fit$interval.type,
    tol = fit$latent$tol.uniroot
  )
  expect_equal(c(sm_mix$lower, sm_mix$upper), bounds, tolerance = 1e-6)
  expect_equal(
    sm_mix$median,
    as.numeric(fit$latent$beta$convolve$quantile(0.5)),
    tolerance = 1e-6
  )
  # Mixture quantiles remain a check against central convolution quantiles.
  qs <- c(0.025, 0.5, 0.975)
  expect_equal(
    as.numeric(fit$latent$beta$convolve$quantile(qs)),
    interactmeta:::interactmeta_gm_quantile(fit$latent$beta$gm, qs),
    tolerance = 0.05
  )
  expect_equal(
    sm_mix$mean,
    interactmeta:::interactmeta_gm_mean(fit$latent$theta$gm) -
      interactmeta:::interactmeta_gm_mean(fit$latent$gamma$gm),
    tolerance = 1e-10
  )

  trials <- fit$latent$trials
  expect_equal(sort(names(fit$latent$m)), sort(trials))
  for (tr in trials) {
    pij <- fit$latent$pi[[tr]]
    sm_m <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$m[[tr]])
    sm_g <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$g[[tr]])
    sm_A <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_A[[tr]])
    sm_B <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_B[[tr]])
    expect_equal(sm_A$mean, sm_m$mean - pij * sm_g$mean, tolerance = 1e-8)
    expect_equal(sm_B$mean, sm_m$mean + (1 - pij) * sm_g$mean, tolerance = 1e-8)
    if (!is.null(fit$latent$mu_A[[tr]]$convolve)) {
      expect_equal(
        fit$latent$mu_A[[tr]]$convolve$quantile(0.5),
        sm_A$median,
        tolerance = 0.08
      )
    }
  }

  # $mu stays the fixed-effect linpred, not the shrunken trial posterior
  mu_mean <- colSums(fit$support$.weight * fit$mu)
  tr1 <- trials[[1L]]
  sm_A1 <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_A[[tr1]])
  expect_gt(abs(mu_mean[[1L]] - sm_A1$mean), 1e-4)

  expect_lt(ov_n * int_n, 10000)
  expect_lt(elapsed, 120)
})

test_that("modifiers() at fixed pi* uses joint (alpha, theta) plus independent gamma", {
  skip_if_not_installed("bayesmeta")

  dat <- cams_latent_example_data()
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
        silent = 2
  ))
  pi_star <- 0.35
  sg <- modifiers(fit, at = pi_star, type = "Overall")
  post <- as.data.frame(fit$support)
  post <- post[order(post$.draw), , drop = FALSE]
  gamma_col <- grep("^b_subgroup", names(post), value = TRUE)
  expect_length(gamma_col, 1L)
  # Stored b_ifrac is beta = theta - gamma. Recover theta on the product
  # grid, form M* from the joint (alpha, theta) row, then add independent
  # gamma: mu_A = M* - pi* gamma, mu_B = M* + (1 - pi*) gamma.
  beta_col <- grep("^b_ifrac", names(post), value = TRUE)
  beta_col <- beta_col[!grepl("ifrac2", beta_col)]
  expect_length(beta_col, 1L)
  theta <- post[[beta_col]] + post[[gamma_col]]
  M_star <- post$b_Intercept + pi_star * theta
  mu_A <- M_star - pi_star * post[[gamma_col]]
  mu_B <- M_star + (1 - pi_star) * post[[gamma_col]]
  dA <- sg$support[sg$support$level == "Younger" & sg$support$type == "Overall", ]
  dB <- sg$support[sg$support$level == "Older" & sg$support$type == "Overall", ]
  dA <- dA[order(dA$.draw), , drop = FALSE]
  dB <- dB[order(dB$.draw), , drop = FALSE]
  expect_equal(unname(dA$effect), unname(mu_A), tolerance = 1e-10)
  expect_equal(unname(dB$effect), unname(mu_B), tolerance = 1e-10)

  expect_error(
    interactmeta:::interactmeta_gm_linear_sum(
      fit$latent$alpha$gm,
      1,
      fit$latent$theta$gm,
      pi_star
    ),
    "independent mixtures"
  )

  ov <- fit$bayesmeta$overall
  nms <- colnames(ov$support$mean)
  i_a <- match("intercept", nms)
  i_t <- match("ifrac", nms)
  w <- as.numeric(ov$support$weight)
  w <- w / sum(w)
  a <- as.numeric(ov$support$mean[, i_a])
  th <- as.numeric(ov$support$mean[, i_t])
  cov_arr <- ov$support$covariance
  va <- pmax(as.numeric(cov_arr[, i_a, i_a]), 0)
  vt <- pmax(as.numeric(cov_arr[, i_t, i_t]), 0)
  cat <- as.numeric(cov_arr[, i_a, i_t])
  mean_M <- a + pi_star * th
  var_node <- va + (pi_star^2) * vt + 2 * pi_star * cat
  mu <- sum(w * mean_M)
  var_joint <- sum(w * (pmax(var_node, 0) + (mean_M - mu)^2))
  mu_a <- sum(w * a)
  mu_th <- sum(w * th)
  V_a <- sum(w * va) + sum(w * (a - mu_a)^2)
  V_th <- sum(w * vt) + sum(w * (th - mu_th)^2)
  Cov <- sum(w * cat) + sum(w * (a - mu_a) * (th - mu_th))
  expect_equal(var_joint, V_a + (pi_star^2) * V_th + 2 * pi_star * Cov, tolerance = 1e-10)

  gm_a <- interactmeta:::interactmeta_gm(a, sqrt(pmax(va, 0)), w)
  gm_th <- interactmeta:::interactmeta_gm(th, sqrt(pmax(vt, 0)), w)
  var_ind <- interactmeta:::interactmeta_gm_var(
    interactmeta:::interactmeta_gm_linear_sum(gm_a, 1, gm_th, pi_star)
  )
  if (abs(Cov) > 1e-8) {
    expect_gt(abs(var_ind - var_joint), 1e-8)
  }
})

test_that("engine-native latent intervals follow bayesmeta post.interval and convolve", {
  skip_if_not_installed("bayesmeta")

  dat <- cams_latent_example_data()
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
        silent = 2
  ))
  expect_equal(fit$interval.type, "central")
  prob <- 0.95
  it <- fit$interval.type

  # Public extractors and native summaries both use fit interval.type.
  int <- interactions(fit, prob = prob)
  pooled <- int$summary[int$summary$type == "Overall", ]
  native_g <- fit$bayesmeta$interaction$post.interval(
    mu.level = prob,
    method = it
  )
  post <- as.data.frame(fit$support)
  g_col <- grep("^b_subgroup", names(post), value = TRUE)[[1L]]
  qi_g <- interactmeta:::comparison_draw_moments(
    post[[g_col]],
    prob = prob,
    weight = post$.weight,
    cond_sd = post$.cond_sd_gamma,
    interval.type = it
  )
  expect_equal(pooled$lci, qi_g$lci, tolerance = 1e-8)
  expect_equal(pooled$uci, qi_g$uci, tolerance = 1e-8)
  expect_true(is.numeric(native_g) && length(native_g) == 2L)

  pars <- interactmeta:::interactmeta_bayesmeta_param_summary(fit, prob = prob)
  ov <- fit$bayesmeta$overall
  which_a <- match(
    interactmeta:::interactmeta_bmr_coef_name(ov, "intercept"),
    colnames(ov$support$mean)
  )
  native_a <- ov$post.interval(
    beta.level = prob,
    which.beta = which_a,
    method = it
  )
  expect_equal(
    unname(pars[c("95% lower", "95% upper"), "Intercept"]),
    unname(as.numeric(native_a)),
    tolerance = 1e-8
  )

  eco <- grep("^Ecological_", colnames(pars), value = TRUE)
  eco_bounds <- interactmeta:::interactmeta_interval_from_distribution(
    fit$latent$beta$convolve$quantile,
    fit$latent$beta$convolve$cdf,
    prob = prob,
    interval.type = it,
    tol = fit$latent$tol.uniroot
  )
  expect_equal(unname(pars[c("95% lower", "95% upper"), eco]), eco_bounds, tolerance = 1e-6)

  blocks <- fit$bayesmeta$trial_blocks
  expect_equal(unname(as.numeric(ov$sigma)), unname(sqrt(blocks$overall$vi)), tolerance = 1e-12)
  expect_equal(
    unname(as.numeric(fit$bayesmeta$interaction$sigma)),
    unname(sqrt(blocks$interaction$vi)),
    tolerance = 1e-12
  )

  expect_error(
    cams(dat, reference = "Younger", silent = 2,
             interval.type = "evidentiary"),
    "derived CAMS sum"
  )
  bad <- dat
  bad$vi[1] <- 0
  expect_error(
    cams(bad, reference = "Younger", silent = 2),
    "strictly positive"
  )
  expect_error(
    cams(dat, reference = "Younger", silent = 2,
             convolve_derived = FALSE),
    "convolve_derived"
  )

  out <- paste(utils::capture.output(print(fit, digits = 2)), collapse = "\n")
  expect_match(out, "central credible intervals", fixed = TRUE)

  # Conditional-mean quantiles are not the published ecological interval.
  beta_col <- grep("^b_ifrac", names(fit$support), value = TRUE)
  beta_col <- beta_col[!grepl("ifrac2", beta_col)]
  expect_length(beta_col, 1L)
  wq <- interactmeta:::interactmeta_weighted_quantile(
    fit$support[[beta_col]], fit$support$.weight, c(0.025, 0.975)
  )
  expect_gt(max(abs(eco_bounds - wq)), 1e-4)
})
