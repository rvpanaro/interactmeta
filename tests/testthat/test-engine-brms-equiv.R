test_that("bayesmeta and brms posteriors agree under identical CAMS priors", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("brms")

  dat <- data.frame(
    trial = rep(paste0("T", 1:4), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 4), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.28, 0.08, 0.22, 0.15, 0.30),
    vi = c(0.04, 0.02, 0.05, 0.03, 0.03, 0.04, 0.02, 0.025)
  )
  p_tau <- cam_prior("half-normal", scale = 0.5)
  p_int <- cam_prior("half-normal", scale = 0.1)
  p_fe <- cam_prior("normal", mean = 0, sd = 10)

  fit_bm <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    tau.prior = p_tau,
    tau.interaction.prior = p_int,
    fixed.prior = p_fe,
    silent = 2
  ))
  fit_br <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "brms",
    tau.prior = p_tau,
    tau.interaction.prior = p_int,
    fixed.prior = p_fe,
    silent = 2,
    engine_args = list(
      chains = 2,
      iter = 2000,
      warmup = 1000,
      cores = 1,
      refresh = 0,
      open_progress = FALSE,
      control = list(adapt_delta = 0.95)
    )
  ))

  expect_equal(fit_bm$engine, "bayesmeta")
  expect_equal(fit_br$engine, "brms")
  expect_equal(fit_bm$method, "bayes")
  expect_equal(fit_br$method, "bayes")
  expect_equal(fit_bm$algorithm, "quadrature")
  expect_equal(fit_br$algorithm, "MCMC")
  expect_equal(fit_bm$priors$tau$scale, fit_br$priors$tau$scale)
  expect_equal(fit_bm$priors$tau_interaction$scale, fit_br$priors$tau_interaction$scale)
  expect_equal(fit_bm$priors$fixed$sd, fit_br$priors$fixed$sd)

  brms_prior <- fit_br$prior
  expect_true(any(brms_prior$class == "sd" & brms_prior$prior == "normal(0, 0.5)"))
  expect_true(any(brms_prior$class == "sd" & brms_prior$prior == "normal(0, 0.1)"))
  expect_true(any(brms_prior$class == "Intercept" & brms_prior$prior == "normal(0, 10)"))

  gamma_col <- grep("^b_subgroup", names(fit_br$draws), value = TRUE)[[1L]]
  beta_col <- grep("^b_ifrac_", names(fit_br$draws), value = TRUE)[[1L]]
  tau_g_col <- grep("^sd_trial__ifrac2", names(fit_br$draws), value = TRUE)[[1L]]

  compare_latent <- function(bm_sm, br_vec, tol_mean, tol_sd, label) {
    expect_equal(bm_sm$mean, mean(br_vec), tolerance = tol_mean, info = paste(label, "mean"))
    expect_equal(bm_sm$sd, stats::sd(br_vec), tolerance = tol_sd, info = paste(label, "sd"))
  }

  # $draws for bayesmeta are conditional means given tau; $latent includes
  # within-tau coefficient variance. FE joints still differ slightly:
  # brms/Stan use independent N(0,10) on (alpha, beta, gamma); blocks use
  # independent N(0,10) on (alpha, theta) and gamma.
  compare_latent(
    interactmeta:::interactmeta_posterior_1d_summary(fit_bm$latent$alpha),
    fit_br$draws$b_Intercept,
    tol_mean = 0.05,
    tol_sd = 0.08,
    "alpha"
  )
  compare_latent(
    interactmeta:::interactmeta_posterior_1d_summary(fit_bm$latent$gamma),
    fit_br$draws[[gamma_col]],
    tol_mean = 0.05,
    tol_sd = 0.08,
    "gamma"
  )
  compare_latent(
    interactmeta:::interactmeta_posterior_1d_summary(fit_bm$latent$beta),
    fit_br$draws[[beta_col]],
    tol_mean = 0.08,
    tol_sd = 0.12,
    "beta"
  )
  compare_latent(
    interactmeta:::interactmeta_discrete_summary(fit_bm$latent$tau_m),
    fit_br$draws$sd_trial__Intercept,
    tol_mean = 0.05,
    tol_sd = 0.05,
    "tau_m"
  )
  compare_latent(
    interactmeta:::interactmeta_discrete_summary(fit_bm$latent$tau_g),
    fit_br$draws[[tau_g_col]],
    tol_mean = 0.04,
    tol_sd = 0.04,
    "tau_g"
  )

  sm_a <- interactmeta:::interactmeta_posterior_1d_summary(fit_bm$latent$alpha)
  sm_g <- interactmeta:::interactmeta_posterior_1d_summary(fit_bm$latent$gamma)
  sg_nonref_br <- fit_br$draws$b_Intercept + fit_br$draws[[gamma_col]]
  expect_equal(sm_a$mean + sm_g$mean, mean(sg_nonref_br), tolerance = 0.08)
})
