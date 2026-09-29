test_that("ecological_slope = FALSE fits on public engines with beta = 0", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  dat <- data.frame(
    trial = rep(c("T1", "T2", "T3"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.10, 0.30, 0.12, 0.28, 0.08, 0.24),
    vi = c(0.04, 0.04, 0.05, 0.03, 0.03, 0.05)
  )

  fit_m <- cam(
    dat,
    reference = "Younger",
    engine = "metafor",
    ecological_slope = FALSE
  )
  expect_false(isTRUE(fit_m$ecological_slope))
  beta_m <- fit_m$draws[[grep("^b_ifrac", names(fit_m$draws), value = TRUE)[[1L]]]]
  expect_equal(as.numeric(beta_m), 0)
  expect_s3_class(subgroups(fit_m), "interactmeta_subgroup_effects")
  expect_s3_class(interactions(fit_m), "interactmeta_interaction_effects")
  expect_s3_class(heterogeneity(fit_m), "interactmeta_heterogeneity")

  fit_b <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    ecological_slope = FALSE,
    silent = 2,
    tau.prior = function(tau) bayesmeta::dhalfnormal(tau, scale = 0.5),
    tau.interaction.prior = function(tau) bayesmeta::dhalfnormal(tau, scale = 0.1)
  ))
  expect_false(isTRUE(fit_b$ecological_slope))
  beta_b <- fit_b$draws[[grep("^b_ifrac", names(fit_b$draws), value = TRUE)[[1L]]]]
  expect_true(all(abs(as.numeric(beta_b)) < 1e-12))
  expect_true(isTRUE(all(fit_b$draws$sd_trial__Intercept >= 0)))
  sg <- subgroups(fit_b)
  expect_s3_class(sg, "interactmeta_subgroup_effects")
  expect_true(nrow(sg$summary) >= 2L)
  het <- heterogeneity(fit_b)
  expect_true(any(grepl("Overall", het$summary$parameter)))
  expect_true(any(grepl("Interaction", het$summary$parameter)))
})

test_that("constrained joint moments recover alpha and gamma", {
  m <- c(0.1, 0.2, 0.15)
  g <- c(0.3, 0.25, 0.35)
  pi <- c(0.4, 0.5, 0.6)
  vm <- rep(0.04, 3)
  vg <- rep(0.05, 3)
  mom <- interactmeta:::interactmeta_joint_alpha_gamma_moments(
    m = m, vm = vm, pi = pi, g = g, vg = vg,
    tau_m = 0, tau_g = 0,
    prior_mean = c(0, 0),
    prior_sd = c(Inf, Inf)
  )
  expect_equal(names(mom$mean), c("alpha", "gamma"))
  expect_true(is.finite(mom$mean[["alpha"]]))
  expect_true(is.finite(mom$mean[["gamma"]]))
  expect_equal(dim(mom$cov), c(2L, 2L))
})
