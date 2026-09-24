test_that("character shorthands normalize to explicit CAMS priors", {
  u <- interactmeta:::normalize_cam_prior("uniform", role = "tau")
  expect_s3_class(u, "cam_prior")
  expect_equal(u$family, "uniform")
  expect_equal(u$lower, 0)
  expect_true(is.infinite(u$upper))
  expect_true(interactmeta:::cam_prior_is_improper_uniform(u))

  s <- interactmeta:::normalize_cam_prior("sqrt", role = "tau")
  expect_equal(s$family, "sqrt")
  expect_true(interactmeta:::cam_prior_is_improper(s))

  p <- cam_prior("half-normal", scale = 0.5)
  expect_equal(format(p), "Half-Normal(0, 0.5)")
  expect_identical(interactmeta:::normalize_cam_prior(p, role = "tau"), p)

  expect_error(
    interactmeta:::normalize_cam_prior("half-normal", role = "tau"),
    "scale"
  )
  expect_error(cam_prior("half-normal", scale = -1), "positive")
})

test_that("interaction prior inherits only when the user supplies tau.prior", {
  dflt <- interactmeta:::cam_priors_resolve()
  expect_equal(dflt$tau$family, "half-normal")
  expect_equal(dflt$tau$scale, 0.5)
  expect_equal(dflt$tau_interaction$family, "half-normal")
  expect_equal(dflt$tau_interaction$scale, 0.1)
  expect_equal(dflt$fixed$family, "normal")
  expect_equal(dflt$fixed$mean, 0)
  expect_equal(dflt$fixed$sd, 10)

  shared <- interactmeta:::cam_priors_resolve(
    tau.prior = cam_prior("half-normal", scale = 0.5),
    user_tau = TRUE
  )
  expect_equal(shared$tau$scale, 0.5)
  expect_equal(shared$tau_interaction$scale, 0.5)

  split <- interactmeta:::cam_priors_resolve(
    tau.prior = cam_prior("half-normal", scale = 0.5),
    tau.interaction.prior = cam_prior("half-normal", scale = 0.1),
    user_tau = TRUE
  )
  expect_equal(split$tau$scale, 0.5)
  expect_equal(split$tau_interaction$scale, 0.1)
})

test_that("bayesmeta translations match the intended densities", {
  skip_if_not_installed("bayesmeta")

  expect_identical(interactmeta:::prior_to_bayesmeta("sqrt"), "sqrt")
  expect_identical(interactmeta:::prior_to_bayesmeta("uniform"), "uniform")

  p_hn <- cam_prior("half-normal", scale = 0.5)
  f_hn <- interactmeta:::prior_to_bayesmeta(p_hn)
  x <- c(0, 0.25, 0.5, 1, 2)
  expect_equal(f_hn(x), bayesmeta::dhalfnormal(x, scale = 0.5), tolerance = 1e-12)
  expect_equal(f_hn(x), 2 * stats::dnorm(x, 0, 0.5), tolerance = 1e-12)

  p_u <- cam_prior("uniform", lower = 0, upper = 2)
  f_u <- interactmeta:::prior_to_bayesmeta(p_u)
  expect_equal(f_u(c(-0.1, 0, 1, 2, 2.1)), c(0, 0.5, 0.5, 0.5, 0), tolerance = 1e-12)

  p_hc <- cam_prior("half-cauchy", scale = 1)
  f_hc <- interactmeta:::prior_to_bayesmeta(p_hc)
  expect_equal(f_hc(x), bayesmeta::dhalfcauchy(x, scale = 1), tolerance = 1e-12)
})

test_that("brms translations match the same mathematical priors", {
  priors <- interactmeta:::cam_priors_resolve(
    tau.prior = cam_prior("half-normal", scale = 0.5),
    tau.interaction.prior = cam_prior("half-normal", scale = 0.1),
    fixed.prior = cam_prior("normal", mean = 0, sd = 10)
  )
  rows <- interactmeta:::prior_to_brms_rows(priors, trial_var = "trial", ifrac2_vars = "ifrac2_1")
  expect_true(any(rows$class == "Intercept" & rows$prior == "normal(0, 10)"))
  expect_true(any(rows$class == "b" & rows$prior == "normal(0, 10)"))
  hn_m <- rows[rows$class == "sd" & rows$coef == "Intercept", , drop = FALSE]
  expect_equal(hn_m$prior, "normal(0, 0.5)")
  expect_equal(hn_m$lb, 0)
  hn_g <- rows[rows$class == "sd" & rows$coef == "ifrac2_1", , drop = FALSE]
  expect_equal(hn_g$prior, "normal(0, 0.1)")
  expect_equal(hn_g$lb, 0)

  u_priors <- interactmeta:::cam_priors_resolve(
    tau.prior = cam_prior("uniform", lower = 0, upper = 2),
    tau.interaction.prior = cam_prior("uniform", lower = 0, upper = 2),
    user_tau = TRUE
  )
  u_rows <- interactmeta:::prior_to_brms_rows(u_priors)
  sd_u <- u_rows[u_rows$class == "sd" & u_rows$coef == "Intercept", , drop = FALSE]
  expect_equal(sd_u$prior, "uniform(0, 2)")
  expect_equal(sd_u$lb, 0)
  expect_equal(sd_u$ub, 2)
})

test_that("improper uniform and sqrt are not silently bounded for brms", {
  expect_error(
    interactmeta:::prior_to_brms_one("uniform", role = "tau"),
    "improper density"
  )
  expect_error(
    interactmeta:::prior_to_brms_one("sqrt", role = "tau"),
    "cannot be represented equivalently"
  )
  expect_error(
    interactmeta:::prior_to_brms_one(cam_prior("uniform"), role = "tau"),
    "improper density"
  )
})

test_that("prior_to_brms() builds a brmsprior with the translated densities", {
  skip_if_not_installed("brms")
  priors <- interactmeta:::cam_priors_resolve()
  bp <- interactmeta:::prior_to_brms(priors, trial_var = "trial", ifrac2_vars = "ifrac2_1")
  expect_s3_class(bp, "brmsprior")
  expect_true(any(bp$class == "sd" & bp$coef == "Intercept" & grepl("normal\\(0, 0.5\\)", bp$prior)))
  expect_true(any(bp$class == "sd" & bp$coef == "ifrac2_1" & grepl("normal\\(0, 0.1\\)", bp$prior)))
  expect_true(any(bp$class == "Intercept" & grepl("normal\\(0, 10\\)", bp$prior)))
})
