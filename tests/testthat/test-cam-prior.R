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

  p <- interactmeta:::cam_prior("half-normal", scale = 0.5)
  expect_equal(format(p), "Half-Normal(0, 0.5)")
  expect_identical(interactmeta:::normalize_cam_prior(p, role = "tau"), p)

  expect_error(
    interactmeta:::normalize_cam_prior("half-normal", role = "tau"),
    "scale"
  )
  expect_error(interactmeta:::cam_prior("half-normal", scale = -1), "positive")
})

test_that("interaction prior inherits only when the user supplies tau.prior", {
  dflt <- interactmeta:::cam_priors_resolve()
  expect_equal(dflt$tau$family, "uniform")
  expect_true(interactmeta:::cam_prior_is_improper_uniform(dflt$tau))
  expect_equal(dflt$tau_interaction$family, "uniform")
  expect_true(interactmeta:::cam_prior_is_improper_uniform(dflt$tau_interaction))
  expect_equal(dflt$fixed$family, "normal")
  expect_equal(dflt$fixed$mean, 0)
  expect_equal(dflt$fixed$sd, 10)

  shared <- interactmeta:::cam_priors_resolve(
    tau.prior = interactmeta:::cam_prior("half-normal", scale = 0.5),
    user_tau = TRUE
  )
  expect_equal(shared$tau$scale, 0.5)
  expect_equal(shared$tau_interaction$scale, 0.5)

  split <- interactmeta:::cam_priors_resolve(
    tau.prior = interactmeta:::cam_prior("half-normal", scale = 0.5),
    tau.interaction.prior = interactmeta:::cam_prior("half-normal", scale = 0.1),
    user_tau = TRUE
  )
  expect_equal(split$tau$scale, 0.5)
  expect_equal(split$tau_interaction$scale, 0.1)
})

test_that("bayesmeta translations match the intended densities", {
  skip_if_not_installed("bayesmeta")

  expect_identical(interactmeta:::prior_to_bayesmeta("sqrt"), "sqrt")
  expect_identical(interactmeta:::prior_to_bayesmeta("uniform"), "uniform")

  p_hn <- interactmeta:::cam_prior("half-normal", scale = 0.5)
  f_hn <- interactmeta:::prior_to_bayesmeta(p_hn)
  x <- c(0, 0.25, 0.5, 1, 2)
  expect_equal(f_hn(x), bayesmeta::dhalfnormal(x, scale = 0.5), tolerance = 1e-12)
  expect_equal(f_hn(x), 2 * stats::dnorm(x, 0, 0.5), tolerance = 1e-12)

  p_u <- interactmeta:::cam_prior("uniform", lower = 0, upper = 2)
  f_u <- interactmeta:::prior_to_bayesmeta(p_u)
  expect_equal(f_u(c(-0.1, 0, 1, 2, 2.1)), c(0, 0.5, 0.5, 0.5, 0), tolerance = 1e-12)

  p_hc <- interactmeta:::cam_prior("half-cauchy", scale = 1)
  f_hc <- interactmeta:::prior_to_bayesmeta(p_hc)
  expect_equal(f_hc(x), bayesmeta::dhalfcauchy(x, scale = 1), tolerance = 1e-12)
})
