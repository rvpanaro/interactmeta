test_that("method = bayes without engine selects bayesmeta", {
  expect_equal(
    interact:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      "bayes"
    ),
    "bayesmeta"
  )
  expect_equal(
    interact:::interactmeta_engine_from_approach(
      c("metafor", "bayesmeta"),
      "bayes"
    ),
    "bayesmeta"
  )
  expect_equal(
    interact:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      NULL,
      method = "Bayesian"
    ),
    "bayesmeta"
  )
  expect_equal(
    interact:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      "REML"
    ),
    "metafor"
  )
  expect_equal(
    interact:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      "ML"
    ),
    "metafor"
  )
  expect_equal(
    interact:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      NULL,
      method = "REML"
    ),
    "metafor"
  )
  expect_equal(
    interact:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      NULL,
      approach = "frequentist"
    ),
    "metafor"
  )
  expect_equal(
    interact:::interactmeta_resolve_cams_method("REML")$backend,
    "metafor"
  )
  expect_equal(
    interact:::interactmeta_resolve_cams_method("Bayesian")$backend,
    "bayesmeta"
  )
})

test_that("public CAMS engines are bayesmeta and metafor", {
  expect_equal(interact:::interactmeta_engines(), c("bayesmeta", "metafor"))
  cap <- interactmeta_capability_matrix()
  expect_true(all(cap$engine[cap$analysis == "cam"] %in% c("bayesmeta", "metafor")))
})

test_that("cams() defaults to bayesmeta DIRECT", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    silent = 2
  ))
  expect_equal(interactmeta_fit_backend(fit), "bayesmeta")
  expect_equal(fit$method, "Bayesian")
  expect_equal(fit$algorithm, "quadrature")
})

test_that("engine / method / algorithm metadata are stored on the fit", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
        silent = 2
  ))
  expect_equal(interactmeta_fit_backend(fit), "bayesmeta")
  expect_equal(fit$method, "Bayesian")
  expect_equal(fit$algorithm, "quadrature")
  expect_false("stanfit" %in% names(fit))
  expect_false("stan_data" %in% names(fit))
  expect_equal(attr(fit, "interactmeta")$backend, "bayesmeta")
  expect_equal(attr(fit, "interactmeta")$method, "Bayesian")
  expect_equal(attr(fit, "interactmeta")$algorithm, "quadrature")
  expect_s3_class(fit$priors, "interactmeta_priors")
  expect_equal(fit$priors$tau$family, "uniform")
  expect_true(interact:::cam_prior_is_improper_uniform(fit$priors$tau))
  expect_true(interact:::cam_prior_is_improper_uniform(fit$priors$tau_interaction))
  expect_match(paste(utils::capture.output(print(fit)), collapse = "\n"), "Bayesian estimates")
})

test_that("shared cam_prior objects are stored identically for bayesmeta", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  p <- interact:::cam_prior("half-normal", scale = 0.5)
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
        tau.prior = p,
    tau.interaction.prior = p,
    silent = 2
  ))
  expect_equal(fit$priors$tau$scale, 0.5)
  expect_equal(fit$priors$tau_interaction$scale, 0.5)
  tau_fun <- fit$bayesmeta$overall$tau.prior
  expect_true(is.function(tau_fun) || identical(tau_fun, "sqrt") || is.character(tau_fun))
  if (is.function(tau_fun)) {
    expect_equal(tau_fun(0.5), bayesmeta::dhalfnormal(0.5, scale = 0.5), tolerance = 1e-10)
  }
})

test_that("omitted method is REML for metafor and bayes for bayesmeta", {
  expect_equal(interact:::interactmeta_default_method("metafor"), "REML")
  expect_equal(interact:::interactmeta_default_method("bayesmeta"), "bayes")
  expect_equal(
    interact:::interactmeta_resolve_method("metafor", NULL),
    "REML"
  )
  expect_equal(
    interact:::interactmeta_resolve_method("bayesmeta", NULL),
    "bayes"
  )
})

test_that("engine = metafor is frequentist and warns on explicit Bayesian priors", {
  skip_if_not_installed("metafor")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    method = "REML",
    silent = 2
  ))
  expect_equal(interactmeta_fit_backend(fit), "metafor")
  expect_equal(fit$method, "REML")
  expect_equal(fit$algorithm, "optimization")
  expect_false("stanfit" %in% names(fit))
  expect_false("stan_data" %in% names(fit))
  expect_s3_class(fit$priors, "interactmeta_priors")
  expect_true(!is.null(fit$metafor$overall))
  expect_true(!is.null(fit$metafor$interaction))

  expect_warning(
    cams(
      dat,
      reference = "Younger",
      method = "REML",
      tau.prior = interact:::cam_prior("half-normal", scale = 0.5),
      silent = 2
    ),
    "frequentist"
  )
})

test_that("metafor CAMS fits a two-trial overall block", {
  skip_if_not_installed("metafor")
  fit <- cams(respire_age, reference = "<65", method = "REML", silent = 2)
  expect_equal(interactmeta_fit_backend(fit), "metafor")
  expect_true(is.finite(fit$support$b_Intercept[[1L]]))
  ifrac_cols <- grep("^b_ifrac_", names(fit$support), value = TRUE)
  expect_true(length(ifrac_cols) >= 1L)
  expect_true(all(vapply(fit$support[ifrac_cols], function(x) is.finite(x[[1L]]), logical(1))))
  expect_equal(as.numeric(fit$metafor$overall$tau2), 0)
  sg <- modifiers(fit, type = "Overall")
  expect_true(all(is.finite(sg$summary$estimate)))
})

test_that("lowercase method = mle is rejected", {
  dat <- cams_engine_example_data()
  expect_error(
    cams(dat, reference = "Younger", method = "mle"),
    "mle"
  )
})

test_that("method = REML selects metafor", {
  skip_if_not_installed("metafor")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cams(dat, reference = "Younger", method = "REML", silent = 2))
  expect_equal(interactmeta_fit_backend(fit), "metafor")
  expect_equal(fit$method, "REML")
})

test_that("engine_args are not confused with model priors", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  expect_error(
    cams(
      dat,
      reference = "Younger",
            tau.prior.overall = "sqrt"
    ),
    "tau.prior"
  )
})

test_that("cams() uses bayesmeta by default with usable accessors", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cams(dat, reference = "Younger", silent = 2))
  expect_s3_class(fit, "interactmeta")
  expect_equal(attr(fit, "interactmeta")$backend, "bayesmeta")
  post <- as.data.frame(fit)
  expect_true("b_Intercept" %in% names(post) || "b_ifrac_Older" %in% names(post))
})

test_that("bayesmeta DIRECT stores orthogonal overall and interaction support", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  fit <- fixture_age_cam_fit()
  expect_s3_class(fit, "interactmeta")
  expect_equal(interactmeta_fit_backend(fit), "bayesmeta")
  expect_true(!is.null(fit$bayesmeta$overall))
  expect_true(!is.null(fit$bayesmeta$interaction))
  expect_true(".weight" %in% names(fit$support))
  expect_true(abs(sum(fit$support$.weight) - 1) < 1e-8)
  ov_n <- length(fit$bayesmeta$overall$support$weight)
  int_n <- nrow(as.data.frame(fit$bayesmeta$interaction$support))
  expect_equal(nrow(fit$support), ov_n * int_n)
  expect_true(!is.null(fit$latent$beta))
  expect_true(!is.null(fit$latent$m))
})
