cams_engine_example_data <- function() {
  data.frame(
    trial = rep(paste0("T", 1:4), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 4), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.28, 0.08, 0.22, 0.15, 0.30),
    vi = c(0.04, 0.02, 0.05, 0.03, 0.03, 0.04, 0.02, 0.025)
  )
}

test_that("method = bayes without engine selects bayesmeta", {
  expect_equal(
    interactmeta:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      "bayes"
    ),
    "bayesmeta"
  )
  expect_equal(
    interactmeta:::interactmeta_engine_from_approach(
      c("metafor", "bayesmeta"),
      "bayes"
    ),
    "bayesmeta"
  )
  expect_equal(
    interactmeta:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      NULL,
      approach = "bayes"
    ),
    "bayesmeta"
  )
  expect_equal(
    interactmeta:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      "REML"
    ),
    "metafor"
  )
  expect_equal(
    interactmeta:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      "ML"
    ),
    "metafor"
  )
  expect_equal(
    interactmeta:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      NULL,
      approach = "REML"
    ),
    "metafor"
  )
  expect_equal(
    interactmeta:::interactmeta_engine_from_approach(
      c("bayesmeta", "metafor"),
      NULL,
      approach = "frequentist"
    ),
    "metafor"
  )
  expect_error(
    interactmeta:::interactmeta_engine_from_approach(
      "bayesmeta",
      "bayes",
      approach = "REML"
    ),
    "metafor"
  )
  expect_error(
    interactmeta:::interactmeta_engine_from_approach(
      "metafor",
      "REML",
      approach = "bayes"
    ),
    "bayesmeta"
  )
})

test_that("public CAMS engines are bayesmeta and metafor", {
  expect_equal(interactmeta:::interactmeta_engines(), c("bayesmeta", "metafor"))
  cap <- interactmeta_capability_matrix()
  expect_true(all(cap$engine[cap$analysis == "cam"] %in% c("bayesmeta", "metafor")))
})

test_that("cam() defaults to bayesmeta DIRECT", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    silent = 2
  ))
  expect_equal(fit$engine, "bayesmeta")
  expect_equal(fit$method, "bayes")
  expect_equal(fit$algorithm, "quadrature")
})

test_that("engine / method / algorithm metadata are stored on the fit", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))
  expect_equal(fit$engine, "bayesmeta")
  expect_equal(fit$method, "bayes")
  expect_equal(fit$algorithm, "quadrature")
  expect_false("stanfit" %in% names(fit))
  expect_false("stan_data" %in% names(fit))
  expect_equal(attr(fit, "interactmeta")$engine, "bayesmeta")
  expect_equal(attr(fit, "interactmeta")$method, "bayes")
  expect_equal(attr(fit, "interactmeta")$algorithm, "quadrature")
  expect_s3_class(fit$priors, "cam_priors")
  expect_equal(fit$priors$tau$family, "uniform")
  expect_true(interactmeta:::cam_prior_is_improper_uniform(fit$priors$tau))
  expect_true(interactmeta:::cam_prior_is_improper_uniform(fit$priors$tau_interaction))
  expect_match(paste(utils::capture.output(print(fit)), collapse = "\n"), "improper uniform")
})

test_that("shared cam_prior objects are stored identically for bayesmeta", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  p <- interactmeta:::cam_prior("half-normal", scale = 0.5)
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
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
  expect_equal(interactmeta:::interactmeta_default_method("metafor"), "REML")
  expect_equal(interactmeta:::interactmeta_default_method("bayesmeta"), "bayes")
  expect_equal(
    interactmeta:::interactmeta_resolve_method("metafor", NULL),
    "REML"
  )
  expect_equal(
    interactmeta:::interactmeta_resolve_method("bayesmeta", NULL),
    "bayes"
  )
})

test_that("engine = metafor is frequentist and warns on explicit Bayesian priors", {
  skip_if_not_installed("metafor")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "metafor",
    silent = 2
  ))
  expect_equal(fit$engine, "metafor")
  expect_equal(fit$method, "REML")
  expect_equal(fit$algorithm, "optimization")
  expect_false("stanfit" %in% names(fit))
  expect_false("stan_data" %in% names(fit))
  expect_s3_class(fit$priors, "cam_priors")
  expect_true(!is.null(fit$metafor$overall))
  expect_true(!is.null(fit$metafor$interaction))

  expect_warning(
    cam(
      dat,
      reference = "Younger",
      engine = "metafor",
      tau.prior = interactmeta:::cam_prior("half-normal", scale = 0.5),
      silent = 2
    ),
    "frequentist"
  )
})

test_that("metafor CAMS fits a two-trial overall block", {
  skip_if_not_installed("metafor")
  fit <- cam(respire_age, reference = "<65", engine = "metafor", silent = 2)
  expect_equal(fit$engine, "metafor")
  expect_true(is.finite(fit$draws$b_Intercept[[1L]]))
  ifrac_cols <- grep("^b_ifrac_", names(fit$draws), value = TRUE)
  expect_true(length(ifrac_cols) >= 1L)
  expect_true(all(vapply(fit$draws[ifrac_cols], function(x) is.finite(x[[1L]]), logical(1))))
  expect_equal(as.numeric(fit$metafor$overall$tau2), 0)
  sg <- subgroups(fit, type = "Pooled")
  expect_true(all(is.finite(sg$summary$estimate)))
})

test_that("cmdstanr and method = mle are rejected", {
  dat <- cams_engine_example_data()
  expect_error(
    cam(dat, reference = "Younger", engine = "cmdstanr"),
    "removed"
  )
  expect_error(
    cam(dat, reference = "Younger", method = "mle"),
    "mle"
  )
})

test_that("method and engine are rejected when incompatible", {
  dat <- cams_engine_example_data()
  expect_error(
    cam(dat, reference = "Younger", engine = "bayesmeta", method = "REML"),
    "method = \"bayes\""
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "brms", method = "ML"),
    "method = \"bayes\""
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "metafor", method = "bayes"),
    "REML"
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "brms", tau.prior = "uniform"),
    "improper density"
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "brms", tau.prior = "sqrt"),
    "cannot be represented equivalently"
  )
})

test_that("engine_args are not confused with model priors", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  expect_error(
    cam(
      dat,
      reference = "Younger",
      engine = "bayesmeta",
      tau.prior.overall = "sqrt"
    ),
    "tau.prior"
  )
})
