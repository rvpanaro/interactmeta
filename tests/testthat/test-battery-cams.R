# Core cams() API: backends, method aliases, removed arguments, modifiers.

test_that("default cams() is Bayesian with bayesmeta backend", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  fit <- suppressWarnings(cams(dat, reference = "Younger", silent = 2))
  fit_b <- suppressWarnings(cams(dat, reference = "Younger", method = "Bayesian", silent = 2))
  fit_bayes <- suppressWarnings(cams(dat, reference = "Younger", method = "bayes", silent = 2))

  expect_s3_class(fit, "interactmeta")
  expect_equal(attr(fit, "interactmeta")$backend, "bayesmeta")
  expect_equal(fit$method, "Bayesian")
  expect_equal(fit$algorithm, "quadrature")
  expect_match(
    paste(capture.output(print(fit)), collapse = "\n"),
    "Bayesian model results"
  )
  expect_equal(
    as.numeric(fit_b$support$.weight),
    as.numeric(fit$support$.weight),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(fit_bayes$support$.weight),
    as.numeric(fit$support$.weight),
    tolerance = 1e-10
  )
})

test_that("REML, ML, and frequentist alias select metafor", {
  skip_if_not_installed("metafor")
  dat <- cams_engine_example_data()
  by_reml <- suppressWarnings(cams(dat, reference = "Younger", method = "REML", silent = 2))
  by_freq <- suppressWarnings(cams(dat, reference = "Younger", method = "frequentist", silent = 2))
  by_ml <- suppressWarnings(cams(dat, reference = "Younger", method = "ML", silent = 2))

  expect_equal(attr(by_reml, "interactmeta")$backend, "metafor")
  expect_equal(by_reml$method, "REML")
  expect_equal(by_freq$method, "REML")
  expect_equal(by_ml$method, "ML")
  expect_true(!is.null(by_reml$metafor$overall))
  expect_true(!is.null(by_reml$metafor$interaction))
})

test_that("removed engine, approach, and engine_args are rejected", {
  dat <- cams_engine_example_data()
  expect_error(
    cams(dat, reference = "Younger", approach = "frequentist", silent = 2),
    "approach"
  )
  expect_error(
    cams(dat, reference = "Younger", engine = "metafor", silent = 2),
    "engine"
  )
  expect_error(
    cams(dat, reference = "Younger", engine_args = list(), silent = 2),
    "engine_args"
  )
})

test_that("tau.prior.overall is not accepted as a model prior name", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_engine_example_data()
  expect_error(
    cams(dat, reference = "Younger", tau.prior.overall = "sqrt", silent = 2),
    "tau.prior"
  )
})

test_that("custom modifier column appears in print header", {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    dose = factor(rep(c("Low", "High"), 3), levels = c("Low", "High")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cams(dat, reference = "Low", modifier = "dose", method = "REML", silent = 2)
  expect_equal(attr(fit, "interactmeta")$modifier, "dose")
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Dose levels:", fixed = TRUE)
  expect_false(grepl("subgroup effects", out, fixed = TRUE))
})

test_that("default fit defers reporting until accessors or print need it", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(es, reference = "Younger", silent = 2))
  expect_false(interactmeta:::interactmeta_reporting_is_ready(fit))
  expect_true(isTRUE(fit$latent$.lazy_trials))
  invisible(print(fit, digits = 2))
  expect_true(interactmeta:::interactmeta_reporting_is_ready(fit))
})

test_that("precompute=TRUE caches identical modifier() results", {
  fit <- fixture_age_cam_fit()
  expect_true(interactmeta:::interactmeta_reporting_is_ready(fit))
  sg1 <- modifiers(fit, type = "Overall")
  sg2 <- modifiers(fit, type = "Overall")
  expect_identical(sg1, sg2)
})