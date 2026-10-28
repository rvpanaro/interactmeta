test_that("cams() defaults skip precompute and use lazy trial latent", {
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(es, reference = "Younger", silent = 2))
  expect_false(interact:::interactmeta_reporting_is_ready(fit))
  expect_true(isTRUE(fit$latent$.lazy_trials))
  expect_null(fit$latent$m[[1L]])
  expect_true(!is.null(fit$latent$beta))
})

test_that("print() lazily builds the parameter matrix", {
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(es, reference = "Younger", silent = 2))
  out <- paste(capture.output(print(fit, digits = 2)), collapse = "\n")
  expect_match(out, "Ecological_Older", fixed = TRUE)
  expect_true(interact:::interactmeta_reporting_is_ready(fit))
})

test_that("shrinkage bundle materializes lazy trial latents", {
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(es, reference = "Younger", silent = 2))
  expect_true(isTRUE(fit$latent$.lazy_trials))
  bundle <- interact:::comparison_shrinkage_rotation_bundle(fit)
  expect_false(isTRUE(bundle$fit$latent$.lazy_trials))
  expect_true(interact:::interactmeta_bayesmeta_trial_latent_ready(bundle$fit))
})

test_that("ensure_bayesmeta_latent() builds trial marginals on demand", {
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(es, reference = "Younger", silent = 2))
  fit <- interact:::interactmeta_ensure_bayesmeta_latent(fit, full_trials = TRUE)
  expect_false(isTRUE(fit$latent$.lazy_trials))
  expect_gt(length(fit$latent$m), 0L)
})

test_that("coarser delta reduces DIRECT support size", {
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit_fine <- suppressWarnings(cams(es, reference = "Younger", silent = 2, delta = 0.01))
  fit_coarse <- suppressWarnings(cams(es, reference = "Younger", silent = 2, delta = 0.02))
  expect_gt(nrow(fit_fine$support), nrow(fit_coarse$support))
})

test_that("tidy() tolerates partial fit$reporting (exponentiate only)", {
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(es, reference = "Younger", silent = 2))
  fit$reporting <- list(exponentiate = TRUE)
  expect_false(interact:::interactmeta_reporting_is_ready(fit))
  out <- suppressWarnings(
    tidy(fit, estimands = "interaction", type = "Overall", exponentiate = TRUE)
  )
  expect_gt(nrow(out), 0L)
})

test_that("fit_args n_support does not error", {
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(
    es,
    reference = "Younger",
    silent = 2,
    fit_args = list(n_support = 40L)
  ))
  expect_s3_class(fit, "interactmeta")
})
