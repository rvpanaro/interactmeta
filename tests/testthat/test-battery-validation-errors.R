# Validation stops and soft-deprecation warnings (fast, minimal fits).

test_that("subgroup_long stops on invalid wide input", {
  expect_error(subgroup_long(1L), "`data` must be a data frame")
  wide <- data.frame(study = "T1", stringsAsFactors = FALSE)
  expect_error(subgroup_long(wide, subgroups = "A", values = "yi"), "Missing wide column")
})

test_that("compare_* guards require correct fit types and arity", {
  skip_if_not_installed("metafor")
  b <- fixture_respire_constructions()
  expect_error(
    interactmeta:::compare_ifrac(b$da, b$cams),
    "compare_ifrac"
  )
  expect_error(compare_heterogeneity(b$da), "at least two")
  expect_error(
    compare(b$da),
    "at least two"
  )
  expect_error(
    interactmeta:::interactmeta_compare_labels(
      list(b$da, b$cams),
      labels = "only_one"
    ),
    "`labels`"
  )
})

test_that("cams REML warns when Bayesian priors are supplied", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  dat <- fixture_cam_tiny()
  attr(dat, "reference") <- "Younger"
  IF_Older <- c(0.3, 0.3, 0.55, 0.55)
  ifrac2_Older <- c(-0.3, 0.7, -0.55, 0.45)
  dat$IF_Older <- IF_Older
  dat$ifrac2_Older <- ifrac2_Older
  expect_warning(
    suppressMessages(cams(
      dat,
      reference = "Younger",
      method = "REML",
      silent = 2,
      tau.prior = function(tau) bayesmeta::dhalfnormal(tau, scale = 0.5)
    )),
    "frequentist"
  )
})

test_that("ifrac stops on invalid reference", {
  skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  expect_error(ifrac(es, reference = 99L), "out of bounds")
})

test_that("modifiers deprecation paths on Bayesian tiny fit", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  expect_warning(modifiers(fit, contribution = 0.5), "contribution")
  expect_error(modifiers(fit, at = 0.5, contribution = 0.5), "only one")
})

test_that("traceplot rejects unavailable type on Bayesian fit", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  expect_error(traceplot(fit, type = "modifier"), "not available")
})