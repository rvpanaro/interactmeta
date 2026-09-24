test_that("engine = bayesmeta fits orthogonal overall + interaction via support grid", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(bayescam(
    es,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))

  expect_s3_class(fit, "interactmeta")
  expect_equal(fit$engine, "bayesmeta")
  expect_equal(attr(fit, "interactmeta")$engine, "bayesmeta")
  expect_true(!is.null(fit$bayesmeta$overall))
  expect_true(!is.null(fit$bayesmeta$interaction))
  expect_true(".weight" %in% names(fit$draws))
  expect_true(abs(sum(fit$draws$.weight) - 1) < 1e-8)
  expect_true(any(grepl("^b_ifrac_", names(fit$draws))))
  expect_true(any(grepl("^b_subgroup", names(fit$draws))))
  expect_true(!is.null(fit$latent$beta))
  expect_true(!is.null(fit$latent$m))

  ov_n <- length(fit$bayesmeta$overall$support$weight)
  int_n <- nrow(as.data.frame(fit$bayesmeta$interaction$support))
  expect_equal(nrow(fit$draws), ov_n * int_n)

  out <- paste(utils::capture.output(print(fit, digits = 2)), collapse = "\n")
  expect_match(out, "overall | interaction partition", fixed = TRUE)
  expect_match(out, "Ecological_", fixed = TRUE)
  expect_match(out, "Interaction_", fixed = TRUE)

  s <- summary(fit, digits = 2)
  expect_s3_class(s, "summary.interactmeta")
  expect_true("tau" %in% colnames(s$parameters))
  expect_true(any(grepl("^Ecological_", colnames(s$parameters))))
  expect_true(any(grepl("^Interaction_", colnames(s$parameters))))
})

test_that("cam() routes to bayesmeta engine", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(bayescam(
    es,
    reference = "Younger",
    silent = 2
  ))
  expect_s3_class(fit, "interactmeta")
  expect_equal(fit$engine, "bayesmeta")
})
