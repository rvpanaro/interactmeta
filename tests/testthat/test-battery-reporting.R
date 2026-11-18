# print(), summary(), and tidy().

test_that("print and summary for Bayesian CAMS use current labels", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_age_cam_fit()

  out <- paste(capture.output(print(fit, digits = 2)), collapse = "\n")
  expect_match(out, "'interactmeta' object", fixed = TRUE)
  expect_match(out, "Bayesian model results", fixed = TRUE)
  expect_match(out, "Ecological_", fixed = TRUE)
  expect_match(out, "method: Bayesian", fixed = TRUE)
  expect_false(grepl("approach:", out, fixed = TRUE))

  s <- summary(fit, digits = 2)
  expect_s3_class(s, "summary.interactmeta")
  expect_true(is.matrix(s$parameters))
  expect_true(any(grepl("^Ecological_", colnames(s$parameters))))

  sout <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(sout, "Summary of 'interactmeta' object", fixed = TRUE)
  expect_match(sout, "heterogeneity:", fixed = TRUE)
})

test_that("tidy.interactmeta returns a tibble with interaction rows", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_age_cam_fit()
  td <- tidy(fit, estimands = "interaction", type = "Overall")
  expect_s3_class(td, "tbl_df")
  expect_gt(nrow(td), 0L)
  expect_true("estimate" %in% names(td))
})

test_that("combined() is separate from removed overall() alias", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  expect_false(exists("overall", envir = asNamespace("interactmeta"), inherits = FALSE))
  cb <- combined(fit, type = "Overall", summary = TRUE)
  expect_s3_class(cb, "interactmeta_combined_effect")
})