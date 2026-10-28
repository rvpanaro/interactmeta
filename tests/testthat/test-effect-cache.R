test_that("cams() precomputes reporting parameters and effect cache", {
  fit <- fixture_age_cam_fit()

  expect_true(interactmeta:::interactmeta_reporting_is_ready(fit))

  sg1 <- modifiers(fit, type = "Overall")
  sg2 <- modifiers(fit, type = "Overall")
  expect_identical(sg2, sg1)

  out <- paste(capture.output(print(fit, digits = 2)), collapse = "\n")
  expect_match(out, "Ecological_Older", fixed = TRUE)
})
