# modifiers(), interactions(), combined(), heterogeneity(), shrinkage().

test_that("effect accessors return expected classes and rows", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()

  sg <- modifiers(fit)
  expect_s3_class(sg, "interactmeta_subgroup_effects")
  expect_equal(nrow(sg$summary), 4L)
  expect_true(all(c("Younger", "Older") %in% sg$summary$level))

  int <- interactions(fit)
  expect_s3_class(int, "interactmeta_interaction_effects")
  expect_true(all(is.finite(int$summary$estimate)))

  cb <- combined(fit, type = "Overall", summary = TRUE)
  expect_s3_class(cb, "interactmeta_combined_effect")
  expect_false(inherits(cb, "interactmeta_overall_effect"))
})

test_that("interaction Overall approximates non-reference minus reference at overall IF", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  int <- interactions(fit, type = "Overall", summary = TRUE)
  sg <- modifiers(fit, type = "Overall", summary = TRUE)
  young <- sg$summary$estimate[
    sg$summary$level == "Younger" & sg$summary$type == "Overall"
  ]
  old <- sg$summary$estimate[
    sg$summary$level == "Older" & sg$summary$type == "Overall"
  ]
  int_est <- int$summary$estimate[int$summary$type == "Overall"]
  expect_equal(int_est, old - young, tolerance = 0.05)
})

test_that("heterogeneity() and shrinkage() work on a cached CAMS fit", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  het <- heterogeneity(fit)
  expect_s3_class(het, "interactmeta_heterogeneity")
  expect_gt(nrow(het$summary), 0L)

  sh <- shrinkage(fit, type = "modifiers", summary = TRUE)
  expect_gt(nrow(sh$summary), 0L)
  expect_equal(unique(sh$support$type), "Shrinkage")
})

test_that("frequentist CAMS exposes Wald-style intervals in print", {
  skip_if_not_installed("metafor")
  fit <- fixture_cam_tiny_fit(method = "REML")
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "REML", fixed = TRUE)
  sg <- modifiers(fit, summary = TRUE)
  expect_true(all(is.finite(sg$summary$estimate)))
})