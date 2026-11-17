# forestplot(), plot(), and themes (slow tier for heavy renders).

test_that("forestplot() returns a forestplot object for CAMS REML", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  fit <- suppressWarnings(cams(respire_age, reference = "<65", method = "REML", silent = 2))
  fp <- forestplot(fit, include_shrinkage = FALSE)
  expect_s3_class(fp, "interactmeta_forestplot")
})

test_that("plot(type = contribution) returns ggplot when ggplot2 is available", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("ggplot2")
  fit <- fixture_cam_tiny_fit()
  p <- plot(fit, type = "contribution", ifrac_var = "IF_Older")
  expect_s3_class(p, "ggplot")
})

test_that("interactmeta_palette() returns named colours", {
  pal <- interactmeta_palette()
  expect_true(is.character(pal))
  expect_false(is.null(names(pal)))
  expect_gt(length(pal), 4L)
})

test_that("full forestplot render with shrinkage matches shrinkage()", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  fit <- fixture_cam_tiny_fit(precompute = TRUE, lazy_latent = FALSE)
  ref <- shrinkage(fit, type = "modifiers", summary = TRUE)$summary
  fp <- forestplot(fit, type = "modifiers", include_shrinkage = TRUE)
  ra <- attr(fp, "panels")[[1L]]
  sh_rows <- ra[ra$row_kind == "shrinkage" & is.finite(ra$yi), , drop = FALSE]
  expect_gt(nrow(sh_rows), 0L)
  expect_equal(nrow(sh_rows), nrow(ref))
})

test_that("forestplot_draw() renders Bayesian CAMS forests", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  fit <- fixture_cam_tiny_fit(precompute = TRUE, lazy_latent = FALSE)
  fp <- forestplot(fit, type = "modifiers", include_shrinkage = TRUE)
  with_null_pdf(expect_invisible(forestplot_draw(fp)))
  fp_int <- forestplot(fit, type = "interaction", include_shrinkage = FALSE)
  with_null_pdf(expect_invisible(forestplot_draw(fp_int)))
})

test_that("forestplot_draw() renders REML CAMS and construction forests", {
  slow()
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  cam <- suppressWarnings(cams(respire_age, reference = "<65", method = "REML", silent = 2))
  fp_cam <- forestplot(cam, include_shrinkage = FALSE)
  with_null_pdf(expect_invisible(forestplot_draw(fp_cam)))
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  fp_da <- forestplot(da_fit, include_shrinkage = FALSE)
  with_null_pdf(expect_invisible(forestplot_draw(fp_da)))
})

test_that("traceplot() draws Bayesian CAMS diagnostics", {
  slow()
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  with_null_pdf(expect_invisible(traceplot(fit, type = "interaction")))
  with_null_pdf(expect_invisible(traceplot(fit, type = "combined", prior = FALSE)))
})

test_that("compare_forestplots() lays out two forests", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("patchwork")
  fit <- fixture_cam_tiny_fit(precompute = TRUE, lazy_latent = FALSE)
  fp_mod <- forestplot(fit, type = "modifiers", include_shrinkage = FALSE)
  fp_int <- forestplot(fit, type = "interaction", include_shrinkage = FALSE)
  with_null_pdf(expect_invisible(compare_forestplots(fp_mod, fp_int, draw = TRUE)))
})

test_that("plot(type = heterogeneity) draws tau marginals", {
  slow()
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  with_null_pdf(expect_invisible(plot(fit, type = "heterogeneity", prior = FALSE)))
})