# forestplot(), plot(), and themes (fast tier; slow draw coverage in test-coverage-sweep.R).

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

test_that("forestplot shrinkage rows match shrinkage() summary", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  fit <- fixture_cam_tiny_fit()
  ref <- shrinkage(fit, type = "modifiers", summary = TRUE)$summary
  fp <- forestplot(fit, type = "modifiers", include_shrinkage = TRUE)
  ra <- attr(fp, "panels")[[1L]]
  sh_rows <- ra[ra$row_kind == "shrinkage" & is.finite(ra$yi), , drop = FALSE]
  expect_gt(nrow(sh_rows), 0L)
  expect_equal(nrow(sh_rows), nrow(ref))
})