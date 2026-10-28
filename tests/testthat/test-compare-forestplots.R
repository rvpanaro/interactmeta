test_that("compare_forestplots() combines two forestplot() results", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("patchwork")

  da_fit <- da(age, reference = "Younger")
  fp_mod <- forestplot(da_fit, include_shrinkage = FALSE)
  fp_int <- forestplot(da_fit, type = "interaction", include_shrinkage = FALSE)

  p <- compare_forestplots(fp_mod, fp_int)
  expect_s3_class(p, "patchwork")

  p_lab <- compare_forestplots(
    fp_mod,
    fp_int,
    labels = c("Modifiers", "Interaction")
  )
  expect_s3_class(p_lab, "patchwork")

  expect_error(
    compare_forestplots(fp_mod, fp_int, labels = c("Only one")),
    "one caption per forest"
  )

  expect_error(
    compare_forestplots(fp_mod),
    "at least two",
    ignore.case = TRUE
  )
  expect_error(
    compare_forestplots(fp_mod, list()),
    "forestplot()"
  )
})

test_that("compare_forestplots(draw = TRUE) prints without error", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("patchwork")

  da_fit <- da(age, reference = "Younger")
  fp_mod <- forestplot(da_fit, include_shrinkage = FALSE)
  fp_int <- forestplot(da_fit, type = "interaction", include_shrinkage = FALSE)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(compare_forestplots(fp_mod, fp_int, draw = TRUE))
})
