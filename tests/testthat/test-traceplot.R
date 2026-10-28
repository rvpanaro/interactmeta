test_that("traceplot() uses native bayesmeta objects and does not refit", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
        reference = "Younger",
    silent = 2
  ))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(traceplot(fit, type = "subgroups"), "not available")

  imports <- parent.env(asNamespace("bayesmeta"))
  mtext_before <- get("mtext", envir = imports, inherits = FALSE)

  called <- 0L
  testthat::with_mocked_bindings(
    {
      objs_int <- traceplot(fit, type = "interaction")
      expect_identical(objs_int[[1L]], fit$bayesmeta$interaction)
    },
    bayesmeta = function(...) {
      called <<- called + 1L
      stop("auxiliary bayesmeta fit is forbidden")
    },
    .package = "bayesmeta"
  )
  expect_equal(called, 0L)
  expect_identical(
    get("mtext", envir = imports, inherits = FALSE),
    mtext_before
  )

  objs_ov <- traceplot(fit, type = "combined")
  expect_length(objs_ov, 1L)
  expect_s3_class(objs_ov[[1L]], "bmr")
  expect_identical(objs_ov[[1L]], fit$bayesmeta$overall)
})

test_that("traceplot() does not substitute bayesmeta on another engine", {
  skip_if_not_installed("metafor")
  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
    method = "REML",
    reference = "Younger",
    silent = 2
  ))
  expect_error(traceplot(fit, type = "interaction"), 'method = "REML"')
  expect_error(traceplot(fit, type = "combined"), "Auxiliary bayesmeta")
})

test_that("traceplot() axis labels depend on type", {
  labs_int <- interact:::comparison_traceplot_axis_labels("interaction")
  expect_equal(labs_int$ylab, "Interaction")
  expect_equal(
    as.character(labs_int$xlab),
    as.character(expression("heterogeneity " * tau[gamma]))
  )

  labs_comb <- interact:::comparison_traceplot_axis_labels("combined")
  expect_equal(labs_comb$ylab, "Effect")
  expect_equal(
    as.character(labs_comb$xlab),
    as.character(expression("heterogeneity " * tau))
  )

  expect_equal(
    interact:::comparison_traceplot_axis_labels(
      "interaction",
      ylab = "custom"
    )$ylab,
    "custom"
  )
})

test_that("comparison_traceplot_default_bmr_X() uses stored trial IF", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
        reference = "Younger",
    silent = 2
  ))
  X <- interact:::comparison_traceplot_default_bmr_X(fit)
  expect_true(is.matrix(X))
  expect_equal(dim(X), c(1L, 2L))
  expect_equal(unname(X[1L, 1L]), 1)
  expect_true(is.finite(X[1L, 2L]))
  expect_true(X[1L, 2L] > 0 && X[1L, 2L] < 1)
})
