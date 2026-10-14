test_that("method = mle is rejected", {
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06),
    n = 20
  )
  expect_error(
    cam(dat, reference = "Younger", method = "mle"),
    "mle"
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "metafor", method = "mle"),
    "mle"
  )
})

test_that("bayescam() uses bayesmeta DIRECT", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.02, 0.05, 0.03),
    n = 10
  )
  fit_bm <- cam(dat, reference = "Younger", engine = "bayesmeta", silent = 2)
  expect_s3_class(fit_bm, "interactmeta")
  expect_equal(fit_bm$engine, "bayesmeta")
  expect_equal(fit_bm$method, "bayes")

  fit <- suppressWarnings(bayescam(
    dat,
    engine = "bayesmeta",
    reference = "Younger",
    silent = 2
  ))
  expect_equal(fit$method, "bayes")
  expect_equal(fit$algorithm, "quadrature")
  expect_false(interact:::comparison_is_mle(fit))
  expect_gt(nrow(fit$draws), 1L)
})
