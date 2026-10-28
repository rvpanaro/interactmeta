test_that("method = mle is rejected", {
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06),
    n = 20
  )
  expect_error(
    cams(dat, reference = "Younger", method = "mle"),
    "mle"
  )
})

test_that("cams() uses bayesmeta DIRECT", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.02, 0.05, 0.03),
    n = 10
  )
  fit_bm <- cams(dat, reference = "Younger", silent = 2)
  expect_s3_class(fit_bm, "interactmeta")
  expect_equal(interactmeta_fit_backend(fit_bm), "bayesmeta")
  expect_equal(fit_bm$method, "Bayesian")

  fit <- suppressWarnings(cams(
    dat,
        reference = "Younger",
    silent = 2
  ))
  expect_equal(fit$method, "Bayesian")
  expect_equal(fit$algorithm, "quadrature")
  expect_gt(nrow(fit$support), 1L)
})
