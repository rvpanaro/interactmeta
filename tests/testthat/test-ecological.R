test_that("ecological_slope is no longer a cams() argument", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("T1", "T2", "T3"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.10, 0.30, 0.12, 0.28, 0.08, 0.24),
    vi = c(0.04, 0.04, 0.05, 0.03, 0.03, 0.05)
  )
  fit <- suppressWarnings(cams(dat, reference = "Younger", silent = 2))
  expect_true(isTRUE(fit$ecological_slope))
  expect_error(
    cams(dat, reference = "Younger", ecological_slope = FALSE, silent = 2),
    "ecological_slope"
  )
  expect_error(
    cams(dat, reference = "Younger", ecological = FALSE, silent = 2),
    "ecological_slope"
  )
})
