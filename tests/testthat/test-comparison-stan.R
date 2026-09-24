test_that("engine = cmdstanr is rejected", {
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "cmdstanr"),
    "cmdstanr"
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "cmdstanr"),
    "removed"
  )
})

test_that("cam() uses the bayesmeta engine by default", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))
  expect_s3_class(fit, "interactmeta")
  expect_equal(attr(fit, "interactmeta")$engine, "bayesmeta")
  post <- as.data.frame(fit)
  expect_true("b_Intercept" %in% names(post) || "b_ifrac_Older" %in% names(post))
  fp <- forestplot(fit)
  expect_s3_class(fp, "interactmeta_forestplot")
})
