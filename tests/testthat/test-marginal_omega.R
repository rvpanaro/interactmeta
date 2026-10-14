test_that("subgroups() defaults to ifrac(fit) draws", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- subgroups(fit, type = "Overall", summary = FALSE)
  expect_equal(me$overall_ifrac_source, "fit")
  expect_equal(length(me$overall_ifrac), length(ifrac(fit)))
})

test_that("subgroups() accepts a custom scalar contribution", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  oi <- 0.35
  me <- subgroups(fit, contribution = oi, type = "Overall", summary = FALSE)
  expect_equal(me$overall_ifrac_source, "user")
  expect_gt(nrow(me$draws), 1L)
  expect_gt(length(unique(me$draws$.draw)), 1L)
})

test_that("combined() accepts explicit contribution from ifrac()", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  oi <- ifrac(fit)
  me <- combined(fit, contribution = oi, type = "Overall", summary = FALSE)
  expect_equal(me$overall_ifrac_source, "user")
})

test_that("max_draws cannot thin a weighted DIRECT support", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cams(dat, reference = "Younger", engine = "bayesmeta", silent = 2)
  n <- nrow(fit$draws)
  expect_gt(n, 5L)
  expect_error(
    subgroups(fit, max_draws = 5L),
    "weighted DIRECT support"
  )
  expect_error(
    predict_cams_by_if(fit, resolution = 3, max_draws = 5L),
    "weighted DIRECT support"
  )
  set.seed(1L)
  want <- sample.int(10L, 4L)
  expect_identical(
    interact:::comparison_marginal_draw_index(10L, max_draws = 4L),
    want
  )
  expect_error(
    interact:::comparison_marginal_draw_index(
      n,
      max_draws = 5L,
      weight = fit$draws$.weight
    ),
    "weighted DIRECT support"
  )
})
