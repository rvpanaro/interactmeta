test_that("interactions() returns contrast summaries at overall IF", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- interactions(fit)
  expect_s3_class(me, "interactmeta_interaction_effects")
  expect_equal(nrow(me$summary), 2L)
  expect_equal(unique(me$summary$interaction), "Older-Younger")
  expect_true(all(is.finite(me$summary$estimate)))
})

test_that("interactions(shrinkage = TRUE) print shows one row per contrast", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- interactions(fit, shrinkage = TRUE)
  expect_equal(nrow(me$summary), 2L)
  once <- interactmeta:::comparison_interaction_shrinkage_summary_once(me$summary)
  expect_equal(nrow(once), 1L)
  expect_equal(unique(once$interaction), "Older-Younger")
})

test_that("interactions(shrinkage = TRUE) returns one row per trial", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- interactions(fit, shrinkage = TRUE)
  expect_true(me$shrinkage)
  expect_equal(unique(me$draws$type), "Shrinkage")
  expect_equal(nrow(me$summary), 2L)
  expect_true(all(c("A", "B") %in% me$summary$trial))
})

test_that("interaction Prediction intervals are wider than Pooled", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- interactions(fit)
  pooled <- me$summary[me$summary$type == "Pooled", , drop = FALSE]
  pred <- me$summary[me$summary$type == "Prediction", , drop = FALSE]
  expect_gt(pred$uci - pred$lci, pooled$uci - pooled$lci)
})

test_that("interaction Pooled matches subgroup difference at overall IF", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  int <- interactions(fit, type = "Pooled", summary = FALSE)
  sg <- subgroups(fit, type = "Pooled", summary = FALSE)
  young <- sg$draws$effect[sg$draws$subgroup == "Younger"]
  old <- sg$draws$effect[sg$draws$subgroup == "Older"]
  int_eff <- int$draws$effect[int$draws$type == "Pooled"]
  expect_equal(length(int_eff), length(young))
  expect_equal(mean(int_eff), mean(old - young), tolerance = 1e-6)
})
