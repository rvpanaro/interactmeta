test_that("new-trial Prediction RE draws are reproducible", {
  post <- data.frame(
    sd_trial__Intercept = rep(0.2, 50),
    sd_trial__ifrac2_Older = rep(0.1, 50)
  )
  fit <- structure(list(), class = "interactmeta")
  attr(fit, "interactmeta") <- list(
    trial_var = "trial",
    ifrac2_var = "ifrac2_Older",
    ifrac2_vars_std = "ifrac2_Older"
  )
  a <- interact:::comparison_new_trial_re_draws(fit, ifrac = 0.4, post = post)
  b <- interact:::comparison_new_trial_re_draws(fit, ifrac = 0.4, post = post)
  expect_equal(a, b)
  c <- interact:::comparison_new_trial_slope_re_draws(fit, post = post)
  d <- interact:::comparison_new_trial_slope_re_draws(fit, post = post)
  expect_equal(c, d)
})

test_that("subgroups() returns summary per subgroup", {
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
  me <- subgroups(fit)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_equal(nrow(me$summary), 4L)
  expect_true(all(c("Younger", "Older") %in% me$summary$subgroup))
  expect_true(all(c("Overall", "Prediction") %in% me$summary$type))
  expect_true(all(is.finite(me$summary$estimate)))
  expect_true(all(is.finite(me$summary$mean)))
})

test_that("subgroups(type = 'Prediction') returns new-trial draws only", {
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
  me <- subgroups(fit, type = "Prediction")
  expect_equal(unique(me$draws$type), "Prediction")
  expect_equal(nrow(me$summary), 2L)
})

test_that("Prediction intervals are wider than Overall in subgroups()", {
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
  me <- subgroups(fit)
  for (sg in c("Younger", "Older")) {
    pooled <- me$summary[me$summary$subgroup == sg & me$summary$type == "Overall", , drop = FALSE]
    pred <- me$summary[me$summary$subgroup == sg & me$summary$type == "Prediction", , drop = FALSE]
    expect_gt(pred$uci - pred$lci, pooled$uci - pooled$lci)
  }
})

test_that("subgroups(weight = interaction) uses interaction-weighted IF", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.01, 0.04, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.6, 0.6),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me_overall <- subgroups(fit, ifrac = "overall", type = "Overall", exponentiate = TRUE, summary = FALSE)
  me_int <- subgroups(
    fit,
    ifrac = "interaction",
    type = "Overall",
    exponentiate = TRUE,
    summary = FALSE
  )
  expect_equal(me_overall$ifrac, "overall")
  expect_equal(me_int$ifrac, "interaction")
  expect_false(isTRUE(all.equal(me_overall$draws$effect, me_int$draws$effect)))
})

test_that("subgroups(exponentiate = TRUE) transforms draws", {
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
  log <- subgroups(fit, type = "Overall", summary = FALSE)
  rat <- subgroups(fit, type = "Overall", exponentiate = TRUE, summary = FALSE)
  expect_true(rat$exponentiate)
  expect_equal(rat$draws$effect, exp(log$draws$effect))
})
