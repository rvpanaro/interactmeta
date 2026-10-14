test_that("combined() returns summary per type", {
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
  me <- combined(fit)
  expect_s3_class(me, "interactmeta_combined_effect")
  expect_equal(nrow(me$summary), 2L)
  expect_true(all(c("Overall", "Prediction") %in% me$summary$type))
  expect_true(all(is.finite(me$summary$estimate)))
  expect_true(all(is.finite(me$summary$mean)))
})

test_that("Prediction intervals are wider than Overall in combined()", {
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
  me <- combined(fit)
  w <- function(type) {
    row <- me$summary[me$summary$type == type, , drop = FALSE]
    row$uci - row$lci
  }
  expect_gt(w("Prediction"), w("Overall"))
})

test_that("combined() differs from reference subgroup pooled effect", {
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
  ov <- combined(fit, type = "Overall", summary = FALSE)
  expect_true(nrow(ov$draws) > 0L)
  expect_true("effect" %in% names(ov$draws))
  sg <- subgroups(fit, type = "Overall", summary = FALSE)
  ref_eff <- sg$draws$effect[sg$draws$subgroup == "Younger"]
  expect_false(isTRUE(all.equal(ov$draws$effect, ref_eff)))
})

test_that("combined(weight = interaction) uses interaction-weighted IF", {
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
  me_overall <- combined(fit, ifrac = "overall", type = "Overall", exponentiate = TRUE, summary = FALSE)
  me_int <- combined(
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

test_that("combined(exponentiate = TRUE) transforms draws", {
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
  log <- combined(fit, type = "Overall", summary = FALSE)
  rat <- combined(fit, type = "Overall", exponentiate = TRUE, summary = FALSE)
  expect_true(rat$exponentiate)
  expect_equal(rat$draws$effect, exp(log$draws$effect))
})
