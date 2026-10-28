# Effect accessors: one cached tiny fit per variant (see helper-fixtures.R).

test_that("new-trial Prediction RE draws are reproducible", {
  post <- data.frame(
    sd_trial__Intercept = rep(0.2, 50),
    sd_trial__ifrac2_Older = rep(0.1, 50)
  )
  fit <- structure(list(), class = "interactmeta")
  attr(fit, "interactmeta") <- list(
    trial_var = "study",
    ifrac2_var = "ifrac2_Older",
    ifrac2_vars_std = "ifrac2_Older"
  )
  a <- interactmeta:::comparison_new_trial_re_draws(fit, ifrac = 0.4, post = post)
  b <- interactmeta:::comparison_new_trial_re_draws(fit, ifrac = 0.4, post = post)
  expect_equal(a, b)
  c <- interactmeta:::comparison_new_trial_slope_re_draws(fit, post = post)
  d <- interactmeta:::comparison_new_trial_slope_re_draws(fit, post = post)
  expect_equal(c, d)
})

test_that("interactions() returns contrast summaries at overall IF", {
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  me <- interactions(fit)
  expect_s3_class(me, "interactmeta_interaction_effects")
  expect_equal(nrow(me$summary), 2L)
  expect_equal(unique(me$summary$interaction), "Older-Younger")
  expect_true(all(is.finite(me$summary$estimate)))
})

test_that("interactions(shrinkage = TRUE) print shows one row per contrast", {
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  me <- interactions(fit, shrinkage = TRUE)
  expect_equal(nrow(me$summary), 2L)
  once <- interactmeta:::comparison_interaction_shrinkage_summary_once(me$summary)
  expect_equal(nrow(once), 1L)
  expect_equal(unique(once$interaction), "Older-Younger")
})

test_that("interactions(shrinkage = TRUE) returns one row per trial", {
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  me <- interactions(fit, shrinkage = TRUE)
  expect_true(me$shrinkage)
  expect_equal(unique(me$support$type), "Shrinkage")
  expect_equal(nrow(me$summary), 2L)
  expect_true(all(c("A", "B") %in% me$summary$trial))
})

test_that("interaction Prediction intervals are wider than Overall", {
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  me <- interactions(fit)
  pooled <- me$summary[me$summary$type == "Overall", , drop = FALSE]
  pred <- me$summary[me$summary$type == "Prediction", , drop = FALSE]
  expect_gt(pred$uci - pred$lci, pooled$uci - pooled$lci)
})

test_that("interaction Overall matches subgroup difference at overall IF", {
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  int <- interactions(fit, type = "Overall", summary = TRUE)
  sg <- modifiers(fit, type = "Overall", summary = TRUE)
  young <- sg$summary$estimate[
    sg$summary$level == "Younger" & sg$summary$type == "Overall"
  ]
  old <- sg$summary$estimate[
    sg$summary$level == "Older" & sg$summary$type == "Overall"
  ]
  int_est <- int$summary$estimate[int$summary$type == "Overall"]
  expect_equal(length(int_est), length(young))
  expect_equal(int_est, old - young, tolerance = 0.01)
})

test_that("modifiers() returns summary per subgroup", {
  fit <- fixture_cam_tiny_fit()
  me <- modifiers(fit)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_equal(nrow(me$summary), 4L)
  expect_true(all(c("Younger", "Older") %in% me$summary$level))
  expect_true(all(c("Overall", "Prediction") %in% me$summary$type))
  expect_true(all(is.finite(me$summary$estimate)))
  expect_true(all(is.finite(me$summary$mean)))
})

test_that("modifiers(type = 'Prediction') returns new-trial draws only", {
  fit <- fixture_cam_tiny_fit()
  me <- modifiers(fit, type = "Prediction")
  expect_equal(unique(me$support$type), "Prediction")
  expect_equal(nrow(me$summary), 2L)
})

test_that("Prediction intervals are wider than Overall in modifiers()", {
  fit <- fixture_cam_tiny_fit()
  me <- modifiers(fit)
  for (sg in c("Younger", "Older")) {
    pooled <- me$summary[me$summary$level == sg & me$summary$type == "Overall", , drop = FALSE]
    pred <- me$summary[me$summary$level == sg & me$summary$type == "Prediction", , drop = FALSE]
    expect_gt(pred$uci - pred$lci, pooled$uci - pooled$lci)
  }
})

test_that("modifiers(weight = interaction) uses interaction-weighted IF", {
  dat <- fixture_cam_tiny(
    ifrac_older = c(0.3, 0.3, 0.6, 0.6),
    ifrac2_older = c(-0.3, 0.7, -0.4, 0.6),
    vi = c(0.01, 0.04, 0.02, 0.02)
  )
  fit <- fixture_cam_tiny_fit(variant = "weighted", dat = dat)
  me_overall <- modifiers(fit, ifrac = "overall", type = "Overall", exponentiate = TRUE, summary = FALSE)
  me_int <- modifiers(
    fit,
    ifrac = "interaction",
    type = "Overall",
    exponentiate = TRUE,
    summary = FALSE
  )
  expect_equal(me_overall$ifrac, "overall")
  expect_equal(me_int$ifrac, "interaction")
  expect_false(isTRUE(all.equal(me_overall$support$effect, me_int$support$effect)))
})

test_that("modifiers(exponentiate = TRUE) transforms draws", {
  fit <- fixture_cam_tiny_fit()
  log <- modifiers(fit, type = "Overall", summary = FALSE)
  rat <- modifiers(fit, type = "Overall", exponentiate = TRUE, summary = FALSE)
  expect_true(rat$exponentiate)
  expect_equal(rat$support$effect, exp(log$support$effect))
})

test_that("modifiers() defaults to ifrac(fit) draws", {
  fit <- fixture_cam_tiny_fit()
  me <- modifiers(fit, type = "Overall", summary = FALSE)
  expect_equal(me$overall_ifrac_source, "fit")
  expect_equal(length(me$overall_ifrac), length(ifrac(fit)))
})

test_that("modifiers() accepts a custom scalar at (reporting IF)", {
  fit <- fixture_cam_tiny_fit()
  oi <- 0.35
  me <- modifiers(fit, at = oi, type = "Overall", summary = FALSE)
  expect_equal(me$overall_ifrac_source, "user")
  expect_gt(nrow(me$support), 1L)
  expect_gt(length(unique(me$support$.draw)), 1L)
})

test_that("ifrac print uses marginal effect columns and single interval label", {
  fit <- fixture_cam_tiny_fit()
  oi_out <- paste(capture.output(print(ifrac(fit))), collapse = "\n")
  expect_match(oi_out, "median", fixed = TRUE)
  expect_match(oi_out, "mean", fixed = TRUE)
  expect_match(oi_out, "sd", fixed = TRUE)
  expect_match(oi_out, "95% CrI", fixed = TRUE)
  expect_false(grepl("Overall 95% CrI", oi_out, fixed = TRUE))
  expect_false(grepl("Prediction", oi_out, fixed = TRUE))
})

test_that("subgroups, combined, and interactions print tables without preamble", {
  fit <- fixture_cam_tiny_fit()
  sg_out <- paste(capture.output(print(modifiers(fit, type = "Overall"))), collapse = "\n")
  expect_match(sg_out, "Overall 95% CrI", fixed = TRUE)
  expect_false(grepl("Reporting IF:", sg_out, fixed = TRUE))
  expect_false(grepl("Subgroup-specific effects", sg_out, fixed = TRUE))

  comb_out <- paste(capture.output(print(combined(fit, type = "Overall"))), collapse = "\n")
  expect_match(comb_out, "combined", fixed = TRUE)
  expect_match(comb_out, "All", fixed = TRUE)
  expect_match(comb_out, "Overall 95% CrI", fixed = TRUE)
  expect_false(grepl("Combined-population treatment effect", comb_out, fixed = TRUE))

  int_out <- paste(capture.output(print(interactions(fit, type = "Overall"))), collapse = "\n")
  expect_match(int_out, "Overall 95% CrI", fixed = TRUE)
  expect_false(grepl("Treatment-by-subgroup interaction", int_out, fixed = TRUE))
})

test_that("contribution= is deprecated in favour of at=", {
  fit <- fixture_cam_tiny_fit()
  expect_warning(
    modifiers(fit, contribution = 0.5, type = "Overall", summary = FALSE),
    "deprecated"
  )
  expect_error(
    modifiers(fit, at = 0.5, contribution = 0.5, type = "Overall", summary = FALSE),
    "only one"
  )
})

test_that("combined() returns summary per type", {
  fit <- fixture_cam_tiny_fit()
  me <- combined(fit)
  expect_s3_class(me, "interactmeta_combined_effect")
  expect_equal(nrow(me$summary), 2L)
  expect_true(all(c("Overall", "Prediction") %in% me$summary$type))
  expect_true(all(is.finite(me$summary$estimate)))
  expect_true(all(is.finite(me$summary$mean)))
})

test_that("Prediction intervals are wider than Overall in combined()", {
  fit <- fixture_cam_tiny_fit()
  me <- combined(fit)
  w <- function(type) {
    row <- me$summary[me$summary$type == type, , drop = FALSE]
    row$uci - row$lci
  }
  expect_gt(w("Prediction"), w("Overall"))
})

test_that("combined() differs from reference subgroup pooled effect", {
  fit <- fixture_cam_tiny_fit()
  ov <- combined(fit, type = "Overall", summary = TRUE)
  sg <- modifiers(fit, type = "Overall", summary = TRUE)
  ref_eff <- sg$summary$estimate[
    sg$summary$level == "Younger" & sg$summary$type == "Overall"
  ]
  expect_false(isTRUE(all.equal(ov$summary$estimate, ref_eff)))
})

test_that("combined(weight = interaction) uses interaction-weighted IF", {
  dat <- fixture_cam_tiny(
    ifrac_older = c(0.3, 0.3, 0.6, 0.6),
    ifrac2_older = c(-0.3, 0.7, -0.4, 0.6),
    vi = c(0.01, 0.04, 0.02, 0.02)
  )
  fit <- fixture_cam_tiny_fit(variant = "combined_weight", dat = dat)
  me_overall <- combined(fit, ifrac = "overall", type = "Overall", exponentiate = TRUE, summary = TRUE)
  me_int <- combined(
    fit,
    ifrac = "interaction",
    type = "Overall",
    exponentiate = TRUE,
    summary = TRUE
  )
  expect_equal(me_overall$ifrac, "overall")
  expect_equal(me_int$ifrac, "interaction")
  expect_false(isTRUE(all.equal(me_overall$summary$estimate, me_int$summary$estimate)))
})

test_that("combined() accepts explicit at from ifrac()", {
  fit <- fixture_cam_tiny_fit()
  oi <- ifrac(fit)
  me <- combined(fit, at = oi, type = "Overall", summary = FALSE)
  expect_equal(me$overall_ifrac_source, "user")
})

test_that("tidy.interactmeta_subgroup_effects() returns summary columns", {
  skip_if_not_installed("generics")
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  me <- modifiers(fit)
  td <- generics::tidy(me)
  expect_s3_class(td, "tbl_df")
  expect_true(all(c("estimate", "std.error", "conf.low", "conf.high") %in% names(td)))
  expect_true(all(c("level", "type") %in% names(td)))
  expect_false("term" %in% names(td))
  expect_false("lci" %in% names(td))
  expect_false("sd" %in% names(td))
  expect_equal(nrow(td), nrow(me$summary))
})

test_that("tidy.interactmeta_interaction_effects() deduplicates shrinkage summaries", {
  skip_if_not_installed("generics")
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  me <- interactions(fit, shrinkage = TRUE)
  td <- generics::tidy(me)
  expect_equal(nrow(td), 1L)
  expect_equal(unique(td$interaction), "Older-Younger")
})

test_that("print methods use biostatistical tables (not tibbles)", {
  skip_if_not_installed("generics")
  fit <- fixture_cam_tiny_fit(variant = "interaction")
  out <- capture.output(print(modifiers(fit)))
  expect_true(any(grepl("level", out, fixed = TRUE)))
  expect_true(any(grepl("Overall 95% CrI", out, fixed = TRUE)))
  expect_true(any(grepl("Prediction", out, fixed = TRUE)))
  expect_false(any(grepl("^\\(", out)))
  expect_false(any(grepl("# A tibble", out, fixed = TRUE)))
  expect_false(any(grepl("^\\s*Prediction\\s+", out)))
  expect_false(any(grepl("^\\s*term\\b", out)))
  ov <- capture.output(print(combined(fit)))
  expect_true(any(grepl("Prediction 95% CrI", ov, fixed = TRUE)))
  expect_false(any(grepl("^\\s*Prediction\\s+", ov)))
  expect_false(any(grepl("^\\s*Overall\\s+", ov)))
})

test_that("max_support cannot thin a weighted DIRECT support", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cams(dat, reference = "Younger", silent = 2)
  n <- nrow(fit$support)
  expect_gt(n, 5L)
  expect_error(modifiers(fit, max_support = 5L), "weighted DIRECT support")
  expect_error(
    predict_cams_by_if(fit, resolution = 3, max_support = 5L),
    "weighted DIRECT support"
  )
  set.seed(1L)
  want <- sample.int(10L, 4L)
  expect_identical(
    interactmeta:::comparison_marginal_draw_index(10L, max_support = 4L),
    want
  )
  expect_error(
    interactmeta:::comparison_marginal_draw_index(
      n,
      max_support = 5L,
      weight = fit$support$.weight
    ),
    "weighted DIRECT support"
  )
})
