# Slow integration smokes: maximize unique line coverage with shared fixtures.
# Behavioural assertions live in test-battery-*.R; this file is for covr breadth.

test_that("coverage: CAMS marginals and reporting", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("ggplot2")
  s <- fixture_coverage_suite()
  fit <- s$bayes_tiny

  for (tp in c("Overall", "Prediction")) {
    suppressMessages(modifiers(fit, type = tp, summary = TRUE))
    suppressMessages(interactions(fit, type = tp, summary = TRUE))
    suppressMessages(combined(fit, type = tp, summary = TRUE))
  }
  suppressMessages(heterogeneity(fit, summary = TRUE))
  suppressMessages(shrinkage(fit, type = "modifiers", summary = TRUE))
  suppressMessages(shrinkage(fit, type = "interaction", summary = TRUE))

  expect_s3_class(tidy(fit), "tbl_df")
  expect_s3_class(summary(fit), "summary.interactmeta")
  expect_output(print(fit))

  with_null_pdf(expect_invisible(plot(fit, type = "contribution", ifrac_var = "IF_Older")))
  ifrac(fit)
  ifrac(s$reml_respire, ifrac = "overall")

  expect_type(native_fit(fit), "list")
  expect_type(observed_effects(fit), "list")
  suppressMessages(standardized_trial_effects(fit, summary = TRUE))
  suppressMessages(shrinkage_components(s$bayes_tiny))
  suppressMessages(shrinkage_diagnostics(s$bayes_tiny))
  suppressMessages(predict_cams_by_if(fit, ifrac_var = "IF_Older", if_grid = c(0.3, 0.5, 0.7)))
  suppressMessages(cam_trial_quantities(fit))
  suppressMessages(cams_trial_quantities(fit))
})

test_that("coverage: forests, traceplot, and heterogeneity plots", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("patchwork")
  s <- fixture_coverage_suite()

  scenarios <- list(
    list(fit = s$bayes_tiny, type = "modifiers", include_shrinkage = TRUE, exponentiate = FALSE),
    list(fit = s$bayes_tiny, type = "interaction", include_shrinkage = FALSE, exponentiate = FALSE),
    list(fit = s$reml_respire, type = "modifiers", include_shrinkage = FALSE, exponentiate = FALSE),
    list(fit = s$da, type = "modifiers", include_shrinkage = FALSE, exponentiate = FALSE),
    list(fit = s$ad, type = "interaction", include_shrinkage = FALSE, exponentiate = FALSE),
    list(fit = s$swada, type = "modifiers", include_shrinkage = FALSE, exponentiate = TRUE)
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fps <- list()
  for (sc in scenarios) {
    args <- sc
    fit <- args$fit
    args$fit <- NULL
    fp <- do.call(forestplot, c(list(x = fit), args))
    expect_s3_class(fp, "interactmeta_forestplot")
    forestplot_height(fp, width = 9.5)
    forestplot_draw(fp)
    if (identical(sc$fit, s$bayes_tiny) && identical(sc$type, "modifiers")) {
      fps$mod <- fp
    }
    if (identical(sc$fit, s$bayes_tiny) && identical(sc$type, "interaction")) {
      fps$int <- fp
    }
  }

  traceplot(s$bayes_tiny, type = "interaction")
  traceplot(s$bayes_tiny, type = "combined", prior = FALSE)
  plot(s$bayes_tiny, type = "heterogeneity", which = "both", prior = FALSE)
  plot(s$reml_respire, type = "heterogeneity", which = "both", prior = FALSE)

  fp_mod <- forestplot(s$bayes_tiny, type = "modifiers", include_shrinkage = FALSE, exponentiate = FALSE)
  fp_int <- forestplot(s$bayes_tiny, type = "interaction", include_shrinkage = FALSE, exponentiate = FALSE)
  compare_forestplots(fp_mod, fp_int, draw = TRUE)
})

test_that("coverage: constructions, compare(), and construction forests", {
  slow()
  skip_if_not_installed("metafor")
  skip_if_not_installed("ggplot2")
  s <- fixture_coverage_suite()
  sens <- s$sens

  tab_int <- compare(sens, estimand = "interaction", type = "Overall")
  tab_mod <- compare(sens, estimand = "modifiers")
  expect_s3_class(tab_int, "interactmeta_compare_interactions")
  expect_s3_class(tab_mod, "interactmeta_compare_subgroups")

  compare_weights(sens$fits$da, sens$fits$ad, sens$fits$swada)
  compare_heterogeneity(sens$fits$da, sens$fits$cams)

  cam_cmp <- compare(s$reml_respire, s$reml_ml, exponentiate = FALSE)
  expect_s3_class(cam_cmp, "interactmeta_compare_interactions")

  ifrac_cmp <- interactmeta:::compare_ifrac(s$reml_respire, s$reml_ml)
  expect_s3_class(ifrac_cmp, "interactmeta_compare_ifrac")

  comb_cmp <- interactmeta:::compare_combined(s$reml_respire, s$reml_ml)
  expect_s3_class(comb_cmp, "interactmeta_compare_combined")

  with_null_pdf({
    p1 <- plot(tab_int)
    expect_s3_class(p1, "ggplot")
    p2 <- plot(ifrac_cmp)
    expect_s3_class(p2, "ggplot")
  })

  conv <- conventional(respire_age, reference = "<65", method = "REML")
  expect_s3_class(conv, "interactmeta_conventional")
  expect_output(print(conv))
  tidy(conv)

  interactmeta_capability_matrix()
  x <- log_hr_from_ci(1.0, 0.8, 1.2, ci_level = 0.95)
  hr_ci_from_log(x$yi, x$sei, ci_level = 0.95)
})

test_that("coverage: misc helpers and wide subgroup_long", {
  slow()
  wide <- data.frame(
    study = "T1",
    yi_Younger = 0.1,
    vi_Younger = 0.04,
    yi_Older = 0.2,
    vi_Older = 0.03,
    stringsAsFactors = FALSE
  )
  long <- subgroup_long(
    wide,
    subgroups = c("Younger", "Older"),
    id_cols = "study",
    values = c("yi", "vi")
  )
  expect_equal(nrow(long), 2L)
  expect_true("subgroup" %in% names(long))

  es <- fixture_age_escalc()
  des <- ifrac(es, reference = "Younger", include_interaction = TRUE)
  expect_gt(nrow(des), nrow(es))

  expect_true(is.character(interactmeta_palette()))
  expect_true(is.character(names(interactmeta_palette())))
  expect_s3_class(interactmeta_theme_report(), "theme")
})

test_that("coverage: expected deprecation and validation errors", {
  slow()
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  da_fit <- da(respire_age, reference = "<65", method = "REML")

  expect_error(modifiers(da_fit, at = 0.5), "CAMS")
  expect_error(traceplot(fit, type = "modifier"), "not available")
  expect_warning(modifiers(fit, contribution = 0.5), "contribution")
  expect_error(modifiers(fit, at = 0.5, contribution = 0.5), "only one")
})