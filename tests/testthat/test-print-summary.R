test_that("print.interactmeta and summary.interactmeta are bayesmeta-style", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  fit <- fixture_age_cam_fit()

  out <- paste(utils::capture.output(print(fit, digits = 2)), collapse = "\n")
  expect_match(out, "'interactmeta' object", fixed = TRUE)
  expect_match(out, "Bayesian estimates", fixed = TRUE)
  expect_match(out, "Interaction_Older", fixed = TRUE)
  expect_false(grepl("ifrac2_Older", out, fixed = TRUE))
  expect_match(out, "Ecological_Older", fixed = TRUE)
  expect_false(grepl("b_subgroupOlder", out, fixed = TRUE))
  expect_match(out, "tau", fixed = TRUE)
  expect_match(out, "central credible intervals", fixed = TRUE)
  expect_match(out, "method: Bayesian", fixed = TRUE)
  expect_false(grepl("approach:", out, fixed = TRUE))
  expect_match(out, "Use summary()", fixed = TRUE)
  expect_false(grepl("combined-population effect", out, fixed = TRUE))
  expect_match(out, "2 Level levels", fixed = TRUE)
  expect_match(out, "Level levels:", fixed = TRUE)
  expect_false(grepl("subgroup effects", out, fixed = TRUE))
  expect_false(grepl("interaction effect", out, fixed = TRUE))
  expect_false(grepl("heterogeneity:", out, fixed = TRUE))

  sg_out <- paste(
    capture.output(print(modifiers(fit, type = "Overall"))),
    collapse = "\n"
  )
  expect_match(sg_out, "Overall 95% CrI", fixed = TRUE)
  expect_false(grepl("Interaction-weighted IF:", sg_out, fixed = TRUE))
  expect_false(grepl("Subgroup-specific effects", sg_out, fixed = TRUE))
  expect_false(grepl("central credible interval", sg_out, fixed = TRUE))
  comb_out <- paste(
    capture.output(print(combined(fit, type = "Overall"))),
    collapse = "\n"
  )
  expect_match(comb_out, "Overall 95% CrI", fixed = TRUE)
  expect_false(grepl("Combined-population treatment effect", comb_out, fixed = TRUE))
  expect_false(grepl("central credible interval", comb_out, fixed = TRUE))

  s <- summary(fit, digits = 2, exponentiate = TRUE)
  expect_s3_class(s, "summary.interactmeta")
  expect_true(is.matrix(s$parameters))
  expect_true(all(c("mode", "median", "mean", "sd", "95% lower", "95% upper") %in%
    rownames(s$parameters)))
  expect_true("tau" %in% colnames(s$parameters))
  expect_true("Ecological_Older" %in% colnames(s$parameters))
  expect_true("Interaction_Older" %in% colnames(s$parameters))
  expect_false("ifrac_Older" %in% colnames(s$parameters))
  expect_false("b_subgroupOlder" %in% colnames(s$parameters))
  expect_true(!is.null(s$data))
  expect_true(ncol(s$data) == 2L)

  sout <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(sout, "Summary of 'interactmeta' object", fixed = TRUE)
  expect_false(grepl("Wald confidence intervals", sout, fixed = TRUE))
  expect_match(sout, "overall via bmr, interaction via bayesmeta", fixed = TRUE)
  expect_match(sout, "Quoted intervals are central credible intervals", fixed = TRUE)
  expect_false(grepl("Use summary()", sout, fixed = TRUE))
  expect_false(grepl("data (", sout, fixed = TRUE))
  expect_false(grepl("\nmodel:", sout, fixed = TRUE))
  expect_match(sout, "heterogeneity:", fixed = TRUE)
  expect_false(grepl("Level-specific effects:", sout, fixed = TRUE))
  expect_false(grepl("combined-population effect:", sout, fixed = TRUE))
  expect_false(grepl("interaction effect:", sout, fixed = TRUE))
  expect_match(sout, "Overall 95% CrI", fixed = TRUE)
  expect_match(sout, "Prediction 95% CrI", fixed = TRUE)
  expect_false(grepl("central credible intervals\\)", sout))
  expect_false(grepl("equal-tailed", sout, fixed = TRUE))
  # Prediction is a column, not a separate type row in the overall block
  expect_false(grepl("\n\\s*Prediction\\s+", sout))
  expect_false(grepl("# A tibble", sout, fixed = TRUE))
  # tidy() remains a raw tibble
  td <- generics::tidy(s$subgroups)
  expect_s3_class(td, "tbl_df")
  expect_true(any(grepl("# A tibble", capture.output(print(td)), fixed = TRUE)))
  expect_true(s$exponentiate)
  expect_gt(s$overall$summary$estimate[[1L]], 0.5)
})

test_that("frequentist parameter summary uses metafor model-results layout", {
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- cams(es, reference = "Younger", method = "REML")
  params <- interact:::interactmeta_metafor_parameters_df(fit)
  mat <- interact:::interactmeta_metafor_param_summary(fit)

  expect_true(all(c("model", "term", "estimate", "se", "zval", "pval", "ci.lb", "ci.ub") %in%
    names(params)))
  expect_equal(
    rownames(mat),
    c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
  )
  expect_true("tau" %in% colnames(mat))
  expect_true("tau_gamma" %in% colnames(mat))
  expect_true("Intercept" %in% colnames(mat))
  expect_true("Ecological_Older" %in% colnames(mat))
  expect_true("Interaction_Older" %in% colnames(mat))

  cs_ov <- coef(summary(fit$metafor$overall))
  expect_equal(
    mat["estimate", "Intercept"],
    as.numeric(cs_ov["intrcpt", "estimate"]),
    tolerance = 1e-10
  )
  cs_int <- coef(summary(fit$metafor$interaction))
  expect_equal(
    mat["estimate", "Interaction_Older"],
    as.numeric(cs_int["intrcpt", "estimate"]),
    tolerance = 1e-10
  )
  lin <- interact:::interactmeta_metafor_lincomb(fit, c(0, 1), -1)
  expect_equal(
    mat["estimate", "Ecological_Older"],
    lin$estimate,
    tolerance = 1e-10
  )

  tau <- params[params$term == "tau" & params$model == "overall", , drop = FALSE]
  expect_true(all(is.finite(tau$ci.lb)))
  expect_true(all(is.finite(tau$ci.ub)))
  expect_equal(
    as.numeric(tau[1L, c("ci.lb", "ci.ub")]),
    as.numeric(as.data.frame(stats::confint(fit$metafor$overall)$random)[
      "tau", c("ci.lb", "ci.ub")
    ])
  )
  tau_i <- params[params$term == "tau" & params$model == "interaction", , drop = FALSE]
  expect_equal(
    as.numeric(tau_i[1L, c("ci.lb", "ci.ub")]),
    as.numeric(as.data.frame(stats::confint(fit$metafor$interaction)$random)[
      "tau", c("ci.lb", "ci.ub")
    ])
  )
})

test_that("print.interactmeta for MLE is no longer a supported engine path", {
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06),
    n = 20
  )
  expect_error(
    cams(dat, reference = "Younger", method = "mle", silent = 2),
    "mle"
  )
})

test_that("print header omits quadrature metadata", {
  meta <- list(
    engine = "bayesmeta",
    method = "bayes",
    inference = "bayesian",
    algorithm = "quadrature",
    mle = FALSE,
    reference = "Younger",
    trials = c("A", "B"),
    subgroups = c("Younger", "Older"),
    n_trials = 2L,
    n_subgroups = 2L,
    n_support = 384L
  )
  out <- paste(
    utils::capture.output(interact:::comparison_fit_print_header(meta)),
    collapse = "\n"
  )
  expect_match(out, "method: Bayesian", fixed = TRUE)
  expect_false(grepl("quadrature", out, ignore.case = TRUE))
  expect_false(grepl("384", out, fixed = TRUE))
})

test_that("comparison_median_qi helpers summarize intervals", {
  set.seed(1)
  x <- rnorm(2000)
  qi <- interact:::comparison_median_qi_vec(x, prob = 0.9)
  expect_equal(qi$estimate, stats::median(x), tolerance = 1e-10)
  expect_true(qi$lci < qi$estimate && qi$uci > qi$estimate)
  expect_equal(qi$prob, 0.9)

  d <- data.frame(
    g = rep(c("a", "b"), each = 500),
    effect = c(rnorm(500, 0), rnorm(500, 1)),
    stringsAsFactors = FALSE
  )
  summ <- interact:::comparison_median_qi_grouped(d, "effect", "g", prob = 0.95)
  expect_equal(nrow(summ), 2L)
  expect_true(all(c("estimate", "mean", "sd", "lci", "uci", "prob") %in% names(summ)))
  expect_true(summ$estimate[summ$g == "b"] > summ$estimate[summ$g == "a"])
})
