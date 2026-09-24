test_that("print.interactmeta and summary.interactmeta are bayesmeta-style", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cam(
    es,
    reference = "Younger",
    engine = "bayesmeta"))

  out <- paste(utils::capture.output(print(fit, digits = 2)), collapse = "\n")
  expect_match(out, "'interactmeta' object", fixed = TRUE)
  expect_match(out, "marginal posterior summary", fixed = TRUE)
  expect_match(out, "Interaction_Older", fixed = TRUE)
  expect_false(grepl("ifrac2_Older", out, fixed = TRUE))
  expect_match(out, "Ecological_Older", fixed = TRUE)
  expect_false(grepl("b_subgroupOlder", out, fixed = TRUE))
  expect_match(out, "tau", fixed = TRUE)
  expect_match(out, "shortest credible intervals", fixed = TRUE)
  expect_match(out, "approach: bayesmeta", fixed = TRUE)
  expect_match(out, "Use summary()", fixed = TRUE)
  expect_false(grepl("combined-population effect", out, fixed = TRUE))
  expect_match(out, "Subgroup:", fixed = TRUE)
  expect_false(grepl("subgroup effects", out, fixed = TRUE))
  expect_false(grepl("interaction effect", out, fixed = TRUE))
  expect_false(grepl("heterogeneity:", out, fixed = TRUE))

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
  expect_false(grepl("Use summary()", sout, fixed = TRUE))
  expect_false(grepl("data (", sout, fixed = TRUE))
  expect_false(grepl("\nmodel:", sout, fixed = TRUE))
  expect_match(sout, "heterogeneity:", fixed = TRUE)
  expect_match(sout, "Subgroup effects", fixed = TRUE)
  expect_match(sout, "combined-population effect", fixed = TRUE)
  expect_match(sout, "interaction effect", fixed = TRUE)
  expect_match(sout, "95% CrI", fixed = TRUE)
  expect_match(sout, "95% Pred", fixed = TRUE)
  # Prediction is a column, not a separate type row in the overall block
  expect_false(grepl("\n\\s*Prediction\\s+", sout))
  expect_false(grepl("# A tibble", sout, fixed = TRUE))
  # tidy() remains a raw tibble
  td <- generics::tidy(s$subgroups)
  expect_s3_class(td, "tbl_df")
  expect_true(any(grepl("# A tibble", capture.output(print(td)), fixed = TRUE)))
  expect_match(sout, "exponentiated", fixed = TRUE)
})

test_that("print.interactmeta for MLE is no longer a supported engine path", {
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06),
    n = 20
  )
  expect_error(
    cam(dat, reference = "Younger", method = "mle", silent = 2),
    "mle"
  )
})

test_that("print header names quadrature nodes, not posterior draws", {
  meta <- list(
    engine = "bayesmeta",
    method = "bayes",
    algorithm = "quadrature",
    mle = FALSE,
    reference = "Younger",
    trials = c("A", "B"),
    subgroups = c("Younger", "Older"),
    n_trials = 2L,
    n_subgroups = 2L,
    n_draws = 384L
  )
  out <- paste(
    utils::capture.output(interactmeta:::comparison_fit_print_header(meta)),
    collapse = "\n"
  )
  expect_match(out, "quadrature nodes: 384", fixed = TRUE)
  expect_false(grepl("posterior draws", out, fixed = TRUE))
  expect_equal(
    interactmeta:::comparison_fit_draws_label(list(engine = "brms", algorithm = "MCMC")),
    "posterior draws"
  )
})

test_that("comparison_median_qi helpers use equal-tail quantiles", {
  set.seed(1)
  x <- rnorm(2000)
  qi <- interactmeta:::comparison_median_qi_vec(x, prob = 0.9)
  expect_equal(qi$estimate, stats::median(x), tolerance = 1e-10)
  expect_true(qi$lci < qi$estimate && qi$uci > qi$estimate)
  expect_equal(qi$prob, 0.9)

  d <- data.frame(
    g = rep(c("a", "b"), each = 500),
    effect = c(rnorm(500, 0), rnorm(500, 1)),
    stringsAsFactors = FALSE
  )
  summ <- interactmeta:::comparison_median_qi_grouped(d, "effect", "g", prob = 0.95)
  expect_equal(nrow(summ), 2L)
  expect_true(all(c("estimate", "mean", "sd", "lci", "uci", "prob") %in% names(summ)))
  expect_true(summ$estimate[summ$g == "b"] > summ$estimate[summ$g == "a"])
})
