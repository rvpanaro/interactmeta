test_that("DA interaction forest uses subgroup-specific heterogeneity footnote", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  da_fit <- da(age, reference = "Younger", method = "REML")
  fp <- forestplot(da_fit, type = "interaction", include_shrinkage = FALSE)
  het <- fp$specs$annotations$heterogeneity
  expect_true("combined" %in% names(het))
  foot <- het[["combined"]]
  expect_match(foot, "tau\\[Younger\\]")
  expect_match(foot, "tau\\[Older\\]")
  expect_false(grepl("tau\\[interaction\\]", foot))
})

test_that("frequentist DA forestplot omits shrinkage rows by default", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  da_fit <- da(age, reference = "Younger", method = "REML")
  fp <- forestplot(da_fit)
  panel <- fp$specs$panels[[1L]]
  expect_false(any(panel$row_kind == "shrinkage", na.rm = TRUE))
})

test_that("AD forest includes interaction Overall and Prediction rows", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  ad_fit <- ad(age, reference = "Younger", method = "REML")
  fp <- forestplot(ad_fit, include_shrinkage = FALSE)
  panel <- fp$specs$panels[[1L]]
  kinds <- panel$row_kind
  expect_true(any(kinds == "overall", na.rm = TRUE))
  expect_true(any(kinds == "prediction", na.rm = TRUE))
  ov <- panel$yi[kinds == "overall"]
  expect_true(any(is.finite(ov)))
  expect_equal(kinds[[1L]], "study")
  expect_true(tail(kinds, 2L)[1L] == "overall")
  expect_true(tail(kinds, 1L) == "prediction")
})

test_that("DA interaction forest Prediction uses DA not AD RE-MA", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("metafor")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    subgroup = rep(c("NIV", "IV"), 2),
    yi = c(0, 0.5, 0.1, 0.8),
    vi = rep(0.05, 4),
    stringsAsFactors = FALSE
  )
  da_fit <- da(dat, reference = "NIV", trial = "study", subgroup = "subgroup", method = "REML")
  ad_fit <- ad(dat, reference = "NIV", trial = "study", subgroup = "subgroup", method = "REML")
  fp_da <- forestplot(da_fit, type = "interaction", include_shrinkage = FALSE)
  fp_ad <- forestplot(ad_fit, include_shrinkage = FALSE)
  ra_da <- fp_da$specs$panels[[1L]]
  ra_ad <- fp_ad$specs$panels[[1L]]
  pred_da <- ra_da$yi[ra_da$row_kind == "prediction"]
  pred_ad <- ra_ad$yi[ra_ad$row_kind == "prediction"]
  expect_length(pred_da, 1L)
  expect_length(pred_ad, 1L)
  pr_ref <- interact:::construction_forest_rma_marginal(
    da_fit$fits$reference,
    0.95,
    "Prediction"
  )
  pr_non <- interact:::construction_forest_rma_marginal(
    da_fit$fits$nonref,
    0.95,
    "Prediction"
  )
  expect_equal(pred_da[[1L]], pr_non$estimate - pr_ref$estimate, tolerance = 1e-10)
  int_rma <- interact:::construction_forest_interaction_rma(da_fit)
  pr_ad <- interact:::construction_forest_rma_marginal(int_rma, 0.95, "Prediction")
  expect_equal(pred_ad[[1L]], pr_ad$estimate, tolerance = 1e-10)
})

test_that("compare_forestplots aligns AD and DA interaction row order", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("patchwork")
  fit_da <- da(age, reference = "Younger", method = "REML")
  fit_ad <- ad(age, reference = "Younger", method = "REML")
  fp_da <- forestplot(fit_da, type = "interaction", include_shrinkage = FALSE)
  fp_ad <- forestplot(fit_ad, include_shrinkage = FALSE)
  expect_equal(
    fp_da$specs$panels[[1L]]$row_kind,
    fp_ad$specs$panels[[1L]]$row_kind
  )
})
