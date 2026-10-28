test_that("SWADA interaction prediction uses common-weight WLS + tau", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("forestplot")

  sw <- swada(respire_age, reference = "<65", method = "REML")
  ad <- ad(respire_age, reference = "<65", method = "REML")

  pr_sw <- interactmeta:::interactmeta_swada_prediction(sw, target = "interaction", prob = 0.95)

  int_rma <- ad$fits$interaction
  pr_metafor <- suppressWarnings(predict(int_rma, level = 95))

  expect_equal(pr_sw$estimate, sw$interaction$estimate, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(pr_sw$se, as.numeric(pr_metafor$se), tolerance = 1e-6)))

  fp_sw <- forestplot(sw, type = "interaction", include_shrinkage = FALSE)
  pred_fp <- fp_sw$specs$panels[[1L]]$yi[
    fp_sw$specs$panels[[1L]]$row_kind == "prediction"
  ]
  expect_length(pred_fp, 1L)
  expect_equal(pred_fp[[1L]], pr_sw$estimate, tolerance = 1e-10)

  int_api <- interactions(sw, type = "Prediction", prob = 0.95)
  expect_equal(nrow(int_api$summary), 1L)
  expect_equal(int_api$summary$type, "Prediction")
  expect_equal(int_api$summary$estimate, pr_sw$estimate, tolerance = 1e-10)
  expect_equal(int_api$summary$lci, pr_sw$ci.lb, tolerance = 1e-10)

  both <- interactions(sw, type = c("Overall", "Prediction"), prob = 0.95)
  expect_equal(nrow(both$summary), 2L)
  expect_setequal(both$summary$type, c("Overall", "Prediction"))

  expect_error(
    interactions(da(respire_age, reference = "<65", method = "REML"), type = "Prediction"),
    "does not support"
  )
})
