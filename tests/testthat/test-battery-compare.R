# compare(), capability matrix, compare_weights().

test_that("capability matrix lists CAMS on bayesmeta and metafor only", {
  cap <- interactmeta_capability_matrix()
  cam <- cap[cap$analysis == "cam", , drop = FALSE]
  expect_true(all(cam$engine %in% c("bayesmeta", "metafor")))
  expect_equal(length(unique(cam$engine)), 2L)
})

test_that("compare() labels DA, AD, SWADA, and CAMS fits", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  ad_fit <- ad(respire_age, reference = "<65", method = "REML")
  sw <- swada(respire_age, reference = "<65", method = "REML")
  cam <- suppressWarnings(cams(respire_age, reference = "<65", method = "REML", silent = 2))

  lab <- suppressWarnings(compare(da_fit, ad_fit, sw, cam))
  expect_s3_class(lab, "interactmeta_compare_interactions")
  expect_equal(nrow(lab), 4L)
  expect_equal(sort(as.character(unique(lab$approach))), c("AD", "CAMS", "DA", "SWADA"))

  w <- compare_weights(da_fit, ad_fit, sw)
  expect_true(all(is.finite(w$da_reference)))
  expect_equal(w$swada_reference, w$swada_nonref, tolerance = 1e-10)
})

test_that("constructions reject CAMS-only marginal arguments", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  expect_error(modifiers(da_fit, at = 0.5), "CAMS")
  expect_error(shrinkage(da_fit), "CAMS")
})