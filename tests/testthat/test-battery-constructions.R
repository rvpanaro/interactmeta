# DA, AD, SWADA, and conventional() bundle.

test_that("da(), ad(), and conventional() agree on subgroup and interaction estimates", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  ad_fit <- ad(respire_age, reference = "<65", method = "REML")
  conv <- conventional(respire_age, reference = "<65", method = "REML")

  expect_s3_class(da_fit, "interactmeta_da")
  expect_s3_class(ad_fit, "interactmeta_ad")
  expect_equal(da_fit$subgroups$estimate, conv$subgroups$estimate, tolerance = 1e-12)
  expect_equal(ad_fit$interaction$estimate, conv$interaction$estimate, tolerance = 1e-12)
})

test_that("DA interaction is difference of subgroup RE means; AD is RE-MA of contrasts", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  ad_fit <- ad(respire_age, reference = "<65", method = "REML")

  mu_a <- da_fit$subgroups$estimate[da_fit$subgroups$label == "<65"]
  mu_b <- da_fit$subgroups$estimate[da_fit$subgroups$label == ">=65"]
  expect_equal(da_fit$interaction$estimate, mu_b - mu_a, tolerance = 1e-12)

  refit <- metafor::rma(
    yi = ad_fit$pairs$g,
    vi = ad_fit$pairs$vg,
    method = "REML"
  )
  expect_equal(as.numeric(ad_fit$interaction$estimate), as.numeric(refit$beta), tolerance = 1e-12)
})

test_that("SWADA collapses DA and AD interaction estimates", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  ad_fit <- ad(respire_age, reference = "<65", method = "REML")
  sw <- swada(respire_age, reference = "<65", method = "REML")

  expect_false(isTRUE(all.equal(
    da_fit$interaction$estimate,
    ad_fit$interaction$estimate,
    tolerance = 1e-8
  )))
  expect_equal(sw$interaction$estimate, ad_fit$interaction$estimate, tolerance = 1e-12)
  expect_true(isTRUE(sw$collapsibility$equal))
})

test_that("CAMS on respire_age fits with REML and readable modifiers", {
  skip_if_not_installed("metafor")
  fit <- suppressWarnings(cams(respire_age, reference = "<65", method = "REML", silent = 2))
  expect_equal(attr(fit, "interactmeta")$backend, "metafor")
  sg <- modifiers(fit, type = "Overall", summary = TRUE)
  expect_true(all(is.finite(sg$summary$estimate)))
})