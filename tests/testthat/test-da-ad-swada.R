test_that("da() and ad() match the conventional() bundle numerically", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65")
  ad_fit <- ad(respire_age, reference = "<65")
  conv <- conventional(respire_age, reference = "<65")

  expect_s3_class(da_fit, "interactmeta_da")
  expect_s3_class(ad_fit, "interactmeta_ad")
  expect_equal(da_fit$subgroups$estimate, conv$subgroups$estimate, tolerance = 1e-12)
  expect_equal(da_fit$interaction$estimate, conv$interaction_da$estimate, tolerance = 1e-12)
  expect_equal(da_fit$interaction$se, conv$interaction_da$se, tolerance = 1e-12)
  expect_equal(ad_fit$interaction$estimate, conv$interaction$estimate, tolerance = 1e-12)
  expect_equal(ad_fit$interaction$se, conv$interaction$se, tolerance = 1e-12)
  expect_equal(ad_fit$interaction$tau2, conv$interaction$tau2, tolerance = 1e-12)
})

test_that("DA is the difference of the two subgroup RE-MA means", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65")
  mu_a <- da_fit$subgroups$estimate[da_fit$subgroups$label == "<65"]
  mu_b <- da_fit$subgroups$estimate[da_fit$subgroups$label == ">=65"]
  expect_equal(
    da_fit$interaction$estimate,
    mu_b - mu_a,
    tolerance = 1e-12
  )
  expect_equal(
    da_fit$interaction$se,
    sqrt(sum(da_fit$subgroups$se^2)),
    tolerance = 1e-12
  )
})

test_that("AD is the RE-MA of g_j = y_B - y_A", {
  skip_if_not_installed("metafor")
  ad_fit <- ad(respire_age, reference = "<65")
  g <- ad_fit$pairs$g
  vg <- ad_fit$pairs$vg
  expect_equal(g, ad_fit$pairs$y_nonref - ad_fit$pairs$y_ref, tolerance = 1e-12)
  expect_equal(vg, ad_fit$pairs$v_ref + ad_fit$pairs$v_nonref, tolerance = 1e-12)
  refit <- metafor::rma(yi = g, vi = vg, method = "REML")
  expect_equal(as.numeric(ad_fit$interaction$estimate), as.numeric(refit$beta), tolerance = 1e-12)
  expect_equal(as.numeric(ad_fit$interaction$se), as.numeric(refit$se), tolerance = 1e-12)
})

test_that("native DA and AD can disagree; SWADA makes them identical", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65")
  ad_fit <- ad(respire_age, reference = "<65")
  sw <- swada(respire_age, reference = "<65")

  expect_false(isTRUE(all.equal(
    da_fit$interaction$estimate,
    ad_fit$interaction$estimate,
    tolerance = 1e-8
  )))
  expect_equal(
    sw$interaction$estimate,
    ad_fit$interaction$estimate,
    tolerance = 1e-12
  )
  expect_true(isTRUE(sw$collapsibility$equal))
  expect_equal(
    sw$collapsibility$mu_nonref_minus_mu_ref,
    sw$interaction$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    sw$collapsibility$mu_nonref_minus_mu_ref,
    sw$subgroups$estimate[sw$subgroups$label == ">=65"] -
      sw$subgroups$estimate[sw$subgroups$label == "<65"],
    tolerance = 1e-12
  )
})

test_that("DA, AD, and SWADA use the documented weight vectors", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65")
  ad_fit <- ad(respire_age, reference = "<65")
  sw <- swada(respire_age, reference = "<65")
  w <- compare_weights(da_fit, ad_fit, sw)

  w_ref <- da_fit$weights$weight_ref / sum(da_fit$weights$weight_ref)
  w_non <- da_fit$weights$weight_nonref / sum(da_fit$weights$weight_nonref)
  w_ad <- ad_fit$weights$weight / sum(ad_fit$weights$weight)
  expect_equal(w$da_reference, w_ref, tolerance = 1e-12)
  expect_equal(w$da_nonref, w_non, tolerance = 1e-12)
  expect_equal(w$ad, w_ad, tolerance = 1e-12)
  expect_equal(w$swada_reference, w$ad, tolerance = 1e-10)
  expect_equal(w$swada_nonref, w$ad, tolerance = 1e-10)
  expect_equal(w$swada_reference, w$swada_nonref, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(w$da_reference, w$da_nonref, tolerance = 1e-8)))
})

test_that("interaction sign is non-reference minus reference", {
  skip_if_not_installed("metafor")
  da_age <- da(respire_age, reference = "<65")
  ad_age <- ad(respire_age, reference = "<65")
  expect_equal(da_age$reference, "<65")
  expect_equal(da_age$nonref, ">=65")
  expect_equal(ad_age$pairs$g, ad_age$pairs$y_nonref - ad_age$pairs$y_ref)

  da_sex <- da(respire_sex, reference = "Male")
  ad_sex <- ad(respire_sex, reference = "Male")
  expect_equal(da_sex$nonref, "Female")
  r1 <- subset(respire_sex, study == "RESPIRE 1")
  g_manual <- as.numeric(
    r1$yi[r1$subgroup == "Female"] - r1$yi[r1$subgroup == "Male"]
  )
  expect_equal(ad_sex$pairs$g[ad_sex$pairs$study == "RESPIRE 1"], g_manual, tolerance = 1e-12)
})

test_that("compare helpers label DA, AD, SWADA, CAMS", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  da_fit <- da(respire_age, reference = "<65")
  ad_fit <- ad(respire_age, reference = "<65")
  sw <- swada(respire_age, reference = "<65")
  cam <- suppressWarnings(cams(
    respire_age,
    reference = "<65",
        silent = 2
  ))
  tab <- suppressWarnings(compare_modifiers(da_fit, sw, cam))
  expect_equal(sort(unique(as.character(tab$approach))), c("CAMS", "DA", "SWADA"))
  ix <- suppressWarnings(compare_interactions(da_fit, ad_fit, sw, cam))
  expect_equal(as.character(ix$approach), c("DA", "AD", "SWADA", "CAMS"))
  expect_equal(ix$contrast[ix$approach == "DA"], ">=65 - <65 (DA)")
  expect_equal(
    unique(as.character(ix$contrast[ix$approach != "DA"])),
    ">=65 - <65"
  )
  expect_equal(ix$estimate[ix$approach == "SWADA"], ix$estimate[ix$approach == "AD"], tolerance = 1e-12)
  # historical conventional() three-argument call still works
  conv <- conventional(respire_age, reference = "<65")
  ix2 <- suppressWarnings(compare_interactions(conv, sw, cam))
  expect_equal(as.character(ix2$approach), c("DA", "AD", "SWADA", "CAMS"))
  expect_equal(ix2$estimate, ix$estimate, tolerance = 1e-10)
})

test_that("exponentiated comparison SEs are delta-method, not blank", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  da_fit <- da(respire_age, reference = "<65")
  ad_fit <- ad(respire_age, reference = "<65")
  sw <- swada(respire_age, reference = "<65")
  cam <- suppressWarnings(cams(
    respire_age,
    reference = "<65",
        silent = 2
  ))
  log_ix <- suppressWarnings(compare_interactions(da_fit, ad_fit, sw, cam, exponentiate = FALSE))
  hr_ix <- suppressWarnings(compare_interactions(da_fit, ad_fit, sw, cam, exponentiate = TRUE))
  wald <- as.character(hr_ix$approach) %in% c("DA", "AD", "SWADA")
  expect_true(all(is.finite(hr_ix$se)))
  expect_equal(
    hr_ix$se[wald],
    log_ix$se[wald] * exp(log_ix$estimate[wald]),
    tolerance = 1e-12
  )
})

test_that("compare_interactions warns when CAMS is not metafor REML", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  da_fit <- da(respire_age, reference = "<65")
  ad_fit <- ad(respire_age, reference = "<65")
  sw <- swada(respire_age, reference = "<65")
  cam <- suppressWarnings(cams(
    respire_age,
    reference = "<65",
        silent = 2
  ))
  expect_warning(
    compare_interactions(da_fit, ad_fit, sw, cam),
    "not like-for-like"
  )
  expect_warning(
    compare_modifiers(da_fit, sw, cam),
    "not like-for-like"
  )
})

test_that("compare_interactions is silent when every fit is metafor REML", {
  skip_if_not_installed("metafor")
  da_fit <- da(age, reference = "Younger")
  ad_fit <- ad(age, reference = "Younger")
  sw <- swada(age, reference = "Younger")
  cam <- cams(age, reference = "Younger", method = "REML", silent = 2)
  expect_equal(interactmeta_fit_backend(da_fit), "metafor")
  expect_equal(interactmeta_fit_backend(cam), "metafor")
  expect_equal(cam$method, "REML")
  expect_silent(compare_interactions(da_fit, ad_fit, sw, cam))
  expect_silent(compare_modifiers(da_fit, sw, cam))
  expect_silent(compare_weights(da_fit, ad_fit, sw))
})
