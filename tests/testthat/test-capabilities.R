test_that("capability matrix lists only implemented combinations", {
  cap <- interactmeta_capability_matrix()
  expect_true(is.data.frame(cap))
  expect_true(all(c(
    "analysis", "engine", "method", "inference",
    "heterogeneity", "posterior", "ifrac", "ifrac_plots"
  ) %in% names(cap)))
  expect_true(all(cap$ifrac[cap$analysis == "cam"]))
  expect_false(any(cap$ifrac[cap$analysis != "cam"]))
  expect_false(any(cap$ifrac_plots[cap$analysis != "cam"]))
  expect_true(any(cap$analysis == "da" & cap$engine == "bayesmeta" & cap$posterior))
  expect_true(any(cap$analysis == "ad" & cap$engine == "bayesmeta" & cap$posterior))
  expect_false(any(cap$analysis == "swada" & cap$engine == "bayesmeta"))
  expect_true(all(unique(cap$engine[cap$analysis == "cam"]) %in% c("bayesmeta", "metafor")))
})

test_that("frequentist DA and AD use only their selected metafor fits", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  ad_fit <- ad(respire_age, reference = "<65", method = "REML")
  expect_equal(da_fit$analysis, "da")
  expect_equal(ad_fit$analysis, "ad")
  expect_equal(interactmeta_fit_backend(da_fit), "metafor")
  expect_equal(interactmeta_fit_backend(ad_fit), "metafor")
  expect_equal(da_fit$inference, "frequentist")
  expect_false(isTRUE(da_fit$capabilities$ifrac))
  expect_false(isTRUE(ad_fit$capabilities$ifrac))
  expect_identical(native_fit(da_fit), da_fit$fits)
  expect_identical(native_fit(ad_fit), ad_fit$fits)
  expect_s3_class(da_fit$fits$reference, "rma")
  expect_s3_class(ad_fit$fits$interaction, "rma")

  int_da <- interactions(da_fit)
  int_ad <- interactions(ad_fit)
  expect_equal(interactmeta_fit_backend(int_da), "metafor")
  expect_equal(interactmeta_fit_backend(int_ad), "metafor")
  expect_equal(int_da$summary$estimate, da_fit$interaction$estimate, tolerance = 1e-12)
  expect_equal(int_ad$summary$estimate, ad_fit$interaction$estimate, tolerance = 1e-12)
  expect_equal(attr(int_da, "interactmeta_origin")$fit_backend, "metafor")
  expect_equal(attr(int_ad, "interactmeta_origin")$fit_backend, "metafor")

  het_da <- heterogeneity(da_fit)
  het_ad <- heterogeneity(ad_fit)
  expect_equal(interactmeta_fit_backend(het_da), "metafor")
  expect_equal(sqrt(max(ad_fit$interaction$tau2, 0)), het_ad$summary$estimate[[1L]], tolerance = 1e-12)
  ci_ref <- stats::confint(da_fit$fits$reference, level = 0.95)
  ci_non <- stats::confint(da_fit$fits$nonref, level = 0.95)
  ci_int <- stats::confint(ad_fit$fits$interaction, level = 0.95)
  expect_equal(
    het_da$summary$lci[het_da$summary$parameter == "tau[<65]"],
    as.numeric(ci_ref$random["tau", "ci.lb"]),
    tolerance = 1e-6
  )
  expect_equal(
    het_da$summary$uci[het_da$summary$parameter == "tau[>=65]"],
    as.numeric(ci_non$random["tau", "ci.ub"]),
    tolerance = 1e-6
  )
  expect_equal(
    het_ad$summary$lci[het_ad$summary$parameter == "tau[interaction]"],
    as.numeric(ci_int$random["tau", "ci.lb"]),
    tolerance = 1e-6
  )
  skip_if_not_installed("forestplot")
  expect_s3_class(
    forestplot(ad_fit, include_shrinkage = FALSE),
    "interactmeta_forestplot"
  )
  fp_ad_reml <- forestplot(ad_fit)
  expect_s3_class(fp_ad_reml, "interactmeta_forestplot")
  expect_false(any(fp_ad_reml$specs$panels[[1L]]$row_kind == "shrinkage", na.rm = TRUE))
  expect_s3_class(
    forestplot(ad_fit, type = "interaction", include_shrinkage = FALSE),
    "interactmeta_forestplot"
  )
  expect_error(forestplot(ad_fit, type = "modifiers"), "AD forests")
  expect_error(
    forestplot(ad_fit, type = "combined"),
    "no longer supported",
    fixed = TRUE
  )
})

test_that("Bayesian DA uses only its selected bayesmeta DA fit", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  da_fit <- da(
    respire_age,
    reference = "<65",
    method = "Bayesian"
  )
  expect_equal(interactmeta_fit_backend(da_fit), "bayesmeta")
  expect_equal(da_fit$method, "Bayesian")
  expect_equal(da_fit$inference, "bayesian")
  expect_true(isTRUE(da_fit$capabilities$posterior))
  expect_false(isTRUE(da_fit$capabilities$ifrac))
  expect_s3_class(da_fit$fits$reference, "bayesmeta")
  expect_s3_class(da_fit$fits$nonref, "bayesmeta")
  expect_true(is.list(da_fit$posterior$interaction))

  gm_ref <- interact:::interactmeta_gm_from_bayesmeta(da_fit$fits$reference)
  gm_non <- interact:::interactmeta_gm_from_bayesmeta(da_fit$fits$nonref)
  gm_da <- interact:::interactmeta_gm_linear_sum(gm_non, 1, gm_ref, -1)
  sm <- interact:::interactmeta_gm_summary(gm_da)
  expect_equal(
    da_fit$interaction$estimate,
    as.numeric(da_fit$posterior$interaction$convolve$quantile(0.5)),
    tolerance = 1e-6
  )
  expect_equal(da_fit$interaction$se, sm$sd, tolerance = 1e-8)
  expect_equal(
    sm$mean,
    interact:::interactmeta_gm_mean(gm_non) - interact:::interactmeta_gm_mean(gm_ref),
    tolerance = 1e-10
  )
  expect_equal(
    sm$sd^2,
    interact:::interactmeta_gm_var(gm_non) + interact:::interactmeta_gm_var(gm_ref),
    tolerance = 1e-8
  )

  int <- interactions(da_fit)
  expect_equal(interactmeta_fit_backend(int), "bayesmeta")
  expect_equal(int$summary$estimate, da_fit$interaction$estimate, tolerance = 1e-12)
  sg <- modifiers(da_fit)
  expect_equal(interactmeta_fit_backend(sg), "bayesmeta")
  het <- heterogeneity(da_fit)
  expect_equal(interactmeta_fit_backend(het), "bayesmeta")
  expect_equal(nrow(het$summary), 3L)
  expect_true(any(het$summary$parameter == "tau[interaction]"))

  expect_error(ifrac(da_fit), "Bayesian posterior")
  expect_error(combined(da_fit), "Bayesian posterior")
  expect_error(plot_contribution(da_fit), "Bayesian posterior")
  skip_if_not_installed("forestplot")
  fp_da <- forestplot(da_fit, include_shrinkage = FALSE)
  expect_s3_class(fp_da, "interactmeta_forestplot")
})

test_that("Bayesian AD uses only its selected bayesmeta AD fit", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  ad_fit <- ad(respire_age, reference = "<65", method = "Bayesian")
  expect_equal(interactmeta_fit_backend(ad_fit), "bayesmeta")
  expect_equal(ad_fit$inference, "bayesian")
  expect_false(isTRUE(ad_fit$capabilities$ifrac))
  expect_s3_class(ad_fit$fits$interaction, "bayesmeta")
  expect_true(is.list(ad_fit$posterior$interaction))

  loc <- interact:::interactmeta_native_location(
    ad_fit$fits$interaction, prob = 0.95, interval.type = "shortest"
  )
  expect_equal(ad_fit$interaction$estimate, loc$median, tolerance = 1e-8)
  expect_equal(ad_fit$interaction$ci.lb, loc$lower, tolerance = 0.05)
  expect_equal(ad_fit$interaction$ci.ub, loc$upper, tolerance = 0.05)
  expect_equal(ad_fit$interaction$tau2, interact:::interactmeta_bayesmeta_tau2_mean(ad_fit$fits$interaction), tolerance = 1e-10)

  int <- interactions(ad_fit)
  expect_equal(interactmeta_fit_backend(int), "bayesmeta")
  expect_equal(int$summary$lci, ad_fit$interaction$ci.lb, tolerance = 1e-12)
  expect_error(modifiers(ad_fit), "interactions\\(fit\\)")
  expect_error(ifrac(ad_fit), "Bayesian posterior")
  expect_error(combined(ad_fit), "Bayesian posterior")
  expect_error(plot_contribution(ad_fit), "Bayesian posterior")
})

test_that("downstream DA/AD accessors do not silently invoke CAMS or the other construction", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  da_fit <- da(
    respire_age,
    reference = "<65",
    method = "Bayesian"
  )
  ad_fit <- ad(respire_age, reference = "<65", method = "Bayesian")
  cam_calls <- 0L
  rma_calls <- 0L
  testthat::with_mocked_bindings(
    {
      invisible(interactions(da_fit))
      invisible(modifiers(da_fit))
      invisible(heterogeneity(da_fit))
      invisible(tidy(da_fit))
      expect_error(ifrac(da_fit), "information-fraction")
      invisible(interactions(ad_fit))
      invisible(heterogeneity(ad_fit))
      expect_error(ifrac(ad_fit), "information-fraction")
    },
    cams = function(...) {
      cam_calls <<- cam_calls + 1L
      stop("must not call CAMS from DA/AD accessors")
    },
    da = function(...) {
      stop("must not refit DA")
    },
    ad = function(...) {
      stop("must not refit AD")
    },
    .package = "interact"
  )
  testthat::with_mocked_bindings(
    {
      invisible(interactions(da_fit))
      invisible(interactions(ad_fit))
    },
    rma = function(...) {
      rma_calls <<- rma_calls + 1L
      stop("metafor::rma must not run downstream of Bayesian DA/AD")
    },
    .package = "metafor"
  )
  expect_equal(cam_calls, 0L)
  expect_equal(rma_calls, 0L)
})

test_that("CAMS accessors do not silently invoke DA or AD", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  fit <- suppressWarnings(cams(
    respire_age,
    reference = "<65",
        silent = 2
  ))
  da_calls <- 0L
  ad_calls <- 0L
  testthat::with_mocked_bindings(
    {
      invisible(interactions(fit))
      invisible(modifiers(fit))
      invisible(combined(fit))
      invisible(ifrac(fit))
      invisible(heterogeneity(fit))
    },
    da = function(...) {
      da_calls <<- da_calls + 1L
      stop("must not call da() from CAMS")
    },
    ad = function(...) {
      ad_calls <<- ad_calls + 1L
      stop("must not call ad() from CAMS")
    },
    .package = "interact"
  )
  expect_equal(da_calls, 0L)
  expect_equal(ad_calls, 0L)
  cap <- interactmeta_capability_matrix()
  cam_cap <- cap[cap$analysis == "cam" & cap$engine == interactmeta_fit_backend(fit), , drop = FALSE]
  expect_true(any(cam_cap$ifrac))
  expect_s3_class(ifrac(fit), "interactmeta_overall_ifrac")
})

test_that("da(approach = 'bayes') uses bayesmeta without naming that engine", {
  skip_if_not_installed("bayesmeta")
  fit <- da(
    respire_age,
    reference = "<65",
    method = "Bayesian"
  )
  expect_equal(interactmeta_fit_backend(fit), "bayesmeta")
  expect_equal(fit$method, "Bayesian")
})

test_that("SWADA rejects Bayesian engines and IF plots still reject Bayesian DA/AD", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  fit_sw <- swada(respire_age, reference = "<65")
  expect_equal(interactmeta_fit_backend(fit_sw), "metafor")
  expect_error(
    swada(respire_age, reference = "<65", method = "Bayesian"),
    "frequentist inference only"
  )
  da_fit <- da(
    respire_age,
    reference = "<65",
    method = "Bayesian"
  )
  expect_error(plot_heterogeneity(da_fit), "information-fraction")
  expect_error(plot_heterogeneity(heterogeneity(da_fit)), "information-fraction")
  p <- plot(da_fit)
  expect_s3_class(p, "ggplot")
})

test_that("DA/AD/SWADA reject CAMS-only marginal arguments", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65", method = "REML")
  ad_fit <- ad(respire_age, reference = "<65", method = "REML")
  sw_fit <- swada(respire_age, reference = "<65")

  expect_error(
    modifiers(da_fit, type = "Prediction"),
    "type = \"Prediction\""
  )
  expect_error(
    interactions(ad_fit, type = "Prediction"),
    "type = \"Prediction\""
  )
  expect_error(modifiers(da_fit, shrinkage = TRUE), "shrinkage = TRUE")
  expect_error(interactions(sw_fit, at = 0.5), "does not accept `at`")
  expect_error(modifiers(da_fit, ifrac = "overall"), "ifrac = \"overall\"")
  expect_error(modifiers(ad_fit, type = "Overall"), "interactions\\(fit\\)")
})
