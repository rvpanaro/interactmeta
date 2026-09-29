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
  expect_false(any(cap$engine %in% c("cmdstanr", "brms")))
  expect_true(all(unique(cap$engine[cap$analysis == "cam"]) %in% c("bayesmeta", "metafor")))
  expect_false(any(cap$analysis %in% c("da", "ad") & cap$engine %in% c("cmdstanr", "brms")))
})

test_that("frequentist DA and AD use only their selected metafor fits", {
  skip_if_not_installed("metafor")
  da_fit <- da(respire_age, reference = "<65", engine = "metafor", method = "REML")
  ad_fit <- ad(respire_age, reference = "<65", engine = "metafor", method = "REML")
  expect_equal(da_fit$analysis, "da")
  expect_equal(ad_fit$analysis, "ad")
  expect_equal(da_fit$engine, "metafor")
  expect_equal(ad_fit$engine, "metafor")
  expect_equal(da_fit$inference, "frequentist")
  expect_false(isTRUE(da_fit$capabilities$ifrac))
  expect_false(isTRUE(ad_fit$capabilities$ifrac))
  expect_identical(native_fit(da_fit), da_fit$fits)
  expect_identical(native_fit(ad_fit), ad_fit$fits)
  expect_s3_class(da_fit$fits$reference, "rma")
  expect_s3_class(ad_fit$fits$interaction, "rma")

  int_da <- interactions(da_fit)
  int_ad <- interactions(ad_fit)
  expect_equal(int_da$engine, "metafor")
  expect_equal(int_ad$engine, "metafor")
  expect_equal(int_da$summary$estimate, da_fit$interaction$estimate, tolerance = 1e-12)
  expect_equal(int_ad$summary$estimate, ad_fit$interaction$estimate, tolerance = 1e-12)
  expect_equal(attr(int_da, "interactmeta_origin")$fit_engine, "metafor")
  expect_equal(attr(int_ad, "interactmeta_origin")$fit_engine, "metafor")

  het_da <- heterogeneity(da_fit)
  het_ad <- heterogeneity(ad_fit)
  expect_equal(het_da$engine, "metafor")
  expect_equal(sqrt(max(ad_fit$interaction$tau2, 0)), het_ad$summary$estimate[[1L]], tolerance = 1e-12)
})

test_that("Bayesian DA uses only its selected bayesmeta DA fit", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  da_fit <- da(
    respire_age,
    reference = "<65",
    engine = "bayesmeta"
  )
  expect_equal(da_fit$engine, "bayesmeta")
  expect_equal(da_fit$method, "bayes")
  expect_equal(da_fit$inference, "bayesian")
  expect_true(isTRUE(da_fit$capabilities$posterior))
  expect_false(isTRUE(da_fit$capabilities$ifrac))
  expect_s3_class(da_fit$fits$reference, "bayesmeta")
  expect_s3_class(da_fit$fits$nonref, "bayesmeta")
  expect_true(is.list(da_fit$posterior$interaction))

  gm_ref <- interactmeta:::interactmeta_gm_from_bayesmeta(da_fit$fits$reference)
  gm_non <- interactmeta:::interactmeta_gm_from_bayesmeta(da_fit$fits$nonref)
  gm_da <- interactmeta:::interactmeta_gm_linear_sum(gm_non, 1, gm_ref, -1)
  sm <- interactmeta:::interactmeta_gm_summary(gm_da)
  expect_equal(
    da_fit$interaction$estimate,
    as.numeric(da_fit$posterior$interaction$convolve$quantile(0.5)),
    tolerance = 1e-6
  )
  expect_equal(da_fit$interaction$se, sm$sd, tolerance = 1e-8)
  expect_equal(
    sm$mean,
    interactmeta:::interactmeta_gm_mean(gm_non) - interactmeta:::interactmeta_gm_mean(gm_ref),
    tolerance = 1e-10
  )
  expect_equal(
    sm$sd^2,
    interactmeta:::interactmeta_gm_var(gm_non) + interactmeta:::interactmeta_gm_var(gm_ref),
    tolerance = 1e-8
  )

  int <- interactions(da_fit)
  expect_equal(int$engine, "bayesmeta")
  expect_equal(int$summary$estimate, da_fit$interaction$estimate, tolerance = 1e-12)
  sg <- subgroups(da_fit)
  expect_equal(sg$engine, "bayesmeta")
  het <- heterogeneity(da_fit)
  expect_equal(het$engine, "bayesmeta")
  expect_equal(nrow(het$summary), 2L)

  expect_error(ifrac(da_fit), "Bayesian posterior")
  expect_error(combined(da_fit), "Bayesian posterior")
  expect_error(plot_contribution(da_fit), "Bayesian posterior")
  expect_error(forestplot(da_fit), "DA")
})

test_that("Bayesian AD uses only its selected bayesmeta AD fit", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  ad_fit <- ad(respire_age, reference = "<65", engine = "bayesmeta")
  expect_equal(ad_fit$engine, "bayesmeta")
  expect_equal(ad_fit$inference, "bayesian")
  expect_false(isTRUE(ad_fit$capabilities$ifrac))
  expect_s3_class(ad_fit$fits$interaction, "bayesmeta")
  expect_true(is.list(ad_fit$posterior$interaction))

  loc <- interactmeta:::interactmeta_native_location(
    ad_fit$fits$interaction, prob = 0.95, interval.type = "shortest"
  )
  expect_equal(ad_fit$interaction$estimate, loc$median, tolerance = 1e-8)
  expect_equal(ad_fit$interaction$ci.lb, loc$lower, tolerance = 1e-8)
  expect_equal(ad_fit$interaction$ci.ub, loc$upper, tolerance = 1e-8)
  expect_equal(ad_fit$interaction$tau2, interactmeta:::interactmeta_bayesmeta_tau2_mean(ad_fit$fits$interaction), tolerance = 1e-10)

  int <- interactions(ad_fit)
  expect_equal(int$engine, "bayesmeta")
  expect_equal(int$summary$lci, ad_fit$interaction$ci.lb, tolerance = 1e-12)
  expect_error(subgroups(ad_fit), "not defined")
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
    engine = "bayesmeta"
  )
  ad_fit <- ad(respire_age, reference = "<65", engine = "bayesmeta")
  cam_calls <- 0L
  rma_calls <- 0L
  testthat::with_mocked_bindings(
    {
      invisible(interactions(da_fit))
      invisible(subgroups(da_fit))
      invisible(heterogeneity(da_fit))
      invisible(tidy(da_fit))
      expect_error(ifrac(da_fit), "information-fraction")
      invisible(interactions(ad_fit))
      invisible(heterogeneity(ad_fit))
      expect_error(ifrac(ad_fit), "information-fraction")
    },
    cam = function(...) {
      cam_calls <<- cam_calls + 1L
      stop("must not call CAMS from DA/AD accessors")
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
    .package = "interactmeta"
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
  fit <- suppressWarnings(cam(
    respire_age,
    reference = "<65",
    engine = "bayesmeta",
    silent = 2
  ))
  da_calls <- 0L
  ad_calls <- 0L
  testthat::with_mocked_bindings(
    {
      invisible(interactions(fit))
      invisible(subgroups(fit))
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
    .package = "interactmeta"
  )
  expect_equal(da_calls, 0L)
  expect_equal(ad_calls, 0L)
  expect_true(isTRUE(fit$capabilities$ifrac))
  expect_s3_class(ifrac(fit), "interactmeta_overall_ifrac")
})

test_that("da(approach = 'bayes') uses bayesmeta without naming that engine", {
  skip_if_not_installed("bayesmeta")
  fit <- da(
    respire_age,
    reference = "<65",
    approach = "bayes"
  )
  expect_equal(fit$engine, "bayesmeta")
  expect_equal(fit$method, "bayes")
})

test_that("SWADA rejects Bayesian engines and IF plots still reject Bayesian DA/AD", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  expect_error(swada(respire_age, reference = "<65", engine = "bayesmeta"), "metafor")
  expect_error(swada(respire_age, reference = "<65", approach = "bayes"), "bayesmeta")
  da_fit <- da(
    respire_age,
    reference = "<65",
    engine = "bayesmeta"
  )
  expect_error(plot_heterogeneity(da_fit), "information-fraction")
  expect_error(plot_heterogeneity(heterogeneity(da_fit)), "information-fraction")
  p <- plot(da_fit)
  expect_s3_class(p, "ggplot")
})

test_that("cmdstanr is not a DA backend", {
  expect_error(
    da(respire_age, reference = "<65", engine = "cmdstanr"),
    "cmdstanr"
  )
  expect_error(
    ad(respire_age, reference = "<65", engine = "brms"),
    "CAMS"
  )
})
