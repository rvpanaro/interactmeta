edss_es <- function() {
  metafor::escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
}

fit_edss_engine <- function(engine, method = NULL, ...) {
  args <- list(
    x = edss_es(),
    reference = "Lower",
    engine = engine,
    silent = 2
  )
  if (!is.null(method)) {
    args$method <- method
  }
  extra <- list(...)
  if (identical(engine, "brms")) {
    if (is.null(extra$chains)) extra$chains <- 1L
    if (is.null(extra$iter)) extra$iter <- 80L
    if (is.null(extra$warmup)) extra$warmup <- 40L
    if (is.null(extra$refresh)) extra$refresh <- 0
  }
  suppressWarnings(do.call(cam, c(args, extra)))
}

expect_origin_engine <- function(obj, engine) {
  org <- interactmeta_origin(obj)
  expect_false(is.null(org))
  expect_equal(org$fit_engine, engine)
}

test_that("fitting engine A does not call engine B (bayesmeta vs metafor)", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  rma_calls <- 0L
  testthat::with_mocked_bindings(
    {
      fit <- fit_edss_engine("bayesmeta")
      expect_equal(fit$engine, "bayesmeta")
      invisible(interactions(fit))
      invisible(subgroups(fit))
      invisible(combined(fit))
      invisible(ifrac(fit))
      invisible(heterogeneity(fit))
      invisible(summary(fit))
    },
    rma = function(...) {
      rma_calls <<- rma_calls + 1L
      stop("metafor::rma must not run for a bayesmeta fit")
    },
    .package = "metafor"
  )
  expect_equal(rma_calls, 0L)
})

test_that("downstream accessors do not call bayesmeta::bayesmeta on a stored fit", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  fit <- fit_edss_engine("bayesmeta")
  called <- 0L
  testthat::with_mocked_bindings(
    {
      invisible(interactions(fit))
      invisible(subgroups(fit))
      invisible(combined(fit))
      invisible(ifrac(fit))
      invisible(heterogeneity(fit))
      invisible(summary(fit))
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      invisible(plot_contribution(fit))
    },
    bayesmeta = function(...) {
      called <<- called + 1L
      stop("bayesmeta::bayesmeta must not run downstream")
    },
    bmr = function(...) {
      called <<- called + 1L
      stop("bayesmeta::bmr must not run downstream")
    },
    .package = "bayesmeta"
  )
  expect_equal(called, 0L)
})

test_that("EDSS end-to-end: bayesmeta quantities share origin and IF", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  fit <- fit_edss_engine("bayesmeta")
  expect_equal(fit$engine, "bayesmeta")
  expect_true(identical(native_fit(fit), fit$bayesmeta))

  int <- interactions(fit, type = "Pooled")
  sg <- subgroups(fit, type = "Pooled")
  cb <- combined(fit, type = "Pooled")
  het <- heterogeneity(fit)
  oi <- ifrac(fit)
  expect_origin_engine(int, "bayesmeta")
  expect_origin_engine(sg, "bayesmeta")
  expect_origin_engine(cb, "bayesmeta")
  expect_origin_engine(het, "bayesmeta")
  expect_origin_engine(oi, "bayesmeta")

  sm <- summary(fit)
  expect_equal(sm$meta$engine, "bayesmeta")
  pooled <- int$summary[int$summary$type == "Pooled", , drop = FALSE]
  expect_true(is.finite(pooled$lci[[1L]]))
  expect_true(pooled$uci[[1L]] > pooled$lci[[1L]])

  lat <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$gamma)
  expect_equal(pooled$estimate[[1L]], lat$median, tolerance = 0.05)
  expect_equal(pooled$sd[[1L]], lat$sd, tolerance = 0.08)

  trial_pi <- attr(oi, "trial_pi")
  stored <- fit$transformed_data$trial_blocks$overall
  expect_equal(
    unname(trial_pi$pi[match(stored$trial, trial_pi$trial)]),
    unname(stored$ifrac),
    tolerance = 1e-10
  )
})

test_that("EDSS end-to-end: metafor derived intervals come from stored rma", {
  skip_if_not_installed("metafor")
  fit <- fit_edss_engine("metafor", method = "REML")
  expect_equal(fit$engine, "metafor")
  expect_true(identical(native_fit(fit), fit$metafor))

  int <- interactions(fit, type = "Pooled")
  expect_origin_engine(int, "metafor")
  pooled <- int$summary[int$summary$type == "Pooled", , drop = FALSE]
  rma <- fit$metafor$interaction
  expect_equal(pooled$estimate[[1L]], as.numeric(rma$beta)[[1L]], tolerance = 1e-8)
  expect_equal(pooled$lci[[1L]], as.numeric(rma$ci.lb)[[1L]], tolerance = 1e-6)
  expect_equal(pooled$uci[[1L]], as.numeric(rma$ci.ub)[[1L]], tolerance = 1e-6)

  cb <- combined(fit, type = "Pooled")
  expect_true(is.finite(cb$summary$lci[[1L]]))
  sg <- subgroups(fit, type = "Pooled")
  expect_true(all(is.finite(sg$summary$lci)))
  expect_origin_engine(cb, "metafor")
  expect_origin_engine(sg, "metafor")
  expect_origin_engine(ifrac(fit), "metafor")
  expect_origin_engine(heterogeneity(fit), "metafor")
})

test_that("EDSS end-to-end: brms quantities stay on that fit when installed", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")
  skip_if_not_installed("metafor")

  fit <- fit_edss_engine("brms", method = "bayes", chains = 1L, iter = 80L, warmup = 40L)
  expect_equal(fit$engine, "brms")
  expect_origin_engine(interactions(fit), "brms")
  expect_origin_engine(subgroups(fit), "brms")
  expect_origin_engine(combined(fit), "brms")
  expect_origin_engine(heterogeneity(fit), "brms")
  expect_origin_engine(ifrac(fit), "brms")
})

test_that("SWADA has no fabricated IF and IF methods reject it", {
  skip_if_not_installed("metafor")
  sw <- swada(age, reference = "Younger")
  expect_s3_class(sw, "interactmeta_swada")
  expect_null(sw$ifrac)
  expect_false("ifrac" %in% names(sw))
  expect_true(!is.null(sw$interaction))
  expect_true(!is.null(sw$weights))
  expect_false(isTRUE(sw$capabilities$ifrac))
  expect_error(ifrac(sw), "information-fraction")
  expect_error(forestplot(sw), "information-fraction")
  expect_error(plot_heterogeneity(sw), "information-fraction")
  expect_error(plot_contribution(sw), "information-fraction")
  expect_error(combined(sw), "information-fraction")
  expect_s3_class(interactions(sw), "interactmeta_construction_interaction")
  expect_s3_class(subgroups(sw), "interactmeta_construction_subgroups")
  p <- plot(sw)
  expect_s3_class(p, "ggplot")

  da_fit <- da(age, reference = "Younger")
  ad_fit <- ad(age, reference = "Younger")
  expect_error(ifrac(da_fit), "DA")
  expect_error(ifrac(ad_fit), "AD")
  expect_error(forestplot(da_fit), "DA")
  expect_error(combined(da_fit), "DA")
  expect_error(subgroups(ad_fit), "not defined")
})

test_that("comparison helpers only combine explicitly supplied fits", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  da_fit <- da(age, reference = "Younger")
  ad_fit <- ad(age, reference = "Younger")
  sw <- swada(age, reference = "Younger")
  fit <- fit_edss_engine("bayesmeta")
  tab <- suppressWarnings(compare_interactions(da_fit, ad_fit, sw, fit))
  expect_true(nrow(tab) >= 3L)
  expect_true(all(c("DA", "AD", "SWADA") %in% as.character(tab$approach)))
})
