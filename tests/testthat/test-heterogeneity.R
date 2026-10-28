test_that("heterogeneity() extracts trial SD parameters", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cams(
    dat,
        reference = "Younger")

  h <- heterogeneity(fit)
  expect_s3_class(h, "interactmeta_heterogeneity")
  expect_equal(nrow(h$summary), 2L)
  expect_true("Combined" %in% h$summary$parameter)
  expect_true("Interaction_Older" %in% h$summary$parameter)
  expect_true(all(h$support$value >= 0, na.rm = TRUE))

  h_int <- heterogeneity(fit, component = "intercept")
  expect_equal(nrow(h_int$summary), 1L)
  expect_equal(h_int$summary$parameter, "Combined")

  tau_med <- as.numeric(fit$bayesmeta$overall$qposterior(tau.p = 0.5))
  expect_equal(
    h$summary$estimate[h$summary$parameter == "Combined"],
    tau_med,
    tolerance = 1e-6
  )
  tau_ci <- fit$bayesmeta$overall$post.interval(
    tau.level = 0.95,
    method = fit$interval.type
  )
  expect_equal(
    unname(c(
      h$summary$lci[h$summary$parameter == "Combined"],
      h$summary$uci[h$summary$parameter == "Combined"]
    )),
    unname(as.numeric(tau_ci)),
    tolerance = 1e-6
  )

  expect_equal(
    h$summary$label[h$summary$type == "slope"],
    paste0("Interaction heterogeneity (", "\u03c4", "_", "\u03b3", ")")
  )

  foot <- heterogeneity_footnote(h, digits = 3)
  expect_match(
    foot,
    paste0("Heterogeneity (", "\u03c4", "): "),
    fixed = TRUE
  )
  expect_false(grepl("% interval", foot, fixed = TRUE))
  expect_match(foot, "\\([0-9.]+, [0-9.]+\\)")
  expect_match(
    foot,
    paste0("Interaction heterogeneity (", "\u03c4", "_", "\u03b3", "):"),
    fixed = TRUE
  )
  expect_false(grepl("Older - Younger", foot, fixed = TRUE))
  expect_true(grepl(" \u00b7 ", foot, fixed = TRUE))
  expect_false(grepl(" | ", foot, fixed = TRUE))
  foot_med <- heterogeneity_footnote(h, include_ci = FALSE, digits = 3)
  expect_false(grepl("\\([0-9.]+, [0-9.]+\\)", foot_med))
})

test_that("interaction heterogeneity footnote omits the contrast label", {
  expect_equal(
    interact:::comparison_heterogeneity_interaction_label("ifrac2_Older", "Younger"),
    paste0("Interaction heterogeneity (", "\u03c4", "_", "\u03b3", ")")
  )
})

test_that("heterogeneity_bayesmeta_tau_prep includes prior density when requested", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  d <- interact:::heterogeneity_bayesmeta_tau_prep(
    fit$bayesmeta$overall,
    prior = TRUE
  )
  expect_false(is.null(d$pdens))
  expect_equal(length(d$pdens), length(d$tau))
  expect_true(any(is.finite(d$pdens)))
})

test_that("plot_heterogeneity() builds a side-by-side ggplot", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  het <- list(
    draws = data.frame(
      parameter = c(rep("Combined", 400), rep("Interaction_Older", 400)),
      type = c(rep("intercept", 400), rep("slope", 400)),
      value = c(abs(stats::rnorm(400, 0, 0.5)), abs(stats::rnorm(400, 0, 0.1))),
      stringsAsFactors = FALSE
    ),
    summary = data.frame()
  )
  class(het) <- c("interactmeta_heterogeneity", "list")

  p <- plot_heterogeneity(het, printplot = FALSE)
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1))))
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1))))
  expect_true(inherits(p$facet, "FacetWrap"))
  expect_equal(p$facet$params$strip.position, "bottom")
  expect_null(p$labels$subtitle)

  p_ov <- plot_heterogeneity(het, which = "overall", printplot = FALSE)
  expect_s3_class(p_ov, "ggplot")
  p_comb <- plot_heterogeneity(het, which = "combined", printplot = FALSE)
  expect_s3_class(p_comb, "ggplot")
  p_int <- plot_heterogeneity(het, which = "interaction", printplot = FALSE)
  expect_s3_class(p_int, "ggplot")
  expect_equal(
    as.character(p_ov$labels$x),
    as.character(expression("heterogeneity " * tau))
  )
  expect_equal(
    as.character(p_int$labels$x),
    as.character(expression("heterogeneity " * tau[gamma]))
  )
  expect_equal(
    p$coordinates$limits$y[[2]],
    max(p_ov$coordinates$limits$y[[2]], p_int$coordinates$limits$y[[2]])
  )

  p2 <- plot(het, printplot = FALSE)
  expect_s3_class(p2, "ggplot")

  expect_error(plot_heterogeneity(list(), printplot = FALSE), "interactmeta_heterogeneity")
  expect_error(
    plot_heterogeneity(het, probs = c(0, 1), printplot = FALSE),
    "probs"
  )

  fit <- fixture_cam_tiny_fit()
  het_fit <- heterogeneity(fit)
  pieces <- interact:::plot_heterogeneity_pieces(het_fit, which = "both")
  prior_df <- interact:::plot_heterogeneity_prior_df(
    pieces = pieces,
    fit = fit,
    prior = TRUE,
    x_hi = 1
  )
  expect_gt(nrow(prior_df), 0L)
  het_gg <- het
  attr(het_gg, "interactmeta") <- list(priors = fit$priors, engine = "ggplot-test")
  p_no_prior <- plot_heterogeneity(het_gg, printplot = FALSE, prior = FALSE)
  p_prior <- plot_heterogeneity(het_gg, printplot = FALSE, prior = TRUE)
  n_line <- function(p) {
    sum(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))
  }
  expect_gt(n_line(p_prior), n_line(p_no_prior))
})

test_that("plot_heterogeneity() uses one bar per discrete tau node", {
  skip_if_not_installed("ggplot2")
  tau <- c(0.05, 0.12, 0.20, 0.35, 0.55, 0.80, 1.10)
  w <- c(0.25, 0.22, 0.18, 0.15, 0.10, 0.07, 0.03)
  # Product-grid repetition like bayesmeta $support
  het <- list(
    draws = data.frame(
      parameter = rep("Combined", length(tau) * 4L),
      type = "intercept",
      value = rep(tau, each = 4L),
      .weight = rep(w / 4, each = 4L),
      stringsAsFactors = FALSE
    ),
    summary = data.frame()
  )
  class(het) <- c("interactmeta_heterogeneity", "list")
  expect_true(interact:::plot_heterogeneity_is_discrete(het$support$value))
  bars <- interact:::plot_heterogeneity_hist(
    het$support$value,
    weight = het$support$.weight,
    bins = 40L
  )
  expect_equal(nrow(bars), length(tau))
  expect_true(all(bars$ymax > 0))
  areas <- (bars$xmax - bars$xmin) * bars$ymax
  expect_equal(sum(areas), 1, tolerance = 1e-8)
})

test_that("plot_heterogeneity() accepts a single MLE draw", {
  skip_if_not_installed("ggplot2")
  het <- list(
    draws = data.frame(
      parameter = c("Combined", "Interaction_Older"),
      type = c("intercept", "slope"),
      value = c(0.22, 0.04),
      stringsAsFactors = FALSE
    ),
    summary = data.frame()
  )
  class(het) <- c("interactmeta_heterogeneity", "list")
  p <- plot_heterogeneity(het, printplot = FALSE)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("plot(het) uses bayesmeta tau plots on the default engine", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")
  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
    reference = "Younger",
        silent = 2
  ))
  het <- heterogeneity(fit)
  expect_false(is.null(attr(het, "interactmeta")))
  expect_equal(
    names(interact:::heterogeneity_native_tau_objects(het)),
    c("tau", "tau_gamma")
  )
  expect_equal(
    as.character(interact:::heterogeneity_tau_xlab("tau")),
    as.character(expression("heterogeneity " * tau))
  )
  expect_equal(
    as.character(interact:::heterogeneity_tau_xlab("tau_gamma")),
    as.character(expression("heterogeneity " * tau[gamma]))
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(het))
  expect_silent(plot(het, which = "interaction", prior = TRUE))
  expect_silent(plot_heterogeneity(het, which = "combined"))
  expect_silent(plot_heterogeneity(het))
  expect_identical(plot_heterogeneity(het, which = "overall"), het)
  objs <- interact:::heterogeneity_native_tau_objects(het)
  d_ov <- interact:::heterogeneity_bayesmeta_tau_prep(objs$tau)
  d_int <- interact:::heterogeneity_bayesmeta_tau_prep(objs$tau_gamma)
  expect_equal(
    interact:::heterogeneity_plot_bayesmeta_ylim(list(d_ov, d_int)),
    c(0, max(d_ov$maxdens, d_int$maxdens))
  )
})

test_that("heterogeneity() attaches metafor profile CIs on tau", {
  skip_if_not_installed("metafor")
  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- cams(
    es,
    reference = "Younger",
    method = "REML",
    silent = 2
  )
  het <- heterogeneity(fit)
  ci_ov <- stats::confint(fit$metafor$overall, level = 0.95)
  ci_int <- stats::confint(fit$metafor$interaction, level = 0.95)
  expect_equal(
    unname(c(
      het$summary$lci[het$summary$parameter == "Combined"],
      het$summary$uci[het$summary$parameter == "Combined"]
    )),
    unname(as.numeric(ci_ov$random["tau", c("ci.lb", "ci.ub")])),
    tolerance = 1e-6
  )
  expect_equal(
    unname(c(
      het$summary$lci[het$summary$parameter == "Interaction_Older"],
      het$summary$uci[het$summary$parameter == "Interaction_Older"]
    )),
    unname(as.numeric(ci_int$random["tau", c("ci.lb", "ci.ub")])),
    tolerance = 1e-6
  )
  foot <- heterogeneity_footnote(het, digits = 3)
  expect_match(foot, "\\([0-9.]+, [0-9.]+\\)")
})

test_that("cutoff shading follows ell(tau) >= cutoff", {
  tau <- c(0, 0.2, 0.5, 0.8, 1)
  ll <- c(-10, -2, 0, -2, -10)
  y_cut <- -3
  polys <- interact:::heterogeneity_metafor_cutoff_polygons(tau, ll, y_cut)
  expect_length(polys, 1L)
  expect_true(all(polys[[1L]]$y >= y_cut - 1e-8))
  expect_true(min(polys[[1L]]$x) <= 0.2 + 1e-6)
  expect_true(max(polys[[1L]]$x) >= 0.8 - 1e-6)
})

test_that("plot_heterogeneity() profiles tau on metafor CAMS", {
  skip_if_not_installed("metafor")
  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- cams(
    es,
    reference = "Younger",
    method = "REML",
    silent = 2
  )
  het <- heterogeneity(fit)
  expect_equal(
    names(interact:::heterogeneity_native_rma_objects(het)),
    c("tau", "tau_gamma")
  )
  d_ov <- interact:::heterogeneity_metafor_tau_prep(
    fit$metafor$overall,
    steps = 8L
  )
  d_int <- interact:::heterogeneity_metafor_tau_prep(
    fit$metafor$interaction,
    steps = 8L
  )
  expect_true(length(d_ov$tau) >= 5L)
  expect_true(length(d_int$ll) >= 5L)
  expect_equal(d_ov$ylab, "restricted log-likelihood")
  expect_equal(
    d_ov$vc,
    sqrt(max(as.numeric(fit$metafor$overall$tau2), 0)),
    tolerance = 1e-8
  )
  expect_equal(d_ov$profile_ci, 0.95)
  expect_length(d_ov$profile_cuts, 1L)
  expect_true(all(d_ov$profile_cuts < d_ov$maxll))
  expect_equal(
    unname(d_ov$profile_cuts[[1L]]),
    d_ov$maxll - 0.5 * stats::qchisq(0.95, df = 1L),
    tolerance = 1e-8
  )
  ylim <- interact:::heterogeneity_plot_metafor_ylim(list(d_ov, d_int))
  expect_equal(ylim[[1L]], min(d_ov$ylim[[1L]], d_int$ylim[[1L]]))
  expect_equal(ylim[[2L]], max(d_ov$ylim[[2L]], d_int$ylim[[2L]]))
  xlim <- interact:::heterogeneity_plot_metafor_xlim(list(d_ov, d_int))
  expect_equal(xlim[[1L]], 0)
  expect_equal(xlim[[2L]], max(d_ov$xlim[[2L]], d_int$xlim[[2L]]))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot_heterogeneity(het), het)
  expect_silent(plot_heterogeneity(het, which = "combined"))
  expect_silent(plot_heterogeneity(het, which = "interaction"))
})
