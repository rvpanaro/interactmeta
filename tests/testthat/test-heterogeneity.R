test_that("heterogeneity() extracts trial SD parameters", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger")

  h <- heterogeneity(fit)
  expect_s3_class(h, "interactmeta_heterogeneity")
  expect_equal(nrow(h$summary), 2L)
  expect_true("Combined" %in% h$summary$parameter)
  expect_true("Interaction_Older" %in% h$summary$parameter)
  expect_true(all(h$draws$value >= 0, na.rm = TRUE))

  h_int <- heterogeneity(fit, component = "intercept")
  expect_equal(nrow(h_int$summary), 1L)
  expect_equal(h_int$summary$parameter, "Combined")

  tau_med <- as.numeric(fit$bayesmeta$overall$qposterior(tau.p = 0.5))
  expect_equal(
    h$summary$estimate[h$summary$parameter == "Combined"],
    tau_med,
    tolerance = 1e-6
  )
  tau_ci <- fit$bayesmeta$overall$post.interval(tau.level = 0.95, method = "shortest")
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
    "Interaction heterogeneity (\u03c4_\u03b3)"
  )

  foot <- heterogeneity_footnote(h, digits = 3)
  expect_match(foot, "Heterogeneity (\u03c4): ", fixed = TRUE)
  expect_false(grepl("% interval", foot, fixed = TRUE))
  expect_match(foot, "\\([0-9.]+, [0-9.]+\\)")
  expect_match(foot, "Interaction heterogeneity (\u03c4_\u03b3):", fixed = TRUE)
  expect_false(grepl("Older - Younger", foot, fixed = TRUE))
  expect_true(grepl(" \u00b7 ", foot, fixed = TRUE))
  expect_false(grepl(" | ", foot, fixed = TRUE))
  foot_med <- heterogeneity_footnote(h, include_ci = FALSE, digits = 3)
  expect_false(grepl("\\([0-9.]+, [0-9.]+\\)", foot_med))
})

test_that("interaction heterogeneity footnote omits the contrast label", {
  expect_equal(
    interact:::comparison_heterogeneity_interaction_label("ifrac2_Older", "Younger"),
    "Interaction heterogeneity (\u03c4_\u03b3)"
  )
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
})

test_that("plot_heterogeneity() uses one bar per discrete tau node", {
  skip_if_not_installed("ggplot2")
  tau <- c(0.05, 0.12, 0.20, 0.35, 0.55, 0.80, 1.10)
  w <- c(0.25, 0.22, 0.18, 0.15, 0.10, 0.07, 0.03)
  # Product-grid repetition like bayesmeta $draws
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
  expect_true(interact:::plot_heterogeneity_is_discrete(het$draws$value))
  bars <- interact:::plot_heterogeneity_hist(
    het$draws$value,
    weight = het$draws$.weight,
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
  fit <- suppressWarnings(cam(
    es,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))
  het <- heterogeneity(fit)
  expect_false(is.null(attr(het, "interactmeta")))
  expect_equal(
    names(interact:::heterogeneity_native_tau_objects(het)),
    c("Combined", "Interaction")
  )
  expect_equal(
    as.character(interact:::heterogeneity_tau_xlab("Combined")),
    as.character(expression("heterogeneity " * tau))
  )
  expect_equal(
    as.character(interact:::heterogeneity_tau_xlab("Interaction")),
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
  d_ov <- interact:::heterogeneity_bayesmeta_tau_prep(objs$Combined)
  d_int <- interact:::heterogeneity_bayesmeta_tau_prep(objs$Interaction)
  expect_equal(
    interact:::heterogeneity_plot_bayesmeta_ylim(list(d_ov, d_int)),
    c(0, max(d_ov$maxdens, d_int$maxdens))
  )
})

test_that("plot_heterogeneity() profiles tau on metafor CAMS", {
  skip_if_not_installed("metafor")
  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- cam(
    es,
    reference = "Younger",
    engine = "metafor",
    silent = 2
  )
  het <- heterogeneity(fit)
  expect_equal(
    names(interact:::heterogeneity_native_rma_objects(het)),
    c("Combined", "Interaction")
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
  ylim <- interact:::heterogeneity_plot_metafor_ylim(list(d_ov, d_int))
  expect_equal(ylim[[1L]], min(d_ov$ylim[[1L]], d_int$ylim[[1L]]))
  expect_equal(ylim[[2L]], max(d_ov$ylim[[2L]], d_int$ylim[[2L]]))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot_heterogeneity(het), het)
  expect_silent(plot_heterogeneity(het, which = "combined"))
  expect_silent(plot_heterogeneity(het, which = "interaction"))
})
