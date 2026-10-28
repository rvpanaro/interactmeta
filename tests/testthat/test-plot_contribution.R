test_that("plot() contribution type returns ggplot", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("ggplot2")
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
  p <- plot(fit, type = "contribution", ifrac_var = "ifrac_Older")
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 4L)
  expect_false(any(vapply(
    p$layers,
    function(ly) inherits(ly$position, "PositionDodge"),
    logical(1)
  )))

  built <- ggplot2::ggplot_build(p)
  point_layer <- which(vapply(
    p$layers,
    function(ly) inherits(ly$geom, "GeomPoint"),
    logical(1)
  ))
  expect_true(length(point_layer) >= 1L)
  xs <- unique(built$data[[point_layer[[1L]]]]$x)
  expect_equal(
    sort(xs),
    sort(100 * unique(fit$data$ifrac_Older)),
    tolerance = 1e-6
  )

  caption <- p$labels$caption
  if (is.language(caption)) {
    caption <- paste(deparse(caption), collapse = "")
  }
  expect_match(caption, "Interaction-weighted information fraction")
  expect_match(caption, "%")
  expect_match(caption, "k = ")
  expect_match(caption, "studies")
  expect_false(grepl("reference:", caption, fixed = TRUE))

  p_no <- plot_contribution(fit, ifrac_var = "ifrac_Older", show_footnote = FALSE)
  expect_null(p_no$labels$caption)

  lty_layer <- vapply(
    p$layers,
    function(ly) "linetype" %in% names(ly$mapping),
    logical(1)
  )
  expect_true(any(lty_layer))
  built_lty <- vapply(
    ggplot2::ggplot_build(p)$plot$scales$scales,
    function(s) identical(s$aesthetics, "linetype"),
    logical(1)
  )
  expect_true(any(built_lty))

  p_no_overall <- plot_contribution(
    fit,
    ifrac_var = "ifrac_Older",
    show_footnote = FALSE,
    show_overall_line = FALSE
  )
  lty_layer_off <- vapply(
    p_no_overall$layers,
    function(ly) "linetype" %in% names(ly$mapping),
    logical(1)
  )
  expect_false(any(lty_layer_off))
})

test_that("attach_observed_effect_interval() exponentiates after the log Wald interval", {
  d <- data.frame(yi = c(-0.5, 0.2), vi = c(0.04, 0.01))
  z <- stats::qnorm(0.975)
  log_int <- interactmeta:::attach_observed_effect_interval(d, prob = 0.95, exponentiate = FALSE)
  expect_equal(log_int$yi, d$yi)
  expect_equal(log_int$yi_lci, d$yi - z * sqrt(d$vi))
  expect_equal(log_int$yi_uci, d$yi + z * sqrt(d$vi))

  ratio <- interactmeta:::attach_observed_effect_interval(d, prob = 0.95, exponentiate = TRUE)
  expect_equal(ratio$yi, exp(d$yi))
  expect_equal(ratio$yi_lci, exp(d$yi - z * sqrt(d$vi)))
  expect_equal(ratio$yi_uci, exp(d$yi + z * sqrt(d$vi)))
})

test_that("plot_contribution(exponentiate = TRUE) exponentiates study points and intervals", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("ggplot2")

  dat <- data.frame(
    study = rep(paste0("T", 1:4), each = 2),
    level = factor(rep(c("Younger", "Older"), 4), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.28, 0.08, 0.22, 0.15, 0.30),
    vi = c(0.04, 0.02, 0.05, 0.03, 0.03, 0.04, 0.02, 0.025)
  )
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    silent = 2
  ))
  p <- plot_contribution(fit, resolution = 11, show_footnote = FALSE, exponentiate = TRUE)
  built <- ggplot2::ggplot_build(p)
  point_i <- which(vapply(
    p$layers,
    function(ly) inherits(ly$geom, "GeomPoint"),
    logical(1)
  ))[[1L]]
  range_i <- which(vapply(
    p$layers,
    function(ly) inherits(ly$geom, "GeomLinerange"),
    logical(1)
  ))[[1L]]
  ribbon_i <- which(vapply(
    p$layers,
    function(ly) inherits(ly$geom, "GeomRibbon"),
    logical(1)
  ))[[1L]]
  pts <- built$data[[point_i]]
  expect_equal(sort(pts$y), sort(exp(dat$yi)), tolerance = 1e-8)
  rng <- built$data[[range_i]]
  expect_true(all(rng$ymin > 0) && all(rng$ymax > 0))
  z <- stats::qnorm(0.975)
  expect_equal(sort(rng$ymin), sort(exp(dat$yi - z * sqrt(dat$vi))), tolerance = 1e-8)
  expect_equal(sort(rng$ymax), sort(exp(dat$yi + z * sqrt(dat$vi))), tolerance = 1e-8)
  rib <- built$data[[ribbon_i]]
  expect_true(all(rib$ymin > 0) && all(rib$ymax > 0))
})
