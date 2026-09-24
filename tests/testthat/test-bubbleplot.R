test_that("detect_ifrac_vars() returns all non-ifrac2 columns", {
  dat <- data.frame(ifrac_a = 1, ifrac2_a = 0, ifrac_b = 2)
  expect_equal(interactmeta:::detect_ifrac_vars(dat), c("ifrac_a", "ifrac_b"))
})

test_that("resolve_bubbleplot_settings() uses comparison fit attribute", {
  fit <- structure(
    list(data = data.frame(
      yi = 1,
      vi = 0.01,
      subgroup = "Older",
      trial = "A",
      ifrac_Older = 0.4
    )),
    class = c("interactmeta", "cam", "list")
  )
  attr(fit, "interactmeta") <- list(
    ifrac_var = "ifrac_Older",
    ifrac2_var = "ifrac2_Older",
    trial_var = "trial",
    subgroup_var = "subgroup",
    effects = "ifrac_Older:subgroup",
    ifrac_data = data.frame(
      trial = "A",
      subgroup = "Older",
      yi = 1,
      vi = 0.01,
      ifrac_Older = 0.4,
      precision = 10
    )
  )
  st <- interactmeta:::resolve_bubbleplot_settings(fit)
  expect_equal(st$ifrac_var, "ifrac_Older")
  expect_equal(st$effects, "ifrac_Older:subgroup")
  expect_equal(st$trial_var, "trial")
})

test_that("resolve_bubbleplot_limits() pads xlim; leaves ylim default", {
  lim <- interactmeta:::resolve_bubbleplot_limits(
    ifrac_vals = c(0.3, 0.5),
    y_vals = c(-1, -0.5),
    xlim = NULL,
    ylim = NULL,
    x_breaks = NULL,
    y_breaks = NULL
  )
  expect_true(lim$xlim[1] < 0.3)
  expect_true(lim$xlim[2] > 0.5)
  expect_true(lim$xlim[2] < 1)
  expect_null(lim$ylim)
  expect_null(lim$y_breaks)

  lim2 <- interactmeta:::resolve_bubbleplot_limits(
    ifrac_vals = c(0.3, 0.5),
    y_vals = c(-1, -0.5),
    xlim = NULL,
    ylim = c(-2, 0),
    x_breaks = NULL,
    y_breaks = NULL
  )
  expect_equal(lim2$ylim, c(-2, 0))
  expect_true(length(lim2$y_breaks) >= 2L)
})

test_that("interactmeta_continuous_breaks() keeps ggplot2 default ticks", {
  expect_s3_class(interactmeta:::interactmeta_continuous_breaks(NULL), "waiver")
  expect_equal(interactmeta:::interactmeta_continuous_breaks(c(-1, 0)), c(-1, 0))
})

test_that("bubbleplot_ifrac_grid() follows observed IF, not always 0-1", {
  g <- interactmeta:::bubbleplot_ifrac_grid(c(0.3, 0.5), resolution = 11)
  expect_equal(length(g), 11L)
  expect_gt(min(g), 0.2)
  expect_lt(max(g), 0.6)
  g0 <- interactmeta:::bubbleplot_ifrac_grid(numeric(), resolution = 5)
  expect_equal(g0, seq(0, 1, length.out = 5))
})

test_that("bubbleplot_default_xlab() matches overall_ifrac style", {
  xlab <- interactmeta:::bubbleplot_default_xlab("ifrac_Older")
  oi_xlab <- interactmeta:::overall_ifrac_default_xlab("Older")
  expect_equal(xlab, oi_xlab)
  expect_equal(xlab, "Information fraction (Older)")
})

test_that("default IF xlab does not leak plotmath tildes", {
  xlab <- interactmeta:::overall_ifrac_default_xlab(">=65")
  expect_type(xlab, "character")
  expect_false(grepl("~", xlab, fixed = TRUE))
  expect_equal(xlab, "Information fraction (>=65)")
  expect_equal(
    interactmeta:::bubbleplot_default_xlab("ifrac_>=65"),
    "Information fraction (>=65)"
  )
})

test_that("interactmeta_as_ggplot_label() drops tildes when parse fails", {
  lab <- interactmeta:::interactmeta_as_ggplot_label("Information~fraction~(>=65)")
  expect_type(lab, "character")
  expect_false(grepl("~", lab, fixed = TRUE))
  expect_equal(lab, "Information fraction (>=65)")
  parsed <- interactmeta:::interactmeta_as_ggplot_label("Information~fraction")
  expect_true(inherits(parsed, "expression") || is.language(parsed))
})

test_that("bubbleplot_ce_estimates() is safe for empty conditional-effects data", {
  expect_equal(interactmeta:::bubbleplot_ce_estimates(data.frame()), numeric())
  expect_equal(
    interactmeta:::bubbleplot_ce_estimates(data.frame(yi = 1)),
    numeric()
  )
})

test_that("resolve_bubbleplot_level_keys() maps display labels to subgroup levels", {
  raw_df <- data.frame(
    subgroup = c("Younger", "Older", "Younger"),
    display = factor(c("Ref", "Other", "Ref"), levels = c("Ref", "Other")),
    yi = c(0, 0, 0)
  )
  keys <- interactmeta:::resolve_bubbleplot_level_keys(
    raw_df = raw_df,
    display_levels = c("Ref", "Other"),
    subgroup_var = "subgroup"
  )
  expect_equal(keys, c("Younger", "Older"))
})

test_that("resolve_ifrac_plot_style() uses red / blue for subgroups", {
  sg <- interactmeta:::subgroup_colour_map(
    c("Younger", "Older"),
    reference = "Younger",
    palette = interactmeta:::interactmeta_palette_name()
  )
  style <- interactmeta:::resolve_ifrac_plot_style(
    display_levels = c("Younger", "Older"),
    level_keys = c("Younger", "Older"),
    reference = "Younger",
    subgroup_labels = NULL,
    legend_math = NULL,
    colors = NULL,
    fills = NULL
  )
  expect_equal(style$colors[["Younger"]], sg$colour[1])
  expect_equal(style$colors[["Older"]], sg$colour[2])
  expect_equal(style$fill_colors[["Younger"]], sg$fill[1])
  expect_equal(style$fill_colors[["Older"]], sg$fill[2])
  expect_false(identical(style$colors[["Younger"]], style$colors[["Older"]]))
  expect_equal(
    interactmeta:::interactmeta_normalize_hex(style$colors[["Younger"]]),
    "#E41A1C"
  )
  expect_equal(
    interactmeta:::interactmeta_normalize_hex(style$colors[["Older"]]),
    "#377EB8"
  )
  expect_equal(style$legend_title, "Subgroup")
})

test_that("subgroup_palette_colours() skips near-white Greys", {
  skip_if_not_installed("paletteer")
  cols <- interactmeta:::subgroup_palette_colours("RColorBrewer::Greys", 2L)
  expect_length(cols, 2L)
  lum <- function(x) mean(grDevices::col2rgb(interactmeta:::interactmeta_normalize_hex(x)))
  expect_lt(lum(cols[[1]]), 200)
  expect_lt(lum(cols[[2]]), 200)
})

test_that("trial_yi_segments() builds two-point trial lines", {
  raw <- data.frame(
    study = rep(c("A", "B"), each = 2),
    ifrac = c(0.3, 0.3, 0.4, 0.4),
    yi = c(0.1, -0.2, 0.05, -0.1),
    display = rep(c("x", "y"), 2)
  )
  seg <- interactmeta:::trial_yi_segments(raw, study_var = "study", ifrac_var = "ifrac")
  expect_equal(nrow(seg), 2L)
  expect_equal(seg$yi_1, c(0.1, 0.05))
})
