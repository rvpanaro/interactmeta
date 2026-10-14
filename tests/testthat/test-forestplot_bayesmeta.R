test_that("as_bayesmeta() builds an object forestplot.bayesmeta accepts", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")

  bma <- as_bayesmeta(
    y = c(-0.2, 0.1, 0.3),
    sigma = c(0.2, 0.15, 0.25),
    labels = c("A", "B", "C"),
    mu = c(median = 0.05, lower = -0.1, upper = 0.2, sd = 0.08),
    prediction = c(median = 0.05, lower = -0.4, upper = 0.5, sd = 0.2),
    shrinkage = data.frame(
      median = c(-0.1, 0.08, 0.2),
      lower = c(-0.3, -0.05, 0.0),
      upper = c(0.1, 0.2, 0.4)
    )
  )
  expect_s3_class(bma, "bayesmeta")
  expect_equal(bma$k, 3L)
  expect_equal(bma$labels, c("A", "B", "C"))

  png(tempfile(fileext = ".png"), width = 7, height = 5, units = "in", res = 72)
  on.exit(dev.off(), add = TRUE)
  expect_silent(
    bayesmeta::forestplot.bayesmeta(bma, plot = TRUE, heterogeneity = FALSE)
  )
})

test_that("forestplot_ensure_nonsummary_row keeps a spacer for v-lines", {
  all_sum <- c(TRUE, TRUE, TRUE)
  mean_v <- c(NA_real_, NA_real_, 0.2)
  out <- interact:::forestplot_ensure_nonsummary_row(all_sum, mean_v)
  expect_equal(out, c(TRUE, FALSE, TRUE))
  already <- c(TRUE, FALSE, TRUE)
  expect_equal(
    interact:::forestplot_ensure_nonsummary_row(already, mean_v),
    already
  )
  header_only <- interact:::forestplot_ensure_nonsummary_row(
    TRUE, NA_real_
  )
  expect_equal(header_only, FALSE)
})

test_that("forestplot_bayesmeta_blocks() pairs study and shrinkage rows", {
  panel <- data.frame(
    yi = c(0.1, 0.05, 0.2, 0.12, 0.15, 0.0),
    ci.lb = c(-0.1, -0.05, 0.0, -0.02, 0.0, -0.3),
    ci.ub = c(0.3, 0.15, 0.4, 0.26, 0.3, 0.3),
    vi = c(0.01, 0.0025, 0.01, 0.0025, 0.005, 0.02),
    row_kind = c("study", "shrinkage", "study", "shrinkage", "overall", "prediction"),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = rep("Younger", 6),
    label = c("S1", "", "S2", "", "mean", "prediction"),
    stringsAsFactors = FALSE
  )
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  expect_length(blocks, 1L)
  expect_equal(nrow(blocks[[1]]$studies), 2L)
  expect_equal(blocks[[1]]$studies$label, c("S1", "S2"))
  expect_equal(nrow(blocks[[1]]$shrinkage), 2L)
  expect_true(!is.null(blocks[[1]]$pooled))
  expect_true(!is.null(blocks[[1]]$prediction))
  expect_length(blocks[[1]]$summaries, 2L)
  expect_true(interact:::forestplot_bayesmeta_block_is_simple(blocks[[1]]))
})

test_that("Combined block keeps both Overall and Prediction summaries", {
  panel <- data.frame(
    yi = c(0.54, 0.50),
    ci.lb = c(0.41, 0.20),
    ci.ub = c(0.71, 0.90),
    vi = c(0.01, 0.04),
    subgroup = c("Overall", "Prediction"),
    row_kind = c("combined", "combined"),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("Combined", "Combined"),
    label = c("Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  expect_length(blocks, 1L)
  expect_equal(blocks[[1]]$heading, "Combined")
  expect_length(blocks[[1]]$summaries, 2L)
  expect_equal(blocks[[1]]$summaries[[1]]$label, "Overall")
  expect_equal(blocks[[1]]$summaries[[2]]$label, "Prediction")
  expect_false(isTRUE(blocks[[1]]$summaries[[1]]$bar))
  expect_true(isTRUE(blocks[[1]]$summaries[[2]]$bar))

  specs <- list(
    panels = list(panel),
    row_labels = labs,
    estimate_col = "yi",
    xlim = list(c(0.1, 2)),
    xticks = list(c(0.25, 0.5, 1, 2)),
    xlab = list("Effect"),
    col_right_heading = list("Effect (95% interval)"),
    base_size = 9,
    pointsize = 4,
    digits = 2,
    exponentiate = TRUE,
    nullval = 1,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = TRUE,
    favours = c("Control", "Experimental"),
    forest_args = list(),
    annotations = list(
      lines = character(),
      ifrac = c(Older = "Information fraction (\u03c0): 0.403 (0.30, 0.50)"),
      heterogeneity = c(
        overall = "Heterogeneity (\u03c4): 0.298 (0.10, 0.60)",
        interaction = "Interaction heterogeneity (\u03c4_\u03b3): 0.050 (0.01, 0.20)"
      )
    ),
    group = "subgroup"
  )
  args <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # No large title; Effect is a column header.
  expect_null(args$title)
  expect_true("Effect (95% interval)" %in% args$labeltext[1, ])
  expect_false("Estimate" %in% args$labeltext[1, ])
  # Both Combined summaries present in labeltext.
  expect_true("Overall" %in% args$labeltext[, 1])
  expect_true("Prediction" %in% args$labeltext[, 1])

  feet <- interact:::forestplot_bayesmeta_footnote_lines(specs)
  expect_length(feet, 1L)
  expect_match(feet[[1L]], "Information fraction (\u03c0)", fixed = TRUE)
  expect_match(feet[[1L]], "Heterogeneity", fixed = TRUE)
  expect_match(feet[[1L]], "Interaction heterogeneity", fixed = TRUE)
})

test_that("forestplot objects get a vertical rule at the graph column", {
  skip_if_not_installed("forestplot")
  panel <- data.frame(
    yi = c(0.2, 0.3),
    ci.lb = c(0.05, 0.1),
    ci.ub = c(0.4, 0.5),
    vi = c(0.01, 0.01),
    row_kind = c("study", "overall"),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("", ""),
    label = c("A", "Overall"),
    stringsAsFactors = FALSE
  )
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  specs <- list(
    panels = list(panel),
    row_labels = labs,
    estimate_col = "yi",
    xlim = list(c(0.07, 0.8)),
    xticks = list(c(0.25, 0.5, 0.75)),
    xlab = list("Information fraction"),
    col_right_heading = list("Information fraction"),
    base_size = 9,
    digits = 2,
    exponentiate = FALSE,
    nullval = 0,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = TRUE,
    favours = FALSE,
    forest_args = list(),
    annotations = list(),
    group = "subgroup"
  )
  args <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  expect_true(is.na(args$zero))
  expect_true(is.finite(args$graph.pos))
  fp <- interact:::forestplot_bayesmeta_build_fp(args)
  extra <- fp$lines$vertical
  expect_true(is.list(extra) && length(extra))
  expect_true(as.character(args$graph.pos) %in% names(extra))
})

test_that("interaction and overall share forestplot::forestplot aesthetics", {
  skip_if_not_installed("forestplot")
  # Interaction-only panel (formerly took the bayesmeta.simple draw path).
  panel_int <- data.frame(
    yi = c(0.1, 0.2, 0.15, 0.18),
    ci.lb = c(-0.1, 0.0, 0.0, -0.1),
    ci.ub = c(0.3, 0.4, 0.3, 0.5),
    vi = c(0.01, 0.01, 0.005, 0.02),
    row_kind = c("study", "study", "overall", "prediction"),
    stringsAsFactors = FALSE
  )
  labs_int <- data.frame(
    heading = rep("Interaction", 4),
    label = c("A", "B", "Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  blocks_int <- interact:::forestplot_bayesmeta_blocks(panel_int, labs_int)
  expect_true(interact:::forestplot_bayesmeta_block_is_simple(blocks_int[[1]]))

  panel_ov <- data.frame(
    yi = c(0.54, 0.50),
    ci.lb = c(0.41, 0.20),
    ci.ub = c(0.71, 0.90),
    vi = c(0.01, 0.04),
    row_kind = c("combined", "combined"),
    stringsAsFactors = FALSE
  )
  labs_ov <- data.frame(
    heading = c("Combined", "Combined"),
    label = c("Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  blocks_ov <- interact:::forestplot_bayesmeta_blocks(panel_ov, labs_ov)
  expect_false(interact:::forestplot_bayesmeta_block_is_simple(blocks_ov[[1]]))

  specs <- list(
    panels = list(panel_int),
    row_labels = labs_int,
    estimate_col = "yi",
    xlim = list(NULL),
    xticks = list(NULL),
    xlab = list("Effect (log scale)"),
    col_right_heading = list("Effect (95% interval)"),
    base_size = 9,
    digits = 2,
    exponentiate = FALSE,
    nullval = 0,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = TRUE,
    favours = NULL,
    forest_args = list(),
    annotations = list(),
    group = "subgroup"
  )
  args_int <- interact:::forestplot_bayesmeta_forestplot_args(blocks_int, specs)
  args_ov <- interact:::forestplot_bayesmeta_forestplot_args(
    blocks_ov,
    modifyList(specs, list(panels = list(panel_ov), row_labels = labs_ov))
  )
  # Same chrome: Effect column header, no large title; summary colours are
  # baked into fn.ci_sum (col$summary is only a fallback).
  expect_equal(args_int$col$summary, args_ov$col$summary)
  expect_equal(args_int$col$summary, interact:::forestplot_overall_summary_colour())
  expect_null(args_int$title)
  expect_null(args_ov$title)
  expect_true("Effect (95% interval)" %in% args_int$labeltext[1, ])
  expect_true("Effect (95% interval)" %in% args_ov$labeltext[1, ])
  expect_true("Overall" %in% args_int$labeltext[, 1])
  expect_true("Prediction" %in% args_int$labeltext[, 1])
  # Interaction vs Combined drawers use different forced colours.
  int_sum <- which(args_int$is.summary & args_int$labeltext[, 1] %in% c("Overall", "Prediction"))
  ov_sum <- which(args_ov$is.summary & args_ov$labeltext[, 1] %in% c("Overall", "Prediction"))
  expect_true(length(int_sum) >= 1L)
  expect_true(length(ov_sum) >= 1L)
  expect_true(is.function(args_int$fn.ci_sum[[int_sum[[1L]]]]))
  expect_true(is.function(args_ov$fn.ci_sum[[ov_sum[[1L]]]]))
  # bayesmeta-compatible axis defaults: no forced clip / xticks.
  expect_null(args_int$clip)
  expect_null(args_ov$clip)
  expect_null(args_int$xticks)
  expect_null(args_ov$xticks)
})

test_that("exponentiated forestplot args drop non-positive CIs and set clip", {
  skip_if_not_installed("forestplot")
  panel <- data.frame(
    yi = c(0.1, 0),
    ci.lb = c(-0.2, -200),
    ci.ub = c(0.4, 0.3),
    vi = c(0.01, 100^2),
    row_kind = c("study", "study"),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("Older", "Older"),
    label = c("A", "B"),
    stringsAsFactors = FALSE
  )
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  lim <- interact:::resolve_forestplot_limits(list(panel), exponentiate = TRUE)
  specs <- list(
    panels = list(panel),
    row_labels = labs,
    estimate_col = "yi",
    xlim = list(lim$xlim),
    xticks = list(lim$xticks),
    xlab = list("Effect"),
    col_right_heading = list("Effect"),
    base_size = 9,
    digits = 2,
    exponentiate = TRUE,
    nullval = 1,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = TRUE,
    favours = NULL,
    forest_args = list(),
    annotations = list(),
    group = "subgroup"
  )
  args <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  expect_true(isTRUE(args$xlog))
  expect_true(all(args$clip > 0))
  est <- as.numeric(args$lower)
  expect_true(all(is.na(est) | est > 0))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fp <- interact:::forestplot_bayesmeta_build_fp(args)
  expect_error(print(fp, new_page = FALSE), NA)
})

test_that("supplied xlim is passed as forestplot clip", {
  skip_if_not_installed("forestplot")
  panel <- data.frame(
    yi = 0.1, ci.lb = -0.1, ci.ub = 0.3, vi = 0.01,
    row_kind = "combined", stringsAsFactors = FALSE
  )
  labs <- data.frame(heading = "Combined", label = "Overall", stringsAsFactors = FALSE)
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  specs <- list(
    panels = list(panel),
    row_labels = labs,
    estimate_col = "yi",
    xlim = list(c(-1, 1)),
    xticks = list(c(-1, 0, 1)),
    xlab = list("x"),
    col_right_heading = list("Effect"),
    base_size = 9,
    digits = 2,
    exponentiate = FALSE,
    nullval = 0,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = FALSE,
    favours = NULL,
    forest_args = list(),
    annotations = list(),
    group = "subgroup"
  )
  args <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  expect_equal(args$clip, c(-1, 1))
  expect_equal(args$xticks, c(-1, 0, 1))
})

test_that("group=trial labels Overall/Prediction with subgroup", {
  skip_if_not_installed("forestplot")
  panel <- data.frame(
    yi = c(0.1, 0.2, 0.15, 0.25, 0.12, 0.18),
    ci.lb = c(-0.1, 0, 0, 0.1, 0, 0),
    ci.ub = c(0.3, 0.4, 0.3, 0.4, 0.3, 0.4),
    vi = 0.01,
    row_kind = c(
      "study", "study",
      "overall", "prediction",
      "overall", "prediction"
    ),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("A", "A", "Younger", "Younger", "Older", "Older"),
    label = c("Younger", "Older", "Overall", "Prediction", "Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  specs <- list(
    panels = list(panel),
    row_labels = labs,
    estimate_col = "yi",
    xlim = list(NULL),
    xticks = list(NULL),
    xlab = list("x"),
    col_right_heading = list("Effect"),
    base_size = 9,
    digits = 2,
    exponentiate = FALSE,
    nullval = 0,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = FALSE,
    favours = NULL,
    forest_args = list(),
    annotations = list(),
    group = "trial"
  )
  args_tr <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  expect_true("Overall (Younger)" %in% args_tr$labeltext[, 1])
  expect_true("Prediction (Younger)" %in% args_tr$labeltext[, 1])
  expect_true("Overall (Older)" %in% args_tr$labeltext[, 1])
  expect_true("Prediction (Older)" %in% args_tr$labeltext[, 1])
  expect_false("Overall" %in% args_tr$labeltext[, 1])
  # Trial heading with studies remains; alone Younger/Older summary headers
  # are omitted (Younger / Older still appear as study labels under A).
  expect_true("A" %in% args_tr$labeltext[, 1])
  expect_equal(sum(args_tr$labeltext[, 1] == "Younger"), 1L)
  expect_equal(sum(args_tr$labeltext[, 1] == "Older"), 1L)

  specs$group <- "subgroup"
  args_sg <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # Summary-only subgroup blocks also use Overall (Younger) and drop the alone
  # heading; study-backed blocks keep a plain heading.
  expect_true("A" %in% args_sg$labeltext[, 1])
  expect_true("Overall (Younger)" %in% args_sg$labeltext[, 1])
  expect_false("Overall" %in% args_sg$labeltext[, 1])
  expect_equal(sum(args_sg$labeltext[, 1] == "Younger"), 1L)

  expect_equal(
    interact:::forestplot_summary_group_label("Overall", "Younger", TRUE),
    "Overall (Younger)"
  )
  expect_equal(
    interact:::forestplot_summary_group_label("Overall", "Combined", TRUE, TRUE),
    "Overall (Combined)"
  )
  expect_equal(
    interact:::forestplot_summary_group_label("Prediction", "Combined", TRUE, TRUE),
    "Prediction (Combined)"
  )
  expect_equal(
    interact:::forestplot_summary_group_label("Overall", "Combined", TRUE, FALSE),
    "Overall"
  )
  expect_equal(
    interact:::forestplot_trial_group_summary_label("Overall", "Younger", "trial"),
    "Overall (Younger)"
  )
  expect_equal(
    interact:::forestplot_trial_group_summary_label("Overall", "Younger", "subgroup"),
    "Overall"
  )
})

test_that("Combined block keeps heading and plain Overall / Prediction labels", {
  skip_if_not_installed("forestplot")
  blocks <- list(list(
    heading = "Combined",
    studies = data.frame(
      yi = numeric(),
      ci.lb = numeric(),
      ci.ub = numeric(),
      label = character(),
      stringsAsFactors = FALSE
    ),
    summaries = list(
      list(median = 0.1, lower = -0.1, upper = 0.3, label = "Overall", bar = FALSE),
      list(median = 0.12, lower = -0.2, upper = 0.4, label = "Prediction", bar = TRUE)
    ),
    pooled = NULL,
    prediction = NULL,
    overall = NULL,
    shrinkage = NULL
  ))
  specs <- list(
    xlab = list("x"),
    col_right_heading = list("Effect"),
    base_size = 9,
    digits = 2,
    exponentiate = FALSE,
    nullval = 0,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = FALSE,
    favours = NULL,
    forest_args = list(),
    annotations = list(),
    group = "subgroup"
  )
  args <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  labs <- args$labeltext[, 1]
  expect_true("Combined" %in% labs)
  expect_true("Overall" %in% labs)
  expect_true("Prediction" %in% labs)
  expect_false("Overall (Combined)" %in% labs)
  expect_false("Prediction (Combined)" %in% labs)
})

test_that("spacer rows use finite zero boxsize", {
  skip_if_not_installed("forestplot")
  panel <- data.frame(
    yi = c(0.1, 0.12, 0.15, 0.2, 0.22, 0.25),
    ci.lb = c(-0.1, -0.1, 0, 0, 0, 0.1),
    ci.ub = c(0.3, 0.3, 0.3, 0.4, 0.4, 0.4),
    vi = 0.01,
    row_kind = c(
      "study", "shrinkage", "overall",
      "study", "shrinkage", "overall"
    ),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("Younger", "Younger", "Younger", "Older", "Older", "Older"),
    label = c("A", "A", "Overall", "A", "A", "Overall"),
    stringsAsFactors = FALSE
  )
  args <- interact:::forestplot_bayesmeta_forestplot_args(
    interact:::forestplot_bayesmeta_blocks(panel, labs),
    list(
      panels = list(panel),
      row_labels = labs,
      estimate_col = "yi",
      xlim = list(NULL),
      xticks = list(NULL),
      xlab = list("x"),
      col_right_heading = list("Effect"),
      base_size = 9,
      digits = 2,
      exponentiate = FALSE,
      nullval = 0,
      show_left = FALSE,
      left_names = character(),
      left_headings = character(),
      show_ci_text = FALSE,
      favours = NULL,
      forest_args = list(),
      annotations = list(),
      group = "subgroup"
    )
  )
  expect_true(is.matrix(args$mean))
  all_na <- apply(args$mean, 1L, function(r) all(is.na(r)))
  expect_true(any(all_na))
  expect_true(all(is.finite(args$boxsize[all_na])))
  expect_true(all(args$boxsize[all_na] == 0))
})

test_that("bayesmeta favours use fp_set_favors under the forest", {
  skip_if_not_installed("forestplot")
  skip_if_not(is.function(get0(
    "fp_set_favors",
    envir = asNamespace("forestplot"),
    inherits = FALSE
  )))
  panel <- data.frame(
    yi = c(0.1, 0.15),
    ci.lb = c(-0.1, 0.0),
    ci.ub = c(0.3, 0.3),
    vi = 0.01,
    row_kind = c("overall", "prediction"),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("Interaction", "Interaction"),
    label = c("Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  specs <- list(
    panels = list(panel),
    row_labels = labs,
    estimate_col = "yi",
    xlim = list(NULL),
    xticks = list(NULL),
    xlab = list("Effect"),
    col_right_heading = list("Effect"),
    base_size = 9,
    digits = 2,
    exponentiate = FALSE,
    nullval = 0,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = FALSE,
    favours = c(
      "Older",
      "Younger"
    ),
    forest_args = list(),
    annotations = list(),
    group = "subgroup"
  )
  args <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # xlab is the axis title only — favours are not space-padded into it.
  expect_equal(args$xlab, "Effect")
  expect_false(grepl("Older", args$xlab, fixed = TRUE))
  expect_false(isTRUE(args$new_page))

  before <- grDevices::dev.cur()
  fp <- interact:::forestplot_bayesmeta_build_fp(args)
  # Building must not leave a new visible device open.
  expect_equal(grDevices::dev.cur(), before)
  fp <- interact:::forestplot_bayesmeta_apply_favours(fp, specs$favours)
  expect_true(is.list(fp$graph_favors))
  expect_equal(
    fp$graph_favors$low,
    "Older"
  )
  expect_equal(
    fp$graph_favors$high,
    "Younger"
  )
  expect_equal(fp$graph_favors$txt_gp$cex, 0.65)
  expect_true(isTRUE(fp$graph_favors$arrows))
})

test_that("group=subgroup inserts a blank row between subgroup blocks", {
  skip_if_not_installed("forestplot")
  panel <- data.frame(
    yi = c(0.1, 0.15, 0.2, 0.25),
    ci.lb = c(-0.1, 0.0, 0.0, 0.1),
    ci.ub = c(0.3, 0.3, 0.4, 0.4),
    vi = 0.01,
    row_kind = c("study", "overall", "study", "overall"),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("Younger", "Younger", "Older", "Older"),
    label = c("A", "Overall", "A", "Overall"),
    stringsAsFactors = FALSE
  )
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  specs <- list(
    panels = list(panel),
    row_labels = labs,
    estimate_col = "yi",
    xlim = list(NULL),
    xticks = list(NULL),
    xlab = list("x"),
    col_right_heading = list("Effect"),
    base_size = 9,
    digits = 2,
    exponentiate = FALSE,
    nullval = 0,
    show_left = FALSE,
    left_names = character(),
    left_headings = character(),
    show_ci_text = FALSE,
    favours = NULL,
    forest_args = list(),
    annotations = list(),
    group = "subgroup"
  )
  args_sg <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # Header + spacer + Younger hdr/study/pooled + between-spacer +
  # Older hdr/study/pooled
  expect_equal(nrow(args_sg$labeltext), 9L)
  expect_equal(args_sg$labeltext[2L, 1], "")
  expect_equal(args_sg$labeltext[3L, 1], "Younger")
  pool_ix <- which(args_sg$labeltext[, 1] == "Overall")
  expect_equal(length(pool_ix), 2L)
  # No blank row immediately before Overall (study label precedes it).
  expect_true(all(args_sg$labeltext[pool_ix - 1L, 1] != ""))
  # Spacer sits between Younger Overall and Older heading.
  yng_pool <- pool_ix[[1L]]
  older_hdr <- which(args_sg$labeltext[, 1] == "Older")[[1L]]
  expect_equal(older_hdr - yng_pool, 2L)
  expect_true(all(is.na(args_sg$mean[yng_pool + 1L])))

  specs$group <- "trial"
  args_tr <- interact:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # Trial nesting with study rows: keep subgroup headers; summaries get
  # Overall (Younger). Blank under column header; no between-subgroup spacers.
  expect_equal(nrow(args_tr$labeltext), 8L)
  expect_equal(args_tr$labeltext[2L, 1], "")
  expect_true("Younger" %in% args_tr$labeltext[, 1])
  pool_tr <- which(grepl("^Overall(\\s*\\(|$)", args_tr$labeltext[, 1]))
  expect_equal(length(pool_tr), 2L)
  expect_true(all(c("Overall (Younger)", "Overall (Older)") %in% args_tr$labeltext[, 1]))

  # Summary-only subgroup blocks: drop alone Younger/Older headers.
  panel_ro <- data.frame(
    yi = c(0.15, 0.2, 0.25, 0.3),
    ci.lb = c(0, 0, 0.1, 0.1),
    ci.ub = c(0.3, 0.4, 0.4, 0.5),
    vi = 0.01,
    row_kind = c("overall", "prediction", "overall", "prediction"),
    stringsAsFactors = FALSE
  )
  labs_ro <- data.frame(
    heading = c("Younger", "Younger", "Older", "Older"),
    label = c("Overall", "Prediction", "Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  args_ro <- interact:::forestplot_bayesmeta_forestplot_args(
    interact:::forestplot_bayesmeta_blocks(panel_ro, labs_ro),
    modifyList(specs, list(group = "subgroup"))
  )
  expect_false("Younger" %in% args_ro$labeltext[, 1])
  expect_false("Older" %in% args_ro$labeltext[, 1])
  expect_true(all(
    c(
      "Overall (Younger)",
      "Prediction (Younger)",
      "Overall (Older)",
      "Prediction (Older)"
    ) %in% args_ro$labeltext[, 1]
  ))

  # Interaction-only: spacer under column header (Interaction leads).
  panel_int <- data.frame(
    yi = c(0.1, 0.15, 0.2),
    ci.lb = c(-0.1, 0.0, 0.0),
    ci.ub = c(0.3, 0.3, 0.4),
    vi = 0.01,
    row_kind = c("study", "overall", "prediction"),
    stringsAsFactors = FALSE
  )
  labs_int <- data.frame(
    heading = rep("Interaction", 3),
    label = c("A", "Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  args_int <- interact:::forestplot_bayesmeta_forestplot_args(
    interact:::forestplot_bayesmeta_blocks(panel_int, labs_int),
    specs
  )
  expect_equal(args_int$labeltext[2L, 1], "")
  int_hdr <- which(args_int$labeltext[, 1] == "Interaction")
  expect_equal(length(int_hdr), 1L)
  expect_equal(int_hdr, 3L)
  pool_int <- which(args_int$labeltext[, 1] == "Overall")
  expect_equal(length(pool_int), 1L)
  expect_equal(args_int$labeltext[pool_int - 1L, 1], "A")

  # Spacer before Combined header as well.
  panel_ov <- data.frame(
    yi = c(0.1, 0.15),
    ci.lb = c(-0.1, 0.0),
    ci.ub = c(0.3, 0.3),
    vi = 0.01,
    row_kind = c("study", "overall"),
    stringsAsFactors = FALSE
  )
  labs_ov <- data.frame(
    heading = rep("Combined", 2),
    label = c("A", "Overall"),
    stringsAsFactors = FALSE
  )
  args_ov <- interact:::forestplot_bayesmeta_forestplot_args(
    interact:::forestplot_bayesmeta_blocks(panel_ov, labs_ov),
    specs
  )
  ov_hdr <- which(args_ov$labeltext[, 1] == "Combined")
  expect_equal(length(ov_hdr), 1L)
  expect_equal(args_ov$labeltext[ov_hdr - 1L, 1], "")

  # Full stack: spacer before Interaction and before Combined after subgroups.
  panel_all <- data.frame(
    yi = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
    ci.lb = c(-0.1, 0.0, 0.0, 0.1, 0.1, 0.15),
    ci.ub = c(0.3, 0.3, 0.4, 0.4, 0.5, 0.55),
    vi = 0.01,
    row_kind = c("study", "overall", "study", "overall", "overall", "overall"),
    stringsAsFactors = FALSE
  )
  labs_all <- data.frame(
    heading = c("Younger", "Younger", "Older", "Older", "Interaction", "Combined"),
    label = c("A", "Overall", "A", "Overall", "Overall", "Overall"),
    stringsAsFactors = FALSE
  )
  args_all <- interact:::forestplot_bayesmeta_forestplot_args(
    interact:::forestplot_bayesmeta_blocks(panel_all, labs_all),
    modifyList(specs, list(group = "subgroup"))
  )
  for (hdr in c("Interaction", "Combined")) {
    ix <- which(args_all$labeltext[, 1] == hdr)[[1L]]
    expect_equal(args_all$labeltext[ix - 1L, 1], "")
  }

  lay_sg <- interact:::forestplot_metafor_subgroup_layout(
    panel, labs, group = "subgroup"
  )
  lay_tr <- interact:::forestplot_metafor_subgroup_layout(
    panel, labs, group = "trial"
  )
  expect_gt(lay_sg$ylim_hi - lay_sg$ylim_lo, lay_tr$ylim_hi - lay_tr$ylim_lo)
})

test_that("print_with_footnotes clears the device when footnotes are absent", {
  skip_if_not_installed("forestplot")

  ns <- asNamespace("forestplot")
  orig <- get("drawForestplotObject", envir = ns)
  n_draw <- 0L
  unlockBinding("drawForestplotObject", ns)
  assign(
    "drawForestplotObject",
    function(obj) {
      n_draw <<- n_draw + 1L
      orig(obj)
    },
    envir = ns
  )
  on.exit(
    {
      assign("drawForestplotObject", orig, envir = ns)
      lockBinding("drawForestplotObject", ns)
    },
    add = TRUE
  )

  pdf(NULL)
  fp <- forestplot::forestplot(
    labeltext = cbind(c("Study", "A"), c("N", "10")),
    mean = c(NA, 0.1),
    lower = c(NA, -0.1),
    upper = c(NA, 0.3),
    is.summary = c(TRUE, FALSE),
    new_page = FALSE,
    xlab = "effect"
  )
  dev.off()

  png(tempfile(fileext = ".png"), width = 600, height = 400)
  on.exit(dev.off(), add = TRUE)

  # Annotated draw, then no-footnote draw: each must call draw once (and the
  # second must start with grid.newpage so it does not ghost the first).
  n_draw <- 0L
  interact:::forestplot_bayesmeta_print_with_footnotes(fp, c("Combined IF: 0.4"))
  expect_equal(n_draw, 1L)

  n_draw <- 0L
  interact:::forestplot_bayesmeta_print_with_footnotes(fp, character())
  expect_equal(n_draw, 1L)
})
