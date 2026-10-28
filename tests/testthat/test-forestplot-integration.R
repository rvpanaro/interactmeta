test_that("forestplot() fit path supports favours and grouping", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.02, 0.05, 0.03),
    n = 10
  )
  fit <- suppressWarnings(cams(
    dat, reference = "Younger",
    ))

  fp <- forestplot(fit, exponentiate = TRUE)
  expect_s3_class(fp, "interactmeta_forestplot")
  expect_identical(fp$engine, "bayesmeta")
  expect_equal(length(fp$specs$panels), 1L)
  expect_equal(
    sort(unique(fp$specs$row_labels$heading)),
    c("All", "Older", "Younger")
  )
  expect_false(any(interactmeta:::forestplot_heading_is_interaction(
    fp$specs$row_labels$heading
  )))
  expect_true(any(interactmeta:::forestplot_heading_is_combined(
    fp$specs$row_labels$heading
  )))
  expect_true(any(attr(fp, "panels")[[1]]$row_kind == "combined", na.rm = TRUE))

  fp_no_comb <- forestplot(fit, include_combined = FALSE)
  expect_false(any(interactmeta:::forestplot_heading_is_combined(
    fp_no_comb$specs$row_labels$heading
  )))

  expect_error(
    forestplot(fit, type = "combined"),
    "no longer supported",
    fixed = TRUE
  )

  fp_sg_only <- forestplot(fit, type = "interaction")
  expect_equal(unique(fp_sg_only$specs$row_labels$heading), "Interaction")

  fp_tr <- forestplot(fit, group = "trial")
  expect_false("Interaction" %in% unique(fp_tr$specs$row_labels$heading))
  expect_true("A" %in% unique(fp_tr$specs$row_labels$heading))
  # Subgroup pooled + prediction use parenthetical labels; alone Younger/Older
  # headers are omitted for those summary-only blocks.
  lay_tr <- interactmeta:::forestplot_metafor_subgroup_layout(
    fp_tr$specs$panels[[1]],
    fp_tr$specs$row_labels,
    heterogeneity = fp_tr$specs$annotations$heterogeneity,
    ifrac = fp_tr$specs$annotations$ifrac
  )
  expect_false(any(c("Younger", "Older") %in% lay_tr$headers$label))
  poly_labs <- vapply(lay_tr$polys, `[[`, character(1), "mlab")
  expect_true(any(grepl("Overall \\(Younger\\)", poly_labs)))
  expect_true(any(grepl("Prediction \\(Younger\\)", poly_labs)))
  expect_true(any(grepl("Overall \\(Older\\)", poly_labs)))
  expect_true(any(grepl("Prediction \\(Older\\)", poly_labs)))

  # Overall sits on the next row after the last study in that subgroup block
  fp_sg <- forestplot(fit, include_shrinkage = FALSE)
  lay_sg <- interactmeta:::forestplot_metafor_subgroup_layout(
    fp_sg$specs$panels[[1]],
    fp_sg$specs$row_labels,
    heterogeneity = fp_sg$specs$annotations$heterogeneity,
    ifrac = fp_sg$specs$annotations$ifrac
  )
  expect_true(all(c("Younger", "Older") %in% lay_sg$headers$label))
  poly_labs_sg <- vapply(lay_sg$polys, `[[`, character(1), "mlab")
  for (lab in c("Overall (Younger)", "Overall (Older)")) {
    hdr <- sub("^Overall \\((.*)\\)$", "\\1", lab)
    hdr_row <- lay_sg$headers$row[match(hdr, lay_sg$headers$label)]
    next_hdrs <- lay_sg$headers$row[lay_sg$headers$row < hdr_row]
    block_lo <- if (length(next_hdrs)) max(next_hdrs) else lay_sg$ylim_lo
    in_block <- lay_sg$study$rows[
      lay_sg$study$rows < hdr_row & lay_sg$study$rows > block_lo
    ]
    pooled_y <- lay_sg$polys[[which(poly_labs_sg == lab)[[1L]]]]$row
    expect_true(length(in_block) >= 1L)
    expect_equal(min(in_block) - pooled_y, 1)
  }

  fp_fav <- forestplot(fit)
  expect_equal(attr(fp_fav, "has_favours"), TRUE)
  expect_equal(
    fp_fav$specs$favours,
    c("Control", "   Experimental")
  )
  # Favours are drawn under the axis (xlab), not reserved as in-plot rows.
  lay_fav <- interactmeta:::forestplot_metafor_subgroup_layout(
    fp_fav$specs$panels[[1]],
    fp_fav$specs$row_labels,
    heterogeneity = fp_fav$specs$annotations$heterogeneity,
    ifrac = fp_fav$specs$annotations$ifrac,
    favours = NULL
  )
  expect_equal(length(lay_fav$favours_rows), 0L)

  fp_dots <- forestplot(
    fit,
        shade = TRUE,
    efac = 1,
    psize = 1.2
  )
  expect_true(isTRUE(fp_dots$specs$forest_args$shade))
  expect_identical(fp_dots$specs$forest_args$efac, 1)
  expect_identical(fp_dots$specs$forest_args$psize, 1.2)
})

test_that("forest.interactmeta is registered on metafor::forest and forwards ...", {
  slow()
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  expect_true(exists("forest.interactmeta", envir = asNamespace("interactmeta"), inherits = FALSE))
  expect_true(is.function(interactmeta::forest.interactmeta))
  body_chr <- paste(deparse(body(interactmeta::forest.interactmeta)), collapse = " ")
  expect_match(body_chr, "forestplot")

  # Dispatch via metafor's generic once the method is registered
  meth <- tryCatch(
    utils::getS3method("forest", "interactmeta", optional = TRUE),
    error = function(e) NULL
  )
  if (!is.null(meth)) {
    expect_identical(body(meth), body(interactmeta::forest.interactmeta))
  }

  # Presentation args in ... still land in forest_args (merge path)
  out <- interactmeta:::forestplot_merge_metafor_args(
    list(x = 1, rows = 1, efac = 0, cex = 0.8),
    list(shade = TRUE, efac = 1, textpos = c(-2, 2), rows = 99L),
    which = "forest"
  )
  expect_true(isTRUE(out$shade))
  expect_identical(out$efac, 1)
  expect_equal(out$textpos, c(-2, 2))
  expect_identical(out$rows, 1) # reserved: user cannot override
})

test_that("forestplot_merge_metafor_args() keeps reserved layout args", {
  slow()
  base <- list(x = 1:2, rows = 3:4, cex = 0.75, efac = 0)
  out <- interactmeta:::forestplot_merge_metafor_args(
    base,
    list(cex = 1, shade = TRUE, rows = 99L, x = 0),
    which = "forest"
  )
  expect_identical(out$x, 1:2)
  expect_identical(out$rows, 3:4)
  expect_identical(out$cex, 1)
  expect_true(isTRUE(out$shade))

  poly <- interactmeta:::forestplot_merge_metafor_args(
    list(x = 0.1, row = 2, mlab = "Overall", col = "black"),
    list(mlab = "nope", col = "red", atransf = exp),
    which = "addpoly"
  )
  expect_identical(poly$mlab, "Overall")
  expect_identical(poly$col, "red")
  expect_identical(poly$atransf, exp)
})

test_that("shrinkage rows omit study name in metafor layout slabs", {
  slow()
  panel <- data.frame(
    yi = c(-0.5, -0.4, -0.2, -0.1),
    ci.lb = c(-0.8, -0.7, -0.5, -0.4),
    ci.ub = c(-0.2, -0.1, 0.1, 0.2),
    vi = 0.04,
    row_kind = c("study", "shrinkage", "study", "shrinkage"),
    colour = "black",
    shape = 15L,
    stringsAsFactors = FALSE
  )
  row_labels <- data.frame(
    heading = c("Younger", "Younger", "Older", "Older"),
    label = c("Trial A", "Trial A", "Trial B", "Trial B"),
    stringsAsFactors = FALSE
  )
  layout <- interactmeta:::forestplot_metafor_subgroup_layout(panel, row_labels)
  expect_equal(
    layout$study$slab,
    c("Trial A", " ", "Trial B", " ")
  )
})

test_that("summary-only layout parks placeholder off diamond rows", {
  slow()
  panel <- data.frame(
    label = c(
      "Overall", "Prediction", "Overall", "Prediction",
      "Overall", "Prediction", "Overall", "Prediction"
    ),
    heading = c(
      "Younger", "Younger", "Older", "Older",
      "Interaction", "Interaction", "All", "All"
    ),
    row_kind = c(rep(c("overall", "prediction"), 3), "combined", "combined"),
    yi = c(-0.78, -0.80, -0.5, -0.6, -0.2, -0.3, -0.7, -0.75),
    ci.lb = c(-1.07, -1.53, -0.8, -1.2, -0.5, -0.9, -0.95, -1.4),
    ci.ub = c(-0.48, -0.06, -0.2, 0.0, 0.1, 0.3, -0.45, -0.1),
    vi = 0.0225,
    colour = "black",
    shape = 15L,
    stringsAsFactors = FALSE
  )
  row_labels <- data.frame(
    label = panel$label,

    heading = panel$heading,
    colour_label = c(rep(NA_character_, 6), "Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  lay <- interactmeta:::forestplot_metafor_subgroup_layout(panel, row_labels)
  poly_rows <- vapply(lay$polys, `[[`, numeric(1), "row")
  expect_true(isTRUE(lay$study$placeholder))
  expect_false(lay$study$rows %in% poly_rows)
  expect_equal(
    sum(lay$headers$label == interactmeta:::forestplot_combined_block_label()),
    1L
  )
})

test_that("forestplot_metafor_textpos() nudges left and right columns", {
  slow()
  xlim <- c(-5, 5)
  expect_equal(interactmeta:::forestplot_metafor_textpos(xlim, right_shift = 0), xlim)
  nudged <- interactmeta:::forestplot_metafor_textpos(xlim, right_shift = 0.12)
  expect_equal(nudged[[1L]], xlim[[1L]])
  expect_lt(nudged[[2L]], xlim[[2L]])
  both <- interactmeta:::forestplot_metafor_textpos(
    xlim,
    right_shift = 0.12,
    left_shift = 0.05
  )
  expect_gt(both[[1L]], xlim[[1L]])
  expect_lt(both[[2L]], xlim[[2L]])
})

test_that("forestplot_metafor_plot_xlim pads for Effect text", {
  slow()
  alim <- c(-1, 1)
  full <- interactmeta:::forestplot_metafor_plot_xlim(
    alim,
    has_ilab = FALSE,
    annotate = TRUE
  )
  expect_lt(full[[1L]], alim[[1L]])
  expect_gt(full[[2L]], alim[[2L]])
  tp <- interactmeta:::forestplot_metafor_textpos(full, right_shift = 0.06)
  expect_lt(tp[[2L]], full[[2L]])
  expect_gt(tp[[2L]], alim[[2L]])
})

test_that("forestplot left columns follow IF / N when present", {
  slow()
  with_cols <- data.frame(
    yi = 0,
    ifrac = "50.0%",
    trial_size = "100",
    stringsAsFactors = FALSE
  )
  expect_equal(
    interactmeta:::forestplot_col_left_names(with_cols),
    c("ifrac", "trial_size")
  )
  expect_equal(
    interactmeta:::forestplot_col_left_headings(with_cols),
    c("IF", "N")
  )
  cleared <- with_cols[, setdiff(names(with_cols), c("ifrac", "trial_size")), , drop = FALSE]
  expect_null(interactmeta:::forestplot_col_left_names(cleared))
  expect_equal(interactmeta:::forestplot_col_left_headings(cleared), character())
})

test_that("forestplot omits N unless n_var is set", {
  slow()
  dat <- data.frame(
    study = c("A", "A"),
    level = c("Younger", "Older"),
    n = c(40, 60),
    stringsAsFactors = FALSE
  )
  expect_null(interactmeta:::forestplot_n_values(dat))
  expect_null(interactmeta:::forestplot_n_values(dat, n_var = FALSE))
  expect_equal(interactmeta:::forestplot_n_values(dat, n_var = "n"), c(40, 60))
  panel <- data.frame(yi = c(0.1, 0.2))
  omitted <- interactmeta:::forestplot_attach_trial_size(panel, dat)
  expect_false("trial_size" %in% names(omitted))
  shown <- interactmeta:::forestplot_attach_trial_size(panel, dat, n_var = "n")
  expect_equal(shown$trial_size, c("100", ""))
})

test_that("estimates_only omits IF and N columns", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
        reference = "Younger"))
  fp_full <- suppressWarnings(forestplot(fit, exponentiate = TRUE))
  fp_est <- suppressWarnings(
    forestplot(fit, exponentiate = TRUE, estimates_only = TRUE)
  )
  expect_equal(
    interactmeta:::forestplot_col_left_names(fp_full$specs$panels[[1L]]),
    "ifrac"
  )
  expect_null(interactmeta:::forestplot_col_left_names(fp_est$specs$panels[[1L]]))
  expect_false("ifrac" %in% names(fp_est$specs$panels[[1L]]))
  expect_false("trial_size" %in% names(fp_est$specs$panels[[1L]]))
  expect_true(all(fp_est$specs$panels[[1L]]$row_kind != "study", na.rm = TRUE))
})

test_that("estimates_only pooled/prediction forest prints", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
    method = "REML",
    reference = "Younger"
  ))
  fp <- suppressWarnings(forestplot(
    fit,
    exponentiate = TRUE,
    estimates_only = TRUE,
    include_shrinkage = FALSE,
    include_combined = FALSE
  ))
  expect_true(all(
    fp$specs$panels[[1L]]$row_kind %in% c("overall", "prediction")
  ))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(print(fp), NA)
})

test_that("forestplot chrome stays stable across type", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cams(
    es,
        reference = "Younger"))

  fps <- lapply(
    c("subgroups", "interaction"),
    function(tp) {
      suppressWarnings(forestplot(fit, type = tp, exponentiate = TRUE))
    }
  )
  base_sizes <- vapply(fps, function(fp) as.numeric(attr(fp, "base_size")), numeric(1))
  expect_true(length(unique(base_sizes)) == 1L)
  expect_equal(base_sizes[[1L]], as.numeric(interactmeta:::interactmeta_forest_base_size()))

  favours <- lapply(fps, function(fp) fp$specs$favours)
  expect_equal(favours[[1L]], c("Control", "   Experimental"))
  expect_equal(favours[[2L]], c("Younger", "   Older"))

  headings <- lapply(fps, function(fp) {
    lh <- fp$specs$left_headings
    if (is.null(lh)) character() else as.character(lh)
  })
  expect_equal(headings[[1L]][[1L]], "IF")
  expect_equal(headings[[2L]][[1L]], "IF (Older)")

  # Auto typography stays off unless requested.
  expect_false(isTRUE(attr(fps[[1L]], "typography_auto")))
})

test_that("forestplot_panel_xlab() defaults by role", {
  slow()
  expect_equal(
    interactmeta:::forestplot_panel_xlab("subgroups", exponentiate = FALSE),
    "Effect (log scale)"
  )
  expect_equal(
    interactmeta:::forestplot_panel_xlab("interaction", exponentiate = TRUE),
    "Interaction"
  )
})

test_that("forest footnotes convert Greek to plotmath", {
  slow()
  line <- paste(
    "Information fraction (\u03c0): 0.40 (0.40, 0.40)",
    paste0("Heterogeneity (", "\u03c4", "): 0.30 (0.13, 0.66)"),
    paste0(
      "Interaction heterogeneity (",
      "\u03c4",
      "_",
      "\u03b3",
      "): 0.10 (0.00, 0.26)"
    ),
    sep = " \u00b7 "
  )
  expr <- interactmeta:::forestplot_footnote_as_plotmath(line)
  d <- paste(deparse(expr), collapse = "")
  expect_true(is.language(expr))
  expect_match(d, "pi", fixed = TRUE)
  expect_match(d, "tau[gamma]", fixed = TRUE)
  expect_match(d, "hat(tau)", fixed = TRUE)
  expect_match(d, "tau", fixed = TRUE)
  expect_false(grepl("\u03c4_y", d, fixed = TRUE))
  expect_false(grepl("_", d, fixed = TRUE))
})


# --- from test-forestplot_bayesmeta.R ---

test_that("as_bayesmeta() builds an object forestplot.bayesmeta accepts", {
  slow()
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
  slow()
  all_sum <- c(TRUE, TRUE, TRUE)
  mean_v <- c(NA_real_, NA_real_, 0.2)
  out <- interactmeta:::forestplot_ensure_nonsummary_row(all_sum, mean_v)
  expect_equal(out, c(TRUE, FALSE, TRUE))
  already <- c(TRUE, FALSE, TRUE)
  expect_equal(
    interactmeta:::forestplot_ensure_nonsummary_row(already, mean_v),
    already
  )
  header_only <- interactmeta:::forestplot_ensure_nonsummary_row(
    TRUE, NA_real_
  )
  expect_equal(header_only, FALSE)
})

test_that("forestplot_bayesmeta_blocks() pairs study and shrinkage rows", {
  slow()
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
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
  expect_length(blocks, 1L)
  expect_equal(nrow(blocks[[1]]$studies), 2L)
  expect_equal(blocks[[1]]$studies$label, c("S1", "S2"))
  expect_equal(nrow(blocks[[1]]$shrinkage), 2L)
  expect_true(!is.null(blocks[[1]]$pooled))
  expect_true(!is.null(blocks[[1]]$prediction))
  expect_length(blocks[[1]]$summaries, 2L)
  expect_true(interactmeta:::forestplot_bayesmeta_block_is_simple(blocks[[1]]))
})

test_that("Combined block keeps both Overall and Prediction summaries", {
  slow()
  panel <- data.frame(
    yi = c(0.54, 0.50),
    ci.lb = c(0.41, 0.20),
    ci.ub = c(0.71, 0.90),
    vi = c(0.01, 0.04),
    level = c("Overall", "Prediction"),
    row_kind = c("combined", "combined"),
    stringsAsFactors = FALSE
  )
  labs <- data.frame(
    heading = c("All", "All"),
    label = c("Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
  expect_length(blocks, 1L)
  expect_equal(blocks[[1]]$heading, "All")
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
      ifrac = character(),
      heterogeneity = c(
        overall = paste0("Heterogeneity (", "\u03c4", "): 0.298 (0.10, 0.60)"),
        interaction = paste0(
          "Interaction heterogeneity (",
          "\u03c4",
          "_",
          "\u03b3",
          "): 0.050 (0.01, 0.20)"
        )
      )
    ),
    group = "subgroup"
  )
  args <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # No large title; Effect is a column header.
  expect_null(args$title)
  expect_true("Effect (95% interval)" %in% args$labeltext[1, ])
  expect_false("Estimate" %in% args$labeltext[1, ])
  # Both Combined summaries present in labeltext.
  expect_true("Overall" %in% args$labeltext[, 1])
  expect_true("Prediction" %in% args$labeltext[, 1])

  feet <- interactmeta:::forestplot_bayesmeta_footnote_lines(specs)
  expect_length(feet, 1L)
  expect_false(grepl("Information fraction (\u03c0)", feet[[1L]], fixed = TRUE))
  expect_match(feet[[1L]], "Heterogeneity", fixed = TRUE)
  expect_match(feet[[1L]], "Interaction heterogeneity", fixed = TRUE)
})

test_that("forestplot_bayesmeta_build_fp does not add a graph-column separator", {
  slow()
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
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
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
  args <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  expect_true(is.na(args$zero))
  expect_true(is.finite(args$graph.pos))
  fp <- interactmeta:::forestplot_bayesmeta_build_fp(args)
  extra <- fp$lines$vertical
  v_name <- paste0("v_", args$graph.pos)
  if (is.list(extra)) {
    expect_false(v_name %in% names(extra))
  }
})

test_that("interaction and overall share forestplot::forestplot aesthetics", {
  slow()
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
  blocks_int <- interactmeta:::forestplot_bayesmeta_blocks(panel_int, labs_int)
  expect_true(interactmeta:::forestplot_bayesmeta_block_is_simple(blocks_int[[1]]))

  panel_ov <- data.frame(
    yi = c(0.54, 0.50),
    ci.lb = c(0.41, 0.20),
    ci.ub = c(0.71, 0.90),
    vi = c(0.01, 0.04),
    row_kind = c("combined", "combined"),
    stringsAsFactors = FALSE
  )
  labs_ov <- data.frame(
    heading = c("All", "All"),
    label = c("Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  blocks_ov <- interactmeta:::forestplot_bayesmeta_blocks(panel_ov, labs_ov)
  expect_false(interactmeta:::forestplot_bayesmeta_block_is_simple(blocks_ov[[1]]))

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
  args_int <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks_int, specs)
  args_ov <- interactmeta:::forestplot_bayesmeta_forestplot_args(
    blocks_ov,
    modifyList(specs, list(panels = list(panel_ov), row_labels = labs_ov))
  )
  # Same chrome: Effect column header, no large title; summary colours are
  # baked into fn.ci_sum (col$summary is only a fallback).
  expect_equal(args_int$col$summary, args_ov$col$summary)
  expect_equal(args_int$col$summary, interactmeta:::forestplot_overall_summary_colour())
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
  slow()
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
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
  lim <- interactmeta:::resolve_forestplot_limits(list(panel), exponentiate = TRUE)
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
  args <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  expect_true(isTRUE(args$xlog))
  expect_true(all(args$clip > 0))
  est <- as.numeric(args$lower)
  expect_true(all(is.na(est) | est > 0))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fp <- interactmeta:::forestplot_bayesmeta_build_fp(args)
  expect_error(print(fp, new_page = FALSE), NA)
})

test_that("supplied xlim is passed as forestplot clip", {
  slow()
  skip_if_not_installed("forestplot")
  panel <- data.frame(
    yi = 0.1, ci.lb = -0.1, ci.ub = 0.3, vi = 0.01,
    row_kind = "combined", stringsAsFactors = FALSE
  )
  labs <- data.frame(heading = "All", label = "Overall", stringsAsFactors = FALSE)
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
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
  args <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  expect_equal(args$clip, c(-1, 1))
  expect_equal(args$xticks, c(-1, 0, 1))
})

test_that("group=trial labels Overall/Prediction with subgroup", {
  slow()
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
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
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
  args_tr <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
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
  args_sg <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # Summary-only subgroup blocks also use Overall (Younger) and drop the alone
  # heading; study-backed blocks keep a plain heading.
  expect_true("A" %in% args_sg$labeltext[, 1])
  expect_true("Overall (Younger)" %in% args_sg$labeltext[, 1])
  expect_false("Overall" %in% args_sg$labeltext[, 1])
  expect_equal(sum(args_sg$labeltext[, 1] == "Younger"), 1L)

  expect_equal(
    interactmeta:::forestplot_summary_group_label("Overall", "Younger", TRUE),
    "Overall (Younger)"
  )
  expect_equal(
    interactmeta:::forestplot_summary_group_label("Overall", "All", TRUE, TRUE),
    "Overall (All)"
  )
  expect_equal(
    interactmeta:::forestplot_summary_group_label("Prediction", "All", TRUE, TRUE),
    "Prediction (All)"
  )
  expect_equal(
    interactmeta:::forestplot_summary_group_label("Overall", "All", TRUE, FALSE),
    "Overall"
  )
  expect_equal(
    interactmeta:::forestplot_trial_group_summary_label("Overall", "Younger", "trial"),
    "Overall (Younger)"
  )
  expect_equal(
    interactmeta:::forestplot_trial_group_summary_label("Overall", "Younger", "subgroup"),
    "Overall"
  )
})

test_that("Combined block keeps heading and plain Overall / Prediction labels", {
  slow()
  skip_if_not_installed("forestplot")
  blocks <- list(list(
    heading = "All",
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
  args <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  labs <- args$labeltext[, 1]
  expect_true("All" %in% labs)
  expect_true("Overall" %in% labs)
  expect_true("Prediction" %in% labs)
  expect_false("Overall (All)" %in% labs)
  expect_false("Prediction (All)" %in% labs)
})

test_that("spacer rows use finite zero boxsize", {
  slow()
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
  args <- interactmeta:::forestplot_bayesmeta_forestplot_args(
    interactmeta:::forestplot_bayesmeta_blocks(panel, labs),
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
  slow()
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
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
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
  args <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
  # xlab is the axis title only — favours are not space-padded into it.
  expect_equal(args$xlab, "Effect")
  expect_false(grepl("Older", args$xlab, fixed = TRUE))
  expect_false(isTRUE(args$new_page))

  before <- grDevices::dev.cur()
  fp <- interactmeta:::forestplot_bayesmeta_build_fp(args)
  # Building must not leave a new visible device open.
  expect_equal(grDevices::dev.cur(), before)
  fp <- interactmeta:::forestplot_bayesmeta_apply_favours(fp, specs$favours)
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
  slow()
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
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
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
  args_sg <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
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
  args_tr <- interactmeta:::forestplot_bayesmeta_forestplot_args(blocks, specs)
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
  args_ro <- interactmeta:::forestplot_bayesmeta_forestplot_args(
    interactmeta:::forestplot_bayesmeta_blocks(panel_ro, labs_ro),
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
  args_int <- interactmeta:::forestplot_bayesmeta_forestplot_args(
    interactmeta:::forestplot_bayesmeta_blocks(panel_int, labs_int),
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
    heading = rep("All", 2),
    label = c("A", "Overall"),
    stringsAsFactors = FALSE
  )
  args_ov <- interactmeta:::forestplot_bayesmeta_forestplot_args(
    interactmeta:::forestplot_bayesmeta_blocks(panel_ov, labs_ov),
    specs
  )
  ov_hdr <- which(args_ov$labeltext[, 1] == "All")
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
    heading = c("Younger", "Younger", "Older", "Older", "Interaction", "All"),
    label = c("A", "Overall", "A", "Overall", "Overall", "Overall"),
    stringsAsFactors = FALSE
  )
  args_all <- interactmeta:::forestplot_bayesmeta_forestplot_args(
    interactmeta:::forestplot_bayesmeta_blocks(panel_all, labs_all),
    modifyList(specs, list(group = "subgroup"))
  )
  for (hdr in c("Interaction", "All")) {
    ix <- which(args_all$labeltext[, 1] == hdr)[[1L]]
    expect_equal(args_all$labeltext[ix - 1L, 1], "")
  }

  lay_sg <- interactmeta:::forestplot_metafor_subgroup_layout(
    panel, labs, group = "subgroup"
  )
  lay_tr <- interactmeta:::forestplot_metafor_subgroup_layout(
    panel, labs, group = "trial"
  )
  expect_gt(lay_sg$ylim_hi - lay_sg$ylim_lo, lay_tr$ylim_hi - lay_tr$ylim_lo)
})

test_that("print_with_footnotes clears the device when footnotes are absent", {
  slow()
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
  interactmeta:::forestplot_bayesmeta_print_with_footnotes(fp, c("Combined IF: 0.4"))
  expect_equal(n_draw, 1L)

  n_draw <- 0L
  interactmeta:::forestplot_bayesmeta_print_with_footnotes(fp, character())
  expect_equal(n_draw, 1L)
})
