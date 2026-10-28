# Panel builder and ordering (fast tier): shared fits via helper-fixtures.

.panel_fit_bm <- NULL
local_panel_fit_bm <- function() {
  if (is.null(.panel_fit_bm)) {
    .panel_fit_bm <<- fixture_cam_tiny_fit(
            dat = data.frame(
        study = rep(c("A", "B"), each = 2),
        level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
        yi = c(0.1, 0.2, 0.15, 0.25),
        vi = c(0.04, 0.02, 0.05, 0.03),
        n = 10
      )
    )
  }
  .panel_fit_bm
}

test_that("comparison_forest_resolve_model_names() warns on unrecognized layers", {
  expect_warning(
    out <- interact:::comparison_forest_resolve_model_names(
      model_names = c("Nope", "Overall"),
      include_shrinkage = FALSE,
      include_overall = TRUE,
      include_prediction = FALSE
    ),
    "Unrecognized `model_names`"
  )
  expect_equal(out, "Overall")
})

test_that("omitted pooled/prediction layers keep NA placeholder rows", {
  overall_sg <- data.frame(
    study = c("0.5", "0.5"),
    level = c("Younger", "Older"),
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.01),
    stringsAsFactors = FALSE
  )
  pred_sg <- overall_sg
  pred_sg$yi <- c(0.15, 0.25)
  empty_int <- overall_sg[0L, , drop = FALSE]

  both <- interact:::comparison_forest_interleave_overall_prediction(
    overall_sg = overall_sg,
    pred_sg = pred_sg,
    overall_int = empty_int,
    pred_int = empty_int,
    reference = "Younger",
    include_overall = TRUE,
    include_prediction = TRUE
  )
  no_pred <- interact:::comparison_forest_interleave_overall_prediction(
    overall_sg = overall_sg,
    pred_sg = pred_sg,
    overall_int = empty_int,
    pred_int = empty_int,
    reference = "Younger",
    include_overall = TRUE,
    include_prediction = FALSE
  )
  expect_equal(nrow(both), nrow(no_pred))
  expect_equal(no_pred$row_kind, rep(c("overall", "prediction"), 2L))
  expect_true(all(is.finite(no_pred$yi[no_pred$row_kind == "overall"])))
  expect_true(all(is.na(no_pred$yi[no_pred$row_kind == "prediction"])))

  no_overall <- interact:::comparison_forest_interleave_overall_prediction(
    overall_sg = overall_sg,
    pred_sg = pred_sg,
    overall_int = empty_int,
    pred_int = empty_int,
    reference = "Younger",
    include_overall = FALSE,
    include_prediction = TRUE
  )
  expect_equal(nrow(both), nrow(no_overall))
  expect_true(all(is.na(no_overall$yi[no_overall$row_kind == "overall"])))
  expect_true(all(is.finite(no_overall$yi[no_overall$row_kind == "prediction"])))
})

test_that("comparison_forest_order_to_template() orders by trial and subgroup", {
  template <- data.frame(
    study = rep(c("B", "A"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(c("Older", "Younger", "Older", "Younger"), levels = c("Younger", "Older")),
    yi = 1:4,
    vi = 0.1,
    stringsAsFactors = FALSE
  )
  out <- interact:::comparison_forest_order_to_template(dat, template, "trial", "subgroup")
  expect_equal(out$trial, c("B", "B", "A", "A"))
  expect_equal(as.character(out), rep(c("Younger", "Older"), 2))
})

test_that("comparison_forest_order_to_template() matches split_ifrac_for_forest order", {
  d <- ifrac(
    yi = c(0.2, 0.1, 0.3, 0.15),
    vi = c(0.01, 0.04, 0.0225, 0.09),
    study = c("B", "B", "A", "A"),
    level = factor(c("Older", "Younger", "Older", "Younger"), levels = c("Younger", "Older")),
    reference = "Younger"
  )
  template <- interact:::split_ifrac_for_forest(d, contrast_sep = "-")
  sg <- interact:::comparison_forest_subgroup_rows(d, "trial", "subgroup", "-")
  out <- interact:::comparison_forest_order_to_template(sg, template, "trial", "subgroup")
  expect_equal(
    paste(out$trial, as.character(out)),
    paste(template$trial, as.character(template))
  )
})

test_that("CAMS forest panels keep log-scale marginal estimates on the ratio scale", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = log(c(1.1, 1.2, 1.15, 1.25)),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Overall"),
    include_shrinkage = FALSE,
    include_overall = TRUE,
    include_prediction = FALSE,
    include_combined = FALSE
  )
  pooled_yi <- panels$resultsA$yi[panels$resultsA$row_kind == "overall"]
  study_yi <- panels$resultsA$yi[panels$obs_row_index]
  expect_true(all(abs(pooled_yi) < 1))
  expect_true(all(abs(study_yi) < 1))
  limits <- interact:::resolve_forestplot_limits(
    list(panels$resultsA, panels$resultsB),
    exponentiate = TRUE
  )
  expect_true(all(limits$xlim[[1]] > 0))
})

test_that("forestplot_render() keeps log-scale panels when exponentiate is TRUE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  built <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Overall"),
    include_shrinkage = FALSE,
    include_overall = TRUE,
    include_prediction = FALSE,
    include_combined = FALSE
  )
  render_args <- c(
    built,
    list(
      prob = 0.95,
      contrast_sep = "-"
    )
  )
  fp_log <- do.call(interact:::forestplot_render, c(render_args, list(exponentiate = FALSE)))
  fp_ratio <- do.call(interact:::forestplot_render, c(render_args, list(exponentiate = TRUE)))
  expect_equal(attr(fp_log, "panels")[[1]]$yi, attr(fp_ratio, "panels")[[1]]$yi)
  limits <- interact:::resolve_forestplot_limits(attr(fp_ratio, "panels"), exponentiate = TRUE)
  expect_true(all(unlist(limits$xlim) > 0))
})

test_that("comparison_forest_overall_subgroup_block() orders subgroups", {
  summary <- data.frame(
    level = factor(c("Older", "Younger"), levels = c("Younger", "Older")),
    type = "Overall",
    mean = c(0.2, 0.1),
    sd = c(0.05, 0.04),
    stringsAsFactors = FALSE
  )
  block <- interact:::comparison_forest_overall_subgroup_block(
    summary,
    "Overall",
    "0.500",
    subgroup_order = c("Younger", "Older")
  )
  expect_equal(as.character(block), c("Younger", "Older"))
})

test_that("comparison_forest_subgroup_order() matches within-trial subgroup order", {
  d <- ifrac(
    yi = c(0.2, 0.1, 0.3, 0.15),
    vi = c(0.01, 0.04, 0.0225, 0.09),
    study = c("B", "B", "A", "A"),
    level = factor(c("Older", "Younger", "Older", "Younger"), levels = c("Younger", "Older")),
    reference = "Younger"
  )
  expect_equal(interact:::comparison_forest_subgroup_order(d, "-"), c("Younger", "Older"))
})

test_that("interleave_comparison_shrinkage_rows() doubles rows with shrinkage estimates", {
  resultsA <- data.frame(
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.01),
    ci.lb = c(0, 0.1),
    ci.ub = c(0.2, 0.3),
    stringsAsFactors = FALSE
  )
  resultsB <- resultsA
  row_labels <- data.frame(
    subheading = c("A", "A"),
    heading = "",
    label = c("Younger", "Older"),
    stringsAsFactors = FALSE
  )
  shrinkage <- list(
    subgroup = data.frame(
      study = c("A", "A"),
      level = c("Younger", "Older"),
      yi = c(0.12, 0.18),
      vi = c(0.008, 0.009),
      stringsAsFactors = FALSE
    ),
    interaction = NULL
  )
  ifrac_disp <- c("80.0%", "20.0%")
  out <- interact:::interleave_comparison_shrinkage_rows(
    resultsA,
    resultsB,
    row_labels,
    shrinkage,
    prob = 0.95,
    ref = "Younger",
    contrast_sep = "-",
    ifrac_display = ifrac_disp
  )
  expect_equal(nrow(out$resultsA), 5L)
  expect_equal(out$resultsA$ifrac, c("80.0%", "", "20.0%", "", ""))
  expect_equal(out$resultsA$yi, c(0.1, 0.12, 0.2, 0.18, NA))
  expect_equal(
    out$row_labels$label,
    c("Younger", "Younger", "Older", "Older", "Older")
  )
  expect_equal(out$obs_row_index, c(1L, 3L))
  expect_equal(out$interaction_shrinkage_row_index, 5L)
  expect_equal(out$subgroup_shrinkage_row_index, c(2L, 4L))
  expect_true(is.finite(out$resultsB$yi[5]))
  expect_true(is.na(out$resultsB$yi[2]))
  styled <- interact:::comparison_forest_apply_shrinkage_colours(
    out$resultsA,
    out$resultsB,
    out$row_labels,
    out$subgroup_shrinkage_row_index,
    out$interaction_shrinkage_row_index,
    "Younger"
  )
  expect_equal(styled$resultsA$fill[2], interact:::comparison_forest_shrinkage_colour())
  expect_equal(styled$resultsA$fill[4], interact:::comparison_forest_shrinkage_colour())
  expect_equal(styled$resultsA$colour[2], interact:::comparison_forest_shrinkage_colour())
})

test_that("comparison_forest_shrinkage returns fitted values", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")

  shrink <- interact:::comparison_forest_shrinkage_data(fit, prob = 0.95)
  expect_equal(nrow(shrink), 4L)
  expect_true(all(is.finite(shrink$yi)))

  sg <- shrinkage(fit, type = "modifiers", summary = TRUE)$summary
  int <- shrinkage(fit, type = "interaction", summary = TRUE)$summary
  ov <- shrinkage(fit, type = "combined", summary = TRUE)$summary
  key_sg <- paste(shrink$trial, shrink, sep = "\r")
  key_fit <- paste(sg$trial, sg, sep = "\r")
  m <- match(key_sg, key_fit)
  expect_false(any(is.na(m)))
  expect_equal(shrink$yi, sg$estimate[m], tolerance = 1e-10)
  expect_equal(shrink$vi, pmax(sg$sd[m]^2, 1e-12), tolerance = 1e-10)

  key_int <- paste(shrink$interaction$trial, shrink$interaction, sep = "\r")
  key_int_fit <- paste(int$trial, int$interaction, sep = "\r")
  mi <- match(key_int, key_int_fit)
  expect_false(any(is.na(mi)))
  expect_equal(shrink$interaction$yi, int$estimate[mi], tolerance = 1e-10)
  expect_equal(shrink$interaction$vi, pmax(int$sd[mi]^2, 1e-12), tolerance = 1e-10)

  expect_equal(sort(shrink$overall$trial), sort(as.character(ov$trial)))
  mo <- match(shrink$overall$trial, as.character(ov$trial))
  expect_equal(shrink$overall$yi, ov$estimate[mo], tolerance = 1e-10)
  expect_equal(shrink$overall$vi, pmax(ov$sd[mo]^2, 1e-12), tolerance = 1e-10)
})

test_that("interleave_comparison_shrinkage_rows() adds interaction shrinkage once per trial", {
  resultsA <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = rep(c("Younger", "Older"), 2),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    ci.lb = c(0, 0.1, 0.05, 0.2),
    ci.ub = c(0.2, 0.3, 0.25, 0.35),
    stringsAsFactors = FALSE
  )
  resultsB <- resultsA
  row_labels <- data.frame(
    subheading = resultsA$trial,
    heading = "",
    label = resultsA,
    stringsAsFactors = FALSE
  )
  shrinkage <- list(
    subgroup = data.frame(
      trial = resultsA$trial,
      subgroup = resultsA,
      yi = resultsA$yi + 0.01,
      vi = resultsA$vi,
      stringsAsFactors = FALSE
    ),
    interaction = data.frame(
      study = c("A", "B"),
      level = c("Older-Younger", "Older-Younger"),
      yi = c(0.11, 0.12),
      vi = c(0.01, 0.01),
      stringsAsFactors = FALSE
    )
  )
  out <- interact:::interleave_comparison_shrinkage_rows(
    resultsA,
    resultsB,
    row_labels,
    shrinkage,
    prob = 0.95,
    ref = "Younger",
    contrast_sep = "-"
  )
  expect_equal(nrow(out$resultsA), 10L)
  int_shr <- out$resultsB[out$interaction_shrinkage_row_index, "yi"]
  expect_equal(sum(is.finite(int_shr)), 2L)
  expect_equal(int_shr[is.finite(int_shr)], c(0.11, 0.12))
  expect_equal(length(out$interaction_shrinkage_row_index), 2L)
  expect_equal(length(out$subgroup_shrinkage_row_index), 4L)
  # Within each trial: obs, sg-shr, (int-shr for non-ref), …
  expect_equal(
    out$row_labels$label,
    c(
      "Younger", "Younger", "Older", "Older", "Older",
      "Younger", "Younger", "Older", "Older", "Older"
    )
  )
  expect_equal(out$obs_row_index, c(1L, 3L, 6L, 8L))
  expect_equal(out$subgroup_shrinkage_row_index, c(2L, 4L, 7L, 9L))
  expect_equal(out$interaction_shrinkage_row_index, c(5L, 10L))
  expect_equal(
    out$row_labels$label[out$subgroup_shrinkage_row_index],
    c("Younger", "Older", "Younger", "Older")
  )
  expect_equal(
    out$row_labels$label[out$interaction_shrinkage_row_index],
    c("Older", "Older")
  )
  sg_shr <- out$resultsA[out$subgroup_shrinkage_row_index, "yi"]
  expect_equal(sum(is.finite(sg_shr)), 4L)
})

test_that("comparison_forest_resolve_model_names() respects include flags", {
  expect_equal(
    interact:::comparison_forest_resolve_model_names(
      model_names = NULL,
      include_shrinkage = TRUE,
      include_overall = FALSE,
      include_prediction = TRUE
    ),
    c("Shrinkage", "Prediction")
  )
  expect_equal(
    interact:::comparison_forest_resolve_model_names(
      model_names = NULL,
      include_shrinkage = FALSE,
      include_overall = FALSE,
      include_prediction = FALSE
    ),
    character(0)
  )
  expect_equal(
    interact:::comparison_forest_resolve_model_names(
      model_names = c("Overall"),
      include_shrinkage = FALSE,
      include_overall = TRUE,
      include_prediction = FALSE
    ),
    "Overall"
  )
  expect_equal(
    interact:::comparison_forest_resolve_model_names(
      model_names = c("Shrinkage", "Overall", "Prediction"),
      include_shrinkage = FALSE,
      include_overall = TRUE,
      include_prediction = TRUE
    ),
    c("Overall", "Prediction")
  )
  expect_equal(
    interact:::comparison_forest_resolve_model_names(
      model_names = c("Overall", "Prediction"),
      include_shrinkage = TRUE,
      include_overall = TRUE,
      include_prediction = TRUE
    ),
    c("Shrinkage", "Overall", "Prediction")
  )
  expect_equal(
    interact:::comparison_forest_resolve_model_names(
      model_names = c("Shrinkage", "Overall"),
      include_shrinkage = TRUE,
      include_overall = TRUE,
      include_prediction = TRUE
    ),
    c("Shrinkage", "Overall", "Prediction")
  )
  expect_equal(
    interact:::comparison_forest_resolve_model_names(
      model_names = NULL,
      include_shrinkage = FALSE,
      include_overall = FALSE,
      include_prediction = TRUE
    ),
    "Prediction"
  )
})

test_that("comparison_forest_build_fit_panels() includes shrinkage when include_shrinkage is TRUE without Shrinkage in model_names", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Overall", "Prediction"),
    include_shrinkage = TRUE,
    include_prediction = TRUE,
  )
  expect_true(any(panels$resultsA$row_kind == "shrinkage", na.rm = TRUE))
  expect_true(panels$include_shrinkage)
})

test_that("comparison_forest_build_fit_panels() omits shrinkage when include_shrinkage is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Overall", "Prediction"),
    include_shrinkage = FALSE,
    include_prediction = TRUE,
  )
  expect_false(any(panels$resultsA$row_kind == "shrinkage", na.rm = TRUE))
  expect_equal(nrow(panels$resultsA), panels$n_study + 6L)
  expect_true(any(panels$resultsA$row_kind == "combined", na.rm = TRUE))
})

test_that("comparison_forest_build_fit_panels() omits combined rows when include_combined is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Overall", "Prediction"),
    include_shrinkage = FALSE,
    include_combined = FALSE,
  )
  expect_false(any(panels$resultsA$row_kind == "combined", na.rm = TRUE))
  expect_equal(nrow(panels$resultsA), panels$n_study + 4L)
})

test_that("comparison_forest_build_fit_panels() keeps NA placeholders when include_prediction is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- cams(dat, reference = "Younger")
  full <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Shrinkage", "Overall", "Prediction"),
    include_shrinkage = TRUE,
    include_overall = TRUE,
    include_prediction = TRUE,
    n_var = "n"
  )
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Shrinkage", "Overall"),
    include_shrinkage = TRUE,
    include_overall = TRUE,
    include_prediction = FALSE,
    n_var = "n"
  )
  expect_equal(nrow(panels$resultsA), nrow(full$resultsA))
  expect_true(any(panels$resultsA$row_kind == "prediction", na.rm = TRUE))
  expect_true(all(is.na(panels$resultsA$yi[panels$resultsA$row_kind == "prediction"])))
  expect_true(any(is.finite(panels$resultsA$yi[panels$resultsA$row_kind == "overall"])))
})

test_that("comparison_forest_build_fit_panels() keeps NA placeholders when include_overall is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- cams(dat, reference = "Younger")
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Shrinkage", "Overall", "Prediction"),
    include_shrinkage = TRUE,
    include_overall = FALSE,
    include_prediction = TRUE,
    n_var = "n"
  )
  expect_true(any(panels$resultsA$row_kind == "overall", na.rm = TRUE))
  expect_true(all(is.na(panels$resultsA$yi[panels$resultsA$row_kind == "overall"])))
  expect_true(any(panels$resultsA$row_kind == "prediction", na.rm = TRUE))
  expect_true(any(panels$resultsA$row_kind == "shrinkage", na.rm = TRUE))
  pred_ix <- which(panels$resultsA$row_kind == "prediction")
  pred_sg_ix <- pred_ix[is.finite(panels$resultsA$yi[pred_ix])]
  expect_true(length(pred_sg_ix) >= 2L)
  expect_true(all(is.finite(panels$resultsA$yi[pred_sg_ix])))
  expect_true(all(is.finite(panels$resultsA$vi[pred_sg_ix])))
})

test_that("comparison_forest_build_fit_panels() includes prediction when include_prediction is TRUE without Prediction in model_names", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = c("Shrinkage", "Overall"),
    include_shrinkage = TRUE,
    include_overall = TRUE,
    include_prediction = TRUE,
  )
  expect_true(any(panels$resultsA$row_kind == "prediction", na.rm = TRUE))
  pred_ix <- which(panels$resultsA$row_kind == "prediction")
  pred_sg_ix <- pred_ix[is.finite(panels$resultsA$yi[pred_ix])]
  expect_true(length(pred_sg_ix) >= 2L)
  expect_true(all(is.finite(panels$resultsA$yi[pred_sg_ix])))
  expect_equal(unique(panels$row_labels$label[pred_sg_ix]), "Prediction")
  expect_true(all(nzchar(panels$row_labels$subheading[pred_ix])))
})

test_that("comparison_forest_build_fit_panels() uses diamonds and grey shrinkage", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  panels <- interact:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_overall = TRUE,
    include_prediction = TRUE,
  )
  expect_true(all(panels$resultsA$row_kind[panels$obs_row_index] == "study"))
  shr_ix <- which(panels$resultsA$row_kind == "shrinkage")
  expect_equal(length(shr_ix), panels$n_study + 2L)
  expect_true(all(nzchar(as.character(panels$row_labels$label[shr_ix]))))
  expect_false(any(interact:::comparison_forest_is_blank_row_label(
    panels$row_labels$label[shr_ix]
  )))
  int_ix <- which(panels$resultsB$row_kind == "shrinkage" & !is.na(panels$resultsB$shape))
  expect_equal(length(int_ix), 2L)
  sg_shr_ix <- shr_ix[!is.na(panels$resultsA$shape[shr_ix])]
  expect_true(all(is.finite(panels$resultsA$yi[sg_shr_ix])))
  expect_true(all(is.finite(panels$resultsA$vi[sg_shr_ix])))
  expect_true(all(is.finite(panels$resultsA$ci.lb[sg_shr_ix])))
  expect_true(all(is.finite(panels$resultsA$ci.ub[sg_shr_ix])))
  expect_true(all(is.finite(panels$resultsB$yi[int_ix])))
  expect_true(all(is.finite(panels$resultsB$vi[int_ix])))
  expect_true(all(is.finite(panels$resultsB$ci.lb[int_ix])))
  expect_true(all(is.finite(panels$resultsB$ci.ub[int_ix])))
  expect_true(all(
    panels$resultsA$fill[sg_shr_ix] == interact:::comparison_forest_shrinkage_colour()
  ))
  overall_ix <- which(panels$resultsA$row_kind == "overall")
  pred_ix <- which(panels$resultsA$row_kind == "prediction")
  combined_ix <- which(panels$resultsA$row_kind == "combined")
  expect_true(length(overall_ix) >= 2L)
  expect_equal(length(pred_ix), 2L)
  expect_equal(length(combined_ix), 2L)
  expect_true(all(
    nzchar(panels$resultsA$ifrac[overall_ix]) &
      grepl("%", panels$resultsA$ifrac[overall_ix], fixed = TRUE)
  ))
  ref <- as.character(interactmeta_reference(fit))
  ov_y <- overall_ix[panels$row_labels$colour_label[overall_ix] == ref]
  ov_o <- overall_ix[panels$row_labels$colour_label[overall_ix] == "Older"]
  if (length(ov_y) && length(ov_o)) {
    pct <- function(x) as.numeric(sub("%.*", "", as.character(x)))
    expect_true(abs(pct(panels$resultsA$ifrac[ov_y[[1L]]]) +
      pct(panels$resultsA$ifrac[ov_o[[1L]]]) - 100) < 0.2)
  }
  expect_true(all(panels$resultsA$ifrac[pred_ix] == ""))
  overall_if <- panels$resultsA$ifrac[
    combined_ix[panels$row_labels$label[combined_ix] == "Overall"]
  ]
  expect_true(length(overall_if) == 1L)
  expect_equal(overall_if, "")
  expect_true(all(
    panels$resultsA$ifrac[
      combined_ix[panels$row_labels$label[combined_ix] == "Prediction"]
    ] == ""
  ))
  expect_true(all(panels$resultsA$shape[combined_ix] == 18L))
  expect_true(all(is.na(panels$resultsB$yi[combined_ix])))
  marg_ix <- c(overall_ix, pred_ix, combined_ix)
  expect_equal(
    unique(panels$row_labels$subheading[marg_ix]),
    c("Younger", "Older", "All")
  )
  expect_equal(
    panels$row_labels$subheading[combined_ix],
    rep(interact:::forestplot_combined_block_label(), 2L)
  )
  expect_equal(panels$row_labels$colour_label[overall_ix], c("Younger", "Older"))
  pred_sg_ix <- pred_ix[is.finite(panels$resultsA$yi[pred_ix])]
  expect_equal(panels$row_labels$colour_label[pred_sg_ix], c("Younger", "Older"))
  # Summary diamonds: blue for all pooled / prediction / overall.
  expect_true(all(panels$resultsA$shape[overall_ix] == 18L))
  pred_int_ix <- pred_ix[is.finite(panels$resultsB$yi[pred_ix])]
  expect_true(all(panels$resultsA$shape[pred_sg_ix] == 18L))
  expect_true(all(panels$resultsB$shape[pred_int_ix] == 18L))
  pooled_cols <- interact:::forestplot_summary_diamond_colours(
    panels$row_labels$colour_label[overall_ix],
    reference = "Younger"
  )
  expect_equal(panels$resultsA$fill[overall_ix], pooled_cols)
  expect_equal(panels$resultsA$colour[overall_ix], pooled_cols)
  pred_cols <- interact:::forestplot_summary_diamond_colours(
    panels$row_labels$colour_label[pred_sg_ix],
    reference = "Younger"
  )
  expect_equal(panels$resultsA$fill[pred_sg_ix], pred_cols)
  expect_equal(panels$resultsA$colour[pred_sg_ix], pred_cols)
  expect_true(all(
    panels$resultsA$fill[combined_ix] == interact:::forestplot_overall_summary_colour()
  ))
  overall_int_ix <- overall_ix[is.finite(panels$resultsB$yi[overall_ix])]
  expect_true(length(overall_int_ix) >= 1L)
  expect_equal(unique(panels$row_labels$label[pred_ix]), "Prediction")
  expect_equal(unique(panels$row_labels$label[overall_ix]), "Overall")
  obs_ix <- panels$obs_row_index
  expect_true(all(panels$resultsA$point_scale[obs_ix] == 1))
  expect_true(all(panels$resultsA$point_scale[shr_ix] == interact:::comparison_forest_shrinkage_point_scale()))
  expect_true(all(panels$resultsB$point_scale[int_ix] == interact:::comparison_forest_shrinkage_point_scale()))
})
