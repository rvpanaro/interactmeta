test_that("forestplot.interactmeta() requires interactmeta fit metadata", {
  skip_if_not_installed("bayesmeta")
  expect_error(
    forestplot.interactmeta(list()),
    "interactmeta|CAMS fit"
  )
  bare <- structure(list(), class = c("interactmeta", "cam", "list"))
  expect_error(
    forestplot.interactmeta(bare),
    "ifrac_data"
  )
})

test_that("comparison_forest_resolve_model_names() warns on unrecognized layers", {
  expect_warning(
    out <- interactmeta:::comparison_forest_resolve_model_names(
      model_names = c("Nope", "Pooled"),
      include_shrinkage = FALSE,
      include_pooled = TRUE,
      include_prediction = FALSE
    ),
    "Unrecognized `model_names`"
  )
  expect_equal(out, "Pooled")
})

test_that("omitted pooled/prediction layers keep NA placeholder rows", {
  pooled_sg <- data.frame(
    trial = c("0.5", "0.5"),
    subgroup = c("Younger", "Older"),
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.01),
    stringsAsFactors = FALSE
  )
  pred_sg <- pooled_sg
  pred_sg$yi <- c(0.15, 0.25)
  empty_int <- pooled_sg[0L, , drop = FALSE]

  both <- interactmeta:::comparison_forest_interleave_pooled_prediction(
    pooled_sg = pooled_sg,
    pred_sg = pred_sg,
    pooled_int = empty_int,
    pred_int = empty_int,
    reference = "Younger",
    include_pooled = TRUE,
    include_prediction = TRUE
  )
  no_pred <- interactmeta:::comparison_forest_interleave_pooled_prediction(
    pooled_sg = pooled_sg,
    pred_sg = pred_sg,
    pooled_int = empty_int,
    pred_int = empty_int,
    reference = "Younger",
    include_pooled = TRUE,
    include_prediction = FALSE
  )
  expect_equal(nrow(both$subgroup), nrow(no_pred$subgroup))
  expect_equal(no_pred$subgroup$row_kind, rep(c("pooled", "prediction"), 2L))
  expect_true(all(is.finite(no_pred$subgroup$yi[no_pred$subgroup$row_kind == "pooled"])))
  expect_true(all(is.na(no_pred$subgroup$yi[no_pred$subgroup$row_kind == "prediction"])))

  no_pooled <- interactmeta:::comparison_forest_interleave_pooled_prediction(
    pooled_sg = pooled_sg,
    pred_sg = pred_sg,
    pooled_int = empty_int,
    pred_int = empty_int,
    reference = "Younger",
    include_pooled = FALSE,
    include_prediction = TRUE
  )
  expect_equal(nrow(both$subgroup), nrow(no_pooled$subgroup))
  expect_true(all(is.na(no_pooled$subgroup$yi[no_pooled$subgroup$row_kind == "pooled"])))
  expect_true(all(is.finite(no_pooled$subgroup$yi[no_pooled$subgroup$row_kind == "prediction"])))
})

test_that("forestplot() validates prob on the comparison path", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger"))
  expect_error(
    forestplot(fit, prob = -0.1),
    "`prob` must be"
  )
  expect_warning(
    forestplot(
      fit,
      model_names = "Nope",
      include_shrinkage = FALSE,
      include_pooled = TRUE,
      include_prediction = FALSE
    ),
    "Unrecognized `model_names`"
  )
})

test_that("comparison_forest_order_to_template() orders by trial and subgroup", {
  template <- data.frame(
    trial = rep(c("B", "A"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(c("Older", "Younger", "Older", "Younger"), levels = c("Younger", "Older")),
    yi = 1:4,
    vi = 0.1,
    stringsAsFactors = FALSE
  )
  out <- interactmeta:::comparison_forest_order_to_template(dat, template, "trial", "subgroup")
  expect_equal(out$trial, c("B", "B", "A", "A"))
  expect_equal(as.character(out$subgroup), rep(c("Younger", "Older"), 2))
})

test_that("comparison_forest_order_to_template() matches split_ifrac_for_forest order", {
  d <- ifrac(
    yi = c(0.2, 0.1, 0.3, 0.15),
    vi = c(0.01, 0.04, 0.0225, 0.09),
    trial = c("B", "B", "A", "A"),
    subgroup = factor(c("Older", "Younger", "Older", "Younger"), levels = c("Younger", "Older")),
    reference = "Younger"
  )
  template <- interactmeta:::split_ifrac_for_forest(d, contrast_sep = "-")$subgroup
  sg <- interactmeta:::comparison_forest_subgroup_rows(d, "trial", "subgroup", "-")
  out <- interactmeta:::comparison_forest_order_to_template(sg, template, "trial", "subgroup")
  expect_equal(
    paste(out$trial, as.character(out$subgroup)),
    paste(template$trial, as.character(template$subgroup))
  )
})

test_that("CAMS forest panels keep log-scale marginal estimates on the ratio scale", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = log(c(1.1, 1.2, 1.15, 1.25)),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Pooled"),
    include_shrinkage = FALSE,
    include_pooled = TRUE,
    include_prediction = FALSE,
    include_overall = FALSE
  )
  pooled_yi <- panels$resultsA$yi[panels$resultsA$row_kind == "pooled"]
  study_yi <- panels$resultsA$yi[panels$obs_row_index]
  expect_true(all(abs(pooled_yi) < 1))
  expect_true(all(abs(study_yi) < 1))
  limits <- interactmeta:::resolve_forestplot_limits(
    list(panels$resultsA, panels$resultsB),
    exponentiate = TRUE
  )
  expect_true(all(limits$xlim[[1]] > 0))
})

test_that("forestplot_render() keeps log-scale panels when exponentiate is TRUE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  built <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Pooled"),
    include_shrinkage = FALSE,
    include_pooled = TRUE,
    include_prediction = FALSE,
    include_overall = FALSE
  )
  render_args <- c(
    built,
    list(
      prob = 0.95,
      contrast_sep = "-"
    )
  )
  fp_log <- do.call(interactmeta:::forestplot_render, c(render_args, list(exponentiate = FALSE)))
  fp_ratio <- do.call(interactmeta:::forestplot_render, c(render_args, list(exponentiate = TRUE)))
  expect_equal(attr(fp_log, "panels")[[1]]$yi, attr(fp_ratio, "panels")[[1]]$yi)
  limits <- interactmeta:::resolve_forestplot_limits(attr(fp_ratio, "panels"), exponentiate = TRUE)
  expect_true(all(unlist(limits$xlim) > 0))
})

test_that("comparison_forest_overall_subgroup_block() orders subgroups", {
  summary <- data.frame(
    subgroup = factor(c("Older", "Younger"), levels = c("Younger", "Older")),
    type = "Pooled",
    mean = c(0.2, 0.1),
    sd = c(0.05, 0.04),
    stringsAsFactors = FALSE
  )
  block <- interactmeta:::comparison_forest_overall_subgroup_block(
    summary,
    "Pooled",
    "0.500",
    subgroup_order = c("Younger", "Older")
  )
  expect_equal(as.character(block$subgroup), c("Younger", "Older"))
})

test_that("comparison_forest_subgroup_order() matches within-trial subgroup order", {
  d <- ifrac(
    yi = c(0.2, 0.1, 0.3, 0.15),
    vi = c(0.01, 0.04, 0.0225, 0.09),
    trial = c("B", "B", "A", "A"),
    subgroup = factor(c("Older", "Younger", "Older", "Younger"), levels = c("Younger", "Older")),
    reference = "Younger"
  )
  expect_equal(interactmeta:::comparison_forest_subgroup_order(d, "-"), c("Younger", "Older"))
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
      trial = c("A", "A"),
      subgroup = c("Younger", "Older"),
      yi = c(0.12, 0.18),
      vi = c(0.008, 0.009),
      stringsAsFactors = FALSE
    ),
    interaction = NULL
  )
  ifrac_disp <- c("80.0%", "20.0%")
  out <- interactmeta:::interleave_comparison_shrinkage_rows(
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
  styled <- interactmeta:::comparison_forest_apply_shrinkage_colours(
    out$resultsA,
    out$resultsB,
    out$row_labels,
    out$subgroup_shrinkage_row_index,
    out$interaction_shrinkage_row_index,
    "Younger"
  )
  expect_equal(styled$resultsA$fill[2], interactmeta:::comparison_forest_shrinkage_colour())
  expect_equal(styled$resultsA$fill[4], interactmeta:::comparison_forest_shrinkage_colour())
  expect_equal(styled$resultsA$colour[2], interactmeta:::comparison_forest_shrinkage_colour())
})

test_that("comparison_forest_shrinkage_data() returns fitted values", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")

  shrink <- interactmeta:::comparison_forest_shrinkage_data(fit, prob = 0.95)
  expect_equal(nrow(shrink$subgroup), 4L)
  expect_true(all(is.finite(shrink$subgroup$yi)))

  sg <- shrinkage(fit, type = "subgroups", summary = TRUE)$summary
  int <- shrinkage(fit, type = "interaction", summary = TRUE)$summary
  ov <- shrinkage(fit, type = "overall", summary = TRUE)$summary
  key_sg <- paste(shrink$subgroup$trial, shrink$subgroup$subgroup, sep = "\r")
  key_fit <- paste(sg$trial, sg$subgroup, sep = "\r")
  m <- match(key_sg, key_fit)
  expect_false(any(is.na(m)))
  expect_equal(shrink$subgroup$yi, sg$estimate[m], tolerance = 1e-10)
  expect_equal(shrink$subgroup$vi, pmax(sg$sd[m]^2, 1e-12), tolerance = 1e-10)

  key_int <- paste(shrink$interaction$trial, shrink$interaction$subgroup, sep = "\r")
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
    trial = rep(c("A", "B"), each = 2),
    subgroup = rep(c("Younger", "Older"), 2),
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
    label = resultsA$subgroup,
    stringsAsFactors = FALSE
  )
  shrinkage <- list(
    subgroup = data.frame(
      trial = resultsA$trial,
      subgroup = resultsA$subgroup,
      yi = resultsA$yi + 0.01,
      vi = resultsA$vi,
      stringsAsFactors = FALSE
    ),
    interaction = data.frame(
      trial = c("A", "B"),
      subgroup = c("Older-Younger", "Older-Younger"),
      yi = c(0.11, 0.12),
      vi = c(0.01, 0.01),
      stringsAsFactors = FALSE
    )
  )
  out <- interactmeta:::interleave_comparison_shrinkage_rows(
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
    interactmeta:::comparison_forest_resolve_model_names(
      model_names = NULL,
      include_shrinkage = TRUE,
      include_pooled = FALSE,
      include_prediction = TRUE
    ),
    c("Shrinkage", "Prediction")
  )
  expect_equal(
    interactmeta:::comparison_forest_resolve_model_names(
      model_names = NULL,
      include_shrinkage = FALSE,
      include_pooled = FALSE,
      include_prediction = FALSE
    ),
    character(0)
  )
  expect_equal(
    interactmeta:::comparison_forest_resolve_model_names(
      model_names = c("Pooled"),
      include_shrinkage = FALSE,
      include_pooled = TRUE,
      include_prediction = FALSE
    ),
    "Pooled"
  )
  expect_equal(
    interactmeta:::comparison_forest_resolve_model_names(
      model_names = c("Shrinkage", "Pooled", "Prediction"),
      include_shrinkage = FALSE,
      include_pooled = TRUE,
      include_prediction = TRUE
    ),
    c("Pooled", "Prediction")
  )
  expect_equal(
    interactmeta:::comparison_forest_resolve_model_names(
      model_names = c("Pooled", "Prediction"),
      include_shrinkage = TRUE,
      include_pooled = TRUE,
      include_prediction = TRUE
    ),
    c("Shrinkage", "Pooled", "Prediction")
  )
  expect_equal(
    interactmeta:::comparison_forest_resolve_model_names(
      model_names = c("Shrinkage", "Pooled"),
      include_shrinkage = TRUE,
      include_pooled = TRUE,
      include_prediction = TRUE
    ),
    c("Shrinkage", "Pooled", "Prediction")
  )
  expect_equal(
    interactmeta:::comparison_forest_resolve_model_names(
      model_names = NULL,
      include_shrinkage = FALSE,
      include_pooled = FALSE,
      include_prediction = TRUE
    ),
    "Prediction"
  )
})

test_that("comparison_forest_build_fit_panels() includes shrinkage when include_shrinkage is TRUE without Shrinkage in model_names", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Pooled", "Prediction"),
    include_shrinkage = TRUE,
    include_prediction = TRUE,
  )
  expect_true(any(panels$resultsA$row_kind == "shrinkage", na.rm = TRUE))
  expect_true(panels$include_shrinkage)
})

test_that("comparison_forest_build_fit_panels() omits shrinkage when include_shrinkage is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Pooled", "Prediction"),
    include_shrinkage = FALSE,
    include_prediction = TRUE,
  )
  expect_false(any(panels$resultsA$row_kind == "shrinkage", na.rm = TRUE))
  expect_equal(nrow(panels$resultsA), panels$n_study + 6L)
  expect_true(any(panels$resultsA$row_kind == "overall", na.rm = TRUE))
})

test_that("comparison_forest_build_fit_panels() omits overall when include_overall is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Pooled", "Prediction"),
    include_shrinkage = FALSE,
    include_overall = FALSE,
  )
  expect_false(any(panels$resultsA$row_kind == "overall", na.rm = TRUE))
  expect_equal(nrow(panels$resultsA), panels$n_study + 4L)
})

test_that("comparison_forest_build_fit_panels() keeps NA placeholders when include_prediction is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  full <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Shrinkage", "Pooled", "Prediction"),
    include_shrinkage = TRUE,
    include_pooled = TRUE,
    include_prediction = TRUE,
    n_var = "n"
  )
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Shrinkage", "Pooled"),
    include_shrinkage = TRUE,
    include_pooled = TRUE,
    include_prediction = FALSE,
    n_var = "n"
  )
  expect_equal(nrow(panels$resultsA), nrow(full$resultsA))
  expect_true(any(panels$resultsA$row_kind == "prediction", na.rm = TRUE))
  expect_true(all(is.na(panels$resultsA$yi[panels$resultsA$row_kind == "prediction"])))
  expect_true(any(is.finite(panels$resultsA$yi[panels$resultsA$row_kind == "pooled"])))
})

test_that("comparison_forest_build_fit_panels() keeps NA placeholders when include_pooled is FALSE", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Shrinkage", "Pooled", "Prediction"),
    include_shrinkage = TRUE,
    include_pooled = FALSE,
    include_prediction = TRUE,
    n_var = "n"
  )
  expect_true(any(panels$resultsA$row_kind == "pooled", na.rm = TRUE))
  expect_true(all(is.na(panels$resultsA$yi[panels$resultsA$row_kind == "pooled"])))
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
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = c("Shrinkage", "Pooled"),
    include_shrinkage = TRUE,
    include_pooled = TRUE,
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
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_pooled = TRUE,
    include_prediction = TRUE,
  )
  expect_true(all(panels$resultsA$row_kind[panels$obs_row_index] == "study"))
  shr_ix <- which(panels$resultsA$row_kind == "shrinkage")
  expect_equal(length(shr_ix), panels$n_study + 2L)
  expect_true(all(nzchar(as.character(panels$row_labels$label[shr_ix]))))
  expect_false(any(interactmeta:::comparison_forest_is_blank_row_label(
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
    panels$resultsA$fill[sg_shr_ix] == interactmeta:::comparison_forest_shrinkage_colour()
  ))
  pooled_ix <- which(panels$resultsA$row_kind == "pooled")
  pred_ix <- which(panels$resultsA$row_kind == "prediction")
  overall_ix <- which(panels$resultsA$row_kind == "overall")
  expect_true(length(pooled_ix) >= 2L)
  expect_equal(length(pred_ix), 2L)
  expect_equal(length(overall_ix), 2L)
  expect_true(all(panels$resultsA$shape[overall_ix] == 18L))
  expect_true(all(is.na(panels$resultsB$yi[overall_ix])))
  marg_ix <- c(pooled_ix, pred_ix, overall_ix)
  expect_equal(
    unique(panels$row_labels$subheading[marg_ix]),
    c("Younger", "Older", "Overall")
  )
  expect_equal(panels$row_labels$subheading[overall_ix], rep("Overall", 2L))
  expect_equal(panels$row_labels$colour_label[pooled_ix], c("Younger", "Older"))
  pred_sg_ix <- pred_ix[is.finite(panels$resultsA$yi[pred_ix])]
  expect_equal(panels$row_labels$colour_label[pred_sg_ix], c("Younger", "Older"))
  # Summary diamonds: blue for all pooled / prediction / overall.
  expect_true(all(panels$resultsA$shape[pooled_ix] == 18L))
  pred_int_ix <- pred_ix[is.finite(panels$resultsB$yi[pred_ix])]
  expect_true(all(panels$resultsA$shape[pred_sg_ix] == 18L))
  expect_true(all(panels$resultsB$shape[pred_int_ix] == 18L))
  pooled_cols <- interactmeta:::forestplot_summary_diamond_colours(
    panels$row_labels$colour_label[pooled_ix],
    reference = "Younger"
  )
  expect_equal(panels$resultsA$fill[pooled_ix], pooled_cols)
  expect_equal(panels$resultsA$colour[pooled_ix], pooled_cols)
  pred_cols <- interactmeta:::forestplot_summary_diamond_colours(
    panels$row_labels$colour_label[pred_sg_ix],
    reference = "Younger"
  )
  expect_equal(panels$resultsA$fill[pred_sg_ix], pred_cols)
  expect_equal(panels$resultsA$colour[pred_sg_ix], pred_cols)
  expect_true(all(
    panels$resultsA$fill[overall_ix] == interactmeta:::forestplot_overall_summary_colour()
  ))
  pooled_int_ix <- pooled_ix[is.finite(panels$resultsB$yi[pooled_ix])]
  expect_true(length(pooled_int_ix) >= 1L)
  expect_equal(unique(panels$row_labels$label[pred_ix]), "Prediction")
  expect_equal(unique(panels$row_labels$label[pooled_ix]), "Pooled")
  obs_ix <- panels$obs_row_index
  expect_true(all(panels$resultsA$point_scale[obs_ix] == 1))
  expect_true(all(panels$resultsA$point_scale[shr_ix] == interactmeta:::comparison_forest_shrinkage_point_scale()))
  expect_true(all(panels$resultsB$point_scale[int_ix] == interactmeta:::comparison_forest_shrinkage_point_scale()))
})

test_that("marginal rows keep subgroup grouping labels for the forest renderer", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger"))
  built <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = NULL,
    include_shrinkage = FALSE,
    include_pooled = TRUE,
    include_prediction = TRUE
  )
  marg_ix <- which(built$resultsA$row_kind %in% c("pooled", "prediction", "overall"))
  expect_true(length(marg_ix) >= 4L)
  older_ix <- marg_ix[built$row_labels$subheading[marg_ix] == "Older"]
  expect_true(length(older_ix) >= 1L)
  fp <- forestplot(fit)
  expect_s3_class(fp, "interactmeta_forestplot")
  expect_true(!is.null(fp$figure) || !is.null(fp$plot))
  expect_identical(fp$engine, "bayesmeta")
})

test_that("forestplot() includes shrinkage overlays by default", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  fp <- forestplot(fit)
  expect_true(any(attr(fp, "panels")[[1]]$row_kind == "shrinkage", na.rm = TRUE))
  fp_off <- forestplot(fit, include_shrinkage = FALSE)
  expect_false(any(attr(fp_off, "panels")[[1]]$row_kind == "shrinkage", na.rm = TRUE))
})

test_that("CAMS forest shrinkage fill is consistent within every trial after render", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.12, 0.22),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.5, 0.5, 0.7, 0.7),
    ifrac2_Older = c(-0.3, 0.7, -0.5, 0.5, -0.7, 0.3)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  fp <- forestplot(fit, include_shrinkage = TRUE)
  panels <- attr(fp, "panels")
  ra <- panels[[1]]
  rb <- panels[[2]]
  shr_ix <- which(ra$row_kind == "shrinkage" & !is.na(ra$shape))
  expect_true(length(shr_ix) >= 6L)
  expect_false(any(ra$fill[shr_ix] %in% c("white", NA_character_)))
  expect_true(all(ra$fill[shr_ix] == interactmeta:::comparison_forest_shrinkage_colour()))
  expect_true(all(ra$colour[shr_ix] == interactmeta:::comparison_forest_shrinkage_colour()))
})

test_that("forest markers are black with gray45 shrinkage", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger"))
  fp <- forestplot(
    fit,
    type = "subgroups",
    include_shrinkage = TRUE
  )
  ra <- attr(fp, "panels")[[1]]
  pooled_ix <- which(ra$row_kind == "pooled")
  pred_ix <- which(ra$row_kind == "prediction" & is.finite(ra$yi))
  study_ix <- which(ra$row_kind == "study")
  shr_ix <- which(ra$row_kind == "shrinkage" & !is.na(ra$shape))
  expect_true(length(pooled_ix) >= 2L)
  expect_true(length(pred_ix) >= 2L)
  rl <- fp$specs$row_labels
  expect_true(all(
    ra$fill[pooled_ix] == interactmeta:::forestplot_summary_diamond_colours(
      rl$heading[pooled_ix],
      reference = "Younger"
    )
  ))
  expect_true(all(
    ra$fill[pred_ix] == interactmeta:::forestplot_summary_diamond_colours(
      rl$heading[pred_ix],
      reference = "Younger"
    )
  ))
  expect_true(all(ra$fill[study_ix] %in% c(interactmeta:::comparison_forest_marker_colour(), "black", "#000000")))
  expect_true(length(shr_ix) >= 1L)
  expect_true(all(ra$fill[shr_ix] == interactmeta:::comparison_forest_shrinkage_colour()))

  fp_int <- forestplot(fit, type = "interaction", include_shrinkage = TRUE)
  ra_int <- attr(fp_int, "panels")[[1]]
  rl_int <- fp_int$specs$row_labels
  int_ix <- which(
    ra_int$row_kind %in% c("pooled", "prediction") &
      rl_int$heading == "Interaction"
  )
  expect_true(length(int_ix) >= 2L)
  expect_true(all(ra_int$fill[int_ix] == interactmeta:::forestplot_interaction_summary_colour()))
})

test_that("include_shrinkage with group=subgroup uses alphabetical trial pairs", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  dat <- data.frame(
    trial = rep(c("C", "A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.12, 0.22),
    vi = 0.01,
    n = 10,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.5, 0.5, 0.7, 0.7),
    ifrac2_Older = c(-0.3, 0.7, -0.5, 0.5, -0.7, 0.3)
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  fp <- forestplot(
    fit,
    type = "subgroups",
    include_shrinkage = TRUE,
    group = "subgroup"
  )
  labs <- fp$specs$row_labels
  rk <- fp$specs$panels[[1]]$row_kind
  yng_ix <- which(labs$heading == "Younger" & rk %in% c("study", "shrinkage"))
  expect_equal(
    labs$label[yng_ix],
    c("A", "A", "B", "B", "C", "C")
  )
  expect_equal(
    rk[yng_ix],
    c("study", "shrinkage", "study", "shrinkage", "study", "shrinkage")
  )

  fp_int <- forestplot(
    fit,
    type = "interaction",
    include_shrinkage = TRUE,
    group = "subgroup"
  )
  labs_int <- fp_int$specs$row_labels
  rk_int <- fp_int$specs$panels[[1]]$row_kind
  int_ix <- which(labs_int$heading == "Interaction" & rk_int %in% c("study", "shrinkage"))
  expect_equal(
    labs_int$label[int_ix],
    c("A", "A", "B", "B", "C", "C")
  )
  expect_equal(
    rk_int[int_ix],
    c("study", "shrinkage", "study", "shrinkage", "study", "shrinkage")
  )
})

test_that("forestplot(fit, exponentiate = TRUE) runs without log-axis warnings", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  w <- character()
  withCallingHandlers(
    fp <- forestplot(fit, exponentiate = TRUE),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  )
  expect_type(fp, "list")
  log_warn <- grep("log|NaN|infinite", w, ignore.case = TRUE)
  expect_length(log_warn, 0)
})

test_that("forestplot(fit) passes weight with shrinkage and exponentiate", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.01, 0.04, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.6, 0.6),
    ifrac2_Older = c(-0.3, 0.7, -0.6, 0.4)
  )
  fit <- cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger")
  fp <- forestplot(
    fit,
    exponentiate = TRUE,
    include_shrinkage = TRUE,
    ifrac = "interaction"
  )
  opts <- attr(fp, "forestplot_fit_options")
  expect_equal(opts$ifrac, "interaction")
  expect_true(opts$exponentiate)
  render <- attr(fp, "forestplot_fit_render")
  expect_equal(render$include_shrinkage, TRUE)
  panels <- attr(fp, "panels")
  expect_true(any(panels[[1]]$row_kind == "shrinkage"))
  expect_true(any(panels[[1]]$row_kind == "pooled"))
  expect_match(attr(fp, "footnote"), "Information fraction \\(\u03c0\\)")

  fp_overall <- forestplot(
    fit,
    exponentiate = TRUE,
    include_shrinkage = TRUE,
    ifrac = "overall"
  )
  expect_match(attr(fp_overall, "footnote"), "Information fraction \\(\u03c0\\)")

  # Pooled / Prediction / Overall diamonds use ifrac(..., ifrac =); panels
  # should differ when overall- vs interaction-weighted ω diverge.
  panels_int <- attr(fp, "panels")[[1]]
  panels_ov <- attr(fp_overall, "panels")[[1]]
  ix_sum <- panels_int$row_kind %in% c("pooled", "prediction", "overall")
  expect_false(isTRUE(all.equal(
    panels_int$yi[ix_sum],
    panels_ov$yi[ix_sum],
    tolerance = 1e-8
  )))
})

test_that("forestplot IF annotation prefix follows weight", {
  expect_equal(
    interactmeta:::forestplot_ifrac_annotation_prefix("overall"),
    "Information fraction (\u03c0)"
  )
  expect_equal(
    interactmeta:::forestplot_ifrac_annotation_prefix("interaction"),
    "Information fraction (\u03c0)"
  )
})

test_that("forestplot(fit) adds heterogeneity and overall-IF annotations", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger")
  fp <- forestplot(
    fit,
        include_shrinkage = FALSE,
    include_pooled = FALSE,
    include_prediction = FALSE,
    include_overall = FALSE
  )
  foot <- attr(fp, "footnote")
  expect_true(is.character(foot) && nzchar(foot))
  expect_match(foot, "Information fraction \\(\u03c0\\)")
  expect_false(grepl("% interval", foot, fixed = TRUE))
  expect_match(foot, "\\([0-9.]+, [0-9.]+\\)")
  expect_match(foot, "Heterogeneity (\u03c4)", fixed = TRUE)
  expect_false(grepl("Overall IF", foot, fixed = TRUE))
  expect_false(grepl("Overall heterogeneity", foot, fixed = TRUE))
  expect_match(foot, "Interaction heterogeneity (\u03c4_\u03b3):", fixed = TRUE)
  expect_false(grepl("Older - Younger", foot, fixed = TRUE))
  expect_true("Older" %in% names(attr(fp, "annotations")$ifrac))
  expect_true("Older" %in% names(attr(fp, "annotations")$subgroup_ifrac))
  expect_true("overall" %in% names(attr(fp, "annotations")$heterogeneity))
  expect_true("interaction" %in% names(attr(fp, "annotations")$heterogeneity))
  expect_equal(length(attr(fp, "annotations")$lines), 0L)

  fp_off <- forestplot(
    fit,
        show_annotations = FALSE,
    include_shrinkage = FALSE,
    include_pooled = FALSE,
    include_prediction = FALSE,
    include_overall = FALSE
  )
  expect_null(attr(fp_off, "footnote"))
})

forestplot_unwrap_row_labels <- function(rl) {
  if (is.list(rl) && !is.data.frame(rl)) {
    rl[[1L]]
  } else {
    rl
  }
}

test_that("all forestplot shrinkage rows match shrinkage() from the fit", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.01, 0.04, 0.01),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger"))

  sg_me <- shrinkage(fit, type = "subgroups", summary = TRUE)$summary
  int_me <- shrinkage(fit, type = "interaction", summary = TRUE)$summary
  ov_me <- shrinkage(fit, type = "overall", summary = TRUE)$summary

  shr <- interactmeta:::comparison_forest_shrinkage_data(
    fit,
    prob = 0.95,
    max_draws = NULL,
    ifrac = "overall"
  )
  # IV mix of subgroup shrinkages differs when subgroup vis differ
  sg <- shr$subgroup
  iv_A <- {
    ii <- which(sg$trial == "A")
    w <- 1 / pmax(sg$vi[ii], 1e-12)
    sum(sg$yi[ii] * w) / sum(w)
  }
  expect_false(isTRUE(all.equal(
    shr$overall$yi[shr$overall$trial == "A"],
    iv_A,
    tolerance = 1e-8
  )))

  # type = "subgroups" (heading = subgroup, label = trial)
  fp_sg <- forestplot(fit, type = "subgroups", include_shrinkage = TRUE)
  ra <- attr(fp_sg, "panels")[[1]]
  rl <- forestplot_unwrap_row_labels(fp_sg$specs$row_labels)
  ix_sg <- which(
    ra$row_kind == "shrinkage" &
      !is.na(ra$shape) &
      is.finite(ra$yi) &
      !interactmeta:::forestplot_heading_is_overall(rl$heading)
  )
  expect_equal(length(ix_sg), 4L)
  for (i in ix_sg) {
    hit <- which(
      as.character(sg_me$trial) == as.character(rl$label[[i]]) &
        as.character(sg_me$subgroup) == as.character(rl$heading[[i]])
    )
    expect_length(hit, 1L)
    expect_equal(ra$yi[[i]], sg_me$estimate[[hit]], tolerance = 1e-8)
    expect_equal(ra$vi[[i]], pmax(sg_me$sd[[hit]]^2, 1e-12), tolerance = 1e-8)
  }

  # type = "interaction"
  fp_int <- forestplot(fit, type = "interaction", include_shrinkage = TRUE)
  rb <- attr(fp_int, "panels")[[2]]
  rl_int <- forestplot_unwrap_row_labels(fp_int$specs$row_labels)
  ix_int <- which(rb$row_kind == "shrinkage" & is.finite(rb$yi))
  expect_equal(length(ix_int), 2L)
  for (i in ix_int) {
    hit <- which(as.character(int_me$trial) == as.character(rl_int$label[[i]]))
    expect_length(hit, 1L)
    expect_equal(rb$yi[[i]], int_me$estimate[[hit]], tolerance = 1e-8)
    expect_equal(rb$vi[[i]], pmax(int_me$sd[[hit]]^2, 1e-12), tolerance = 1e-8)
  }

  # type = "overall"
  fp_ov <- forestplot(fit, type = "overall", include_shrinkage = TRUE)
  ra_ov <- attr(fp_ov, "panels")[[1]]
  rl_ov <- forestplot_unwrap_row_labels(fp_ov$specs$row_labels)
  ix_ov <- which(
    ra_ov$row_kind == "shrinkage" &
      interactmeta:::forestplot_heading_is_overall(rl_ov$heading)
  )
  expect_equal(length(ix_ov), 2L)
  for (i in ix_ov) {
    tid <- as.character(rl_ov$label[[i]])
    hit <- which(as.character(ov_me$trial) == tid)
    expect_length(hit, 1L)
    expect_equal(ra_ov$yi[[i]], ov_me$estimate[[hit]], tolerance = 1e-8)
    expect_equal(ra_ov$vi[[i]], pmax(ov_me$sd[[hit]]^2, 1e-12), tolerance = 1e-8)
  }

  # group = "trial" still uses subgroup shrinkage from the fit
  fp_tr <- forestplot(
    fit,
    type = "subgroups",
    group = "trial",
    include_shrinkage = TRUE
  )
  ra_tr <- attr(fp_tr, "panels")[[1]]
  rl_tr <- forestplot_unwrap_row_labels(fp_tr$specs$row_labels)
  ix_tr <- which(
    ra_tr$row_kind == "shrinkage" &
      !is.na(ra_tr$shape) &
      is.finite(ra_tr$yi) &
      !interactmeta:::forestplot_heading_is_overall(rl_tr$heading)
  )
  expect_equal(length(ix_tr), 4L)
  for (i in ix_tr) {
    hit <- which(
      as.character(sg_me$trial) == as.character(rl_tr$heading[[i]]) &
        as.character(sg_me$subgroup) == as.character(rl_tr$label[[i]])
    )
    expect_length(hit, 1L)
    expect_equal(ra_tr$yi[[i]], sg_me$estimate[[hit]], tolerance = 1e-8)
  }
})
