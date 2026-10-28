test_that("forestplot.interactmeta() requires interactmeta fit metadata", {
  skip_if_not_installed("bayesmeta")
  expect_error(
    forestplot.interactmeta(list()),
    "interactmeta|CAMS fit"
  )
  bare <- structure(list(), class = c("interactmeta", "list"))
  expect_error(
    forestplot.interactmeta(bare),
    "ifrac_data"
  )
})

test_that("forestplot() validates prob on the comparison path", {
  slow()
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
  fit <- suppressWarnings(cams(
    dat,
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
      include_overall = TRUE,
      include_prediction = FALSE
    ),
    "Unrecognized `model_names`"
  )
})

test_that("marginal rows keep subgroup grouping labels for the forest renderer", {
  slow()
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
  fit <- suppressWarnings(cams(
    dat,
        reference = "Younger"))
  built <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = NULL,
    include_shrinkage = FALSE,
    include_overall = TRUE,
    include_prediction = TRUE
  )
  marg_ix <- which(built$resultsA$row_kind %in% c("overall", "prediction", "combined"))
  expect_true(length(marg_ix) >= 4L)
  older_ix <- marg_ix[built$row_labels$subheading[marg_ix] == "Older"]
  expect_true(length(older_ix) >= 1L)
  fp <- forestplot(fit)
  expect_s3_class(fp, "interactmeta_forestplot")
  expect_true(!is.null(fp$figure) || !is.null(fp$plot))
  expect_identical(fp$engine, "bayesmeta")
})

test_that("forestplot() includes shrinkage overlays by default", {
  slow()
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
  fp <- forestplot(fit)
  expect_true(any(attr(fp, "panels")[[1]]$row_kind == "shrinkage", na.rm = TRUE))
  fp_off <- forestplot(fit, include_shrinkage = FALSE)
  expect_false(any(attr(fp_off, "panels")[[1]]$row_kind == "shrinkage", na.rm = TRUE))
})

test_that("CAMS forest shrinkage fill is consistent within every trial after render", {
  slow()
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.12, 0.22),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.5, 0.5, 0.7, 0.7),
    ifrac2_Older = c(-0.3, 0.7, -0.5, 0.5, -0.7, 0.3)
  )
  fit <- cams(dat, reference = "Younger")
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
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
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
  fit <- suppressWarnings(cams(
    dat,
        reference = "Younger"))
  fp <- forestplot(
    fit,
    type = "modifiers",
    include_shrinkage = TRUE
  )
  ra <- attr(fp, "panels")[[1]]
  overall_ix <- which(ra$row_kind == "overall")
  pred_ix <- which(ra$row_kind == "prediction" & is.finite(ra$yi))
  study_ix <- which(ra$row_kind == "study")
  shr_ix <- which(ra$row_kind == "shrinkage" & !is.na(ra$shape))
  expect_true(length(overall_ix) >= 2L)
  expect_true(length(pred_ix) >= 2L)
  rl <- fp$specs$row_labels
  expect_true(all(
    ra$fill[overall_ix] == interactmeta:::forestplot_summary_diamond_colours(
      rl$heading[overall_ix],
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
    ra_int$row_kind %in% c("overall", "prediction") &
      rl_int$heading == "Interaction"
  )
  expect_true(length(int_ix) >= 2L)
  expect_true(all(ra_int$fill[int_ix] == interactmeta:::forestplot_interaction_summary_colour()))
})

test_that("include_shrinkage with group=subgroup uses alphabetical trial pairs", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  dat <- data.frame(
    study = rep(c("C", "A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.12, 0.22),
    vi = 0.01,
    n = 10,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.5, 0.5, 0.7, 0.7),
    ifrac2_Older = c(-0.3, 0.7, -0.5, 0.5, -0.7, 0.3)
  )
  fit <- suppressWarnings(cams(
    dat, reference = "Younger",
    ))
  fp <- forestplot(
    fit,
    type = "modifiers",
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
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    n = 10
  )
  fit <- suppressWarnings(cams(
    dat, reference = "Younger",
    ))
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
  slow()
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.01, 0.04, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.6, 0.6),
    ifrac2_Older = c(-0.3, 0.7, -0.6, 0.4)
  )
  fit <- cams(
    dat,
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
  expect_true(any(panels[[1]]$row_kind == "overall"))
  expect_false(grepl("Information fraction \\(\u03c0\\)", attr(fp, "footnote")))

  fp_overall <- forestplot(
    fit,
    exponentiate = TRUE,
    include_shrinkage = TRUE,
    ifrac = "overall"
  )
  expect_false(grepl("Information fraction \\(\u03c0\\)", attr(fp_overall, "footnote")))

  # Overall / Prediction / Combined diamonds use ifrac(..., ifrac =); panels
  # should differ when overall- vs interaction-weighted ω diverge.
  panels_int <- attr(fp, "panels")[[1]]
  panels_ov <- attr(fp_overall, "panels")[[1]]
  ix_sum <- panels_int$row_kind %in% c("overall", "prediction", "combined")
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
  slow()
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
  fit <- cams(
    dat,
        reference = "Younger")
  fp <- forestplot(
    fit,
        include_shrinkage = FALSE,
    include_overall = FALSE,
    include_prediction = FALSE,
    include_combined = FALSE
  )
  foot <- attr(fp, "footnote")
  expect_true(is.character(foot) && nzchar(foot))
  expect_false(grepl("Information fraction \\(\u03c0\\)", foot))
  expect_match(foot, paste0("Heterogeneity (", "\u03c4", ")"), fixed = TRUE)
  expect_false(grepl("Combined IF", foot, fixed = TRUE))
  expect_false(grepl("Combined heterogeneity", foot, fixed = TRUE))
  expect_match(
    foot,
    paste0("Interaction heterogeneity (", "\u03c4", "_", "\u03b3", "):"),
    fixed = TRUE
  )
  expect_false(grepl("Older - Younger", foot, fixed = TRUE))
  expect_length(attr(fp, "annotations")$ifrac, 0L)
  expect_true("combined" %in% names(attr(fp, "annotations")$heterogeneity))
  expect_true("interaction" %in% names(attr(fp, "annotations")$heterogeneity))
  expect_equal(length(attr(fp, "annotations")$lines), 0L)

  fp_off <- forestplot(
    fit,
        show_annotations = FALSE,
    include_shrinkage = FALSE,
    include_overall = FALSE,
    include_prediction = FALSE,
    include_combined = FALSE
  )
  expect_null(attr(fp_off, "footnote"))
})

test_that("all forestplot shrinkage rows match shrinkage() from the fit", {
  slow()
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")

  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.01, 0.04, 0.01),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45),
    n = 10
  )
  fit <- suppressWarnings(cams(
    dat,
        reference = "Younger"))

  sg_me <- shrinkage(fit, type = "modifiers", summary = TRUE)$summary
  int_me <- shrinkage(fit, type = "interaction", summary = TRUE)$summary
  ov_me <- shrinkage(fit, type = "combined", summary = TRUE)$summary

  shr <- interactmeta:::comparison_forest_shrinkage_data(
    fit,
    prob = 0.95,
    max_support = NULL,
    ifrac = "overall"
  )
  # IV mix of subgroup shrinkages differs when subgroup vis differ
  sg <- shr
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

  # type = "modifiers" (heading = subgroup, label = trial)
  fp_sg <- forestplot(fit, type = "modifiers", include_shrinkage = TRUE)
  ra <- attr(fp_sg, "panels")[[1]]
  rl <- forestplot_unwrap_row_labels(fp_sg$specs$row_labels)
  ix_sg <- which(
    ra$row_kind == "shrinkage" &
      !is.na(ra$shape) &
      is.finite(ra$yi) &
      !interactmeta:::forestplot_heading_is_combined(rl$heading)
  )
  expect_equal(length(ix_sg), 4L)
  for (i in ix_sg) {
    hit <- which(
      as.character(sg_me$trial) == as.character(rl$label[[i]]) &
        as.character(sg_me) == as.character(rl$heading[[i]])
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

  expect_error(
    forestplot(fit, type = "combined"),
    "no longer supported",
    fixed = TRUE
  )

  # group = "trial" still uses subgroup shrinkage from the fit
  fp_tr <- forestplot(
    fit,
    type = "modifiers",
    group = "trial",
    include_shrinkage = TRUE
  )
  ra_tr <- attr(fp_tr, "panels")[[1]]
  rl_tr <- forestplot_unwrap_row_labels(fp_tr$specs$row_labels)
  ix_tr <- which(
    ra_tr$row_kind == "shrinkage" &
      !is.na(ra_tr$shape) &
      is.finite(ra_tr$yi) &
      !interactmeta:::forestplot_heading_is_combined(rl_tr$heading)
  )
  expect_equal(length(ix_tr), 4L)
  for (i in ix_tr) {
    hit <- which(
      as.character(sg_me$trial) == as.character(rl_tr$heading[[i]]) &
        as.character(sg_me) == as.character(rl_tr$label[[i]])
    )
    expect_length(hit, 1L)
    expect_equal(ra_tr$yi[[i]], sg_me$estimate[[hit]], tolerance = 1e-8)
  }
})
