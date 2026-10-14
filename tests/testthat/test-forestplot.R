test_that("forestplot_compose_row_label() groups consecutive headings", {
  heading <- c("ADVANCE", "ADVANCE", "AFFIRM", "AFFIRM")
  label <- c("Younger", "Older", "Younger", "Older")
  flat <- interact:::forestplot_compose_row_label(heading, label)
  expect_equal(flat, paste(heading, label, sep = " / "))

  grouped <- interact:::forestplot_compose_row_label(
    heading, label, group_by_heading = TRUE
  )
  expect_match(grouped[[1]], "^ADVANCE\\s+Younger$")
  expect_match(grouped[[2]], "^\\s+Older$")

  sg_first <- interact:::forestplot_compose_row_label(
    c("Younger", "Younger", "Older", "Older"),
    c("ADVANCE", "AFFIRM", "ADVANCE", "AFFIRM"),
    group_by_heading = TRUE
  )
  expect_match(sg_first[[1]], "^Younger\\s+ADVANCE$")
  expect_match(sg_first[[2]], "^\\s+AFFIRM$")
})

test_that("forestplot_group_headings() blanks repeats within a block", {
  expect_equal(
    interact:::forestplot_group_headings(c("A", "A", "B", "B", "B")),
    c("A", "", "B", "", "")
  )
  expect_equal(interact:::forestplot_group_headings(character()), character())
  expect_equal(interact:::forestplot_group_headings("Solo"), "Solo")
})

test_that("forestplot() exposes CAMS formals for completion", {
  fmls <- names(formals(forestplot))
  expect_true(all(
    c("type", "include_shrinkage", "include_overall", "include_prediction",
      "include_combined", "estimates_only", "exponentiate") %in% fmls
  ))
  expect_false("printplot" %in% fmls)
  expect_true(isTRUE(eval(formals(forestplot)$include_shrinkage)))
  expect_equal(
    eval(formals(forestplot)$type),
    c("subgroups", "interaction", "combined")
  )
})

test_that("forestplot() forwards non-CAMS objects to forestplot::forestplot", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  e <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  # escalc → bayesmeta::forestplot.escalc (should not hit the CAMS path)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(forestplot(e))
})

test_that("resolve_forestplot_label_col() defaults and aliases", {
  expect_equal(interact:::resolve_forestplot_label_col(), "trial")
  expect_equal(interact:::resolve_forestplot_label_col(label = "study"), "study")
  expect_equal(interact:::resolve_forestplot_label_col(trial = "study"), "study")
  expect_null(interact:::resolve_forestplot_subgroup_col())
  expect_equal(
    interact:::resolve_forestplot_subgroup_col(NULL, data = edss),
    "subgroup"
  )
  expect_null(interact:::resolve_forestplot_subgroup_col(FALSE, data = edss))
  expect_equal(interact:::resolve_forestplot_subgroup_col("arm"), "arm")
})

test_that("resolve_forestplot_favours() defaults", {
  expect_equal(
    interact:::resolve_forestplot_favours(NULL),
    c("Control", "   Experimental")
  )
  expect_equal(
    interact:::resolve_forestplot_favours(NULL, type = "subgroups"),
    c("Control", "   Experimental")
  )
  expect_null(interact:::resolve_forestplot_favours(FALSE))
  expect_equal(
    interact:::resolve_forestplot_favours(c("A", "B")),
    c("A", "B")
  )
  expect_error(interact:::resolve_forestplot_favours("only-one"))
  # Interaction defaults name reference (left) and non-reference (right).
  expect_equal(
    interact:::resolve_forestplot_favours(
      NULL,
      type = "interaction",
      reference = "Younger",
      nonref_subgroup = "Older"
    ),
    c("Younger", "   Older")
  )
  expect_equal(
    interact:::forestplot_interaction_favours("Younger", "Older"),
    c("Younger", "Older")
  )
  # Fall back to Control / Experimental if subgroup labels are missing.
  expect_equal(
    interact:::resolve_forestplot_favours(
      NULL,
      type = "interaction",
      reference = "Younger",
      nonref_subgroup = ""
    ),
    c("Control", "   Experimental")
  )
  expect_equal(
    interact:::resolve_forestplot_favours(
      NULL,
      type = "combined",
      reference = "Younger",
      nonref_subgroup = ""
    ),
    c("Control", "   Experimental")
  )
})

test_that("layout reserves Favours rows under Interaction and Combined prediction", {
  panel <- data.frame(
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.08, 0.12, 0.1),
    vi = 0.01,
    ci.lb = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.08, 0.12, 0.1) - 0.2,
    ci.ub = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.08, 0.12, 0.1) + 0.2,
    shape = c(15L, 15L, 18L, 18L, 18L, 18L, 18L, 18L),
    colour = "black",
    row_kind = c(
      "study", "study", "overall", "prediction",
      "overall", "prediction", "combined", "combined"
    ),
    stringsAsFactors = FALSE
  )
  row_labels <- data.frame(
    heading = c(
      "Younger", "Older", "Younger", "Younger",
      "Interaction", "Interaction", "Combined", "Combined"
    ),
    label = c(
      "A", "B", "Overall", "Prediction",
      "Overall", "Prediction", "Overall", "Prediction"
    ),
    stringsAsFactors = FALSE
  )
  lay <- interact:::forestplot_metafor_subgroup_layout(
    panel,
    row_labels,
    favours = c("Control", "Experimental")
  )
  expect_equal(length(lay$favours_rows), 2L)
  expect_true(all(is.finite(lay$favours_rows)))

  poly_rows <- vapply(lay$polys, `[[`, numeric(1), "row")
  poly_labs <- vapply(lay$polys, `[[`, character(1), "mlab")
  int_pred <- poly_rows[poly_labs == "Prediction (Interaction)"]
  ov_pred <- poly_rows[poly_labs == "Combined: Prediction"]
  expect_equal(length(int_pred), 1L)
  expect_equal(length(ov_pred), 1L)
  # Favours sit just below the prediction diamond (after the 0.25 gap).
  expect_equal(lay$favours_rows[[1L]], int_pred[[1L]] - 1.25)
  expect_equal(lay$favours_rows[[2L]], ov_pred[[1L]] - 1.25)

  lay_off <- interact:::forestplot_metafor_subgroup_layout(
    panel,
    row_labels,
    favours = NULL
  )
  expect_equal(length(lay_off$favours_rows), 0L)

  expect_equal(
    interact:::forestplot_metafor_xlab("Effect", NULL),
    "Effect"
  )
  expect_equal(
    interact:::forestplot_metafor_xlab("", c("L", "R")),
    ""
  )
  expect_equal(
    interact:::forestplot_metafor_xlab("Effect", c("L", "R")),
    "Effect"
  )
})

test_that("layout places Favours under last subgroup prediction when Int/Combined absent", {
  panel <- data.frame(
    yi = c(0.1, 0.15, 0.25, 0.2, 0.3),
    vi = 0.01,
    ci.lb = c(0.1, 0.15, 0.25, 0.2, 0.3) - 0.2,
    ci.ub = c(0.1, 0.15, 0.25, 0.2, 0.3) + 0.2,
    shape = 15L,
    colour = "black",
    row_kind = c("study", "overall", "prediction", "overall", "prediction"),
    stringsAsFactors = FALSE
  )
  row_labels <- data.frame(
    heading = c("Younger", "Younger", "Younger", "Older", "Older"),
    label = c("A", "Overall", "Prediction", "Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  lay <- interact:::forestplot_metafor_subgroup_layout(
    panel,
    row_labels,
    favours = c("Control", "Experimental")
  )
  expect_equal(length(lay$favours_rows), 1L)
  poly_rows <- vapply(lay$polys, `[[`, numeric(1), "row")
  poly_labs <- vapply(lay$polys, `[[`, character(1), "mlab")
  older_pred <- poly_rows[poly_labs == "Prediction (Older)"]
  expect_equal(length(older_pred), 1L)
  expect_equal(lay$favours_rows[[1L]], older_pred[[1L]] - 1.25)
})

test_that("forestplot_apply_default_colours() uses blue summaries", {
  ra <- data.frame(
    yi = c(0.1, 0.2),
    shape = c(15L, 18L),
    colour = c("red", "red"),
    fill = c("red", "red"),
    row_kind = c("study", "overall"),
    stringsAsFactors = FALSE
  )
  rb <- ra
  rb$shape[[2L]] <- NA_integer_
  st <- interact:::forestplot_apply_default_colours(ra, rb)
  expect_equal(st$resultsA$shape[[1L]], 15L)
  expect_equal(
    st$resultsA$colour[[1L]],
    interact:::comparison_forest_marker_colour()
  )
  expect_equal(
    st$resultsA$colour[[2L]],
    interact:::forestplot_default_summary_colour()
  )
  expect_equal(
    st$resultsA$fill[[2L]],
    interact:::forestplot_default_summary_colour()
  )
})

test_that("forestplot_style_colours() uses default summary blue", {
  cols <- interact:::forestplot_style_colours()
  expect_equal(cols$summary, interact:::forestplot_default_summary_colour())
  expect_equal(cols$marker, interact:::comparison_forest_marker_colour())
})

test_that("forestplot_default_efac() is modest diamond height", {
  expect_equal(interact:::forestplot_default_efac(), c(0, 1.5))
  expect_equal(interact:::forestplot_default_poly_efac(), 1.5)
})

test_that("forestplot_default_poly_fonts() uses bold faces", {
  fonts <- interact:::forestplot_default_poly_fonts(family = "")
  expect_equal(unname(fonts), c(2L, 2L))
  expect_equal(names(fonts), c("", ""))
  fonts2 <- interact:::forestplot_default_poly_fonts(family = "sans")
  expect_equal(unname(fonts2), c(2L, 2L))
  expect_equal(names(fonts2), c("sans", "sans"))
})

test_that("comparison_forest_apply_effect_styles() colours by heading", {
  ra <- data.frame(
    yi = c(0.1, 0.2, 0.15, 0.05),
    shape = 18L,
    colour = "black",
    fill = "black",
    bold = FALSE,
    row_kind = c("overall", "overall", "overall", "combined"),
    stringsAsFactors = FALSE
  )
  rb <- ra
  rb$yi <- NA_real_
  rb$shape <- NA_integer_
  rl <- data.frame(
    heading = c("Younger", "Older", "Interaction", "Combined"),
    label = c("Overall", "Overall", "Overall", "Overall")
  )
  out <- interact:::comparison_forest_apply_effect_styles(
    ra, rb, rl, ref = "Younger", marker = interact:::forestplot_default_summary_colour(), bold_summary = TRUE
  )
  expect_equal(out$resultsA$colour[[1L]], interact:::forestplot_reference_summary_colour())
  expect_equal(out$resultsA$colour[[2L]], interact:::forestplot_nonref_summary_colour())
  expect_equal(out$resultsA$colour[[3L]], interact:::forestplot_interaction_summary_colour())
  expect_equal(out$resultsA$colour[[4L]], interact:::forestplot_overall_summary_colour())
  expect_true(all(out$resultsA$bold))
})

test_that("forestplot_summary_diamond_colour() uses blue for all roles", {
  blue <- interact:::forestplot_default_summary_colour()
  expect_equal(interact:::forestplot_summary_diamond_colour("Combined"), blue)
  expect_equal(interact:::forestplot_summary_diamond_colour("Interaction"), blue)
  expect_equal(
    interact:::forestplot_summary_diamond_colour(
      "Younger",
      reference = "Younger",
      subgroup_levels = c("Younger", "Older")
    ),
    blue
  )
  expect_equal(
    interact:::forestplot_summary_diamond_colour(
      "Older",
      reference = "Younger",
      subgroup_levels = c("Younger", "Older")
    ),
    blue
  )
  expect_equal(
    interact:::forestplot_summary_diamond_colour("ADVANCE", reference = "Younger"),
    blue
  )
})

test_that("summary diamond colours match for subgroup and trial nesting keys", {
  ref <- "Younger"
  expect_equal(
    interact:::forestplot_summary_diamond_colours(
      c("Younger", "Older", "Interaction", "Combined"),
      reference = ref
    ),
    interact:::forestplot_summary_diamond_colours(
      c("Younger", "Older", "Interaction", "Combined"),
      reference = ref,
      subgroup_levels = c("A", "B", "Younger", "Older")
    )
  )
})

test_that("forestplot_panel_df() attaches ifrac_display", {
  d <- data.frame(yi = 0.1, vi = 0.01, stringsAsFactors = FALSE)
  out <- interact:::forestplot_panel_df(d, prob = 0.95, ifrac_display = "50.0%")
  expect_equal(out$ifrac, "50.0%")
})

test_that("forestplot_stack_panels() adds overall study and shrinkage", {
  resultsA <- data.frame(
    yi = c(0.1, 0.12, 0.3, 0.28, 0.2, 0.25),
    vi = c(0.04, 0.04, 0.04, 0.04, 0.01, 0.01),
    ci.lb = c(-0.1, -0.1, 0.1, 0.1, 0, 0),
    ci.ub = c(0.3, 0.3, 0.5, 0.5, 0.4, 0.5),
    shape = c(15L, 15L, 15L, 15L, 18L, 18L),
    colour = "black",
    bold = FALSE,
    row_kind = c(
      "study", "shrinkage", "study", "shrinkage", "combined", "combined"
    ),
    stringsAsFactors = FALSE
  )
  resultsB <- resultsA
  resultsB$yi[c(1L, 2L, 3L, 4L)] <- c(0.05, 0.06, NA, NA)
  resultsB$shape[c(3L, 4L)] <- NA_integer_
  row_labels <- data.frame(
    heading = c("Younger", "Younger", "Older", "Older", "Combined", "Combined"),
    label = c("A", "A", "A", "A", "Overall", "Prediction"),
    colour_label = c("Younger", "Younger", "Older", "Older", "Combined", "Combined"),
    stringsAsFactors = FALSE
  )
  stacked <- interact:::forestplot_stack_panels(
    resultsA, resultsB, row_labels, ref = "Younger"
  )
  ov <- interact:::forestplot_heading_is_combined(stacked$row_labels$heading)
  expect_true(any(stacked$resultsA$row_kind[ov] == "study"))
  expect_true(any(stacked$resultsA$row_kind[ov] == "shrinkage"))
  expect_true(any(stacked$resultsA$row_kind[ov] == "combined"))

  filtered <- interact:::forestplot_filter_type(
    stacked$resultsA,
    stacked$resultsB,
    stacked$row_labels,
    type = "combined"
  )
  expect_true(all(interact:::forestplot_heading_is_combined(filtered$row_labels$heading)))
  expect_true(any(filtered$resultsA$row_kind == "study"))
  expect_true(any(filtered$resultsA$row_kind == "shrinkage"))

  # Prefer trial-level overall shrinkage from shrinkage(type = "combined")
  # over IV-combining subgroup shrinkage rows.
  overall_shrinkage <- data.frame(
    trial = "A",
    yi = 0.99,
    vi = 0.05,
    stringsAsFactors = FALSE
  )
  stacked_ov <- interact:::forestplot_stack_panels(
    resultsA,
    resultsB,
    row_labels,
    ref = "Younger",
    overall_shrinkage = overall_shrinkage
  )
  ov2 <- interact:::forestplot_heading_is_combined(stacked_ov$row_labels$heading)
  shr_ix <- which(ov2 & stacked_ov$resultsA$row_kind == "shrinkage")
  expect_equal(length(shr_ix), 1L)
  expect_equal(stacked_ov$resultsA$yi[[shr_ix]], 0.99)
  expect_equal(stacked_ov$resultsA$vi[[shr_ix]], 0.05)
})

test_that("forestplot_nonref_ifrac_display() uses the non-reference IF", {
  expect_equal(
    interact:::forestplot_nonref_ifrac_display(
      c("60.00%", "40.00%"),
      c("Younger", "Older"),
      ref = "Younger"
    ),
    "40.00%"
  )
  expect_equal(
    interact:::forestplot_nonref_subgroup(c("Younger", "Older"), ref = "Younger"),
    "Older"
  )
  expect_equal(
    interact:::forestplot_heading_with_nonref("Combined", "Older"),
    "Combined (Older)"
  )
  expect_equal(
    interact:::forestplot_col_left_headings(
      data.frame(ifrac = "40%", trial_size = "10"),
      nonref_subgroup = "Older"
    ),
    c("IF", "N")
  )
  expect_equal(
    interact:::forestplot_col_left_headings(
      data.frame(ifrac = "40%", trial_size = "10")
    ),
    c("IF", "N")
  )
  expect_equal(
    interact:::forestplot_nonref_ifrac_display(
      c("55.00%"),
      c("High"),
      ref = "Low"
    ),
    "55.00%"
  )
  expect_equal(
    interact:::forestplot_nonref_ifrac_display(character(), character(), ref = "A"),
    ""
  )
})

test_that("forestplot_restrict_ifrac_to_subgroup() keeps only the quoted subgroup", {
  ra <- data.frame(
    ifrac = c("60%", "40%", "60%", "40%", "40%"),
    row_kind = c("study", "study", "overall", "overall", "study"),
    stringsAsFactors = FALSE
  )
  rl <- data.frame(
    heading = c("Younger", "Older", "Younger", "Older", "Interaction"),
    label = c("A", "A", "Overall", "Overall", "A"),
    colour_label = c("Younger", "Older", "Younger", "Older", "Interaction"),
    stringsAsFactors = FALSE
  )
  out <- interact:::forestplot_restrict_ifrac_to_subgroup(ra, rl, "Older")
  expect_equal(out$ifrac, c("", "40%", "", "40%", "40%"))
  # No restriction when subgroup is empty (type = "subgroups").
  out_all <- interact:::forestplot_restrict_ifrac_to_subgroup(ra, rl, "")
  expect_equal(out_all$ifrac, ra$ifrac)
})

test_that("forestplot_stack_panels() Interaction and Combined keep IF and study size", {
  resultsA <- data.frame(
    yi = c(0.1, 0.3, 0.15, 0.05, 0.2, 0.25),
    vi = c(0.04, 0.04, 0.01, 0.02, 0.01, 0.01),
    ci.lb = c(-0.1, 0.1, 0.05, -0.3, 0, 0),
    ci.ub = c(0.3, 0.5, 0.25, 0.4, 0.4, 0.5),
    shape = c(15L, 15L, 18L, 18L, 18L, 18L),
    colour = "black",
    bold = FALSE,
    ifrac = c("60.00%", "40.00%", "", "", "", ""),
    trial_size = c("200", "", "", "", "", ""),
    row_kind = c(
      "study", "study", "overall", "prediction", "combined", "combined"
    ),
    stringsAsFactors = FALSE
  )
  resultsB <- resultsA
  resultsB$yi <- c(0.05, NA, 0.08, 0.02, NA, NA)
  resultsB$shape <- c(15L, NA, 18L, 18L, NA, NA)
  row_labels <- data.frame(
    heading = c("Younger", "Older", "Younger", "Younger", "Combined", "Combined"),
    label = c("A", "A", "Overall", "Prediction", "Overall", "Prediction"),
    colour_label = c("Younger", "Older", "Younger", "Younger", "Combined", "Combined"),
    stringsAsFactors = FALSE
  )
  stacked <- interact:::forestplot_stack_panels(
    resultsA, resultsB, row_labels, ref = "Younger"
  )

  int <- interact:::forestplot_heading_is_interaction(stacked$row_labels$heading) &
    stacked$resultsA$row_kind == "study"
  ov <- interact:::forestplot_heading_is_combined(stacked$row_labels$heading) &
    stacked$resultsA$row_kind == "study"
  expect_true(any(int))
  expect_true(any(ov))
  expect_true(all(nzchar(stacked$resultsA$ifrac[int])))
  expect_true(all(nzchar(stacked$resultsA$trial_size[int])))
  expect_true(all(nzchar(stacked$resultsA$ifrac[ov])))
  expect_true(all(nzchar(stacked$resultsA$trial_size[ov])))
  # IF is the non-reference percentage; group name is on the IF column header.
  expect_equal(unique(stacked$resultsA$ifrac[int]), "40.00%")
  expect_equal(unique(stacked$resultsA$ifrac[ov]), "40.00%")
  expect_equal(unique(stacked$row_labels$heading[int]), "Interaction")
  expect_equal(unique(stacked$row_labels$heading[ov]), "Combined")
  expect_equal(
    interact:::forestplot_resolve_nonref_subgroup(
      stacked$row_labels, stacked$resultsA, ref = "Younger"
    ),
    "Older"
  )
  expect_equal(
    interact:::forestplot_col_left_headings(
      stacked$resultsA,
      nonref_subgroup = "Older"
    )[[1L]],
    "IF"
  )
  expect_equal(unique(stacked$resultsA$trial_size[int]), "200")
  expect_equal(unique(stacked$resultsA$trial_size[ov]), "200")

  # Columns remain available for type = interaction / overall filters;
  # IF header wording stays fixed for stable column width.
  nonref <- interact:::forestplot_resolve_nonref_subgroup(
    stacked$row_labels, stacked$resultsA, ref = "Younger"
  )
  for (tp in c("interaction", "combined")) {
    filtered <- interact:::forestplot_filter_type(
      stacked$resultsA, stacked$resultsB, stacked$row_labels, type = tp
    )
    expect_true("ifrac" %in% names(filtered$resultsA))
    expect_true("trial_size" %in% names(filtered$resultsA))
    expect_equal(
      interact:::forestplot_col_left_names(filtered$resultsA),
      c("ifrac", "trial_size")
    )
    expect_equal(
      interact:::forestplot_col_left_headings(
        filtered$resultsA,
        nonref_subgroup = nonref
      )[[1L]],
      "IF"
    )
    # Filtered panels alone cannot recover the label (no Younger/Older rows).
    expect_equal(
      interact:::forestplot_resolve_nonref_subgroup(
        filtered$row_labels, filtered$resultsA, ref = "Younger"
      ),
      ""
    )
  }
})

test_that("forestplot_stack_panels() Interaction matches Combined aesthetics", {
  resultsA <- data.frame(
    yi = c(0.1, 0.12, 0.3, 0.28, 0.15, 0.05),
    vi = c(0.04, 0.04, 0.04, 0.04, 0.01, 0.02),
    ci.lb = c(-0.1, -0.1, 0.1, 0.1, 0.05, -0.3),
    ci.ub = c(0.3, 0.3, 0.5, 0.5, 0.25, 0.4),
    shape = c(22L, 15L, 22L, 15L, 18L, 18L),
    colour = "black",
    bold = FALSE,
    row_kind = c(
      "study", "shrinkage", "study", "shrinkage", "overall", "prediction"
    ),
    stringsAsFactors = FALSE
  )
  resultsB <- resultsA
  # Interaction contrasts / summaries live on panel B (open-circle shapes historically).
  resultsB$shape <- c(21L, 16L, 21L, 16L, 18L, 18L)
  resultsB$yi[c(3L, 4L)] <- NA_real_
  row_labels <- data.frame(
    heading = c("Younger", "Younger", "Older", "Older", "Younger", "Younger"),
    label = c("A", "A", "A", "A", "Younger", "Younger"),
    colour_label = c("Younger", "Younger", "Older", "Older", "Younger", "Younger"),
    stringsAsFactors = FALSE
  )
  stacked <- interact:::forestplot_stack_panels(
    resultsA, resultsB, row_labels, ref = "Younger"
  )
  int <- interact:::forestplot_heading_is_interaction(stacked$row_labels$heading)
  expect_true(any(int))
  expect_true(any(stacked$resultsA$row_kind[int] == "study"))
  expect_true(any(stacked$resultsA$row_kind[int] == "overall"))
  expect_true(any(stacked$resultsA$row_kind[int] == "prediction"))
  # Filled squares like Combined studies — not open circles (21).
  expect_true(all(stacked$resultsA$shape[int & stacked$resultsA$row_kind == "study"] == 15L))
  expect_equal(
    stacked$row_labels$label[int & stacked$resultsA$row_kind == "overall"],
    "Overall"
  )
  expect_equal(
    stacked$row_labels$label[int & stacked$resultsA$row_kind == "prediction"],
    "Prediction"
  )

  # Same summary chrome as Combined via forestplot::forestplot args.
  panel <- stacked$resultsA[int, , drop = FALSE]
  labs <- stacked$row_labels[int, , drop = FALSE]
  blocks <- interact:::forestplot_bayesmeta_blocks(panel, labs)
  expect_false(isTRUE(blocks[[1]]$summaries[[1]]$bar))
  expect_true(isTRUE(blocks[[1]]$summaries[[2]]$bar))
  expect_equal(blocks[[1]]$summaries[[1]]$label, "Overall")
  expect_equal(blocks[[1]]$summaries[[2]]$label, "Prediction")
})

test_that("forestplot_filter_type() keeps interaction study and shrinkage", {
  resultsA <- data.frame(
    yi = c(0.1, 0.12, 0.2, 0.22),
    shape = 15L,
    row_kind = c("study", "shrinkage", "overall", "prediction"),
    stringsAsFactors = FALSE
  )
  resultsB <- resultsA
  row_labels <- data.frame(
    heading = rep("Interaction", 4L),
    label = c("A", "A", "Overall", "Prediction"),
    stringsAsFactors = FALSE
  )
  out <- interact:::forestplot_filter_type(
    resultsA, resultsB, row_labels, type = "interaction"
  )
  expect_equal(nrow(out$resultsA), 4L)
  expect_true(all(c("study", "shrinkage") %in% out$resultsA$row_kind))
})

test_that("forestplot_filter_type() keeps the requested block", {
  resultsA <- data.frame(
    yi = c(-0.5, -0.4, -0.35, -0.2, -0.15, -0.7, -0.65),
    row_kind = c(
      "study", "overall", "prediction",
      "study", "overall",
      "combined", "combined"
    ),
    bold = c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  resultsB <- resultsA
  row_labels <- data.frame(
    heading = c(
      "Younger", "Younger", "Younger",
      "Interaction", "Interaction",
      "Combined", "Combined"
    ),
    label = c(
      "A", "Overall", "Prediction",
      "A", "Overall",
      "Overall", "Prediction"
    ),
    stringsAsFactors = FALSE
  )

  sg <- interact:::forestplot_filter_type(
    resultsA, resultsB, row_labels, type = "subgroups"
  )
  expect_equal(unique(sg$row_labels$heading), "Younger")
  expect_equal(nrow(sg$resultsA), 3L)
  expect_true(all(c("Overall", "Prediction") %in%
    sg$row_labels$label[sg$row_labels$heading == "Younger"]))
  expect_false(any(interact:::forestplot_heading_is_interaction(sg$row_labels$heading)))
  expect_false(any(interact:::forestplot_heading_is_combined(sg$row_labels$heading)))

  int <- interact:::forestplot_filter_type(
    resultsA, resultsB, row_labels, type = "interaction"
  )
  expect_equal(unique(int$row_labels$heading), "Interaction")
  expect_equal(nrow(int$resultsA), 2L)

  ov <- interact:::forestplot_filter_type(
    resultsA, resultsB, row_labels, type = "combined"
  )
  expect_equal(unique(ov$row_labels$heading), "Combined")
  expect_equal(nrow(ov$resultsA), 2L)
  expect_equal(ov$row_labels$label, c("Overall", "Prediction"))
})

test_that("forestplot_annotation_lines() flattens IF and heterogeneity", {
  ann <- list(
    lines = "Extra note",
    ifrac = c(Younger = "IF Younger", Older = "IF Older"),
    subgroup_ifrac = character(),
    heterogeneity = c(
      combined = "Het combined",
      interaction = "Het interaction"
    )
  )
  lines <- interact:::forestplot_annotation_lines(ann)
  expect_equal(
    lines,
    paste(
      c("IF Younger", "IF Older", "Het combined", "Het interaction", "Extra note"),
      collapse = " \u00b7 "
    )
  )
  expect_equal(interact:::forestplot_annotation_lines(NULL), character())

  rows <- interact:::forestplot_annotation_footnote_rows(ann)
  expect_length(rows, 1L)
  expect_equal(
    rows[[1L]],
    c("IF Younger", "IF Older", "Het combined", "Het interaction", "Extra note")
  )
})

test_that("forestplot_filter_annotations_type() scopes notes", {
  ann <- list(
    lines = character(),
    ifrac = c(Younger = "IF Younger", Older = "IF Older"),
    subgroup_ifrac = character(),
    heterogeneity = c(
      interaction = "Interaction het",
      combined = "Combined het"
    )
  )
  sg <- interact:::forestplot_filter_annotations_type(ann, "subgroups")
  expect_true(all(c("Younger", "Older") %in% names(sg$ifrac)))
  expect_true(all(c("interaction", "combined") %in% names(sg$heterogeneity)))

  int <- interact:::forestplot_filter_annotations_type(ann, "interaction")
  expect_equal(length(int$ifrac), 0L)
  expect_equal(names(int$heterogeneity), "interaction")

  ov <- interact:::forestplot_filter_annotations_type(ann, "combined")
  expect_equal(names(ov$heterogeneity), "combined")
})

test_that("forestplot_row_information_fraction() omits interaction / contrast rows", {
  d <- data.frame(
    trial = c("A", "A", "A", "B", "B", "B"),
    subgroup = c("g1", "g2", "g1-g2", "g1", "g2", "g1-g2"),
    yi = c(0.2, 0.1, 0.1, 0.3, 0.15, 0.15),
    vi = c(0.01, 0.04, 0.05, 0.0225, 0.09, 0.1),
    stringsAsFactors = FALSE
  )
  out <- interact:::forestplot_row_information_fraction(d, contrast_sep = "-")
  expect_equal(length(out), 6L)
  expect_true(all(nzchar(out[c(1L, 2L, 4L, 5L)])))
  expect_equal(out[c(3L, 6L)], c("", ""))
})

test_that("resolve_forestplot_limits() keeps ratio-scale xlim positive when exponentiate", {
  panels <- list(
    data.frame(yi = -0.12, ci.lb = -0.25, ci.ub = 0.05),
    data.frame(yi = 0, ci.lb = -0.05, ci.ub = 0.1)
  )
  lim <- interact:::resolve_forestplot_limits(panels, exponentiate = TRUE)
  expect_true(all(unlist(lim$xlim) > 0))

  # Imputed missing-subgroup CIs (vi = big_var) exp() underflow to 0.
  huge <- list(
    data.frame(yi = 0, ci.lb = -200, ci.ub = 0.3, vi = 100^2)
  )
  lim_huge <- interact:::resolve_forestplot_limits(huge, exponentiate = TRUE)
  expect_true(all(unlist(lim_huge$xlim) > 0))
  expect_gte(min(unlist(lim_huge$xlim)), 0.05)
})

test_that("forestplot_nice_axis() uses even 0.5-spaced ticks on the additive scale", {
  nice <- interact:::forestplot_nice_axis(c(-1.72, 1.41), exponentiate = FALSE)
  # Symmetric about 0 after outward rounding (half-width 2).
  expect_equal(nice$xlim, c(-2, 2))
  expect_equal(nice$xticks, seq(-2, 2, by = 0.5))
  # No irregular endpoint like -1.8 as a labelled tick.
  expect_false(-1.8 %in% nice$xticks)
})

test_that("forestplot_nice_axis() uses powers of 2 on the ratio scale", {
  nice <- interact:::forestplot_nice_axis(c(0.18, 3.9), exponentiate = TRUE)
  expect_true(1 %in% nice$xticks)
  expect_equal(nice$xticks, c(0.125, 0.25, 0.5, 1, 2, 4, 8))
  expect_equal(nice$xlim, c(0.125, 8))
  # No irregular labels like 0.18 / 0.71.
  expect_false(any(abs(nice$xticks - 0.18) < 1e-8))
  expect_false(any(abs(nice$xticks - 0.71) < 1e-8))
})

test_that("forestplot_nice_axis() keeps 1 as the centred ratio tick", {
  nice <- interact:::forestplot_nice_axis(c(0.4, 1.8), exponentiate = TRUE)
  expect_equal(nice$xticks, c(0.25, 0.5, 1, 2, 4))
  expect_equal(mean(log(nice$xlim)), log(1), tolerance = 1e-12)
})

test_that("explicit xlim and xticks are not recentred onto the null", {
  lim <- interact:::forestplot_apply_nice_axis(
    c(0.07, 0.8),
    xticks = c(0.25, 0.5, 0.75),
    exponentiate = FALSE,
    nullval = 0
  )
  expect_equal(lim$xlim, c(0.07, 0.8))
  expect_equal(lim$xticks, c(0.25, 0.5, 0.75))
})

test_that("resolve_forestplot_limits() rounds outward onto an even tick grid", {
  panels <- list(
    data.frame(
      yi = c(-1.5, 1.2),
      ci.lb = c(-1.72, 0.9),
      ci.ub = c(-1.1, 1.41)
    )
  )
  lim <- interact:::resolve_forestplot_limits(panels, exponentiate = FALSE)
  ticks <- lim$xticks
  xl <- lim$xlim
  expect_equal(diff(ticks), rep(0.5, length(ticks) - 1L))
  expect_true(0 %in% ticks)
  expect_lte(xl[[1L]], -1.72)
  expect_gte(xl[[2L]], 1.41)
})

test_that("forestplot_default_xlab() depends on exponentiate", {
  expect_equal(
    interact:::forestplot_default_xlab(FALSE),
    "Effect (log scale)"
  )
  expect_equal(
    interact:::forestplot_default_xlab(TRUE),
    "Effect"
  )
})

test_that("split_ifrac_for_forest builds contrast rows when absent", {
  d <- ifrac(
    yi = c(0.2, 0.1, 0.3, 0.15),
    vi = c(0.01, 0.04, 0.0225, 0.09),
    trial = c("A", "A", "B", "B"),
    subgroup = c("g1", "g2", "g1", "g2"),
    reference = "g2"
  )
  panels <- interact:::split_ifrac_for_forest(d, contrast_sep = "-")
  expect_equal(nrow(panels$subgroup), 4L)
  expect_equal(nrow(panels$interaction), 2L)
  expect_equal(as.character(panels$interaction$subgroup), c("g1-g2", "g1-g2"))
})

test_that("align_interaction_to_subgroups matches nrow to subgroup rows", {
  sub <- data.frame(
    trial = c("A", "A"),
    subgroup = c("g1", "g2"),
    yi = c(0.2, 0.1),
    vi = c(0.01, 0.04),
    stringsAsFactors = FALSE
  )
  int <- data.frame(
    trial = "A",
    subgroup = "g1-g2",
    yi = 0.1,
    vi = 0.05,
    stringsAsFactors = FALSE
  )
  aligned <- interact:::align_interaction_to_subgroups(sub, int, reference = "g2")
  expect_equal(nrow(aligned$interaction), 2L)
  expect_equal(aligned$interaction$yi, c(0.1, NA_real_))
})

test_that("order_forest_rows() supports subgroup-first nesting", {
  d <- data.frame(
    trial = c("B", "A", "B", "A"),
    subgroup = c("Older", "Older", "Younger", "Younger"),
    yi = 1:4,
    vi = 0.01,
    stringsAsFactors = FALSE
  )
  sg <- interact:::order_forest_rows(d, group = "subgroup")
  expect_equal(as.character(sg$subgroup), c("Older", "Older", "Younger", "Younger"))
  # Trials are alphabetical within each subgroup block
  expect_equal(as.character(sg$trial), c("A", "B", "A", "B"))
  tr <- interact:::order_forest_rows(d, group = "trial")
  expect_equal(as.character(tr$trial), c("B", "B", "A", "A"))
})

test_that("forestplot() fit path supports favours and grouping", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.02, 0.05, 0.03),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))

  fp <- forestplot(fit, exponentiate = TRUE)
  expect_s3_class(fp, "interactmeta_forestplot")
  expect_identical(fp$engine, "bayesmeta")
  expect_equal(length(fp$specs$panels), 1L)
  expect_equal(
    unique(fp$specs$row_labels$heading),
    c("Younger", "Older")
  )
  expect_false(any(interact:::forestplot_heading_is_interaction(
    fp$specs$row_labels$heading
  )))
  expect_false(any(interact:::forestplot_heading_is_combined(
    fp$specs$row_labels$heading
  )))

  fp_ov <- forestplot(fit, type = "combined")
  expect_equal(unique(fp_ov$specs$row_labels$heading), "Combined")
  expect_equal(attr(fp_ov, "type"), "combined")

  fp_sg_only <- forestplot(fit, type = "interaction")
  expect_equal(unique(fp_sg_only$specs$row_labels$heading), "Interaction")

  fp_tr <- forestplot(fit, group = "trial")
  expect_false("Interaction" %in% unique(fp_tr$specs$row_labels$heading))
  expect_true("A" %in% unique(fp_tr$specs$row_labels$heading))
  # Subgroup pooled + prediction use parenthetical labels; alone Younger/Older
  # headers are omitted for those summary-only blocks.
  lay_tr <- interact:::forestplot_metafor_subgroup_layout(
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
  lay_sg <- interact:::forestplot_metafor_subgroup_layout(
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
  lay_fav <- interact:::forestplot_metafor_subgroup_layout(
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
  skip_if_not_installed("bayesmeta"); skip_if_not_installed("forestplot")
  expect_true(exists("forest.interactmeta", envir = asNamespace("interact"), inherits = FALSE))
  expect_true(is.function(interact::forest.interactmeta))
  body_chr <- paste(deparse(body(interact::forest.interactmeta)), collapse = " ")
  expect_match(body_chr, "forestplot")

  # Dispatch via metafor's generic once the method is registered
  meth <- tryCatch(
    utils::getS3method("forest", "interactmeta", optional = TRUE),
    error = function(e) NULL
  )
  if (!is.null(meth)) {
    expect_identical(body(meth), body(interact::forest.interactmeta))
  }

  # Presentation args in ... still land in forest_args (merge path)
  out <- interact:::forestplot_merge_metafor_args(
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
  base <- list(x = 1:2, rows = 3:4, cex = 0.75, efac = 0)
  out <- interact:::forestplot_merge_metafor_args(
    base,
    list(cex = 1, shade = TRUE, rows = 99L, x = 0),
    which = "forest"
  )
  expect_identical(out$x, 1:2)
  expect_identical(out$rows, 3:4)
  expect_identical(out$cex, 1)
  expect_true(isTRUE(out$shade))

  poly <- interact:::forestplot_merge_metafor_args(
    list(x = 0.1, row = 2, mlab = "Overall", col = "black"),
    list(mlab = "nope", col = "red", atransf = exp),
    which = "addpoly"
  )
  expect_identical(poly$mlab, "Overall")
  expect_identical(poly$col, "red")
  expect_identical(poly$atransf, exp)
})

test_that("shrinkage rows omit study name in metafor layout slabs", {
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
  layout <- interact:::forestplot_metafor_subgroup_layout(panel, row_labels)
  expect_equal(
    layout$study$slab,
    c("Trial A", " ", "Trial B", " ")
  )
})

test_that("summary-only layout parks placeholder off diamond rows", {
  panel <- data.frame(
    label = c(
      "Overall", "Prediction", "Overall", "Prediction",
      "Overall", "Prediction", "Overall", "Prediction"
    ),
    heading = c(
      "Younger", "Younger", "Older", "Older",
      "Interaction", "Interaction", "Combined", "Combined"
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
  lay <- interact:::forestplot_metafor_subgroup_layout(panel, row_labels)
  poly_rows <- vapply(lay$polys, `[[`, numeric(1), "row")
  expect_true(isTRUE(lay$study$placeholder))
  expect_false(lay$study$rows %in% poly_rows)
  expect_equal(sum(lay$headers$label == "Combined"), 1L)
})

test_that("forestplot_metafor_textpos() nudges left and right columns", {
  xlim <- c(-5, 5)
  expect_equal(interact:::forestplot_metafor_textpos(xlim, right_shift = 0), xlim)
  nudged <- interact:::forestplot_metafor_textpos(xlim, right_shift = 0.12)
  expect_equal(nudged[[1L]], xlim[[1L]])
  expect_lt(nudged[[2L]], xlim[[2L]])
  both <- interact:::forestplot_metafor_textpos(
    xlim,
    right_shift = 0.12,
    left_shift = 0.05
  )
  expect_gt(both[[1L]], xlim[[1L]])
  expect_lt(both[[2L]], xlim[[2L]])
})

test_that("forestplot_metafor_plot_xlim pads for Effect text", {
  alim <- c(-1, 1)
  full <- interact:::forestplot_metafor_plot_xlim(
    alim,
    has_ilab = FALSE,
    annotate = TRUE
  )
  expect_lt(full[[1L]], alim[[1L]])
  expect_gt(full[[2L]], alim[[2L]])
  tp <- interact:::forestplot_metafor_textpos(full, right_shift = 0.06)
  expect_lt(tp[[2L]], full[[2L]])
  expect_gt(tp[[2L]], alim[[2L]])
})

test_that("forestplot left columns follow IF / N when present", {
  with_cols <- data.frame(
    yi = 0,
    ifrac = "50.0%",
    trial_size = "100",
    stringsAsFactors = FALSE
  )
  expect_equal(
    interact:::forestplot_col_left_names(with_cols),
    c("ifrac", "trial_size")
  )
  expect_equal(
    interact:::forestplot_col_left_headings(with_cols),
    c("IF", "N")
  )
  cleared <- with_cols[, setdiff(names(with_cols), c("ifrac", "trial_size")), drop = FALSE]
  expect_null(interact:::forestplot_col_left_names(cleared))
  expect_equal(interact:::forestplot_col_left_headings(cleared), character())
})

test_that("forestplot omits N unless n_var is set", {
  dat <- data.frame(
    trial = c("A", "A"),
    subgroup = c("Younger", "Older"),
    n = c(40, 60),
    stringsAsFactors = FALSE
  )
  expect_null(interact:::forestplot_n_values(dat))
  expect_null(interact:::forestplot_n_values(dat, n_var = FALSE))
  expect_equal(interact:::forestplot_n_values(dat, n_var = "n"), c(40, 60))
  panel <- data.frame(yi = c(0.1, 0.2))
  omitted <- interact:::forestplot_attach_trial_size(panel, dat)
  expect_false("trial_size" %in% names(omitted))
  shown <- interact:::forestplot_attach_trial_size(panel, dat, n_var = "n")
  expect_equal(shown$trial_size, c("100", ""))
})

test_that("estimates_only omits IF and N columns", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cam(
    es,
    engine = "bayesmeta",
    reference = "Younger"))
  fp_full <- suppressWarnings(forestplot(fit, exponentiate = TRUE))
  fp_est <- suppressWarnings(
    forestplot(fit, exponentiate = TRUE, estimates_only = TRUE)
  )
  expect_equal(
    interact:::forestplot_col_left_names(fp_full$specs$panels[[1L]]),
    "ifrac"
  )
  expect_null(interact:::forestplot_col_left_names(fp_est$specs$panels[[1L]]))
  expect_false("ifrac" %in% names(fp_est$specs$panels[[1L]]))
  expect_false("trial_size" %in% names(fp_est$specs$panels[[1L]]))
  expect_true(all(fp_est$specs$panels[[1L]]$row_kind != "study", na.rm = TRUE))
})

test_that("estimates_only pooled/prediction forest prints", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cam(
    es,
    engine = "metafor",
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
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  fit <- suppressWarnings(cam(
    es,
    engine = "bayesmeta",
    reference = "Younger"))

  fps <- lapply(
    c("subgroups", "interaction", "combined"),
    function(tp) {
      suppressWarnings(forestplot(fit, type = tp, exponentiate = TRUE))
    }
  )
  base_sizes <- vapply(fps, function(fp) as.numeric(attr(fp, "base_size")), numeric(1))
  expect_true(length(unique(base_sizes)) == 1L)
  expect_equal(base_sizes[[1L]], as.numeric(interact:::interactmeta_forest_base_size()))

  favours <- lapply(fps, function(fp) fp$specs$favours)
  expect_equal(favours[[1L]], c("Control", "   Experimental"))
  expect_equal(favours[[2L]], c("Younger", "   Older"))
  expect_equal(favours[[3L]], c("Control", "   Experimental"))

  headings <- lapply(fps, function(fp) {
    lh <- fp$specs$left_headings
    if (is.null(lh)) character() else as.character(lh)
  })
  expect_true(all(vapply(headings, function(h) {
    length(h) >= 1L && identical(h[[1L]], "IF")
  }, logical(1))))

  # Auto typography stays off unless requested.
  expect_false(isTRUE(attr(fps[[1L]], "typography_auto")))
})

test_that("forestplot_panel_xlab() defaults by role", {
  expect_equal(
    interact:::forestplot_panel_xlab("subgroups", exponentiate = FALSE),
    "Effect (log scale)"
  )
  expect_equal(
    interact:::forestplot_panel_xlab("interaction", exponentiate = TRUE),
    "Interaction"
  )
})

test_that("forest footnotes convert Greek to plotmath", {
  line <- paste(
    "Information fraction (\u03c0): 0.40 (0.40, 0.40)",
    "Heterogeneity (\u03c4): 0.30 (0.13, 0.66)",
    "Interaction heterogeneity (\u03c4_\u03b3): 0.10 (0.00, 0.26)",
    sep = " \u00b7 "
  )
  expr <- interact:::forestplot_footnote_as_plotmath(line)
  d <- paste(deparse(expr), collapse = "")
  expect_true(is.language(expr))
  expect_match(d, "pi", fixed = TRUE)
  expect_match(d, "tau[gamma]", fixed = TRUE)
  expect_match(d, "tau", fixed = TRUE)
  expect_false(grepl("\u03c4_y", d, fixed = TRUE))
  expect_false(grepl("_", d, fixed = TRUE))
})
