test_that("forestplot_compose_row_label() groups consecutive headings", {
  heading <- c("ADVANCE", "ADVANCE", "AFFIRM", "AFFIRM")
  label <- c("Younger", "Older", "Younger", "Older")
  flat <- interactmeta:::forestplot_compose_row_label(heading, label)
  expect_equal(flat, paste(heading, label, sep = " / "))

  grouped <- interactmeta:::forestplot_compose_row_label(
    heading, label, group_by_heading = TRUE
  )
  expect_match(grouped[[1]], "^ADVANCE\\s+Younger$")
  expect_match(grouped[[2]], "^\\s+Older$")

  sg_first <- interactmeta:::forestplot_compose_row_label(
    c("Younger", "Younger", "Older", "Older"),
    c("ADVANCE", "AFFIRM", "ADVANCE", "AFFIRM"),
    group_by_heading = TRUE
  )
  expect_match(sg_first[[1]], "^Younger\\s+ADVANCE$")
  expect_match(sg_first[[2]], "^\\s+AFFIRM$")
})

test_that("forestplot_group_headings() blanks repeats within a block", {
  expect_equal(
    interactmeta:::forestplot_group_headings(c("A", "A", "B", "B", "B")),
    c("A", "", "B", "", "")
  )
  expect_equal(interactmeta:::forestplot_group_headings(character()), character())
  expect_equal(interactmeta:::forestplot_group_headings("Solo"), "Solo")
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
    c("modifiers", "interaction", "subgroups")
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
  expect_equal(interactmeta:::resolve_forestplot_label_col(), "trial")
  expect_equal(interactmeta:::resolve_forestplot_label_col(label = "study"), "study")
  expect_equal(interactmeta:::resolve_forestplot_label_col(trial = "study"), "study")
  expect_null(interactmeta:::resolve_forestplot_subgroup_col())
  expect_equal(
    interactmeta:::resolve_forestplot_subgroup_col(NULL, data = edss),
    "subgroup"
  )
  expect_null(interactmeta:::resolve_forestplot_subgroup_col(FALSE, data = edss))
  expect_equal(interactmeta:::resolve_forestplot_subgroup_col("arm"), "arm")
})

test_that("resolve_forestplot_favours() defaults", {
  expect_equal(
    interactmeta:::resolve_forestplot_favours(NULL),
    c("Control", "   Experimental")
  )
  expect_equal(
    interactmeta:::resolve_forestplot_favours(NULL, type = "modifiers"),
    c("Control", "   Experimental")
  )
  expect_null(interactmeta:::resolve_forestplot_favours(FALSE))
  expect_equal(
    interactmeta:::resolve_forestplot_favours(c("A", "B")),
    c("A", "B")
  )
  expect_error(interactmeta:::resolve_forestplot_favours("only-one"))
  # Interaction defaults name reference (left) and non-reference (right).
  expect_equal(
    interactmeta:::resolve_forestplot_favours(
      NULL,
      type = "interaction",
      reference = "Younger",
      nonref_subgroup = "Older"
    ),
    c("Younger", "   Older")
  )
  expect_equal(
    interactmeta:::forestplot_interaction_favours("Younger", "Older"),
    c("Younger", "Older")
  )
  # Fall back to Control / Experimental if subgroup labels are missing.
  expect_equal(
    interactmeta:::resolve_forestplot_favours(
      NULL,
      type = "interaction",
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
      "Interaction", "Interaction", "All", "All"
    ),
    label = c(
      "A", "B", "Overall", "Prediction",
      "Overall", "Prediction", "Overall", "Prediction"
    ),
    stringsAsFactors = FALSE
  )
  lay <- interactmeta:::forestplot_metafor_subgroup_layout(
    panel,
    row_labels,
    favours = c("Control", "Experimental")
  )
  expect_equal(length(lay$favours_rows), 2L)
  expect_true(all(is.finite(lay$favours_rows)))

  poly_rows <- vapply(lay$polys, `[[`, numeric(1), "row")
  poly_labs <- vapply(lay$polys, `[[`, character(1), "mlab")
  int_pred <- poly_rows[poly_labs == "Prediction (Interaction)"]
  ov_pred <- poly_rows[poly_labs == "All: Prediction"]
  expect_equal(length(int_pred), 1L)
  expect_equal(length(ov_pred), 1L)
  # Favours sit just below the prediction diamond (after the 0.25 gap).
  expect_equal(lay$favours_rows[[1L]], int_pred[[1L]] - 1.25)
  expect_equal(lay$favours_rows[[2L]], ov_pred[[1L]] - 1.25)

  lay_off <- interactmeta:::forestplot_metafor_subgroup_layout(
    panel,
    row_labels,
    favours = NULL
  )
  expect_equal(length(lay_off$favours_rows), 0L)

  expect_equal(
    interactmeta:::forestplot_metafor_xlab("Effect", NULL),
    "Effect"
  )
  expect_equal(
    interactmeta:::forestplot_metafor_xlab("", c("L", "R")),
    ""
  )
  expect_equal(
    interactmeta:::forestplot_metafor_xlab("Effect", c("L", "R")),
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
  lay <- interactmeta:::forestplot_metafor_subgroup_layout(
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
  st <- interactmeta:::forestplot_apply_default_colours(ra, rb)
  expect_equal(st$resultsA$shape[[1L]], 15L)
  expect_equal(
    st$resultsA$colour[[1L]],
    interactmeta:::comparison_forest_marker_colour()
  )
  expect_equal(
    st$resultsA$colour[[2L]],
    interactmeta:::forestplot_default_summary_colour()
  )
  expect_equal(
    st$resultsA$fill[[2L]],
    interactmeta:::forestplot_default_summary_colour()
  )
})

test_that("forestplot_style_colours() uses default summary blue", {
  cols <- interactmeta:::forestplot_style_colours()
  expect_equal(cols$summary, interactmeta:::forestplot_default_summary_colour())
  expect_equal(cols$marker, interactmeta:::comparison_forest_marker_colour())
})

test_that("forestplot_default_efac() is modest diamond height", {
  expect_equal(interactmeta:::forestplot_default_efac(), c(0, 1.5))
  expect_equal(interactmeta:::forestplot_default_poly_efac(), 1.5)
})

test_that("forestplot_default_poly_fonts() uses bold faces", {
  fonts <- interactmeta:::forestplot_default_poly_fonts(family = "")
  expect_equal(unname(fonts), c(2L, 2L))
  expect_equal(names(fonts), c("", ""))
  fonts2 <- interactmeta:::forestplot_default_poly_fonts(family = "sans")
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
    heading = c("Younger", "Older", "Interaction", "All"),
    label = c("Overall", "Overall", "Overall", "Overall")
  )
  out <- interactmeta:::comparison_forest_apply_effect_styles(
    ra, rb, rl, ref = "Younger", marker = interactmeta:::forestplot_default_summary_colour(), bold_summary = TRUE
  )
  expect_equal(out$resultsA$colour[[1L]], interactmeta:::forestplot_reference_summary_colour())
  expect_equal(out$resultsA$colour[[2L]], interactmeta:::forestplot_nonref_summary_colour())
  expect_equal(out$resultsA$colour[[3L]], interactmeta:::forestplot_interaction_summary_colour())
  expect_equal(out$resultsA$colour[[4L]], interactmeta:::forestplot_overall_summary_colour())
  expect_true(all(out$resultsA$bold))
})

test_that("forestplot_summary_diamond_colour() uses blue for all roles", {
  blue <- interactmeta:::forestplot_default_summary_colour()
  expect_equal(interactmeta:::forestplot_summary_diamond_colour("All"), blue)
  expect_equal(interactmeta:::forestplot_summary_diamond_colour("Interaction"), blue)
  expect_equal(
    interactmeta:::forestplot_summary_diamond_colour(
      "Younger",
      reference = "Younger",
      subgroup_levels = c("Younger", "Older")
    ),
    blue
  )
  expect_equal(
    interactmeta:::forestplot_summary_diamond_colour(
      "Older",
      reference = "Younger",
      subgroup_levels = c("Younger", "Older")
    ),
    blue
  )
  expect_equal(
    interactmeta:::forestplot_summary_diamond_colour("ADVANCE", reference = "Younger"),
    blue
  )
})

test_that("summary diamond colours match for subgroup and trial nesting keys", {
  ref <- "Younger"
  expect_equal(
    interactmeta:::forestplot_summary_diamond_colours(
      c("Younger", "Older", "Interaction", "All"),
      reference = ref
    ),
    interactmeta:::forestplot_summary_diamond_colours(
      c("Younger", "Older", "Interaction", "All"),
      reference = ref,
      subgroup_levels = c("A", "B", "Younger", "Older")
    )
  )
})

test_that("forestplot_panel_df() attaches ifrac_display", {
  d <- data.frame(yi = 0.1, vi = 0.01, stringsAsFactors = FALSE)
  out <- interactmeta:::forestplot_panel_df(d, prob = 0.95, ifrac_display = "50.0%")
  expect_equal(out$ifrac, "50.0%")
})

test_that("forestplot_stack_panels() adds Combined summaries without trial rows", {
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
    heading = c("Younger", "Younger", "Older", "Older", "All", "All"),
    label = c("A", "A", "A", "A", "Overall", "Prediction"),
    colour_label = c("Younger", "Younger", "Older", "Older", "All", "All"),
    stringsAsFactors = FALSE
  )
  stacked <- interactmeta:::forestplot_stack_panels(
    resultsA, resultsB, row_labels, ref = "Younger"
  )
  ov <- interactmeta:::forestplot_heading_is_combined(stacked$row_labels$heading)
  expect_false(any(stacked$resultsA$row_kind[ov] == "study"))
  expect_false(any(stacked$resultsA$row_kind[ov] == "shrinkage"))
  expect_true(all(stacked$resultsA$row_kind[ov] == "combined"))
  expect_equal(
    sort(unique(stacked$row_labels$label[ov])),
    c("Overall", "Prediction")
  )

  all_only <- interactmeta:::forestplot_filter_type(
    stacked$resultsA,
    stacked$resultsB,
    stacked$row_labels,
    type = "modifiers",
    include_combined = TRUE
  )
  expect_true(any(interactmeta:::forestplot_heading_is_combined(all_only$row_labels$heading)))
  expect_true(any(all_only$resultsA$row_kind == "study"))
})

test_that("forestplot_nonref_ifrac_display() uses the non-reference IF", {
  expect_equal(
    interactmeta:::forestplot_nonref_ifrac_display(
      c("60.00%", "40.00%"),
      c("Younger", "Older"),
      ref = "Younger"
    ),
    "40.00%"
  )
  expect_equal(
    interactmeta:::forestplot_nonref_subgroup(c("Younger", "Older"), ref = "Younger"),
    "Older"
  )
  expect_equal(
    interactmeta:::forestplot_heading_with_nonref("All", "Older"),
    "All (Older)"
  )
  expect_equal(
    interactmeta:::forestplot_col_left_headings(
      data.frame(ifrac = "40%", trial_size = "10"),
      nonref_subgroup = "Older"
    ),
    c("IF (Older)", "N")
  )
  expect_equal(
    interactmeta:::forestplot_col_left_headings(
      data.frame(ifrac = "40%", trial_size = "10")
    ),
    c("IF", "N")
  )
  expect_equal(
    interactmeta:::forestplot_nonref_ifrac_display(
      c("55.00%"),
      c("High"),
      ref = "Low"
    ),
    "55.00%"
  )
  expect_equal(
    interactmeta:::forestplot_nonref_ifrac_display(character(), character(), ref = "A"),
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
  out <- interactmeta:::forestplot_restrict_ifrac_to_subgroup(ra, rl, "Older")
  expect_equal(out$ifrac, c("", "40%", "", "40%", "40%"))
  # No restriction when subgroup is empty (type = "modifiers").
  out_all <- interactmeta:::forestplot_restrict_ifrac_to_subgroup(ra, rl, "")
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
    heading = c("Younger", "Older", "Younger", "Younger", "All", "All"),
    label = c("A", "A", "Overall", "Prediction", "Overall", "Prediction"),
    colour_label = c("Younger", "Older", "Younger", "Younger", "All", "All"),
    stringsAsFactors = FALSE
  )
  stacked <- interactmeta:::forestplot_stack_panels(
    resultsA, resultsB, row_labels, ref = "Younger"
  )

  int <- interactmeta:::forestplot_heading_is_interaction(stacked$row_labels$heading) &
    stacked$resultsA$row_kind == "study"
  comb <- interactmeta:::forestplot_heading_is_combined(stacked$row_labels$heading)
  expect_true(any(int))
  expect_true(any(comb))
  expect_false(any(comb & stacked$resultsA$row_kind == "study"))
  expect_true(all(stacked$resultsA$row_kind[comb] == "combined"))
  expect_true(all(nzchar(stacked$resultsA$ifrac[int])))
  expect_true(all(nzchar(stacked$resultsA$trial_size[int])))
  # IF is the non-reference percentage; group name is on the IF column header.
  expect_equal(unique(stacked$resultsA$ifrac[int]), "40.00%")
  expect_equal(unique(stacked$row_labels$heading[int]), "Interaction")
  expect_equal(unique(stacked$row_labels$heading[comb]), "All")
  expect_equal(
    interactmeta:::forestplot_resolve_nonref_subgroup(
      stacked$row_labels, stacked$resultsA, ref = "Younger"
    ),
    "Older"
  )
  expect_equal(
    interactmeta:::forestplot_col_left_headings(
      stacked$resultsA,
      nonref_subgroup = "Older"
    )[[1L]],
    "IF (Older)"
  )
  expect_equal(unique(stacked$resultsA$trial_size[int]), "200")

  # Columns remain available for type = interaction / subgroups filters.
  nonref <- interactmeta:::forestplot_resolve_nonref_subgroup(
    stacked$row_labels, stacked$resultsA, ref = "Younger"
  )
  for (tp in c("interaction", "modifiers")) {
    filtered <- interactmeta:::forestplot_filter_type(
      stacked$resultsA, stacked$resultsB, stacked$row_labels, type = tp
    )
    expect_true("ifrac" %in% names(filtered$resultsA))
    expect_true("trial_size" %in% names(filtered$resultsA))
    expect_equal(
      interactmeta:::forestplot_col_left_names(filtered$resultsA),
      c("ifrac", "trial_size")
    )
    nr_hdr <- if (identical(tp, "interaction")) nonref else ""
    expect_equal(
      interactmeta:::forestplot_col_left_headings(
        filtered$resultsA,
        nonref_subgroup = nr_hdr
      )[[1L]],
      if (identical(tp, "interaction")) "IF (Older)" else "IF"
    )
    if (identical(tp, "interaction")) {
      # Interaction-only panels have no Younger/Older study rows.
      expect_equal(
        interactmeta:::forestplot_resolve_nonref_subgroup(
          filtered$row_labels, filtered$resultsA, ref = "Younger"
        ),
        ""
      )
    }
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
  stacked <- interactmeta:::forestplot_stack_panels(
    resultsA, resultsB, row_labels, ref = "Younger"
  )
  int <- interactmeta:::forestplot_heading_is_interaction(stacked$row_labels$heading)
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
  blocks <- interactmeta:::forestplot_bayesmeta_blocks(panel, labs)
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
  out <- interactmeta:::forestplot_filter_type(
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
      "All", "All"
    ),
    label = c(
      "A", "Overall", "Prediction",
      "A", "Overall",
      "Overall", "Prediction"
    ),
    stringsAsFactors = FALSE
  )

  sg <- interactmeta:::forestplot_filter_type(
    resultsA, resultsB, row_labels, type = "modifiers"
  )
  expect_equal(
    sort(unique(sg$row_labels$heading)),
    c("All", "Younger")
  )
  expect_equal(nrow(sg$resultsA), 5L)
  expect_true(all(c("Overall", "Prediction") %in%
    sg$row_labels$label[sg$row_labels$heading == "Younger"]))
  expect_false(any(interactmeta:::forestplot_heading_is_interaction(sg$row_labels$heading)))
  expect_true(any(interactmeta:::forestplot_heading_is_combined(sg$row_labels$heading)))

  sg_no_comb <- interactmeta:::forestplot_filter_type(
    resultsA, resultsB, row_labels,
    type = "modifiers",
    include_combined = FALSE
  )
  expect_equal(unique(sg_no_comb$row_labels$heading), "Younger")
  expect_false(any(interactmeta:::forestplot_heading_is_combined(sg_no_comb$row_labels$heading)))

  int <- interactmeta:::forestplot_filter_type(
    resultsA, resultsB, row_labels, type = "interaction"
  )
  expect_equal(unique(int$row_labels$heading), "Interaction")
  expect_equal(nrow(int$resultsA), 2L)

  expect_error(
    interactmeta:::forestplot_filter_type(
      resultsA, resultsB, row_labels, type = "combined"
    ),
    "no longer supported",
    fixed = TRUE
  )
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
  lines <- interactmeta:::forestplot_annotation_lines(ann)
  expect_equal(
    lines,
    paste(
      c("IF Younger", "IF Older", "Het combined", "Het interaction", "Extra note"),
      collapse = " \u00b7 "
    )
  )
  expect_equal(interactmeta:::forestplot_annotation_lines(NULL), character())

  rows <- interactmeta:::forestplot_annotation_footnote_rows(ann)
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
  sg <- interactmeta:::forestplot_filter_annotations_type(ann, "modifiers")
  expect_true(all(c("Younger", "Older") %in% names(sg$ifrac)))
  expect_true(all(c("interaction", "combined") %in% names(sg$heterogeneity)))

  int <- interactmeta:::forestplot_filter_annotations_type(ann, "interaction")
  expect_equal(length(int$ifrac), 0L)
  expect_equal(names(int$heterogeneity), "interaction")

  expect_error(
    interactmeta:::forestplot_filter_annotations_type(ann, "combined"),
    "no longer supported",
    fixed = TRUE
  )
})

test_that("forestplot_row_information_fraction() omits interaction / contrast rows", {
  d <- data.frame(
    study = c("A", "A", "A", "B", "B", "B"),
    level = c("g1", "g2", "g1-g2", "g1", "g2", "g1-g2"),
    yi = c(0.2, 0.1, 0.1, 0.3, 0.15, 0.15),
    vi = c(0.01, 0.04, 0.05, 0.0225, 0.09, 0.1),
    stringsAsFactors = FALSE
  )
  out <- interactmeta:::forestplot_row_information_fraction(d, contrast_sep = "-")
  expect_equal(length(out), 6L)
  expect_true(all(nzchar(out[c(1L, 2L, 4L, 5L)])))
  expect_equal(out[c(3L, 6L)], c("", ""))
})

test_that("resolve_forestplot_limits() keeps ratio-scale xlim positive when exponentiate", {
  panels <- list(
    data.frame(yi = -0.12, ci.lb = -0.25, ci.ub = 0.05),
    data.frame(yi = 0, ci.lb = -0.05, ci.ub = 0.1)
  )
  lim <- interactmeta:::resolve_forestplot_limits(panels, exponentiate = TRUE)
  expect_true(all(unlist(lim$xlim) > 0))

  # Imputed missing-subgroup CIs (vi = big_var) exp() underflow to 0.
  huge <- list(
    data.frame(yi = 0, ci.lb = -200, ci.ub = 0.3, vi = 100^2)
  )
  lim_huge <- interactmeta:::resolve_forestplot_limits(huge, exponentiate = TRUE)
  expect_true(all(unlist(lim_huge$xlim) > 0))
  expect_gte(min(unlist(lim_huge$xlim)), 0.05)
})

test_that("forestplot_nice_axis() uses even 0.5-spaced ticks on the additive scale", {
  nice <- interactmeta:::forestplot_nice_axis(c(-1.72, 1.41), exponentiate = FALSE)
  # Symmetric about 0 after outward rounding (half-width 2).
  expect_equal(nice$xlim, c(-2, 2))
  expect_equal(nice$xticks, seq(-2, 2, by = 0.5))
  # No irregular endpoint like -1.8 as a labelled tick.
  expect_false(-1.8 %in% nice$xticks)
})

test_that("forestplot_nice_axis() uses powers of 2 on the ratio scale", {
  nice <- interactmeta:::forestplot_nice_axis(c(0.18, 3.9), exponentiate = TRUE)
  expect_true(1 %in% nice$xticks)
  expect_equal(nice$xticks, c(0.125, 0.25, 0.5, 1, 2, 4, 8))
  expect_equal(nice$xlim, c(0.125, 8))
  # No irregular labels like 0.18 / 0.71.
  expect_false(any(abs(nice$xticks - 0.18) < 1e-8))
  expect_false(any(abs(nice$xticks - 0.71) < 1e-8))
})

test_that("forestplot_nice_axis() keeps 1 as the centred ratio tick", {
  nice <- interactmeta:::forestplot_nice_axis(c(0.4, 1.8), exponentiate = TRUE)
  expect_equal(nice$xticks, c(0.25, 0.5, 1, 2, 4))
  expect_equal(mean(log(nice$xlim)), log(1), tolerance = 1e-12)
})

test_that("explicit xlim and xticks are not recentred onto the null", {
  lim <- interactmeta:::forestplot_apply_nice_axis(
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
  lim <- interactmeta:::resolve_forestplot_limits(panels, exponentiate = FALSE)
  ticks <- lim$xticks
  xl <- lim$xlim
  expect_equal(diff(ticks), rep(0.5, length(ticks) - 1L))
  expect_true(0 %in% ticks)
  expect_lte(xl[[1L]], -1.72)
  expect_gte(xl[[2L]], 1.41)
})

test_that("forestplot_default_xlab() depends on exponentiate", {
  expect_equal(
    interactmeta:::forestplot_default_xlab(FALSE),
    "Effect (log scale)"
  )
  expect_equal(
    interactmeta:::forestplot_default_xlab(TRUE),
    "Effect"
  )
})

test_that("split_ifrac_for_forest builds contrast rows when absent", {
  d <- ifrac(
    yi = c(0.2, 0.1, 0.3, 0.15),
    vi = c(0.01, 0.04, 0.0225, 0.09),
    study = c("A", "A", "B", "B"),
    level = c("g1", "g2", "g1", "g2"),
    reference = "g2"
  )
  panels <- interactmeta:::split_ifrac_for_forest(d, contrast_sep = "-")
  expect_equal(nrow(panels), 4L)
  expect_equal(nrow(panels$interaction), 2L)
  expect_equal(as.character(panels$interaction), c("g1-g2", "g1-g2"))
})

test_that("align_interaction_to_subgroups matches nrow to subgroup rows", {
  sub <- data.frame(
    study = c("A", "A"),
    level = c("g1", "g2"),
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
  aligned <- interactmeta:::align_interaction_to_subgroups(sub, int, reference = "g2")
  expect_equal(nrow(aligned$interaction), 2L)
  expect_equal(aligned$interaction$yi, c(0.1, NA_real_))
})

test_that("order_forest_rows() supports subgroup-first nesting", {
  d <- data.frame(
    study = c("B", "A", "B", "A"),
    level = c("Older", "Older", "Younger", "Younger"),
    yi = 1:4,
    vi = 0.01,
    stringsAsFactors = FALSE
  )
  sg <- interactmeta:::order_forest_rows(d, group = "subgroup")
  expect_equal(as.character(sg), c("Older", "Older", "Younger", "Younger"))
  # Trials are alphabetical within each subgroup block
  expect_equal(as.character(sg$trial), c("A", "B", "A", "B"))
  tr <- interactmeta:::order_forest_rows(d, group = "trial")
  expect_equal(as.character(tr$trial), c("B", "B", "A", "A"))
})
test_that("forestplot_figure_inches scales with rows and footnotes", {
  fp <- structure(
    list(
      method = "REML",
      specs = list(
        panels = list(data.frame(yi = seq_len(10L))),
        row_labels = data.frame(
          heading = rep(c("A", "B"), each = 5L),
          label = as.character(seq_len(10L)),
          stringsAsFactors = FALSE
        ),
        annotations = list(
          lines = character(),
          ifrac = character(),
          heterogeneity = character()
        )
      ),
      figure = structure(list(method = "REML"), class = "interactmeta_forestplot_metafor")
    ),
    class = "interactmeta_forestplot",
    panels = list(data.frame(yi = seq_len(10L)))
  )
  dims <- interactmeta:::forestplot_figure_inches(fp)
  # Full forests: ~0.28 in per layout row + ~2 in chrome.
  expect_match(dims$formula, "0\\.28\\*")
  expect_gte(dims$width, 12)
  expect_gte(
    dims$width,
    interactmeta:::forestplot_default_width_height_ratio() * dims$height - 1e-9
  )
  h <- forestplot_height(fp, width = 8)
  expect_equal(h$width, 8)
  expect_gte(h$height, 3)
  expect_lt(h$height, 0.9 * 12.5)
  expect_true(is.finite(h$n_rows) && h$n_rows >= 10L)

  # A caller-chosen width must not inherit height inflated from the default
  # 12.5 in recommendation (that made vignette forests too tall).
  h10 <- forestplot_height(fp, width = 10)
  expect_equal(h10$width, 10)
  expect_equal(h10$height, h$height)

  # Annotation footnotes are counted for sizing.
  fp$specs$annotations <- list(
    lines = character(),
    ifrac = c(Older = "Information fraction (\u03c0): 0.40 (0.30, 0.50)"),
    heterogeneity = c(
      overall = paste0("Heterogeneity (", "\u03c4", "): 0.10 (0.01, 0.30)"),
      interaction = paste0(
        "Interaction heterogeneity (",
        "\u03c4",
        "_",
        "\u03b3",
        "): 0.05 (0.01, 0.20)"
      )
    )
  )
  dims_ann <- interactmeta:::forestplot_figure_inches(fp)
  expect_gte(dims_ann$n_footnotes, 1L)
  expect_gt(dims_ann$n_footnotes, dims$n_footnotes)
})

test_that("forestplot_height scales for a large multi-subgroup forest", {
  # ~4 blocks × 7 studies + pooled/prediction ≈ age-like default forest.
  n <- 36L
  headings <- rep(c("Younger", "Older", "Interaction", "All"), each = 9L)
  fp <- structure(
    list(
      method = "REML",
      specs = list(
        panels = list(data.frame(
          yi = seq_len(n),
          ci.lb = seq_len(n) - 0.2,
          ci.ub = seq_len(n) + 0.2,
          row_kind = rep(c(rep("study", 7L), "overall", "prediction"), 4L)
        )),
        row_labels = data.frame(
          heading = headings,
          label = paste0("row", seq_len(n)),
          stringsAsFactors = FALSE
        ),
        estimate_col = "yi"
      ),
      figure = structure(list(method = "REML"), class = "interactmeta_forestplot_metafor")
    ),
    class = "interactmeta_forestplot",
    panels = list(data.frame(
      yi = seq_len(n),
      ci.lb = seq_len(n) - 0.2,
      ci.ub = seq_len(n) + 0.2,
      row_kind = rep(c(rep("study", 7L), "overall", "prediction"), 4L)
    ))
  )
  dims <- forestplot_height(fp, width = 10)
  expect_equal(dims$width, 10)
  expect_gte(dims$n_rows, n)
  expect_true(!is.null(dims$formula))
  expect_gte(dims$height, 0.28 * (dims$n_rows + 1L) + 2)
  expect_lte(dims$height, 24)
  # Fixed fig.height = 6 would be far too short for this forest.
  expect_gt(dims$height, 8)

  auto <- interactmeta:::forestplot_figure_inches(fp)
  expect_gte(
    auto$width,
    interactmeta:::forestplot_default_width_height_ratio() * auto$height - 1e-9
  )
  flat <- interactmeta:::forestplot_figure_inches(fp, width_height = FALSE)
  expect_lte(flat$width, auto$width + 1e-9)
})
test_that("forestplot_metafor_cex() shrinks on a short device", {
  tall <- interactmeta:::forestplot_metafor_cex(9, n_layout_rows = 20, device_height_in = 10)
  short <- interactmeta:::forestplot_metafor_cex(9, n_layout_rows = 20, device_height_in = 3)
  expect_lt(short, tall)
  expect_lte(short, 0.75)
})

test_that("forestplot_display_metrics() scales with device and rows", {
  m_tall <- interactmeta:::forestplot_display_metrics(
    n_rows = 20,
    device_width_in = 12,
    device_height_in = 10,
    has_favours = TRUE
  )
  m_short <- interactmeta:::forestplot_display_metrics(
    n_rows = 20,
    device_width_in = 12,
    device_height_in = 4,
    has_favours = TRUE
  )
  expect_lt(m_short$base_size, m_tall$base_size)
  expect_true(is.finite(m_tall$panel.height))
  expect_true(is.finite(m_short$panel.height))
})

test_that("forestplot_n_layout_rows counts headers", {
  fp <- structure(
    list(
      method = "REML",
      specs = list(
        panels = list(data.frame(yi = 1:4)),
        row_labels = data.frame(
          heading = c("A", "A", "B", "B"),
          label = letters[1:4],
          stringsAsFactors = FALSE
        )
      )
    ),
    class = "interactmeta_forestplot",
    panels = list(data.frame(yi = 1:4))
  )
  # 4 content + 2 unique headings + 1 table header
  expect_equal(interactmeta:::forestplot_n_layout_rows(fp), 7L)
})
