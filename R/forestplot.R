#' Forest plot for a fitted CAMS model
#'
#' Draws a forest for a [cam()] fit using [bayesmeta::forestplot.bayesmeta]
#' (single-block) or [forestplot::forestplot()] with the same conventions
#' (multi-block): observed study estimates, optional grey **shrinkage**
#' overlays, and **mean** / **prediction** / **Overall** summaries. An
#' **information fraction** column shows each row's share of trial precision.
#'
#' For non-[interactmeta] input (e.g. [metafor::escalc]), arguments are forwarded
#' to [forestplot::forestplot()], which dispatches to
#' [bayesmeta::forestplot.escalc] when **bayesmeta** is loaded (interactmeta loads it
#' on attach).
#'
#' By default (`type = "subgroups"`) the plot shows **subgroup** blocks
#' (each with **Pooled** and **Prediction**), with overall-IF notes under the
#' x-axis. Use `type = "interaction"` or `type = "overall"` / `"combined"` for those blocks.
#'
#' @param x A CAMS fit from [cam()] (`interactmeta`), or another object handled by
#'   [forestplot::forestplot()] (e.g. `escalc`).
#' @param n_var Column name for per-study patient counts shown in an **N**
#'   column. Default `FALSE` (omit). Set to `"n"` or another column name to
#'   show sizes. When `"n"` is missing but `treat.total` and `cont.total` are
#'   present (2-by-2 event-count data), study size is taken as their sum.
#' @param contrast_sep Separator between subgroup labels on interaction rows
#'   (default `"-"`).
#' @param model_names Optional character vector selecting marginal layers among
#'   `"Shrinkage"`, `"Pooled"`, and `"Prediction"` (intersected with the
#'   `include_*` flags). Prefer `include_shrinkage` / `include_pooled` /
#'   `include_prediction` for new code. When `NULL` or the default `"Model 1"`,
#'   layers follow those flags alone.
#' @param include_shrinkage,include_pooled,include_prediction,include_overall
#'   Include shrinkage, pooled, prediction, and overall marginal rows.
#'   Shrinkage, pooled, prediction, and overall default to `TRUE`. Grey
#'   shrinkage intervals appear alongside quoted estimates (bayesmeta style)
#'   in the active `type` block. Set `include_shrinkage = FALSE` to omit them.
#'   Switching off pooled or prediction keeps blank placeholder rows so the
#'   figure height and spacing stay close to the full plot.
#' @param type Which forest block(s) to draw. One of:
#'   - `"subgroups"` (default): subgroup blocks only (each with **Pooled** and
#'     **Prediction** when `include_pooled` / `include_prediction` are `TRUE`)
#'   - `"interaction"`: Interaction block only (Pooled + Prediction)
#'   - `"overall"` / `"combined"`: Combined-population block only (Pooled + Prediction)
#'
#'   Other features (favours, annotations, `include_shrinkage`, …) apply
#'   within the selected type.
#' @param ifrac How the overall information fraction ω is pooled via [ifrac()]
#'   for **annotations** and for **Pooled / Prediction / Overall** diamonds:
#'   `"interaction"` (default; τ_γ and subgroup sampling variances) or
#'   `"overall"` (τ and design SEs). Shrinkage greys are the trial-specific
#'   CAMS reconstructions from [shrinkage()] (shrunken \eqn{m_j,g_j} inverted
#'   at each trial's own \eqn{\pi_j}); they are not population fitted values.
#'   Observed study rows remain at each trial's own π_j.
#' @param xlab Character; x-axis title. Default is
#'   `"Effect (log scale)"` when `exponentiate = FALSE`, or
#'   `"Effect"` when `exponentiate = TRUE` (ratio-scale ticks via
#'   `xlog = TRUE`).
#' @param group How to nest study rows. `"subgroup"` (default) groups by
#'   subgroup level with trial names as row labels; `"trial"` groups by trial
#'   with subgroup names as row labels.
#' @param show_annotations If `TRUE` (default), annotate under the x-axis with
#'   overall information-fraction quantiles and heterogeneity credibility
#'   intervals.
#' @param annotation_digits Decimal places in annotation text (default `2`).
#' @param favours Optional length-2 character labels placed under the x-axis
#'   (left / right of the null). Default is `c("Control", "Experimental")`
#'   for `type = "subgroups"` / `"overall"`. For `type = "interaction"`, the
#'   default names the reference and non-reference subgroups (same left/right
#'   roles). Set to `FALSE` to suppress.
#' @param show_ci_text Logical; show estimate and CI text. Default `TRUE`.
#' @param exponentiate Logical; if `TRUE`, show estimates on the ratio scale.
#' @param col_right_heading Column heading for the effect column.
#' @param prob Coverage probability for confidence intervals (default `0.95`).
#' @param study_fill, study_colour Optional overrides for study-row colours.
#' @param base_size, pointsize, digits, stroke, fill,
#'   summary_fill, summary_colour, plot.margin, which_panel
#'   Plotting options passed through to the bayesmeta / forestplot renderer.
#' @param estimates_only If `TRUE`, drop observed study rows (keep shrinkage /
#'   summary diamonds) and omit the IF and study-size columns.
#' @param xlim,xticks Optional axis limits and tick positions on the displayed
#'   effect scale. When both are `NULL` (default) and `exponentiate = FALSE`,
#'   limits follow [bayesmeta::forestplot.bayesmeta] / [forestplot::forestplot]
#'   defaults (`clip = c(-Inf, Inf)`, automatic ticks). When
#'   `exponentiate = TRUE`, a strictly positive clip and power-of-two ticks
#'   are chosen automatically (forestplot's log axis cannot include 0). When
#'   `xlim` is supplied, values are rounded outward onto a nice grid (and
#'   `xticks` are filled if omitted).
#' @param plot.margin,panel.height Unused (kept for compatibility).
#' @param grouping Optional vector of subheadings for row labels.
#' @param reference,label,subgroup,trial,contrast,big_var,include_interaction,pooled
#'   Ignored for CAMS fits; retained for back-compatibility.
#' @param ... Further arguments passed through to [forestplot::forestplot()]
#'   when the plot is drawn (same escape hatch as
#'   [bayesmeta::forestplot.bayesmeta]'s `...`). For non-CAMS `x`, unused CAMS
#'   formals are ignored and `...` is forwarded to [forestplot::forestplot()].
#'
#' @return For a CAMS fit, a `interactmeta_forestplot` object with attributes
#'   `"panels"`, `"panel.height"` (when set), and `"base_size"`. Drawing is
#'   deferred to [print.interactmeta_forestplot()] / [forestplot_draw()] (so
#'   `fp <- forestplot(fit)` does not draw; a top-level `forestplot(fit)` does).
#'   Otherwise whatever [forestplot::forestplot()] returns.
#'
#' @section Aesthetics vs content:
#' Marker colours, blue summary diamonds / prediction bars, headers,
#' favours, IF column headings, and default typography (`base_size` /
#' `pointsize`) are fixed CAMS chrome (same for every `type`). Only which
#' blocks appear changes with `type`. Pass `typography_auto = TRUE` in `...`
#' to re-fit fonts to the device and row count. Axis scaling and other
#' [forestplot::forestplot()] knobs default to bayesmeta-compatible
#' behaviour unless you override them (`xlim`, `xticks`, or `...`).
#'
#' @section Related entry points:
#' - [bayesmeta::forestplot.escalc] — study-only forest from [escalc()]
#' - [bayesmeta::forestplot.bayesmeta] — fitted bayesmeta object
#' - [metafor::forest] on a CAMS fit — [forest.interactmeta] delegates here
#' - [forestplot::forestplot] — underlying S3 generic (`forestplot::forestplot(fit)`
#'   still dispatches to the CAMS method)
#'
#' @examples
#' \dontrun{
#' e <- escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
#' forestplot(e)                          # → forestplot.escalc
#' fit <- cam(e, reference = "Lower")
#' forestplot(fit, exponentiate = TRUE, type = "subgroups")
#' forest(fit, exponentiate = TRUE)       # same via metafor::forest
#' }
#' @rdname forestplot
#' @export
forestplot <- function(x,
                       reference = NULL,
                       label = NULL,
                       modifier = NULL,
                       subgroup = NULL,
                       trial = NULL,
                       n_var = FALSE,
                       contrast = "all",
                       contrast_sep = "-",
                       big_var = 100^2,
                       include_interaction = FALSE,
                       pooled = NULL,
                       model_names = "Model 1",
                       include_shrinkage = TRUE,
                       include_pooled = TRUE,
                       include_prediction = TRUE,
                       include_overall = TRUE,
                       ifrac = c("interaction", "overall"),
                       xlab = NULL,
                       group = "subgroup",
                       show_annotations = TRUE,
                       annotation_digits = 2,
                       favours = NULL,
                       show_ci_text = NULL,
                       exponentiate = FALSE,
                       col_right_heading = NULL,
                       stroke = interactmeta_forest_stroke(),
                       digits = 2,
                       fill = "white",
                       study_fill = NULL,
                       study_colour = NULL,
                       summary_fill = NULL,
                       summary_colour = NULL,
                       pointsize = interactmeta_forest_pointsize(),
                       base_size = interactmeta_forest_base_size(),
                       prob = 0.95,
                       plot.margin = interactmeta_forestplot_margin(),
                       xlim = NULL,
                       xticks = NULL,
                       grouping = NULL,
                       panel.height = NULL,
                       estimates_only = FALSE,
                       type = c("subgroups", "interaction", "overall", "combined"),
                       which_panel = NULL,
                       ...) {
  if (is.character(modifier) && length(modifier) == 1L && nzchar(modifier)) {
    subgroup <- modifier
  }
  interactmeta_abort_non_cams(x, "forestplot()")
  if (!inherits(x, "interactmeta")) {
    return(forestplot::forestplot(x, ...))
  }

  check_forestplot_deps()
  # Keep typography fixed across options / row counts unless the user opts in
  # via `typography_auto = TRUE`.
  dots <- list(...)
  typography_auto <- isTRUE(dots$typography_auto)
  dots$typography_auto <- NULL
  panel_height_auto <- missing(panel.height) || is.null(panel.height)
  group <- match.arg(group, c("subgroup", "trial"))
  type <- comparison_resolve_effect_type(match.arg(type))
  if (is.null(show_ci_text)) {
    show_ci_text <- TRUE
  }
  show_ci_text <- isTRUE(show_ci_text)
  sg_pair <- forestplot_fit_reference_nonref(x)
  favours <- resolve_forestplot_favours(
    favours = favours,
    exponentiate = exponentiate,
    type = type,
    reference = sg_pair$reference,
    nonref_subgroup = sg_pair$nonref
  )
  validate_forestplot_args(
    prob = prob,
    which_panel = which_panel,
    xlim = xlim,
    xticks = xticks,
    fill = fill,
    n_var = n_var
  )

  comparison_args <- list(
    fit = x,
    prob = prob,
    contrast_sep = contrast_sep,
    include_shrinkage = include_shrinkage,
    include_pooled = include_pooled,
    include_prediction = include_prediction,
    include_overall = include_overall,
    exponentiate = exponentiate,
    ifrac = match.arg(ifrac),
    n_var = n_var,
    group = group,
    show_annotations = show_annotations,
    annotation_digits = annotation_digits,
    favours = favours,
    show_ci_text = show_ci_text,
    typography_auto = typography_auto,
    xlab = xlab,
    col_right_heading = col_right_heading,
    stroke = stroke,
    digits = digits,
    fill = fill,
    study_fill = study_fill,
    study_colour = study_colour,
    summary_fill = summary_fill,
    summary_colour = summary_colour,
    pointsize = pointsize,
    base_size = base_size,
    plot.margin = plot.margin,
    xlim = xlim,
    xticks = xticks,
    grouping = grouping,
    estimates_only = estimates_only,
    type = type,
    which_panel = which_panel
  )
  if (!identical(model_names, "Model 1")) {
    comparison_args$model_names <- model_names
  }
  if (!panel_height_auto) {
    comparison_args$panel.height <- panel.height
  }
  if (!typography_auto) {
    comparison_args$base_size <- base_size
    comparison_args$pointsize <- pointsize
  }

  # Legacy printplot / quiet: drawing is via print() / forestplot_draw() only.
  # `dots` was captured above (with `typography_auto` already stripped).
  dots$printplot <- NULL
  dots$quiet <- NULL
  for (nm in names(comparison_args)) {
    dots[[nm]] <- NULL
  }
  do.call(comparison_forestplot_from_fit, c(comparison_args, dots))
}

#' S3 method for [forestplot::forestplot] on CAMS fits
#'
#' Thin adapter so `forestplot::forestplot(fit, ...)` reaches the same
#' implementation as [forestplot()]. Prefer [forestplot()] for argument
#' completion in the IDE.
#'
#' @inheritParams forestplot
#' @rdname forestplot
#' @exportS3Method forestplot::forestplot interactmeta
forestplot.interactmeta <- function(x, ...) {
  check_comparison_fit(x)
  do.call(forestplot, c(list(x = x), list(...)))
}

#' Forest plot method for CAMS fits ([metafor::forest] generic)
#'
#' Metafor-idiomatic entry point for a [cam()] fit. Delegates to
#' [forestplot()] with the same arguments so CAMS options (`type`,
#' `group`, …) share one implementation.
#'
#' @param x A CAMS fit from [cam()].
#' @param ... Passed to [forestplot()] (CAMS options and/or forestplot
#'   presentation args).
#'
#' @return A `interactmeta_forestplot` object (see [forestplot()]).
#'
#' @examples
#' \dontrun{
#' fit <- cam(escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss),
#'             reference = "Lower")
#' forest(fit, exponentiate = TRUE)
#' }
#'
#' @seealso [forestplot()], [metafor::forest()]
#' @rdname forestplot
#' @export
#' @exportS3Method metafor::forest interactmeta
forest.interactmeta <- function(x, ...) {
  forestplot(x, ...)
}

#' Render forest panels (shared by fit overlays)
#' @noRd
forestplot_render <- function(resultsA,
                                         resultsB,
                                         row_labels,
                                         n_study,
                                         ref,
                                         prob,
                                         contrast_sep,
                                         xlab = NULL,
                                         col_right_heading = NULL,
                                         summary_fill = NULL,
                                         summary_colour = NULL,
                                         fill = "white",
                                         study_fill = NULL,
                                         study_colour = NULL,
                                         grouping = NULL,
                                         estimates_only = FALSE,
                                         type = c("subgroups", "interaction", "overall", "combined"),
                                         which_panel = NULL,
                                         xlim = NULL,
                                         xticks = NULL,
                                         base_size = interactmeta_forest_base_size(),
                                         pointsize = interactmeta_forest_pointsize(),
                                         digits = 2,
                                         stroke = interactmeta_forest_stroke(),
                                         plot.margin = interactmeta_forestplot_margin(),
                                         panel.height = NULL,
                                         exponentiate = FALSE,
                                         favours = NULL,
                                         show_ci_text = FALSE,
                                         typography_auto = FALSE,
                                         obs_row_index = NULL,
                                         include_shrinkage = TRUE,
                                         include_prediction = TRUE,
                                         annotations = NULL,
                                         group = "subgroup",
                                         ...) {
  dots <- list(...)
  if ("layout" %in% names(dots)) {
    # Former side-by-side layout removed; always stack vertically.
    dots$layout <- NULL
  }
  if ("group" %in% names(dots)) {
    group <- dots$group
    dots$group <- NULL
  }
  group <- match.arg(as.character(group), c("subgroup", "trial"))
  if ("annotations" %in% names(dots)) {
    annotations <- dots$annotations
    dots$annotations <- NULL
  }
  if ("footnote" %in% names(dots)) {
    # Prefer structured annotations; keep a plain footnote string if supplied.
    if (is.null(annotations) && !is.null(dots$footnote) && nzchar(as.character(dots$footnote)[[1L]])) {
      annotations <- list(
        lines = as.character(dots$footnote),
        subgroup_ifrac = character()
      )
    }
    dots$footnote <- NULL
  }
  if ("show_heterogeneity_footnote" %in% names(dots)) {
    dots$show_heterogeneity_footnote <- NULL
  }
  if ("heterogeneity_digits" %in% names(dots)) {
    dots$heterogeneity_digits <- NULL
  }
  if ("panel.height" %in% names(dots)) {
    panel.height <- dots$panel.height
    dots$panel.height <- NULL
  }
  if ("include_shrinkage" %in% names(dots)) {
    include_shrinkage <- dots$include_shrinkage
    dots <- dots[!names(dots) %in% "include_shrinkage"]
  }
  overall_shrinkage <- NULL
  if ("overall_shrinkage" %in% names(dots)) {
    overall_shrinkage <- dots$overall_shrinkage
    dots <- dots[!names(dots) %in% "overall_shrinkage"]
  }
  if (is.null(overall_shrinkage)) {
    overall_shrinkage <- attr(resultsA, "overall_shrinkage", exact = TRUE)
  }
  if ("include_prediction" %in% names(dots)) {
    include_prediction <- dots$include_prediction
    dots <- dots[!names(dots) %in% "include_prediction"]
  }
  if ("include_pooled" %in% names(dots)) {
    dots <- dots[!names(dots) %in% "include_pooled"]
  }
  if ("include_prediction" %in% names(dots)) {
    dots <- dots[!names(dots) %in% "include_prediction"]
  }
  if ("estimates_only" %in% names(dots)) {
    estimates_only <- dots$estimates_only
    dots <- dots[!names(dots) %in% "estimates_only"]
  }
  if ("type" %in% names(dots)) {
    type <- dots$type
    dots <- dots[!names(dots) %in% "type"]
  }
  if ("result" %in% names(dots)) {
    warning("`result` is deprecated; use `type` instead.", call. = FALSE)
    type <- dots$result
    dots <- dots[!names(dots) %in% "result"]
  }
  estimates_only <- isTRUE(estimates_only)
  type <- comparison_resolve_effect_type(as.character(type)[[1L]])
  # Legacy printplot / quiet ignored — draw via print() / forestplot_draw().
  dots$printplot <- NULL
  dots$quiet <- NULL
  if ("exponentiate" %in% names(dots)) {
    exponentiate <- dots$exponentiate
    dots <- dots[!names(dots) %in% "exponentiate"]
  }
  if ("style" %in% names(dots)) {
    # Appearance is fixed to CAMS default chrome; ignore legacy style= args.
    dots <- dots[!names(dots) %in% "style"]
  }
  if ("favours" %in% names(dots)) {
    favours <- dots$favours
    dots <- dots[!names(dots) %in% "favours"]
  }
  if ("show_ci_text" %in% names(dots)) {
    show_ci_text <- dots$show_ci_text
    dots <- dots[!names(dots) %in% "show_ci_text"]
  }
  if ("typography_auto" %in% names(dots)) {
    typography_auto <- dots$typography_auto
    dots <- dots[!names(dots) %in% "typography_auto"]
  }
  exponentiate <- isTRUE(exponentiate)
  if (is.null(show_ci_text)) {
    show_ci_text <- TRUE
  }
  show_ci_text <- isTRUE(show_ci_text)
  # Favours defaults (esp. interaction subgroup labels) are resolved after
  # stacking once `nonref_subgroup` is known.
  typography_auto <- isTRUE(typography_auto)

  render_snapshot <- list(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    n_study = n_study,
    ref = ref,
    prob = prob,
    contrast_sep = contrast_sep,
    xlab = xlab,
    col_right_heading = col_right_heading,
    summary_fill = summary_fill,
    summary_colour = summary_colour,
    fill = fill,
    study_fill = study_fill,
    study_colour = study_colour,
    grouping = grouping,
    estimates_only = estimates_only,
    which_panel = which_panel,
    xlim = xlim,
    xticks = xticks,
    base_size = base_size,
    pointsize = pointsize,
    digits = digits,
    stroke = stroke,
    plot.margin = plot.margin,
    panel.height = panel.height,
    exponentiate = exponentiate,
    favours = favours,
    show_ci_text = show_ci_text,
    typography_auto = typography_auto,
    obs_row_index = obs_row_index,
    include_shrinkage = include_shrinkage,
    include_prediction = include_prediction,
    group = group
  )
  if (length(dots)) {
    render_snapshot <- c(render_snapshot, dots)
  }

  sum_f <- if (!is.null(summary_fill)) summary_fill else "white"
  sum_c <- if (!is.null(summary_colour)) summary_colour else interactmeta_accent_colour()

  if (!"shape" %in% names(resultsA)) {
    resultsA$shape <- ifelse(resultsA$bold, 18, 22)
    resultsB$shape <- ifelse(resultsB$bold, 18, 21)
  }

  xlab_user_supplied <- !(
    is.null(xlab) || !length(xlab) || all(!nzchar(as.character(xlab)))
  )
  if (!xlab_user_supplied) {
    xlab <- forestplot_default_xlab(exponentiate)
  }
  if (is.null(col_right_heading) || !length(col_right_heading)) {
    col_right_heading <- paste0("Effect\n(", 100 * prob, "% interval)")
  }

  row_labels$heading <- row_labels$subheading
  row_labels$subheading <- NA_character_
  row_labels$variable <- paste0(
    apply(row_labels, 1, paste0, collapse = ""),
    "\r",
    seq_len(nrow(row_labels))
  )
  resultsA$variable <- row_labels$variable
  resultsB$variable <- row_labels$variable
  bold <- is.na(row_labels$heading)
  study_rows <- if (!is.null(obs_row_index)) obs_row_index else seq_len(n_study)

  if (!is.null(summary_fill)) {
    resultsB$fill[study_rows] <- resultsA$fill[study_rows]
  } else if (is.list(fill)) {
    if (length(fill) != 2) {
      stop("`fill` must be a list of length 2 when provided as a list.", call. = FALSE)
    }
    resultsA$fill <- fill[[1]]
    resultsB$fill <- fill[[2]]
  } else if (!is.null(fill)) {
    resultsA$fill[study_rows] <- fill
    resultsB$fill[study_rows] <- fill
  }

  if (is.null(study_fill) && is.null(study_colour)) {
    # Studies are forest_study (black); shrinkage is forest_shrinkage (applied below).
    marker <- comparison_forest_marker_colour()
    ix_obs <- study_rows
    resultsA$fill[ix_obs] <- marker
    resultsA$colour[ix_obs] <- marker
    resultsB$fill[ix_obs] <- marker
    resultsB$colour[ix_obs] <- marker
    resultsB$colour[is.na(resultsB$yi)] <- NA_character_
    resultsB$fill[is.na(resultsB$yi)] <- NA_character_
  } else {
    study_f <- if (!is.null(study_fill)) {
      study_fill
    } else {
      comparison_forest_marker_colour()
    }
    study_c <- if (!is.null(study_colour)) {
      study_colour
    } else {
      comparison_forest_marker_colour()
    }
    ix_obs <- study_rows
    resultsA$fill[ix_obs] <- study_f
    resultsA$colour[ix_obs] <- study_c
    resultsB$fill[ix_obs] <- study_f
    resultsB$colour[ix_obs] <- study_c
    resultsB$colour[is.na(resultsB$yi)] <- NA_character_
  }

  layer_on <- interactmeta:::comparison_forest_layer_on
  has_shrinkage_rows <- layer_on(include_shrinkage) &&
    "row_kind" %in% names(resultsA) &&
    any(resultsA$row_kind == "shrinkage", na.rm = TRUE)
  if (has_shrinkage_rows && "colour_label" %in% names(row_labels)) {
    ix_sg_shr <- which(resultsA$row_kind == "shrinkage" & !is.na(resultsA$shape))
    ix_int_shr <- which(resultsB$row_kind == "shrinkage" & !is.na(resultsB$shape))
    styled_shr <- comparison_forest_apply_shrinkage_colours(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      subgroup_shrinkage_row_index = ix_sg_shr,
      interaction_shrinkage_row_index = ix_int_shr,
      ref = ref
    )
    resultsA <- styled_shr$resultsA
    resultsB <- styled_shr$resultsB
  }

  # Default chrome: black study squares. Summary diamond colours are applied
  # after stacking once block headings (Younger / Older / Interaction / Overall)
  # are final.
  styled <- forestplot_apply_default_colours(resultsA, resultsB)
  resultsA <- styled$resultsA
  resultsB <- styled$resultsB

  stacked <- forestplot_stack_panels(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    ref = ref,
    overall_shrinkage = overall_shrinkage
  )
  resultsA <- stacked$resultsA
  resultsB <- stacked$resultsB
  row_labels <- stacked$row_labels
  n_study <- stacked$n_study

  has_marginal_rows <- "row_kind" %in% names(resultsA) &&
    any(resultsA$row_kind %in% c("pooled", "prediction", "overall"), na.rm = TRUE)
  if (has_marginal_rows) {
    style_cols <- forestplot_style_colours()
    styled_marg <- comparison_forest_apply_effect_styles(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      ref = ref,
      contrast_sep = contrast_sep,
      marker = style_cols$summary,
      bold_summary = TRUE
    )
    resultsA <- styled_marg$resultsA
    resultsB <- styled_marg$resultsB
  }

  # Resolve before type filter: interaction / overall-only panels drop the
  # Younger/Older study rows that carry the non-reference label for the IF header.
  nonref_subgroup <- forestplot_resolve_nonref_subgroup(
    row_labels = row_labels,
    resultsA = resultsA,
    ref = ref
  )
  # Subgroup-only forests show each row's own IF — no non-ref header suffix.
  if (identical(type, "subgroups")) {
    nonref_subgroup <- ""
  }
  favours_user <- favours
  favours <- resolve_forestplot_favours(
    favours = favours_user,
    exponentiate = exponentiate,
    type = type,
    reference = ref,
    nonref_subgroup = nonref_subgroup
  )
  render_snapshot$favours <- favours

  filtered <- forestplot_filter_type(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    type = type
  )
  resultsA <- filtered$resultsA
  resultsB <- filtered$resultsB
  row_labels <- filtered$row_labels
  # IF cells for the non-reference subgroup (header is fixed "IF").
  resultsA <- forestplot_restrict_ifrac_to_subgroup(
    resultsA = resultsA,
    row_labels = row_labels,
    subgroup = nonref_subgroup
  )
  n_study <- sum(resultsA$row_kind == "study", na.rm = TRUE)
  annotations <- forestplot_filter_annotations_type(annotations, type)
  obs_row_index <- which(resultsA$row_kind == "study")
  bold <- if ("bold" %in% names(resultsA)) {
    as.logical(resultsA$bold)
  } else {
    rep(FALSE, nrow(resultsA))
  }
  bold[is.na(bold)] <- FALSE

  if ("colour_label" %in% names(row_labels)) {
    row_labels$colour_label <- NULL
  }

  if (!is.null(grouping)) {
    row_labels$subheading <- grouping
  }

  if (estimates_only) {
    if ("row_kind" %in% names(resultsA)) {
      idx <- resultsA$row_kind != "study"
    } else {
      idx <- resultsB$bold
      idx[is.na(idx)] <- FALSE
    }
    row_labels <- row_labels[idx, , drop = FALSE]
    resultsA <- resultsA[idx, , drop = FALSE]
    resultsB <- resultsB[idx, , drop = FALSE]
    bold <- bold[idx]
    # Study-level IF / N columns are meaningless without study rows.
    drop_left <- intersect(c("ifrac", "trial_size"), names(resultsA))
    if (length(drop_left)) {
      resultsA <- resultsA[, setdiff(names(resultsA), drop_left), drop = FALSE]
    }
    drop_left_b <- intersect(c("ifrac", "trial_size"), names(resultsB))
    if (length(drop_left_b)) {
      resultsB <- resultsB[, setdiff(names(resultsB), drop_left_b), drop = FALSE]
    }
  }

  resultsA$stroke <- stroke
  resultsB$stroke <- stroke

  if ("point_scale" %in% names(resultsA)) {
    resultsA$point_scale[is.na(resultsA$point_scale)] <- 1
    if ("point_scale" %in% names(resultsB)) {
      resultsB$point_scale[is.na(resultsB$point_scale)] <- 1
    } else {
      resultsB$point_scale <- resultsA$point_scale
    }
  }
  if ("plot_yi" %in% names(resultsA)) {
    resultsA$est_plot <- ifelse(is.na(resultsA$yi), resultsA$plot_yi, resultsA$yi)
    resultsB$est_plot <- ifelse(is.na(resultsB$yi), resultsB$plot_yi, resultsB$yi)
  }

  # One panel per forest type.
  panels_to_use <- list(resultsA)
  panel_roles <- type
  row_labels_panels <- row_labels
  if (is.list(col_right_heading)) {
    col_right_heading <- list(col_right_heading[[1L]])
  } else {
    col_right_heading <- list(col_right_heading)
  }
  n_panels <- length(panels_to_use)
  # Role-specific axis title unless the user supplied `xlab`.
  xlab_for_roles <- if (isTRUE(xlab_user_supplied)) xlab else NULL
  xlab <- lapply(panel_roles, function(role) {
    forestplot_panel_xlab(role, exponentiate = exponentiate, xlab = xlab_for_roles)
  })
  has_any_favours <- !is.null(favours)

  panels_to_use <- forestplot_drop_left_on_secondary(panels_to_use)
  metrics <- forestplot_apply_display_metrics(
    panels = panels_to_use,
    base_size = base_size,
    pointsize = pointsize,
    panel.height = panel.height,
    typography_auto = typography_auto,
    has_favours = has_any_favours,
    show_ci_text = show_ci_text
  )
  base_size <- metrics$base_size
  pointsize <- metrics$pointsize
  panel.height <- metrics$panel.height
  render_snapshot$base_size <- base_size
  render_snapshot$pointsize <- pointsize
  render_snapshot$panel.height <- panel.height
  render_snapshot$typography_auto <- typography_auto
  render_snapshot$show_ci_text <- show_ci_text
  render_snapshot$type <- type

  forest_nullval <- if (isTRUE(exponentiate)) 1 else 0

  if (isTRUE(exponentiate) && is.null(xlim)) {
    # forestplot xlog cannot use the default clip = c(-Inf, Inf). log(0)
    # CIs (e.g. imputed missing subgroups with huge vi) make getTicks()
    # fail with "exponential scale that goes <= 0, bottom = '0'".
    limits <- resolve_forestplot_limits(
      panels = panels_to_use,
      xlim = NULL,
      xticks = xticks,
      exponentiate = TRUE,
      nullval = forest_nullval
    )
    xlim <- limits$xlim
    xticks <- limits$xticks
  } else if (!is.null(xlim)) {
    limits <- resolve_forestplot_limits(
      panels = panels_to_use,
      xlim = xlim,
      xticks = xticks,
      exponentiate = exponentiate,
      nullval = forest_nullval
    )
    xlim <- limits$xlim
    xticks <- limits$xticks
  }
  # xlim NULL + xticks set (additive scale): pass ticks through; leave clip
  # open (bayesmeta style). Ratio-scale plots always set clip above.

  estimate_col <- if ("est_plot" %in% names(resultsA)) {
    "est_plot"
  } else if ("plot_yi" %in% names(resultsA)) {
    "plot_yi"
  } else {
    "yi"
  }

  fp <- forestplot_render_interactmeta_panels(
    panels = panels_to_use,
    row_labels = row_labels_panels,
    xlim = xlim,
    xticks = xticks,
    xlab = xlab,
    col_right_heading = col_right_heading,
    estimate_col = estimate_col,
    base_size = base_size,
    pointsize = pointsize,
    digits = digits,
    exponentiate = exponentiate,
    nullval = forest_nullval,
    plot.margin = plot.margin,
    panel.height = panel.height,
    show_left = TRUE,
    left_names = forestplot_col_left_names(resultsA),
    left_headings = forestplot_col_left_headings(
      resultsA,
      nonref_subgroup = nonref_subgroup
    ),
    show_ci_text = show_ci_text,
    favours = favours,
    forest_args = dots,
    annotations = annotations,
    group = group,
    reference = ref,
    panel_roles = panel_roles
  )

  attr(fp, "panels") <- list(resultsA, resultsB)
  attr(fp, "panel.height") <- panel.height
  attr(fp, "base_size") <- base_size
  attr(fp, "forestplot_render") <- render_snapshot
  attr(fp, "typography_auto") <- typography_auto
  attr(fp, "display_inches") <- metrics$device_inches
  attr(fp, "has_favours") <- !is.null(favours)
  attr(fp, "show_ci_text") <- show_ci_text
  attr(fp, "type") <- type
  if (!is.null(annotations) &&
        (length(annotations$lines) ||
           length(annotations$ifrac) ||
           length(annotations$heterogeneity))) {
    attr(fp, "footnote") <- paste(
      c(
        annotations$lines,
        unname(annotations$ifrac),
        unname(annotations$heterogeneity)
      ),
      collapse = "\n"
    )
    attr(fp, "annotations") <- annotations
  }
  fp <- forestplot_add_class(fp)

  fp
}

#' Line width for forest CI / null / diamond geoms (fraction of `base_size`)
#' @noRd
forestplot_base_line_size <- function(base_size = interactmeta_forest_base_size()) {
  base_size / 16
}

#' Draw a `forestplot()` result on the active graphics device.
#' @noRd
forestplot_print_result <- function(fp) {
  suppressWarnings(
    suppressMessages({
      forestplot_maybe_print(fp)
    })
  )
  invisible(fp)
}

#' Drop left-side annotation columns from secondary forest panels.
#' @noRd
forestplot_drop_left_on_secondary <- function(panels) {
  drop <- intersect(c("ifrac", "trial_size"), names(panels[[1L]]))
  if (length(panels) < 2L || !length(drop)) {
    return(panels)
  }
  panels[-1L] <- lapply(panels[-1L], function(p) {
    p[, setdiff(names(p), drop), drop = FALSE]
  })
  panels
}

#' First panel value from a specs field (scalar or length-1 list)
#' @noRd
forestplot_specs_panel_value <- function(x, i = 1L) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.list(x) && !is.data.frame(x) && !is.character(x) && length(x)) {
    return(x[[min(i, length(x))]])
  }
  x
}

check_forestplot_deps <- function() {
  if (!requireNamespace("bayesmeta", quietly = TRUE)) {
    stop(
      "Package `bayesmeta` is required for `forestplot()`. ",
      "Install it with install.packages(\"bayesmeta\").",
      call. = FALSE
    )
  }
  if (!requireNamespace("forestplot", quietly = TRUE)) {
    stop(
      "Package `forestplot` is required for `forestplot()` ",
      "(bayesmeta's plotting backend). ",
      "Install it with install.packages(\"forestplot\").",
      call. = FALSE
    )
  }
}

#' Per-row patient counts for forest study-size labels
#'
#' Resolves `n_var` to a numeric vector aligned with `dat`, or `NULL` to omit
#' labels. Default `"n"` also accepts escalc-style `treat.total` + `cont.total`
#' (summed) when `n` is absent.
#' @noRd
forestplot_n_values <- function(dat, n_var = FALSE) {
  if (identical(n_var, FALSE)) {
    return(NULL)
  }
  if (!is.character(n_var) || length(n_var) != 1L || !nzchar(n_var)) {
    stop("`n_var` must be `FALSE` or a single column name.", call. = FALSE)
  }
  if (n_var %in% names(dat)) {
    return(as.numeric(dat[[n_var]]))
  }
  if (identical(n_var, "n") &&
        all(c("treat.total", "cont.total") %in% names(dat))) {
    return(as.numeric(dat[["treat.total"]]) + as.numeric(dat[["cont.total"]]))
  }
  if (identical(n_var, "n")) {
    warning(
      "Column `n` not found; trial size labels are omitted. ",
      "Pass `n_var` with the patient-count column name, or `n_var = FALSE` to silence.",
      call. = FALSE
    )
    return(NULL)
  }
  stop(
    sprintf("`n_var` column `%s` is not in the forest plot data.", n_var),
    call. = FALSE
  )
}

#' Resolve patient-count column name for forest plot trial labels
#' @noRd
forestplot_resolve_n_var <- function(dat, n_var = FALSE) {
  if (identical(n_var, FALSE)) {
    return(NULL)
  }
  if (!is.character(n_var) || length(n_var) != 1L || !nzchar(n_var)) {
    stop("`n_var` must be `FALSE` or a single column name.", call. = FALSE)
  }
  if (n_var %in% names(dat)) {
    return(n_var)
  }
  # Arm totals are handled in forestplot_n_values(); treat as available.
  if (identical(n_var, "n") &&
        all(c("treat.total", "cont.total") %in% names(dat))) {
    return("n")
  }
  if (identical(n_var, "n")) {
    warning(
      "Column `n` not found; trial size labels are omitted. ",
      "Pass `n_var` with the patient-count column name, or `n_var = FALSE` to silence.",
      call. = FALSE
    )
    return(NULL)
  }
  stop(
    sprintf("`n_var` column `%s` is not in the forest plot data.", n_var),
    call. = FALSE
  )
}

#' Total patients per trial from subgroup rows
#' @noRd
forestplot_trial_n_lookup <- function(dat,
                                      trial_col = "trial",
                                      subgroup_col = "subgroup",
                                      n_var = FALSE,
                                      contrast_sep = "-") {
  n_vals_all <- forestplot_n_values(dat, n_var)
  if (is.null(n_vals_all)) {
    return(NULL)
  }
  sg <- dat
  n_vals <- n_vals_all
  if (subgroup_col %in% names(sg)) {
    keep <- is_ifrac_component_row(sg[[subgroup_col]], contrast_sep = contrast_sep)
    sg <- sg[keep, , drop = FALSE]
    n_vals <- n_vals_all[keep]
  }
  if (!nrow(sg) || !(trial_col %in% names(sg))) {
    return(NULL)
  }
  if (!any(is.finite(n_vals))) {
    return(NULL)
  }
  tr <- as.character(sg[[trial_col]])
  split_n <- split(n_vals, tr)
  stats::setNames(
    vapply(split_n, function(x) sum(x, na.rm = TRUE), numeric(1L)),
    names(split_n)
  )
}

#' Per-row trial size labels (once on the first row of each trial)
#' @noRd
forestplot_row_trial_size <- function(df,
                                      trial_col = "trial",
                                      n_lookup = NULL) {
  n <- nrow(df)
  out <- rep("", n)
  if (is.null(n_lookup) || !length(n_lookup) || !(trial_col %in% names(df))) {
    return(out)
  }
  tr <- as.character(df[[trial_col]])
  n_vals <- unname(n_lookup[tr])
  ok <- is.finite(n_vals) & n_vals > 0
  if (!any(ok)) {
    return(out)
  }
  first_row <- !duplicated(tr)
  ix <- which(ok & first_row)
  out[ix] <- format(n_vals[ix], big.mark = ",", trim = TRUE, scientific = FALSE)
  out
}

#' Add `trial_size` column to a forest panel when patient counts are available
#' @noRd
forestplot_attach_trial_size <- function(panel_df,
                                         source_df,
                                         trial_col = "trial",
                                         subgroup_col = "subgroup",
                                         n_var = FALSE,
                                         contrast_sep = "-") {
  if (!nrow(panel_df) || nrow(panel_df) != nrow(source_df)) {
    stop("`panel_df` and `source_df` must have the same number of rows.", call. = FALSE)
  }
  n_lookup <- forestplot_trial_n_lookup(
    source_df,
    trial_col = trial_col,
    subgroup_col = subgroup_col,
    n_var = n_var,
    contrast_sep = contrast_sep
  )
  if (is.null(n_lookup)) {
    return(panel_df)
  }
  panel_df$trial_size <- forestplot_row_trial_size(
    source_df,
    trial_col = trial_col,
    n_lookup = n_lookup
  )
  panel_df
}

#' @noRd
forestplot_col_left_names <- function(panel_df) {
  cols <- character(0)
  if ("ifrac" %in% names(panel_df)) {
    cols <- c(cols, "ifrac")
  }
  if ("trial_size" %in% names(panel_df)) {
    cols <- c(cols, "trial_size")
  }
  if (!length(cols)) {
    NULL
  } else {
    cols
  }
}

#' @noRd
forestplot_col_left_headings <- function(panel_df, nonref_subgroup = NULL) {
  # Keep headings fixed across type / subgroup so column widths stay stable.
  # `nonref_subgroup` is retained for call compatibility.
  force(nonref_subgroup)
  heads <- character(0)
  if ("ifrac" %in% names(panel_df)) {
    heads <- c(heads, "IF")
  }
  if ("trial_size" %in% names(panel_df)) {
    heads <- c(heads, "N")
  }
  heads
}

#' Blank IF cells that are not for the named subgroup
#'
#' When `subgroup` is set (interaction / overall forests), keep IF only on
#' rows belonging to that subgroup; Interaction / Overall study rows already
#' carry the non-ref IF and are retained. The IF column header itself stays
#' fixed as `"IF"`.
#' @noRd
forestplot_restrict_ifrac_to_subgroup <- function(resultsA,
                                                  row_labels,
                                                  subgroup) {
  sg <- as.character(subgroup %||% "")[[1L]]
  if (!nzchar(sg) || is.null(resultsA) || !"ifrac" %in% names(resultsA)) {
    return(resultsA)
  }
  n <- nrow(resultsA)
  if (!n || is.null(row_labels) || nrow(row_labels) != n) {
    return(resultsA)
  }
  colour <- if ("colour_label" %in% names(row_labels)) {
    as.character(row_labels$colour_label)
  } else {
    rep("", n)
  }
  heading <- if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    rep("", n)
  }
  colour[is.na(colour)] <- ""
  heading[is.na(heading)] <- ""
  ifrac <- as.character(resultsA$ifrac)
  ifrac[is.na(ifrac)] <- ""

  keep <- colour == sg | heading == sg
  # Interaction / Overall rows already store the quoted non-ref IF.
  is_int_ov <- forestplot_heading_is_interaction(heading) |
    forestplot_heading_is_overall(heading)
  keep <- keep | (is_int_ov & nzchar(ifrac))

  resultsA$ifrac[!keep] <- ""
  resultsA
}

#' Non-reference subgroup label from forest row labels (for IF column header)
#' @noRd
forestplot_resolve_nonref_subgroup <- function(row_labels,
                                               resultsA,
                                               ref = NULL) {
  if (is.null(row_labels) || !nrow(row_labels)) {
    return("")
  }
  cl <- if ("colour_label" %in% names(row_labels)) {
    as.character(row_labels$colour_label)
  } else if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    return("")
  }
  cl[is.na(cl)] <- ""
  if (!is.null(resultsA) && "row_kind" %in% names(resultsA)) {
    rk <- as.character(resultsA$row_kind)
    cl <- cl[rk == "study" | is.na(rk) | !nzchar(rk)]
  }
  # Drop Interaction / Overall block tags so only subgroup levels remain.
  cl <- cl[
    nzchar(cl) &
      !forestplot_heading_is_interaction(cl) &
      !forestplot_heading_is_overall(cl)
  ]
  forestplot_nonref_subgroup(cl, ref = ref)
}

validate_forestplot_data <- function(data) {
  if (!is.data.frame(data)) {
    stop("`x` must be a data frame or `escalc` / `ifrac` output.", call. = FALSE)
  }
  if (!nrow(data)) {
    stop("`data` must have at least one row for the forest plot.", call. = FALSE)
  }
  req <- c("trial", "subgroup", "yi", "vi")
  missing <- setdiff(req, names(data))
  if (length(missing)) {
    stop(
      "`data` is missing required columns: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  data
}

#' Validate top-level forest plot arguments
#' @noRd
validate_forestplot_args <- function(prob = 0.95,
                                     which_panel = NULL,
                                     xlim = NULL,
                                     xticks = NULL,
                                     fill = "white",
                                     n_var = FALSE) {
  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) ||
        prob <= 0 || prob > 1) {
    stop("`prob` must be a single finite number in (0, 1].", call. = FALSE)
  }
  if (!is.null(which_panel)) {
    ok <- length(which_panel) == 1L &&
      is.numeric(which_panel) &&
      !is.na(which_panel) &&
      as.integer(which_panel) %in% c(1L, 2L) &&
      identical(as.numeric(which_panel), as.numeric(as.integer(which_panel)))
    if (!ok) {
      stop("`which_panel` must be `NULL`, `1`, or `2`.", call. = FALSE)
    }
  }
  validate_forestplot_axis_arg(xlim, arg = "xlim")
  validate_forestplot_xticks_arg(xticks)
  if (is.list(fill) && !is.data.frame(fill) && length(fill) != 2L) {
    stop("`fill` must be a list of length 2 when provided as a list.", call. = FALSE)
  }
  if (!(identical(n_var, FALSE) ||
          (is.character(n_var) && length(n_var) == 1L && nzchar(n_var)))) {
    stop("`n_var` must be `FALSE` or a single column name.", call. = FALSE)
  }
  invisible(NULL)
}

#' Validate `xlim` shapes for forest plots
#' @noRd
validate_forestplot_axis_arg <- function(x, arg = "xlim") {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  check_vec <- function(v, label) {
    if (!is.numeric(v) || length(v) != 2L || any(!is.finite(v))) {
      stop(
        sprintf("`%s` must be a length-2 finite numeric vector%s.", arg, label),
        call. = FALSE
      )
    }
  }
  if (is.list(x)) {
    if (!length(x)) {
      stop(sprintf("`%s` must be a non-empty list when provided as a list.", arg), call. = FALSE)
    }
    for (i in seq_along(x)) {
      check_vec(x[[i]], sprintf(" (element %d)", i))
    }
    return(invisible(NULL))
  }
  check_vec(x, "")
  invisible(NULL)
}

#' Validate `xticks` shapes for forest plots (any length ≥ 1)
#' @noRd
validate_forestplot_xticks_arg <- function(x) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  check_ticks <- function(v, label) {
    if (!is.numeric(v) || !length(v) || any(!is.finite(v))) {
      stop(
        sprintf("`xticks` must be a finite numeric vector%s.", label),
        call. = FALSE
      )
    }
  }
  if (is.list(x)) {
    if (!length(x)) {
      stop("`xticks` must be a non-empty list when provided as a list.", call. = FALSE)
    }
    for (i in seq_along(x)) {
      check_ticks(x[[i]], sprintf(" (element %d)", i))
    }
    return(invisible(NULL))
  }
  check_ticks(x, "")
  invisible(NULL)
}

#' Accept only known forest-plot input classes
#' @noRd
validate_forestplot_x <- function(x) {
  if (inherits(x, "interactmeta")) {
    return(invisible(NULL))
  }
  if (is.data.frame(x) || inherits(x, "escalc")) {
    return(invisible(NULL))
  }
  stop(
    "`x` must be a data frame, `escalc` / `ifrac()` output, or a CAMS fit from `cam()`.",
    call. = FALSE
  )
}

#' Resolve the study / group label column for [forestplot()]
#' @noRd
resolve_forestplot_label_col <- function(label = NULL, trial = NULL) {
  if (!is.null(label) && !is.null(trial) && !identical(as.character(label), as.character(trial))) {
    stop(
      "Pass only one of `label` or `trial` (soft-deprecated alias for `label`).",
      call. = FALSE
    )
  }
  col <- if (!is.null(label)) {
    label
  } else if (!is.null(trial)) {
    trial
  } else {
    "trial"
  }
  if (!is.character(col) || length(col) != 1L || !nzchar(col)) {
    stop("`label` must be a single non-empty column name.", call. = FALSE)
  }
  col
}

#' Resolve the subgroup column for [forestplot()]
#'
#' Returns `NULL` for univariate forests (`subgroup = FALSE`, or `NULL` when
#' `data` has no `"subgroup"` column).
#' @noRd
resolve_forestplot_subgroup_col <- function(subgroup = NULL, data = NULL) {
  if (isFALSE(subgroup)) {
    return(NULL)
  }
  if (is.null(subgroup)) {
    nms <- tryCatch(names(data), error = function(e) NULL)
    if (!is.null(nms) && "subgroup" %in% nms) {
      return("subgroup")
    }
    return(NULL)
  }
  if (!is.character(subgroup) || length(subgroup) != 1L || !nzchar(subgroup)) {
    stop(
      "`subgroup` must be a single non-empty column name, NULL, or FALSE.",
      call. = FALSE
    )
  }
  subgroup
}

#' Prepare one-row-per-estimate data for a univariate forest
#' @noRd
forestplot_prepare_univariate_data <- function(x, label_col, contrast_sep = "-") {
  if (is_ifrac_output(x)) {
    d <- tibble::as_tibble(x)
    if ("subgroup" %in% names(d)) {
      d <- d[
        !is_ifrac_contrast_row(d$subgroup, contrast_sep = contrast_sep),
        ,
        drop = FALSE
      ]
    }
  } else {
    d <- tryCatch(
      as.data.frame(x, stringsAsFactors = FALSE),
      error = function(e) {
        stop(
          "`x` must be a data frame, `escalc`, or `ifrac` output.",
          call. = FALSE
        )
      }
    )
  }

  if (!is.data.frame(d) || !nrow(d)) {
    stop("`x` must be a non-empty data frame for the forest plot.", call. = FALSE)
  }
  if (!label_col %in% names(d)) {
    stop(
      "Label column `", label_col, "` not found in `x`.",
      call. = FALSE
    )
  }
  missing <- setdiff(c("yi", "vi"), names(d))
  if (length(missing)) {
    stop(
      "`x` is missing required columns: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  trial <- as.character(d[[label_col]])
  sg_lab <- if ("subgroup" %in% names(d)) {
    as.character(d$subgroup)
  } else {
    rep("", nrow(d))
  }
  sg_lab[is.na(sg_lab)] <- ""

  out <- data.frame(
    trial = trial,
    subgroup = sg_lab,
    yi = as.numeric(d$yi),
    vi = as.numeric(d$vi),
    stringsAsFactors = FALSE
  )
  # Carry through size / extra columns (n, custom n_var, etc.)
  skip <- c(label_col, "trial", "subgroup", "yi", "vi")
  for (nm in setdiff(names(d), skip)) {
    if (!nm %in% names(out)) {
      out[[nm]] <- d[[nm]]
    }
  }
  out
}

#' Univariate (single-panel) forest plot — no subgroup interaction panel
#' @noRd
forestplot_univariate <- function(x,
                                  label_col,
                                  n_var = FALSE,
                                  pooled = NULL,
                                  model_names = "Model 1",
                                  xlab = NULL,
                                  favours = NULL,
                                  show_ci_text = FALSE,
                                  typography_auto = FALSE,
                                  exponentiate = FALSE,
                                  col_right_heading = NULL,
                                  stroke = interactmeta_forest_stroke(),
                                  digits = 2,
                                  fill = "white",
                                  study_fill = NULL,
                                  study_colour = NULL,
                                  summary_fill = NULL,
                                  summary_colour = NULL,
                                  pointsize = interactmeta_forest_pointsize(),
                                  base_size = interactmeta_forest_base_size(),
                                  prob = 0.95,
                                  plot.margin = interactmeta_forestplot_margin(),
                                  xlim = NULL,
                                  xticks = NULL,
                                  grouping = NULL,
                                  panel.height = NULL,
                                  estimates_only = FALSE,
                                  ...) {
  dat <- forestplot_prepare_univariate_data(x, label_col = label_col)
  if (!nrow(dat)) {
    stop("No rows found for the univariate forest plot.", call. = FALSE)
  }

  # Information fraction across all rows (study-level weights)
  ifrac_disp <- forestplot_univariate_ifrac_display(dat)

  resultsA <- forestplot_panel_df(
    dat,
    prob = prob,
    ifrac_display = ifrac_disp
  )
  resultsA <- forestplot_attach_univariate_size(
    resultsA,
    dat,
    n_var = n_var
  )

  # Dummy interaction panel (unused; which_panel = 1)
  resultsB <- resultsA
  resultsB$yi <- NA_real_
  resultsB$vi <- NA_real_
  resultsB$ci.lb <- NA_real_
  resultsB$ci.ub <- NA_real_
  if ("est_plot" %in% names(resultsB)) {
    resultsB$est_plot <- NA_real_
  }
  if ("plot_yi" %in% names(resultsB)) {
    resultsB$plot_yi <- NA_real_
  }

  n_study <- nrow(resultsA)
  # Study name as the row label; keep subgroup text only when it disambiguates
  # multiple rows that share a trial id
  trial_chr <- as.character(dat$trial)
  sg_chr <- as.character(dat$subgroup)
  dup_trial <- ave(trial_chr, trial_chr, FUN = length) > 1L
  row_lab <- ifelse(dup_trial & nzchar(sg_chr), sg_chr, "")
  row_labels <- data.frame(
    subheading = trial_chr,
    heading = "",
    label = row_lab,
    stringsAsFactors = FALSE
  )

  sum_f <- if (!is.null(summary_fill)) summary_fill else "white"
  sum_c <- if (!is.null(summary_colour)) summary_colour else interactmeta_accent_colour()

  resultsA$fill <- if (!is.null(study_fill)) study_fill else "white"
  resultsA$colour <- if (!is.null(study_colour)) {
    study_colour
  } else {
    comparison_forest_marker_colour()
  }
  resultsA$bold <- FALSE
  resultsB$fill <- NA_character_
  resultsB$colour <- NA_character_
  resultsB$bold <- FALSE
  resultsA$shape <- 22L
  resultsB$shape <- NA_integer_

  if (!is.null(pooled)) {
    pooled <- validate_pooled(pooled, contrast_sep = "-")
    extra <- append_pooled_rows(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      pooled_subgroup = pooled$subgroup,
      pooled_interaction = NULL,
      model_names = model_names,
      prob = prob,
      study_labels = row_labels$label
    )
    resultsA <- extra$resultsA
    resultsB <- extra$resultsB
    row_labels <- extra$row_labels
    is_summary <- seq_len(nrow(resultsA)) > n_study
    resultsA$fill[is_summary] <- sum_f
    resultsA$colour[is_summary] <- sum_c
    resultsA$bold[is_summary] <- TRUE
    resultsA$shape[is_summary] <- 18L
    resultsB$fill[is_summary] <- NA_character_
    resultsB$colour[is_summary] <- NA_character_
    resultsB$bold[is_summary] <- TRUE
    resultsB$shape[is_summary] <- NA_integer_
  }

  if (is.null(col_right_heading) || !length(col_right_heading)) {
    col_right_heading <- paste0("Effect\n(", 100 * prob, "% interval)")
  }
  if (is.null(xlab) || !length(xlab) || all(!nzchar(as.character(xlab)))) {
    xlab <- forestplot_default_xlab(exponentiate)
  }

  forestplot_render(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    n_study = n_study,
    ref = NA_character_,
    prob = prob,
    contrast_sep = "-",
    xlab = xlab,
    col_right_heading = col_right_heading,
    summary_fill = summary_fill,
    summary_colour = summary_colour,
    fill = fill,
    study_fill = if (!is.null(study_fill)) study_fill else "white",
    study_colour = if (!is.null(study_colour)) {
      study_colour
    } else {
      comparison_forest_marker_colour()
    },
    grouping = grouping,
    estimates_only = estimates_only,
    which_panel = 1L,
    xlim = xlim,
    xticks = xticks,
    base_size = base_size,
    pointsize = pointsize,
    digits = digits,
    stroke = stroke,
    plot.margin = plot.margin,
    panel.height = panel.height,
    exponentiate = exponentiate,
    favours = favours,
    show_ci_text = show_ci_text,
    typography_auto = typography_auto,
    obs_row_index = NULL,
    ...
  )
}

#' Study-level information-fraction display for univariate forests
#' @noRd
forestplot_univariate_ifrac_display <- function(dat) {
  w <- 1 / dat$vi
  w[!is.finite(w) | w < 0] <- NA_real_
  total <- sum(w, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_character_, nrow(dat)))
  }
  pct <- 100 * w / total
  ifelse(
    is.finite(pct),
    format(round(pct, 1), nsmall = 1, trim = TRUE),
    NA_character_
  )
}

#' Attach optional N column for univariate forests
#' @noRd
forestplot_attach_univariate_size <- function(panel_df, dat, n_var = FALSE) {
  n_vals <- forestplot_n_values(dat, n_var)
  if (is.null(n_vals)) {
    return(panel_df)
  }
  panel_df$trial_size <- ifelse(
    is.finite(n_vals),
    as.character(as.integer(round(n_vals))),
    ""
  )
  panel_df
}

split_ifrac_for_forest <- function(data, contrast_sep, group = "subgroup") {
  group <- match.arg(group, c("subgroup", "trial"))
  is_contrast <- is_ifrac_contrast_row(data$subgroup, contrast_sep = contrast_sep)
  subgroup <- data[!is_contrast, , drop = FALSE]
  interaction <- data[is_contrast, , drop = FALSE]

  if (nrow(interaction) == 0L && nrow(subgroup) > 0L) {
    interaction <- build_contrast_rows_from_subgroups(
      subgroup,
      contrast_sep = contrast_sep,
      reference = resolve_forest_reference(subgroup, data)
    )
  }

  subgroup <- order_forest_rows(subgroup, group = group)
  interaction <- order_forest_rows(
    interaction,
    group = group,
    by_trial = unique(subgroup$trial)
  )

  list(subgroup = subgroup, interaction = interaction)
}

resolve_forest_reference <- function(subgroup, data, reference = NULL) {
  # Prefer an explicit / stored reference; never silently invent one for CAMS.
  if (!is.null(reference) && nzchar(as.character(reference)[[1L]])) {
    return(interactmeta_reference(
      if (is.data.frame(data)) data else subgroup,
      reference = reference,
      subgroup_var = "subgroup"
    ))
  }
  if (is.data.frame(data) && "reference" %in% names(data)) {
    return(interactmeta_reference(data, subgroup_var = "subgroup"))
  }
  if (is.data.frame(subgroup) && "reference" %in% names(subgroup)) {
    return(interactmeta_reference(subgroup, subgroup_var = "subgroup"))
  }
  stop(
    "Cannot resolve forest reference: supply `reference` or a `reference` column.",
    call. = FALSE
  )
}

build_contrast_rows_from_subgroups <- function(subgroup, contrast_sep, reference) {
  ref <- reference[[1]]
  trials <- unique(subgroup$trial)
  others <- setdiff(sort(unique(as.character(subgroup$subgroup))), ref)
  if (!length(others)) {
    return(subgroup[0, , drop = FALSE])
  }

  wide <- subgroup |>
    dplyr::select("trial", "subgroup", "yi", "vi") |>
    tidyr::pivot_wider(
      names_from = "subgroup",
      values_from = c("yi", "vi"),
      names_glue = "{.value}_{subgroup}"
    )

  out <- lapply(trials, function(tr) {
    row <- wide[wide$trial == tr, , drop = FALSE]
    if (!nrow(row)) {
      return(NULL)
    }
    dplyr::bind_rows(lapply(others, function(lab) {
      yi1 <- row[[paste0("yi_", lab)]]
      yi2 <- row[[paste0("yi_", ref)]]
      vi1 <- row[[paste0("vi_", lab)]]
      vi2 <- row[[paste0("vi_", ref)]]
      dplyr::tibble(
        trial = tr,
        subgroup = paste0(lab, contrast_sep, ref),
        yi = yi1 - yi2,
        vi = vi1 + vi2
      )
    }))
  })
  dplyr::bind_rows(out)
}

align_interaction_to_subgroups <- function(subgroup, interaction, reference = NULL) {
  if (!nrow(subgroup)) {
    return(list(subgroup = subgroup, interaction = interaction))
  }
  if (!nrow(interaction)) {
    blank <- subgroup
    blank$yi <- NA_real_
    blank$vi <- NA_real_
    return(list(subgroup = subgroup, interaction = blank))
  }

  int_one <- interaction |>
    dplyr::group_by(.data$trial) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select("trial", "yi", "vi")

  joined <- dplyr::left_join(
    subgroup |> dplyr::select("trial", "subgroup"),
    int_one,
    by = "trial"
  )

  if (!is.null(reference) && length(reference) == 1L && nzchar(reference)) {
    is_ref <- as.character(joined$subgroup) == reference
    joined$yi[is_ref] <- NA_real_
    joined$vi[is_ref] <- NA_real_
  }

  list(subgroup = subgroup, interaction = joined)
}

order_forest_rows <- function(df,
                              by_trial = NULL,
                              group = "subgroup") {
  if (!nrow(df)) {
    return(df)
  }
  group <- match.arg(group, c("subgroup", "trial"))
  # Subgroup nesting uses alphabetical trial labels within each block (classic
  # metafor forest ordering). Trial nesting keeps first-appearance block order
  # unless `by_trial` is supplied.
  trial_ord <- if (is.null(by_trial)) {
    tr <- unique(as.character(df$trial))
    if (identical(group, "subgroup")) {
      sort(tr)
    } else {
      tr
    }
  } else {
    by_trial
  }
  df$trial <- factor(as.character(df$trial), levels = trial_ord)
  sg_was_factor <- is.factor(df$subgroup)
  sg_levels <- if (sg_was_factor) levels(df$subgroup) else unique(as.character(df$subgroup))
  df$subgroup <- factor(as.character(df$subgroup), levels = sg_levels)

  if (identical(group, "subgroup")) {
    df <- df[order(df$subgroup, df$trial), , drop = FALSE]
  } else {
    df <- df[order(df$trial, df$subgroup), , drop = FALSE]
  }
  df$trial <- as.character(df$trial)
  if (sg_was_factor) {
    df$subgroup <- factor(as.character(df$subgroup), levels = sg_levels)
  } else {
    df$subgroup <- as.character(df$subgroup)
  }
  df
}

#' Subgroup identity used for forest colouring (independent of display nesting)
#' @noRd
forestplot_row_colour_label <- function(row_lab) {
  if ("colour_label" %in% names(row_lab)) {
    cl <- as.character(row_lab$colour_label[[1L]])
    if (!is.na(cl) && nzchar(cl)) {
      return(cl)
    }
  }
  as.character(row_lab$label[[1L]])
}

#' Non-reference subgroup component used for Interaction / Overall rows
#'
#' @return A list with `ifrac` (plain percentage string) and `subgroup`
#'   (non-reference label for block headers).
#' @noRd
forestplot_nonref_component <- function(ifrac_vals,
                                        subgroup_vals,
                                        ref = NULL) {
  ifrac_vals <- as.character(ifrac_vals)
  subgroup_vals <- as.character(subgroup_vals)
  ifrac_vals[is.na(ifrac_vals)] <- ""
  subgroup_vals[is.na(subgroup_vals)] <- ""
  empty <- list(ifrac = "", subgroup = "")
  n <- length(ifrac_vals)
  if (!n) {
    return(empty)
  }
  ref_chr <- if (is.null(ref) || !length(ref) ||
                  !nzchar(as.character(ref)[[1L]])) {
    NA_character_
  } else {
    as.character(ref)[[1L]]
  }
  if (!is.na(ref_chr)) {
    idx <- which(
      nzchar(subgroup_vals) &
        subgroup_vals != ref_chr &
        nzchar(ifrac_vals)
    )
    if (!length(idx)) {
      idx <- which(
        nzchar(subgroup_vals) & subgroup_vals != ref_chr
      )
    }
    if (!length(idx)) {
      idx <- which(nzchar(ifrac_vals))
    }
  } else {
    idx <- which(nzchar(ifrac_vals) | nzchar(subgroup_vals))
  }
  if (!length(idx)) {
    return(empty)
  }
  i <- idx[[1L]]
  list(ifrac = ifrac_vals[[i]], subgroup = subgroup_vals[[i]])
}

#' Plain non-reference IF percentage for trial-level rows
#' @noRd
forestplot_nonref_ifrac_display <- function(ifrac_vals,
                                            subgroup_vals,
                                            ref = NULL) {
  forestplot_nonref_component(ifrac_vals, subgroup_vals, ref = ref)$ifrac
}

#' Non-reference subgroup label (for Interaction / Overall headers)
#' @noRd
forestplot_nonref_subgroup <- function(subgroup_vals, ref = NULL) {
  forestplot_nonref_component(
    ifrac_vals = rep("", length(subgroup_vals)),
    subgroup_vals = subgroup_vals,
    ref = ref
  )$subgroup
}

#' Block heading with optional non-reference subgroup, e.g. `"Overall (Older)"`
#' @noRd
forestplot_heading_with_nonref <- function(base, subgroup) {
  base <- as.character(base)[[1L]]
  subgroup <- as.character(subgroup %||% "")[[1L]]
  if (!nzchar(subgroup)) {
    base
  } else {
    paste0(base, " (", subgroup, ")")
  }
}

#' Whether a forest heading is an Interaction block (possibly with subgroup)
#' @noRd
forestplot_heading_is_interaction <- function(heading) {
  grepl("^Interaction(\\s*\\(|$)", as.character(heading))
}

#' Whether a forest heading is an Overall block (possibly with subgroup)
#' @noRd
forestplot_heading_is_overall <- function(heading) {
  grepl("^Overall(\\s*\\(|$)", as.character(heading))
}

#' Stack subgroup, interaction, and overall into one vertical forest panel
#'
#' Order: subgroup blocks (study / shrinkage / pooled / prediction), then an
#' Interaction block (one row per trial + optional shrinkage + summaries), then
#' Overall (trial-level overall study / optional shrinkage + overall diamonds).
#'
#' @noRd
forestplot_stack_panels <- function(resultsA,
                                    resultsB,
                                    row_labels,
                                    ref = NULL,
                                    overall_shrinkage = NULL) {
  n <- nrow(resultsA)
  if (!identical(nrow(resultsB), n) || !identical(nrow(row_labels), n)) {
    stop("Panel rows and labels must have the same length.", call. = FALSE)
  }
  if (is.null(overall_shrinkage)) {
    overall_shrinkage <- attr(resultsA, "overall_shrinkage", exact = TRUE)
  }

  rk <- if ("row_kind" %in% names(resultsA)) {
    as.character(resultsA$row_kind)
  } else {
    rep("study", n)
  }
  rk[is.na(rk) | !nzchar(rk)] <- "study"

  heading <- if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    as.character(row_labels$subheading)
  }
  lab <- as.character(row_labels$label)
  colour_lab <- if ("colour_label" %in% names(row_labels)) {
    as.character(row_labels$colour_label)
  } else {
    lab
  }
  heading[is.na(heading)] <- ""
  lab[is.na(lab)] <- ""
  colour_lab[is.na(colour_lab)] <- ""

  study_ix <- which(rk == "study")
  nest_by_subgroup <- TRUE
  if (length(study_ix)) {
    nest_by_subgroup <- mean(
      colour_lab[study_ix] == heading[study_ix],
      na.rm = TRUE
    ) >= 0.5
  }
  trial_id <- if (isTRUE(nest_by_subgroup)) lab else heading

  blank_from <- function(template, source = NULL) {
    out <- template[1L, , drop = FALSE]
    for (col in names(out)) {
      if (is.character(out[[col]])) {
        out[[col]] <- ""
      } else if (is.logical(out[[col]])) {
        out[[col]] <- FALSE
      } else {
        out[[col]] <- NA
      }
    }
    if (!is.null(source)) {
      shared <- intersect(names(out), names(source))
      out[1L, shared] <- source[1L, shared, drop = FALSE]
    }
    if ("point_scale" %in% names(out) &&
          (is.na(out$point_scale[[1L]]) || !length(out$point_scale))) {
      out$point_scale <- 1
    }
    out
  }

  parts_panel <- list()
  parts_lab <- list()
  push_row <- function(panel_row, lab_row) {
    parts_panel[[length(parts_panel) + 1L]] <<- panel_row
    parts_lab[[length(parts_lab) + 1L]] <<- lab_row
  }

  # Interaction-only shrinkage lives on panel B (A shape is NA).
  is_int_shr <- rk == "shrinkage" &
    is.na(resultsA$shape) &
    !is.na(resultsB$shape)

  sg_ix <- which(
    rk %in% c("study", "shrinkage", "pooled", "prediction") & !is_int_shr
  )
  for (i in sg_ix) {
    # Keep NA pooled/prediction placeholders so switched-off layers reserve
    # the same vertical space as the full plot.
    push_row(resultsA[i, , drop = FALSE], row_labels[i, , drop = FALSE])
  }

  int_header <- "Interaction"

  # One interaction study row per trial, immediately followed by its shrinkage
  # row when present. Under subgroup nesting, trials are alphabetical.
  int_study <- which(rk == "study" & is.finite(resultsB$yi))
  int_shr_ix <- which(is_int_shr & is.finite(resultsB$yi))

  first_index_by_trial <- function(idx) {
    out <- integer(0)
    names <- character(0)
    for (i in idx) {
      tid <- trial_id[[i]]
      if (!nzchar(tid) || tid %in% names) {
        next
      }
      names <- c(names, tid)
      out <- c(out, i)
    }
    stats::setNames(out, names)
  }

  study_by_trial <- first_index_by_trial(int_study)
  shr_by_trial <- first_index_by_trial(int_shr_ix)
  trial_order <- unique(c(names(study_by_trial), names(shr_by_trial)))
  if (isTRUE(nest_by_subgroup)) {
    trial_order <- sort(trial_order)
  }

  for (tid in trial_order) {
    if (tid %in% names(study_by_trial)) {
      i <- study_by_trial[[tid]]
      row <- blank_from(resultsA, resultsB[i, , drop = FALSE])
      # Non-reference subgroup IF (with group label) + study size.
      if ("ifrac" %in% names(row)) {
        tr_ix <- which(
          rk == "study" & !is_int_shr & trial_id == tid & is.finite(resultsA$yi)
        )
        row$ifrac <- forestplot_nonref_ifrac_display(
          ifrac_vals = if ("ifrac" %in% names(resultsA)) {
            resultsA$ifrac[tr_ix]
          } else {
            character(0)
          },
          subgroup_vals = colour_lab[tr_ix],
          ref = ref
        )
      }
      if ("trial_size" %in% names(resultsA)) {
        sizes <- resultsA$trial_size[trial_id == tid]
        sizes <- sizes[nzchar(as.character(sizes))]
        row$trial_size <- if (length(sizes)) sizes[[1L]] else ""
      }
      row$row_kind <- "study"
      # Match Overall block aesthetics: filled study squares (not open circles).
      row$shape <- 15L
      row$bold <- FALSE
      rl <- row_labels[i, , drop = FALSE]
      rl$heading <- int_header
      if ("subheading" %in% names(rl)) {
        rl$subheading <- NA_character_
      }
      rl$label <- tid
      if ("colour_label" %in% names(rl)) {
        rl$colour_label <- int_header
      }
      push_row(row, rl)
    }

    if (tid %in% names(shr_by_trial)) {
      i <- shr_by_trial[[tid]]
      row <- blank_from(resultsA, resultsB[i, , drop = FALSE])
      if ("ifrac" %in% names(row)) {
        row$ifrac <- ""
      }
      if ("trial_size" %in% names(row)) {
        row$trial_size <- ""
      }
      row$row_kind <- "shrinkage"
      # Match Overall shrinkage marker (filled square; grey fill applied earlier).
      row$shape <- 15L
      if ("point_scale" %in% names(resultsB)) {
        row$point_scale <- resultsB$point_scale[[i]]
      }
      rl <- row_labels[i, , drop = FALSE]
      rl$heading <- int_header
      if ("subheading" %in% names(rl)) {
        rl$subheading <- NA_character_
      }
      rl$label <- tid
      if ("colour_label" %in% names(rl)) {
        rl$colour_label <- int_header
      }
      push_row(row, rl)
    }
  }

  for (kind in c("pooled", "prediction")) {
    ii_finite <- which(rk == kind & is.finite(resultsB$yi))
    ii_any <- which(rk == kind)
    if (!length(ii_any)) {
      next
    }
    i <- if (length(ii_finite)) ii_finite[[1L]] else ii_any[[1L]]
    row <- blank_from(resultsA, resultsB[i, , drop = FALSE])
    if ("ifrac" %in% names(row)) {
      row$ifrac <- ""
    }
    if ("trial_size" %in% names(row)) {
      row$trial_size <- ""
    }
    row$row_kind <- kind
    if (is.finite(resultsB$yi[[i]])) {
      row$shape <- comparison_forest_diamond_shape()
      row$bold <- TRUE
    } else {
      # Placeholder for a switched-off interaction layer.
      row$yi <- NA_real_
      row$vi <- NA_real_
      if ("ci.lb" %in% names(row)) {
        row$ci.lb <- NA_real_
      }
      if ("ci.ub" %in% names(row)) {
        row$ci.ub <- NA_real_
      }
      row$shape <- NA_integer_
      row$bold <- FALSE
    }
    rl <- row_labels[i, , drop = FALSE]
    rl$heading <- int_header
    if ("subheading" %in% names(rl)) {
      rl$subheading <- NA_character_
    }
    # Same labels as the Overall block (Pooled diamond, Prediction bar).
    rl$label <- if (identical(kind, "pooled")) "Pooled" else "Prediction"
    if ("colour_label" %in% names(rl)) {
      rl$colour_label <- int_header
    }
    push_row(row, rl)
  }

  ov_header <- "Overall"
  # Include NA Overall placeholders so omitted pooled/prediction keep space.
  ov_ix <- which(rk == "overall")
  ov_finite <- which(rk == "overall" & is.finite(resultsA$yi))

  # Trial-level overall study (+ shrinkage) under Overall, parallel to the
  # Interaction block: precision-weighted within-trial combination of subgroup
  # endpoints (and of subgroup shrinkage when present).
  if (length(ov_finite)) {
    sg_study_ix <- which(
      rk == "study" & is.finite(resultsA$yi) & is.finite(resultsA$vi) & !is_int_shr
    )
    sg_shr_ix <- which(
      rk == "shrinkage" &
        is.finite(resultsA$yi) &
        is.finite(resultsA$vi) &
        !is_int_shr &
        !is.na(resultsA$shape)
    )
    ov_trials <- unique(trial_id[sg_study_ix])
    ov_trials <- ov_trials[nzchar(ov_trials)]
    if (isTRUE(nest_by_subgroup)) {
      ov_trials <- sort(ov_trials)
    }

    for (tid in ov_trials) {
      st_ii <- sg_study_ix[trial_id[sg_study_ix] == tid]
      if (!length(st_ii)) {
        next
      }
      yi_s <- resultsA$yi[st_ii]
      vi_s <- pmax(resultsA$vi[st_ii], 1e-12)
      w <- 1 / vi_s
      ov_yi <- sum(yi_s * w) / sum(w)
      ov_vi <- 1 / sum(w)

      row <- blank_from(resultsA, resultsA[st_ii[[1L]], , drop = FALSE])
      row$yi <- ov_yi
      row$vi <- ov_vi
      if ("ci.lb" %in% names(row) && "ci.ub" %in% names(row)) {
        half <- stats::qnorm(0.975) * sqrt(ov_vi)
        row$ci.lb <- ov_yi - half
        row$ci.ub <- ov_yi + half
      }
      # Non-reference subgroup IF (with group label) + study size.
      if ("ifrac" %in% names(row)) {
        row$ifrac <- forestplot_nonref_ifrac_display(
          ifrac_vals = if ("ifrac" %in% names(resultsA)) {
            resultsA$ifrac[st_ii]
          } else {
            character(0)
          },
          subgroup_vals = colour_lab[st_ii],
          ref = ref
        )
      }
      if ("trial_size" %in% names(resultsA)) {
        sizes <- resultsA$trial_size[st_ii]
        sizes <- sizes[nzchar(as.character(sizes))]
        row$trial_size <- if (length(sizes)) sizes[[1L]] else ""
      }
      row$row_kind <- "study"
      row$shape <- if (!is.na(resultsA$shape[[st_ii[[1L]]]])) {
        as.integer(resultsA$shape[[st_ii[[1L]]]])
      } else {
        15L
      }
      row$bold <- FALSE
      rl <- row_labels[st_ii[[1L]], , drop = FALSE]
      rl$heading <- ov_header
      if ("subheading" %in% names(rl)) {
        rl$subheading <- NA_character_
      }
      rl$label <- tid
      if ("colour_label" %in% names(rl)) {
        rl$colour_label <- ov_header
      }
      push_row(row, rl)

      shr_ii <- sg_shr_ix[trial_id[sg_shr_ix] == tid]
      ov_shr_hit <- NULL
      if (!is.null(overall_shrinkage) && nrow(overall_shrinkage)) {
        hit <- overall_shrinkage[
          as.character(overall_shrinkage$trial) == tid,
          ,
          drop = FALSE
        ]
        if (nrow(hit) && is.finite(hit$yi[[1L]]) && is.finite(hit$vi[[1L]])) {
          ov_shr_hit <- hit[1L, , drop = FALSE]
        }
      }
      if (is.null(ov_shr_hit) && !length(shr_ii)) {
        next
      }
      if (!is.null(ov_shr_hit)) {
        shr_yi <- as.numeric(ov_shr_hit$yi[[1L]])
        shr_vi <- pmax(as.numeric(ov_shr_hit$vi[[1L]]), 1e-12)
        shr_template_ix <- if (length(shr_ii)) shr_ii[[1L]] else st_ii[[1L]]
      } else {
        yi_h <- resultsA$yi[shr_ii]
        vi_h <- pmax(resultsA$vi[shr_ii], 1e-12)
        wh <- 1 / vi_h
        shr_yi <- sum(yi_h * wh) / sum(wh)
        shr_vi <- 1 / sum(wh)
        shr_template_ix <- shr_ii[[1L]]
      }

      row <- blank_from(resultsA, resultsA[shr_template_ix, , drop = FALSE])
      row$yi <- shr_yi
      row$vi <- shr_vi
      if ("ci.lb" %in% names(row) && "ci.ub" %in% names(row)) {
        half <- stats::qnorm(0.975) * sqrt(shr_vi)
        row$ci.lb <- shr_yi - half
        row$ci.ub <- shr_yi + half
      }
      if ("ifrac" %in% names(row)) {
        row$ifrac <- ""
      }
      if ("trial_size" %in% names(row)) {
        row$trial_size <- ""
      }
      row$row_kind <- "shrinkage"
      row$shape <- if (!is.na(resultsA$shape[[shr_template_ix]])) {
        as.integer(resultsA$shape[[shr_template_ix]])
      } else {
        15L
      }
      if ("point_scale" %in% names(resultsA)) {
        row$point_scale <- resultsA$point_scale[[shr_template_ix]]
      }
      row$bold <- FALSE
      rl <- row_labels[shr_template_ix, , drop = FALSE]
      rl$heading <- ov_header
      if ("subheading" %in% names(rl)) {
        rl$subheading <- NA_character_
      }
      rl$label <- tid
      if ("colour_label" %in% names(rl)) {
        rl$colour_label <- ov_header
      }
      push_row(row, rl)
    }
  }

  for (i in ov_ix) {
    rl <- row_labels[i, , drop = FALSE]
    rl$heading <- ov_header
    if ("subheading" %in% names(rl)) {
      rl$subheading <- NA_character_
    }
    push_row(resultsA[i, , drop = FALSE], rl)
  }

  if (!length(parts_panel)) {
    stop("Stacked forest panel is empty.", call. = FALSE)
  }

  # Align columns before rbind (some rows add `row_kind` / extras later).
  all_panel_cols <- unique(unlist(lapply(parts_panel, names), use.names = FALSE))
  parts_panel <- lapply(parts_panel, function(p) {
    missing <- setdiff(all_panel_cols, names(p))
    for (nm in missing) {
      p[[nm]] <- NA
    }
    p[all_panel_cols]
  })
  all_lab_cols <- unique(unlist(lapply(parts_lab, names), use.names = FALSE))
  parts_lab <- lapply(parts_lab, function(p) {
    missing <- setdiff(all_lab_cols, names(p))
    for (nm in missing) {
      p[[nm]] <- NA
    }
    p[all_lab_cols]
  })

  panel <- do.call(rbind, parts_panel)
  labels <- do.call(rbind, parts_lab)
  rownames(panel) <- NULL
  rownames(labels) <- NULL

  labels$variable <- paste0(
    apply(labels, 1L, paste0, collapse = ""),
    "\r",
    seq_len(nrow(labels))
  )
  panel$variable <- labels$variable

  list(
    resultsA = panel,
    resultsB = panel,
    row_labels = labels,
    n_study = sum(panel$row_kind == "study", na.rm = TRUE)
  )
}

#' Default x-axis title for a forest panel role
#' @noRd
forestplot_panel_xlab <- function(role, exponentiate = FALSE, xlab = NULL) {
  if (!is.null(xlab) && !identical(xlab, FALSE)) {
    lab <- as.character(xlab)[[1L]]
    if (nzchar(lab)) {
      return(lab)
    }
  }
  if (identical(role, "interaction")) {
    if (isTRUE(exponentiate)) {
      "Interaction"
    } else {
      "Interaction (log scale)"
    }
  } else {
    if (isTRUE(exponentiate)) "Effect" else "Effect (log scale)"
  }
}

#' Plotmath expression for a forest footnote so Greek renders
#'
#' Footnotes are stored as UTF-8 (`π`, `τ`, `τ_γ`) for tests and print.
#' grid / base graphics draw them with plotmath (`pi`, `tau`, `tau[gamma]`)
#' so the Greek letters do not depend on the device font.
#' @noRd
forestplot_footnote_as_plotmath <- function(txt) {
  txt <- paste(enc2utf8(as.character(txt)), collapse = " \u00b7 ")
  txt <- txt[[1L]]
  if (!nzchar(txt) || is.na(txt)) {
    return("")
  }
  tok_tg <- "<<TG>>"
  tok_tau <- "<<TAU>>"
  tok_pi <- "<<PI>>"
  tok_dot <- "<<DOT>>"
  txt <- gsub("\u03c4_\u03b3", tok_tg, txt, fixed = TRUE)
  txt <- gsub("\u03c4_y", tok_tg, txt, fixed = TRUE)
  txt <- gsub("tau_gamma", tok_tg, txt, fixed = TRUE)
  txt <- gsub("\u03c4", tok_tau, txt, fixed = TRUE)
  txt <- gsub("\u03c0", tok_pi, txt, fixed = TRUE)
  txt <- gsub(" \u00b7 ", tok_dot, txt, fixed = TRUE)
  txt <- gsub("\u00b7", tok_dot, txt, fixed = TRUE)
  txt <- gsub("\\", "\\\\", txt, fixed = TRUE)
  txt <- gsub("\"", "\\\"", txt, fixed = TRUE)
  txt <- gsub(tok_tg, '" * tau[gamma] * "', txt, fixed = TRUE)
  txt <- gsub(tok_tau, '" * tau * "', txt, fixed = TRUE)
  txt <- gsub(tok_pi, '" * pi * "', txt, fixed = TRUE)
  txt <- gsub(tok_dot, '" %.% "', txt, fixed = TRUE)
  expr <- tryCatch(
    parse(text = paste0('"', txt, '"'))[[1L]],
    error = function(e) NULL
  )
  if (is.null(expr)) {
    paste(as.character(txt), collapse = "")
  } else {
    expr
  }
}

#' Keep only the requested forest type block(s) after stacking
#' @noRd
forestplot_filter_type <- function(resultsA,
                                       resultsB,
                                       row_labels,
                                       type = c("subgroups", "interaction", "overall", "combined")) {
  type <- comparison_resolve_effect_type(match.arg(type))
  n <- nrow(resultsA)
  if (!n) {
    stop("Forest panel is empty after stacking.", call. = FALSE)
  }
  heading <- if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    rep("", n)
  }
  heading[is.na(heading)] <- ""
  # Forward-fill blank headings within a block (layout may leave some empty).
  last <- ""
  for (i in seq_len(n)) {
    if (nzchar(heading[[i]])) {
      last <- heading[[i]]
    } else if (nzchar(last)) {
      heading[[i]] <- last
    }
  }

  idx <- if (identical(type, "subgroups")) {
    # Subgroup blocks only; Overall is type = "overall".
    !forestplot_heading_is_interaction(heading) &
      !forestplot_heading_is_overall(heading)
  } else if (identical(type, "interaction")) {
    forestplot_heading_is_interaction(heading)
  } else {
    forestplot_heading_is_overall(heading)
  }
  if (!any(idx)) {
    stop(
      sprintf(
        "`type = \"%s\"` produced an empty forest. ",
        type
      ),
      "Check `include_pooled` / `include_prediction` / `include_overall` ",
      "and that the fit has rows for this type.",
      call. = FALSE
    )
  }
  list(
    resultsA = resultsA[idx, , drop = FALSE],
    resultsB = resultsB[idx, , drop = FALSE],
    row_labels = row_labels[idx, , drop = FALSE]
  )
}

#' Drop annotation notes that do not belong to the active type
#' @noRd
forestplot_filter_annotations_type <- function(annotations,
                                                  type = c("subgroups", "interaction", "overall", "combined")) {
  type <- comparison_resolve_effect_type(match.arg(type))
  annotations <- forestplot_normalize_annotations(annotations)
  # `subgroups` keeps the full footnote stack (overall IF + both heterogeneity
  # notes). Interaction / overall types keep only their slice.
  if (identical(type, "subgroups")) {
    return(annotations)
  }
  het <- annotations$heterogeneity
  ifrac <- annotations$ifrac
  if (identical(type, "interaction")) {
    ifrac <- character()
    het <- het[intersect(names(het), "interaction")]
  } else {
    # overall: keep overall het; IF notes are subgroup-specific — keep as
    # trailing unmatched notes (layout places them after remaining blocks).
    het <- het[intersect(names(het), "overall")]
  }
  annotations$heterogeneity <- het
  annotations$ifrac <- ifrac
  annotations
}

#' Row labels for study forest panels
#'
#' @param group `"subgroup"`: heading = subgroup, label = trial (metafor-style).
#'   `"trial"`: heading = trial, label = subgroup.
#' @noRd
forestplot_study_row_labels <- function(trial, subgroup, group = "subgroup") {
  group <- match.arg(group, c("subgroup", "trial"))
  trial <- as.character(trial)
  subgroup <- as.character(subgroup)
  if (identical(group, "subgroup")) {
    data.frame(
      subheading = subgroup,
      heading = "",
      label = trial,
      colour_label = subgroup,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      subheading = trial,
      heading = "",
      label = subgroup,
      colour_label = subgroup,
      stringsAsFactors = FALSE
    )
  }
}

#' Plotmath expression for a forest footnote so Greek renders
#'
#' Footnotes are stored as UTF-8 (`π`, `τ`, `τ_γ`) for tests and print.
#' grid / base graphics draw them with plotmath (`pi`, `tau`, `tau[gamma]`)
#' so the Greek letters do not depend on the device font.
#' @noRd
forestplot_footnote_as_plotmath <- function(txt) {
  txt <- paste(enc2utf8(as.character(txt)), collapse = " \u00b7 ")
  txt <- txt[[1L]]
  if (!nzchar(txt) || is.na(txt)) {
    return("")
  }
  tok_tg <- "<<TG>>"
  tok_tau <- "<<TAU>>"
  tok_pi <- "<<PI>>"
  tok_dot <- "<<DOT>>"
  txt <- gsub("\u03c4_\u03b3", tok_tg, txt, fixed = TRUE)
  txt <- gsub("\u03c4_y", tok_tg, txt, fixed = TRUE)
  txt <- gsub("tau_gamma", tok_tg, txt, fixed = TRUE)
  txt <- gsub("\u03c4", tok_tau, txt, fixed = TRUE)
  txt <- gsub("\u03c0", tok_pi, txt, fixed = TRUE)
  txt <- gsub(" \u00b7 ", tok_dot, txt, fixed = TRUE)
  txt <- gsub("\u00b7", tok_dot, txt, fixed = TRUE)
  txt <- gsub("\\", "\\\\", txt, fixed = TRUE)
  txt <- gsub("\"", "\\\"", txt, fixed = TRUE)
  txt <- gsub(tok_tg, '" * tau[gamma] * "', txt, fixed = TRUE)
  txt <- gsub(tok_tau, '" * tau * "', txt, fixed = TRUE)
  txt <- gsub(tok_pi, '" * pi * "', txt, fixed = TRUE)
  txt <- gsub(tok_dot, '" %.% "', txt, fixed = TRUE)
  expr <- tryCatch(
    parse(text = paste0('"', txt, '"'))[[1L]],
    error = function(e) NULL
  )
  if (is.null(expr)) {
    paste(as.character(txt), collapse = "")
  } else {
    expr
  }
}

#' Default x-axis title for forest plots
#'
#' CAMS stores effects on the estimation scale (often log OR / log HR). When
#' `exponentiate = FALSE` the axis is that scale; when `TRUE`, tick labels are
#' shown on the ratio scale via metafor's `atransf = exp`.
#' @noRd
forestplot_default_xlab <- function(exponentiate = FALSE) {
  if (isTRUE(exponentiate)) {
    "Effect"
  } else {
    "Effect (log scale)"
  }
}

#' Evenly spaced "round" tick positions for a forest axis
#'
#' Rounds `xlim` outward onto a regular step grid (default 0.5 on the additive
#' scale), then expands symmetrically about `nullval` so the reference line
#' sits in the centre of the forest band. Ratio-scale axes do the same on the
#' log scale.
#' @noRd
forestplot_nice_axis <- function(xlim,
                                 exponentiate = FALSE,
                                 nullval = NULL,
                                 step = NULL) {
  xlim <- as.numeric(xlim)
  if (length(xlim) != 2L || any(!is.finite(xlim))) {
    stop("`xlim` must be a length-2 finite numeric vector.", call. = FALSE)
  }
  if (xlim[[1L]] > xlim[[2L]]) {
    xlim <- rev(xlim)
  }
  if (is.null(nullval)) {
    nullval <- if (isTRUE(exponentiate)) 1 else 0
  }
  nullval <- as.numeric(nullval)[1L]

  if (isTRUE(exponentiate)) {
    return(forestplot_nice_axis_ratio(xlim, nullval = nullval))
  }

  span <- diff(xlim)
  if (is.null(step) || !is.finite(step) || step <= 0) {
    # Prefer 0.5 for typical forest spans so ticks stay visually even and
    # readable (avoid pretty()'s irregular mixes like -1.8, -1, 0, 0.5).
    if (span <= 4) {
      step <- 0.5
    } else if (span <= 8) {
      step <- 1
    } else {
      candidates <- c(1, 2, 5, 10)
      n_ticks <- span / candidates
      step <- candidates[[which.min(abs(n_ticks - 7))]]
    }
  }

  lo <- floor(xlim[[1L]] / step) * step
  hi <- ceiling(xlim[[2L]] / step) * step
  if (is.finite(nullval)) {
    # Centre the null line: equal extent on both sides, snapped to the step.
    half <- max(nullval - lo, hi - nullval)
    half <- ceiling(half / step) * step
    lo <- nullval - half
    hi <- nullval + half
  }
  # Guard against floating-point drift (e.g. 0.1 + 0.2).
  ticks <- round(seq(lo, hi, by = step) / step) * step
  ticks <- unique(as.numeric(ticks))
  list(xlim = c(lo, hi), xticks = ticks)
}

#' @noRd
forestplot_nice_axis_ratio <- function(xlim, nullval = 1) {
  xlim <- pmax(as.numeric(xlim), 1e-8)
  nullval <- max(as.numeric(nullval)[1L], 1e-8)
  # Powers of 2 about the null (usually 1): …, 0.25, 0.5, 1, 2, 4, …
  # so ticks are log-even and the reference line has an explicit label.
  log2_lo <- log2(xlim[[1L]] / nullval)
  log2_hi <- log2(xlim[[2L]] / nullval)
  k_lo <- floor(log2_lo)
  k_hi <- ceiling(log2_hi)
  # Centre the null on a symmetric power-of-2 grid.
  k_half <- max(abs(k_lo), abs(k_hi), 1L)
  ks <- seq.int(-k_half, k_half)
  ticks <- as.numeric(nullval * (2^ks))
  list(xlim = range(ticks), xticks = ticks)
}

#' Apply nice axis rounding to one xlim / xticks pair
#' @noRd
forestplot_apply_nice_axis <- function(xlim,
                                       xticks = NULL,
                                       exponentiate = FALSE,
                                       nullval = NULL) {
  # Honour an explicit tick set with the supplied limits. Recentering onto
  # the null line (e.g. MAP IF forests with xlim = c(0.07, 0.8)) would pull
  # 0 inside the clip and draw a second vertical rule next to the graph
  # separator.
  if (is.numeric(xticks) && length(xticks)) {
    return(list(xlim = as.numeric(xlim), xticks = as.numeric(xticks)))
  }
  nice <- forestplot_nice_axis(
    xlim,
    exponentiate = exponentiate,
    nullval = nullval
  )
  list(xlim = nice$xlim, xticks = nice$xticks)
}

#' Data-driven x-axis limits for forest plot panels
#' @noRd
resolve_forestplot_limits <- function(panels,
                                      xlim = NULL,
                                      xticks = NULL,
                                      exponentiate = FALSE,
                                      nullval = NULL) {
  if (is.null(nullval)) {
    nullval <- if (isTRUE(exponentiate)) 1 else 0
  }
  if (!is.null(xlim)) {
    if (is.list(xlim)) {
      out_xlim <- vector("list", length(xlim))
      out_ticks <- vector("list", length(xlim))
      ticks_in <- if (is.list(xticks)) xticks else rep(list(xticks), length(xlim))
      for (i in seq_along(xlim)) {
        nice <- forestplot_apply_nice_axis(
          xlim[[i]],
          xticks = ticks_in[[i]],
          exponentiate = exponentiate,
          nullval = nullval
        )
        out_xlim[[i]] <- nice$xlim
        out_ticks[[i]] <- nice$xticks
      }
      return(list(xlim = out_xlim, xticks = out_ticks))
    }
    nice <- forestplot_apply_nice_axis(
      xlim,
      xticks = xticks,
      exponentiate = exponentiate,
      nullval = nullval
    )
    return(list(xlim = nice$xlim, xticks = nice$xticks))
  }

  panel_xlim <- lapply(panels, function(df) {
    vals <- c(df$yi, df$ci.lb, df$ci.ub)
    if ("est_plot" %in% names(df)) {
      vals <- c(vals, df$est_plot)
    } else if ("plot_yi" %in% names(df)) {
      vals <- c(vals, df$plot_yi)
    }
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      return(if (isTRUE(exponentiate)) c(0.25, 4) else c(-1, 1))
    }
    if (isTRUE(exponentiate)) {
      vals <- exp(vals)
      # Drop underflow (exp of huge negative CIs is 0) so imputed missing
      # rows with vi = big_var do not collapse the ratio axis to 0.
      vals <- vals[is.finite(vals) & vals > 1e-6]
      if (!length(vals)) {
        return(c(0.25, 4))
      }
      vals <- c(vals, nullval)
      rng <- range(vals, na.rm = TRUE)
      pad <- 0.08
      return(c(
        max(rng[[1]] * (1 - pad), 0.05),
        rng[[2]] * (1 + pad)
      ))
    }
    vals <- c(vals, nullval)
    rng <- range(vals, na.rm = TRUE)
    pad <- diff(rng) * 0.08
    if (!is.finite(pad) || pad <= 0) {
      pad <- max(0.05 * max(abs(rng), 1e-6), 1e-6)
    }
    c(rng[[1]] - pad, rng[[2]] + pad)
  })

  xlim <- if (length(panel_xlim) == 1L) {
    panel_xlim[[1]]
  } else {
    panel_xlim
  }

  if (is.list(xlim)) {
    ticks_in <- if (is.null(xticks)) {
      rep(list(NULL), length(xlim))
    } else if (is.list(xticks)) {
      xticks
    } else {
      rep(list(xticks), length(xlim))
    }
    out_xlim <- vector("list", length(xlim))
    out_ticks <- vector("list", length(xlim))
    for (i in seq_along(xlim)) {
      nice <- forestplot_apply_nice_axis(
        xlim[[i]],
        xticks = ticks_in[[i]],
        exponentiate = exponentiate,
        nullval = nullval
      )
      out_xlim[[i]] <- nice$xlim
      out_ticks[[i]] <- nice$xticks
    }
    return(list(xlim = out_xlim, xticks = out_ticks))
  }

  nice <- forestplot_apply_nice_axis(
    xlim,
    xticks = xticks,
    exponentiate = exponentiate,
    nullval = nullval
  )
  list(xlim = nice$xlim, xticks = nice$xticks)
}

#' @noRd
forestplot_row_precision <- function(df) {
  n <- nrow(df)
  if (!n) {
    return(numeric(0))
  }
  if ("precision" %in% names(df)) {
    prec <- as.numeric(df$precision)
    if ("vi" %in% names(df)) {
      miss <- !is.finite(prec)
      prec[miss] <- 1 / df$vi[miss]
    }
    return(prec)
  }
  if ("vi" %in% names(df)) {
    return(1 / df$vi)
  }
  rep(NA_real_, n)
}

#' Per-row information fraction: precision / sum(precision) within trial.
#'
#' This is each **subgroup's own** information share for the forest left
#' column. On the reference row it equals \eqn{1-\pi_j}; on the non-reference
#' row it equals the CAMS design \eqn{\pi_j}=`ifrac_<nonref>`. Shrinkage and
#' the model always use the non-reference \eqn{\pi_j}, not the reference share.
#' @noRd
forestplot_row_information_fraction <- function(df,
                                                contrast_sep = "-",
                                                digits = 1) {
  n <- nrow(df)
  if (!n || !"trial" %in% names(df)) {
    return(rep("", n))
  }

  prec <- forestplot_row_precision(df)
  is_component <- rep(TRUE, n)
  if ("subgroup" %in% names(df)) {
    is_component <- is_ifrac_component_row(df$subgroup, contrast_sep = contrast_sep)
  }

  out <- rep("", n)
  for (tr in unique(df$trial)) {
    in_trial <- df$trial == tr
    denom <- sum(prec[in_trial & is_component], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) {
      next
    }
    # Subgroup components only — never print IF on interaction / contrast rows.
    idx <- which(in_trial & is_component)
    vals <- prec[idx] / denom
    out[idx] <- forestplot_format_ifrac(vals, digits = digits)
  }
  out
}

#' @noRd
forestplot_format_ifrac <- function(x, digits = 1) {
  out <- rep("", length(x))
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  pct <- 100 * x[ok]
  out[ok] <- paste0(format(round(pct, digits), nsmall = digits, trim = TRUE), "%")
  out
}

#' @noRd
forestplot_interleave_ifrac_display <- function(ifrac_display) {
  if (!length(ifrac_display)) {
    return(character(0))
  }
  out <- rep("", 2L * length(ifrac_display))
  out[seq(1L, length(out), by = 2L)] <- ifrac_display
  out
}

forestplot_panel_df <- function(df,
                                prob,
                                ifrac_display = NULL) {
  if (!nrow(df)) {
    return(data.frame(
      yi = numeric(0),
      vi = numeric(0),
      ci.lb = numeric(0),
      ci.ub = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  # Prefer precomputed posterior QI (shrinkage greys) over normal approx
  has_qi <- all(c("ci.lb", "ci.ub") %in% names(df)) &&
    any(is.finite(df$ci.lb) & is.finite(df$ci.ub))
  if (has_qi) {
    out <- data.frame(
      yi = as.numeric(df$yi),
      vi = as.numeric(df$vi),
      ci.lb = as.numeric(df$ci.lb),
      ci.ub = as.numeric(df$ci.ub),
      stringsAsFactors = FALSE
    )
    out[is.na(df$yi), c("yi", "ci.lb", "ci.ub")] <- NA
  } else {
    esc <- suppressWarnings(
      suppressMessages(metafor::escalc(yi = df$yi, vi = df$vi))
    )
    summ <- suppressWarnings(
      suppressMessages(summary(esc, level = 100 * prob, transf = identity))
    )
    out <- data.frame(
      yi = summ$yi,
      vi = df$vi,
      ci.lb = summ$ci.lb,
      ci.ub = summ$ci.ub,
      stringsAsFactors = FALSE
    )
    bad_ub <- !is.na(summ$ci.ub) & summ$ci.ub > 195
    out[bad_ub, c("yi", "ci.lb", "ci.ub")] <- NA
    out[is.na(df$yi), c("yi", "ci.lb", "ci.ub")] <- NA
  }

  if (!is.null(ifrac_display)) {
    if (length(ifrac_display) != nrow(out)) {
      stop("`ifrac_display` must have one value per forest row.", call. = FALSE)
    }
    out$ifrac <- ifrac_display
    out$ifrac[is.na(out$yi)] <- ""
  }

  out
}

validate_pooled <- function(pooled, contrast_sep) {
  if (!is.list(pooled) || is.data.frame(pooled)) {
    stop(
      "`pooled` must be a named list with elements `subgroup` and/or `interaction`.",
      call. = FALSE
    )
  }
  nms <- names(pooled)
  if (is.null(nms) || !any(nzchar(nms))) {
    stop(
      "`pooled` must be a named list with elements `subgroup` and/or `interaction`.",
      call. = FALSE
    )
  }
  allowed <- c("subgroup", "interaction")
  unknown <- setdiff(nms[nzchar(nms)], allowed)
  if (length(unknown)) {
    stop(
      "`pooled` has unknown element(s): ",
      paste(unknown, collapse = ", "),
      ". Use `subgroup` and/or `interaction` only.",
      call. = FALSE
    )
  }
  usable <- intersect(allowed, nms)
  usable <- usable[vapply(usable, function(nm) !is.null(pooled[[nm]]), logical(1L))]
  if (!length(usable)) {
    stop(
      "`pooled` must include at least one of `subgroup` or `interaction`.",
      call. = FALSE
    )
  }
  for (nm in usable) {
    validate_forestplot_data(pooled[[nm]])
  }
  pooled
}

append_pooled_rows <- function(resultsA,
                               resultsB,
                               row_labels,
                               pooled_subgroup,
                               pooled_interaction,
                               model_names,
                               prob,
                               study_labels) {
  if (is.null(pooled_subgroup)) {
    return(list(resultsA = resultsA, resultsB = resultsB, row_labels = row_labels))
  }

  num_subgroups <- length(unique(as.character(pooled_subgroup$subgroup)))
  k <- nrow(pooled_subgroup)
  n <- nrow(row_labels)
  additional_rows <- num_subgroups * length(model_names)

  pooled_sg <- forestplot_panel_df(pooled_subgroup, prob = prob)
  rownames(pooled_sg) <- paste(pooled_subgroup$trial, pooled_subgroup$subgroup, sep = "|")

  extra_labels <- data.frame(
    subheading = rep(model_names, each = num_subgroups),
    heading = sub("\\|.*", "", rownames(pooled_sg)),
    label = sub(".*\\|", "", rownames(pooled_sg)),
    stringsAsFactors = FALSE
  )

  row_labels <- rbind(
    row_labels,
    extra_labels[seq_len(additional_rows), , drop = FALSE]
  )

  pad <- function(df, n_add) {
    if (!n_add) {
      return(df)
    }
    blank <- df[1, , drop = FALSE]
    for (col in names(blank)) {
      if (col %in% c("ifrac", "trial_size")) {
        blank[[col]] <- ""
      } else if (is.character(blank[[col]])) {
        blank[[col]] <- ""
      } else if (is.logical(blank[[col]])) {
        blank[[col]] <- FALSE
      } else {
        blank[[col]] <- NA
      }
    }
    rbind(df, blank[rep(1L, n_add), , drop = FALSE])
  }

  resultsA <- pad(resultsA, additional_rows)
  resultsB <- pad(resultsB, additional_rows)
  resultsA[(n + 1):(n + additional_rows), ] <- pooled_sg[seq_len(additional_rows), ]

  if (!is.null(pooled_interaction) && nrow(pooled_interaction)) {
    rows_to_update <- which(extra_labels$subheading %in% model_names)
    if (length(rows_to_update)) {
      ref_idx <- 1L
      include_indices <- seq(rows_to_update[1], rows_to_update[k], num_subgroups) + ref_idx - 1L
      int_summ <- forestplot_panel_df(pooled_interaction, prob = prob)
      n_int <- length(include_indices)
      if (nrow(int_summ) >= 2 * n_int) {
        resultsB[include_indices, c("yi", "ci.lb", "ci.ub")] <-
          int_summ[seq(2, 2 * n_int, 2), c("yi", "ci.lb", "ci.ub")]
      }
    }
  } else {
    resultsB[(n + 1):(n + additional_rows), ] <- resultsA[(n + 1):(n + additional_rows), ]
  }

  list(resultsA = resultsA, resultsB = resultsB, row_labels = row_labels)
}
