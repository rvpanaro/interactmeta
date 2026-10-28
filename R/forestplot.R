# ---- forestplot.R ----
#' Forest plot for a fitted model
#'
#' Draws a forest for a [cams()] fit using [bayesmeta::forestplot.bayesmeta]
#' (single-block) or [forestplot::forestplot()] with the same conventions
#' (multi-block): observed study estimates, optional grey **shrinkage**
#' overlays, and **Overall** / **Prediction** / **Combined** summaries. An
#' **information fraction** column shows each row's share of trial precision.
#'
#' For non-[interactmeta] input (e.g. [metafor::escalc]), arguments are forwarded
#' to [forestplot::forestplot()], which dispatches to
#' [bayesmeta::forestplot.escalc] when **bayesmeta** is loaded (interactmeta loads it
#' on attach).
#'
#' By default (`type = "modifiers"`) the plot shows **modifier-level** panels
#' (each with **Overall** and **Prediction**), with overall-IF notes under the
#' x-axis. Use `type = "interaction"` for the interaction block only.
#' `type = "subgroups"` is deprecated (same as `"modifiers"`).
#'
#' @param x A CAMS fit from [cams()] (`interactmeta`), a construction fit from
#'   [da()], [ad()], [swada()], or [conventional()], or another object handled by
#'   [forestplot::forestplot()] (e.g. `escalc`).
#' @param n_var Column name for per-study patient counts shown in an **N**
#'   column. Default `FALSE` (omit). Set to `"n"` or another column name to
#'   show sizes. When `"n"` is missing but `treat.total` and `cont.total` are
#'   present (2-by-2 event-count data), study size is taken as their sum.
#' @param contrast_sep Separator between subgroup labels on interaction rows
#'   (default `"-"`).
#' @param model_names Optional character vector selecting marginal layers among
#'   `"Shrinkage"`, `"Overall"`, and `"Prediction"` (intersected with the
#'   `include_*` flags). Prefer `include_shrinkage` / `include_overall` /
#'   `include_prediction` for new code. When `NULL` or the default `"Model 1"`,
#'   layers follow those flags alone.
#' @param include_shrinkage,include_overall,include_prediction,include_combined
#'   Include shrinkage, overall-effect, prediction, and combined-population rows.
#'   All four layers default to `TRUE`. Grey
#'   shrinkage intervals appear alongside quoted estimates (bayesmeta style)
#'   in the active `type` block. Shrinkage is included by default for Bayesian
#'   fits only; frequentist (`REML` / `ML`) forests omit shrinkage rows.
#'   Set `include_shrinkage = FALSE` to omit them on Bayesian fits.
#'   For `type = "modifiers"`, `include_combined` controls a trailing
#'   **All** section (**Overall** / **Prediction** for [combined()] at
#'   the reporting IF), not the subgroup panels. Switching off overall or
#'   prediction keeps blank placeholder rows so the figure height and spacing
#'   stay close to the full plot.
#' @param type Which forest block(s) to draw. One of:
#'   - `"modifiers"` (default): one section per modifier level (reference and
#'     non-reference), each with **Overall** and **Prediction** when enabled;
#'     optional trailing **All** summaries when `include_combined = TRUE`
#'   - `"interaction"`: Interaction block only (Overall + Prediction)
#'
#'   Other features (favours, annotations, `include_shrinkage`, …) apply
#'   within the selected type.
#' @param ifrac How the overall information fraction ω is pooled via [ifrac()]
#'   for **annotations** and for **Overall / Prediction / Combined** diamonds:
#'   `"interaction"` (default; τ_γ and subgroup sampling variances) or
#'   `"overall"` (τ and design SEs). Shrinkage greys are the trial-specific
#'   CAMS reconstructions from [shrinkage()] (shrunken \eqn{m_j,g_j} inverted
#'   at each trial's own \eqn{\pi_j}); they are not population fitted values.
#'   Observed trial rows remain at each trial's own π_j.
#' @param xlab Character; x-axis title. Default is
#'   `"Effect (log scale)"` when `exponentiate = FALSE`, or
#'   `"Effect"` when `exponentiate = TRUE` (ratio-scale ticks via
#'   `xlog = TRUE`).
#' @param group How to nest trial rows. `"subgroup"` (default) groups by
#'   subgroup level with trial names as row labels; `"trial"` groups by trial
#'   with subgroup names as row labels.
#' @param show_annotations If `TRUE` (default), annotate under the x-axis with
#'   overall information-fraction quantiles and heterogeneity credibility
#'   intervals.
#' @param annotation_digits Decimal places in annotation text (default `2`).
#' @param favours Optional length-2 character labels placed under the x-axis
#'   (left / right of the null). Default is `c("Control", "Experimental")`
#'   for `type = "modifiers"`. For `type = "interaction"`, the
#'   default names the reference and non-reference subgroups (same left/right
#'   roles), and the IF column header is `IF (non-reference subgroup)`.
#'   Set to `FALSE` to suppress favours.
#' @param show_ci_text Logical; show estimate and CI text. Default `TRUE`.
#' @param exponentiate Logical; if `TRUE`, show estimates on the ratio scale.
#' @param col_right_heading Column heading for the effect column.
#' @param prob Coverage probability for confidence intervals (default `0.95`).
#' @param study_fill,study_colour Optional overrides for study-row colours.
#' @param base_size,pointsize,digits,stroke,fill,summary_fill,summary_colour,which_panel
#'   Plotting options passed through to the bayesmeta / forestplot renderer.
#' @param estimates_only If `TRUE`, drop observed trial rows (keep shrinkage /
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
#' @param reference,label,modifier,subgroup,trial,contrast,big_var,pooled
#'   Ignored for CAMS fits; retained for back-compatibility.
#' @param ... Further arguments passed through to [forestplot::forestplot()]
#'   when the plot is drawn (same escape hatch as
#'   [bayesmeta::forestplot.bayesmeta]'s `...`). For non-CAMS `x`, unused CAMS
#'   formals are ignored and `...` is forwarded to [forestplot::forestplot()].
#'
#' Construction fits use the same layout (study rows, shrinkage, Overall /
#' Prediction, optional Combined) without CAMS information fractions.
#' [ad()] uses `type = "interaction"` by default; only the interaction block
#' is available (`type = "modifiers"` errors).
#'
#' @return For a CAMS or construction fit, a `interactmeta_forestplot` object with attributes
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
#' fit <- cams(e, reference = "Lower")
#' forestplot(fit, exponentiate = TRUE, type = "modifiers")
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
                       pooled = NULL,
                       model_names = "Model 1",
                       include_shrinkage = TRUE,
                       include_overall = TRUE,
                       include_prediction = TRUE,
                       include_combined = TRUE,
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
                       type = c("modifiers", "interaction", "subgroups"),
                       which_panel = NULL,
                       ...) {
  if (is.character(modifier) && length(modifier) == 1L && nzchar(modifier)) {
    subgroup <- modifier
  }
  if (inherits(x, "interactmeta") ||
        interactmeta_is_non_cams_analysis(x) ||
        inherits(x, interactmeta_non_cams_classes())) {
    include_shrinkage <- forestplot_resolve_include_shrinkage(
      fit = x,
      include_shrinkage = include_shrinkage,
      call = match.call(expand.dots = FALSE)
    )
  }
  if (interactmeta_is_non_cams_analysis(x) ||
        inherits(x, interactmeta_non_cams_classes())) {
    check_forestplot_deps()
    dots <- list(...)
    typography_auto <- isTRUE(dots$typography_auto)
    dots$typography_auto <- NULL
    panel_height_auto <- missing(panel.height) || is.null(panel.height)
    group <- match.arg(group, c("subgroup", "trial"))
    type <- forestplot_resolve_fit_type(
      fit = x,
      type = type,
      call = match.call(expand.dots = FALSE)
    )
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
    construction_args <- list(
      fit = x,
      prob = prob,
      contrast_sep = contrast_sep,
      include_shrinkage = include_shrinkage,
      include_overall = include_overall,
      include_prediction = include_prediction,
      include_combined = include_combined,
      exponentiate = exponentiate,
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
      type = if (identical(type, "subgroups")) "modifiers" else type,
      which_panel = which_panel
    )
    if (!identical(model_names, "Model 1")) {
      construction_args$model_names <- model_names
    }
    if (!panel_height_auto) {
      construction_args$panel.height <- panel.height
    }
    if (!typography_auto) {
      construction_args$base_size <- base_size
      construction_args$pointsize <- pointsize
    }
    dots$printplot <- NULL
    dots$quiet <- NULL
    for (nm in names(construction_args)) {
      dots[[nm]] <- NULL
    }
    return(do.call(construction_forestplot_from_fit, c(construction_args, dots)))
  }
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
  fcall <- match.call(expand.dots = FALSE)
  type_user <- if (is.null(fcall$type)) {
    "modifiers"
  } else {
    forestplot_parse_panel_type(type)
  }
  if (identical(type_user, "subgroups")) {
    .Deprecated(msg = 'type = "subgroups" is deprecated; use type = "modifiers".')
  }
  type <- forestplot_resolve_type(type_user)
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
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
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
    type = if (identical(type, "subgroups")) "modifiers" else type,
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
#' @rdname forestplot
#' @exportS3Method forestplot::forestplot interactmeta
forestplot.interactmeta <- function(x, ...) {
  check_comparison_fit(x)
  do.call(forestplot, c(list(x = x), list(...)))
}

#' @rdname forestplot
#' @exportS3Method forestplot::forestplot interactmeta_da
forestplot.interactmeta_da <- function(x, ...) {
  forestplot(x, ...)
}

#' @rdname forestplot
#' @exportS3Method forestplot::forestplot interactmeta_ad
forestplot.interactmeta_ad <- function(x, ...) {
  forestplot(x, ...)
}

#' @rdname forestplot
#' @exportS3Method forestplot::forestplot interactmeta_swada
forestplot.interactmeta_swada <- function(x, ...) {
  forestplot(x, ...)
}

#' @rdname forestplot
#' @exportS3Method forestplot::forestplot interactmeta_conventional
forestplot.interactmeta_conventional <- function(x, ...) {
  forestplot(x, ...)
}

#' @rdname forestplot
#' @exportS3Method metafor::forest interactmeta_da
forest.interactmeta_da <- function(x, ...) {
  forestplot(x, ...)
}

#' @rdname forestplot
#' @exportS3Method metafor::forest interactmeta_ad
forest.interactmeta_ad <- function(x, ...) {
  forestplot(x, ...)
}

#' @rdname forestplot
#' @exportS3Method metafor::forest interactmeta_swada
forest.interactmeta_swada <- function(x, ...) {
  forestplot(x, ...)
}

#' @rdname forestplot
#' @exportS3Method metafor::forest interactmeta_conventional
forest.interactmeta_conventional <- function(x, ...) {
  forestplot(x, ...)
}

#' Forest plot method for CAMS fits ([metafor::forest] generic)
#'
#' Metafor-idiomatic entry point for a [cams()] fit. Delegates to
#' [forestplot()] with the same arguments so CAMS options (`type`,
#' `group`, …) share one implementation.
#'
#' @param x A CAMS fit from [cams()].
#' @param ... Passed to [forestplot()] (CAMS options and/or forestplot
#'   presentation args).
#'
#' @return A `interactmeta_forestplot` object (see [forestplot()]).
#'
#' @examples
#' \dontrun{
#' fit <- cams(escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss),
#'             reference = "Lower")
#' forest(fit, exponentiate = TRUE)
#' }
#'
#' @seealso [forestplot()], [metafor::forest()]
#' @rdname forestplot
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
                                         type = c("modifiers", "interaction", "subgroups"),
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
                                         include_combined = TRUE,
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
  if ("include_overall" %in% names(dots)) {
    dots <- dots[!names(dots) %in% "include_overall"]
  }
  if ("include_combined" %in% names(dots)) {
    include_combined <- dots$include_combined
    dots <- dots[!names(dots) %in% "include_combined"]
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
  type <- forestplot_resolve_type(as.character(type)[[1L]])
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
    include_combined = include_combined,
    group = group
  )
  if (length(dots)) {
    render_snapshot <- c(render_snapshot, dots)
  }

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

  layer_on <- comparison_forest_layer_on
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
  # after stacking once block headings (Younger / Older / Interaction / Combined)
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
    any(resultsA$row_kind %in% c("overall", "prediction", "combined"), na.rm = TRUE)
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
    type = type,
    include_combined = include_combined
  )
  resultsA <- filtered$resultsA
  resultsB <- filtered$resultsB
  row_labels <- filtered$row_labels
  # IF cells for the non-reference subgroup when `type = "interaction"`.
  resultsA <- forestplot_restrict_ifrac_to_subgroup(
    resultsA = resultsA,
    row_labels = row_labels,
    subgroup = nonref_subgroup
  )
  n_study <- sum(resultsA$row_kind == "study", na.rm = TRUE)
  fit_analysis <- dots$fit_analysis
  if ("fit_analysis" %in% names(dots)) {
    dots$fit_analysis <- NULL
  }
  annotations <- forestplot_filter_annotations_type(
    annotations,
    type = type,
    fit_analysis = fit_analysis
  )
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

#' Total patients per trial from subgroup rows
#' @noRd
forestplot_trial_n_lookup <- function(dat,
                                      trial_col = NULL,
                                      subgroup_col = NULL,
                                      n_var = FALSE,
                                      contrast_sep = "-") {
  trial_col <- trial_col %||% interactmeta_resolve_study_col(dat)
  subgroup_col <- subgroup_col %||% interactmeta_resolve_level_col(dat)
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

#' IF column header (`IF` or `IF (non-reference subgroup)`)
#' @noRd
forestplot_if_column_heading <- function(nonref_subgroup = NULL) {
  sg <- as.character(nonref_subgroup %||% "")[[1L]]
  if (!nzchar(sg)) {
    "IF"
  } else {
    paste0("IF (", sg, ")")
  }
}

#' @noRd
forestplot_col_left_headings <- function(panel_df, nonref_subgroup = NULL) {
  heads <- character(0)
  if ("ifrac" %in% names(panel_df)) {
    heads <- c(heads, forestplot_if_column_heading(nonref_subgroup))
  }
  if ("trial_size" %in% names(panel_df)) {
    heads <- c(heads, "N")
  }
  heads
}

#' Blank IF cells that are not for the named subgroup
#'
#' When `subgroup` is set (interaction / overall forests), keep IF only on
#' rows belonging to that subgroup; Interaction / Combined study rows already
#' carry the non-ref IF and are retained. With `type = "interaction"`, the
#' column header is `IF (non-reference subgroup)`; subgroup forests use `"IF"`.
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
  # Interaction / Combined rows already store the quoted non-ref IF.
  is_int_ov <- forestplot_heading_is_interaction(heading) |
    forestplot_heading_is_combined(heading)
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
  # Drop Interaction / Combined block tags so only subgroup levels remain.
  cl <- cl[
    nzchar(cl) &
      !forestplot_heading_is_interaction(cl) &
      !forestplot_heading_is_combined(cl)
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
  study_col <- interactmeta_resolve_study_col(data)
  level_col <- interactmeta_resolve_level_col(data)
  req <- c(study_col, level_col, "yi", "vi")
  missing <- setdiff(req, names(data))
  if (length(missing)) {
    stop(
      "`data` is missing required columns: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  # Internal forest code expects `trial` / `subgroup` (metafor-era names).
  if (!"trial" %in% names(data)) {
    data$trial <- data[[study_col]]
  }
  if (!"subgroup" %in% names(data)) {
    data$subgroup <- data[[level_col]]
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
    if (!is.null(nms) && length(nms)) {
      col <- interactmeta_resolve_level_col(data, fallback = NULL)
      if (col %in% nms) {
        return(col)
      }
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

split_ifrac_for_forest <- function(data, contrast_sep, group = "subgroup") {
  group <- match.arg(group, c("subgroup", "trial"))
  level_col <- interactmeta_resolve_level_col(data)
  is_contrast <- is_ifrac_contrast_row(data[[level_col]], contrast_sep = contrast_sep)
  subgroup <- data[!is_contrast, , drop = FALSE]
  interaction <- data[is_contrast, , drop = FALSE]

  if (nrow(interaction) == 0L && nrow(subgroup) > 0L) {
    interaction <- build_contrast_rows_from_subgroups(
      subgroup,
      contrast_sep = contrast_sep,
      reference = resolve_forest_reference(subgroup, data)
    )
  }

  study_col <- interactmeta_resolve_study_col(subgroup)
  subgroup <- order_forest_rows(subgroup, group = group)
  interaction <- order_forest_rows(
    interaction,
    group = group,
    by_trial = unique(as.character(subgroup[[study_col]]))
  )

  list(subgroup = subgroup, interaction = interaction)
}

resolve_forest_reference <- function(subgroup, data, reference = NULL) {
  # Prefer an explicit / stored reference; never silently invent one for CAMS.
  if (!is.null(reference) && nzchar(as.character(reference)[[1L]])) {
    return(interactmeta_reference(
      if (is.data.frame(data)) data else subgroup,
      reference = reference,
      subgroup_var = "level"
    ))
  }
  if (is.data.frame(data) && "reference" %in% names(data)) {
    return(interactmeta_reference(data, subgroup_var = "level"))
  }
  if (is.data.frame(subgroup) && "reference" %in% names(subgroup)) {
    return(interactmeta_reference(subgroup, subgroup_var = "level"))
  }
  stop(
    "Cannot resolve forest reference: supply `reference` or a `reference` column.",
    call. = FALSE
  )
}

build_contrast_rows_from_subgroups <- function(subgroup, contrast_sep, reference) {
  ref <- reference[[1]]
  study_col <- interactmeta_resolve_study_col(subgroup)
  level_col <- interactmeta_resolve_level_col(subgroup)
  trials <- unique(subgroup[[study_col]])
  others <- setdiff(sort(unique(as.character(subgroup[[level_col]]))), ref)
  if (!length(others)) {
    return(subgroup[0, , drop = FALSE])
  }

  wide <- subgroup |>
    dplyr::select(dplyr::all_of(study_col), dplyr::all_of(level_col), "yi", "vi") |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(level_col),
      values_from = c("yi", "vi"),
      names_glue = sprintf("{.value}_{%s}", level_col)
    )

  out <- lapply(trials, function(tr) {
    row <- wide[wide[[study_col]] == tr, , drop = FALSE]
    if (!nrow(row)) {
      return(NULL)
    }
    dplyr::bind_rows(lapply(others, function(lab) {
      yi1 <- row[[paste0("yi_", lab)]]
      yi2 <- row[[paste0("yi_", ref)]]
      vi1 <- row[[paste0("vi_", lab)]]
      vi2 <- row[[paste0("vi_", ref)]]
      row_out <- dplyr::tibble(
        yi = yi1 - yi2,
        vi = vi1 + vi2
      )
      row_out[[study_col]] <- tr
      row_out[[level_col]] <- paste0(lab, contrast_sep, ref)
      row_out
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

  trial_col <- if ("trial" %in% names(subgroup)) {
    "trial"
  } else {
    interactmeta_resolve_study_col(subgroup)
  }
  level_col <- if ("subgroup" %in% names(subgroup)) {
    "subgroup"
  } else {
    interactmeta_resolve_level_col(subgroup)
  }
  if (!trial_col %in% names(interaction)) {
    interaction$trial <- interaction[[interactmeta_resolve_study_col(interaction)]]
  }

  int_one <- interaction |>
    dplyr::group_by(.data[[trial_col]]) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(c(trial_col, "yi", "vi")))

  joined <- dplyr::left_join(
    subgroup |> dplyr::select(dplyr::all_of(c(trial_col, level_col))),
    int_one,
    by = trial_col
  )
  names(joined)[names(joined) == level_col] <- "subgroup"
  if (!identical(trial_col, "trial")) {
    joined$trial <- joined[[trial_col]]
  }

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
  trial_col <- if ("trial" %in% names(df)) {
    "trial"
  } else {
    interactmeta_resolve_study_col(df)
  }
  level_col <- if ("subgroup" %in% names(df)) {
    "subgroup"
  } else {
    interactmeta_resolve_level_col(df)
  }
  if (!trial_col %in% names(df)) {
    stop("Forest data is missing a study / trial column.", call. = FALSE)
  }
  if (!level_col %in% names(df)) {
    stop("Forest data is missing a level / subgroup column.", call. = FALSE)
  }
  # Subgroup nesting uses alphabetical trial labels within each block (classic
  # metafor forest ordering). Trial nesting keeps first-appearance block order
  # unless `by_trial` is supplied.
  trial_ord <- if (is.null(by_trial)) {
    tr <- unique(as.character(df[[trial_col]]))
    if (identical(group, "subgroup")) {
      sort(tr)
    } else {
      tr
    }
  } else {
    by_trial
  }
  df[[trial_col]] <- factor(as.character(df[[trial_col]]), levels = trial_ord)
  sg_was_factor <- is.factor(df[[level_col]])
  sg_levels <- if (sg_was_factor) {
    levels(df[[level_col]])
  } else {
    unique(as.character(df[[level_col]]))
  }
  df[[level_col]] <- factor(as.character(df[[level_col]]), levels = sg_levels)

  if (identical(group, "subgroup")) {
    df <- df[order(df[[level_col]], df[[trial_col]]), , drop = FALSE]
  } else {
    df <- df[order(df[[trial_col]], df[[level_col]]), , drop = FALSE]
  }
  df[[trial_col]] <- as.character(df[[trial_col]])
  if (sg_was_factor) {
    df[[level_col]] <- factor(as.character(df[[level_col]]), levels = sg_levels)
  } else {
    df[[level_col]] <- as.character(df[[level_col]])
  }
  if (!identical(trial_col, "trial") && !"trial" %in% names(df)) {
    df$trial <- df[[trial_col]]
  }
  if (!identical(level_col, "subgroup") && !"subgroup" %in% names(df)) {
    df$subgroup <- df[[level_col]]
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

#' Non-reference subgroup component used for Interaction / Combined rows
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

#' Non-reference subgroup label (for Interaction / Combined headers)
#' @noRd
forestplot_nonref_subgroup <- function(subgroup_vals, ref = NULL) {
  forestplot_nonref_component(
    ifrac_vals = rep("", length(subgroup_vals)),
    subgroup_vals = subgroup_vals,
    ref = ref
  )$subgroup
}

#' Forest block heading for the combined population ([combined()] estimand)
#' @noRd
forestplot_combined_block_label <- function() {
  "All"
}

#' Block heading with optional non-reference subgroup, e.g. `"All (Older)"`
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

#' Whether a forest heading is a Combined block (possibly with subgroup)
#' @noRd
forestplot_heading_is_combined <- function(heading) {
  lab <- forestplot_combined_block_label()
  h <- as.character(heading)
  grepl(paste0("^(", lab, "|Combined)(\\s*\\(|$)"), h)
}

#' Stack subgroup, interaction, and combined effects into one vertical forest panel
#'
#' Order: subgroup blocks (study / shrinkage / overall / prediction), then an
#' Interaction block (one row per trial + optional shrinkage + summaries), then
#' Combined (**Overall** / **Prediction** for [combined()] at the reporting IF
#' only; no trial-level rows in that block).
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
    rk %in% c("study", "shrinkage", "overall", "prediction") & !is_int_shr
  )
  if ("colour_label" %in% names(row_labels)) {
    cl <- as.character(row_labels$colour_label)
    int_marg <- rk %in% c("overall", "prediction") &
      !is.na(cl) &
      cl == "Interaction"
    sg_ix <- sg_ix[!int_marg[sg_ix]]
  }
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
      # Match Combined block aesthetics: filled study squares (not open circles).
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
      # Match Combined shrinkage marker (filled square; grey fill applied earlier).
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

  for (kind in c("overall", "prediction")) {
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
    # Same labels as the Combined block (Overall diamond, Prediction bar).
    rl$label <- if (identical(kind, "overall")) "Overall" else "Prediction"
    if ("colour_label" %in% names(rl)) {
      rl$colour_label <- int_header
    }
    push_row(row, rl)
  }

  ov_header <- forestplot_combined_block_label()
  # Include NA Combined placeholders so omitted pooled/prediction keep space.
  ov_ix <- which(rk == "combined")

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
  tok_tau_hat <- "<<TAUHAT>>"
  tok_tg_hat <- "<<TGHAT>>"
  tok_pi <- "<<PI>>"
  tok_dot <- "<<DOT>>"
  txt <- gsub(
    paste0("at fitted interaction heterogeneity ", "\u03c4", "_", "\u03b3"),
    tok_tg_hat,
    txt,
    fixed = TRUE
  )
  txt <- gsub("at fitted interaction heterogeneity tau_gamma", tok_tg_hat, txt, fixed = TRUE)
  txt <- gsub(paste0("at fitted heterogeneity ", "\u03c4"), tok_tau_hat, txt, fixed = TRUE)
  txt <- gsub("at fitted heterogeneity tau", tok_tau_hat, txt, fixed = TRUE)
  txt <- gsub(paste0("at fitted ", "\u03c4"), tok_tau_hat, txt, fixed = TRUE)
  txt <- gsub("at fitted tau", tok_tau_hat, txt, fixed = TRUE)
  txt <- gsub("\u03c4_\u03b3", tok_tg, txt, fixed = TRUE)
  txt <- gsub("\u03c4_y", tok_tg, txt, fixed = TRUE)
  txt <- gsub("tau_gamma", tok_tg, txt, fixed = TRUE)
  txt <- gsub("\u03c4", tok_tau, txt, fixed = TRUE)
  txt <- gsub("\u03c0", tok_pi, txt, fixed = TRUE)
  txt <- gsub(" \u00b7 ", tok_dot, txt, fixed = TRUE)
  txt <- gsub("\u00b7", tok_dot, txt, fixed = TRUE)
  txt <- gsub("\\", "\\\\", txt, fixed = TRUE)
  txt <- gsub("\"", "\\\"", txt, fixed = TRUE)
  txt <- gsub(tok_tg_hat, '" * hat(tau[gamma]) * "', txt, fixed = TRUE)
  txt <- gsub(tok_tau_hat, '" * hat(tau) * "', txt, fixed = TRUE)
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

#' Default / effective `include_shrinkage` for forest plots
#'
#' Shrinkage rows are Bayesian-only; frequentist fits never show them.
#' @noRd
forestplot_resolve_include_shrinkage <- function(fit,
                                                   include_shrinkage,
                                                   call) {
  if (!interactmeta_is_bayesian_fit(fit)) {
    return(FALSE)
  }
  if (is.null(call$include_shrinkage)) {
    return(TRUE)
  }
  isTRUE(include_shrinkage)
}

#' Resolve `forestplot()` `type` for construction / CAMS fits
#'
#' AD forests only support the interaction block; when `type` is omitted,
#' `"interaction"` is used. Explicit `type = "modifiers"` still errors.
#' @noRd
forestplot_resolve_fit_type <- function(
    fit,
    type = c("modifiers", "interaction", "subgroups"),
    call = NULL) {
  if (is.null(call)) {
    call <- match.call(expand.dots = FALSE)
  }
  cap <- interactmeta_fit_capabilities(fit)
  subgroup_panels <- isTRUE(cap$subgroups)
  if (!subgroup_panels && is.null(call$type)) {
    return("interaction")
  }
  type_user <- if (is.null(call$type)) {
    "modifiers"
  } else {
    forestplot_parse_panel_type(type)
  }
  if (identical(type_user, "subgroups")) {
    .Deprecated(msg = 'type = "subgroups" is deprecated; use type = "modifiers".')
  }
  type_resolved <- forestplot_resolve_type(type_user)
  forestplot_assert_panel_type_for_fit(
    fit = fit,
    type_resolved = type_resolved,
    type_user = type_user,
    subgroup_panels = subgroup_panels
  )
  if (identical(type_resolved, "subgroups")) {
    "modifiers"
  } else {
    type_resolved
  }
}

#' Error when `forestplot()` `type` is not implemented for this analysis
#' @noRd
forestplot_assert_panel_type_for_fit <- function(fit,
                                                 type_resolved,
                                                 type_user,
                                                 subgroup_panels = NULL) {
  if (is.null(subgroup_panels)) {
    subgroup_panels <- isTRUE(interactmeta_fit_capabilities(fit)$subgroups)
  }
  if (subgroup_panels || identical(type_resolved, "interaction")) {
    return(invisible(type_resolved))
  }
  analysis <- interactmeta_fit_analysis(fit)
  acr <- interactmeta_analysis_acronym(analysis)
  stop(
    acr,
    " forests support only `type = \"interaction\"` ",
    "(the interaction block). ",
    "`type = \"modifiers\"` draws subgroup modifier panels and the optional ",
    "**All** section (`include_combined`); those are not defined for ",
    acr,
    ". You passed `type = ", shQuote(type_user), "`.",
    call. = FALSE
  )
}

#' Parse `forestplot()` `type` (reject marginal-layer names and unknown values)
#' @noRd
forestplot_parse_panel_type <- function(type) {
  if (missing(type) || !length(type)) {
    return("modifiers")
  }
  raw <- as.character(type)[[1L]]
  if (identical(raw, "subgroup")) {
    raw <- "modifiers"
  }
  if (identical(raw, "combined")) {
    stop(
      "`type = \"combined\"` is no longer supported in `forestplot()`. ",
      "Use `type = \"modifiers\"` with `include_combined = TRUE` (default) ",
      "for the combined-population **All** section, or ",
      "`include_combined = FALSE` to omit it.",
      call. = FALSE
    )
  }
  allowed <- c("modifiers", "interaction", "subgroups")
  if (!raw %in% allowed) {
    stop(
      "`type` must be `\"modifiers\"` or `\"interaction\"`. ",
      "You passed ", shQuote(raw), ".",
      call. = FALSE
    )
  }
  raw
}

#' Resolve `forestplot()` `type` (subgroups or interaction only)
#' @noRd
forestplot_resolve_type <- function(type) {
  type <- forestplot_parse_panel_type(type)
  if (identical(type, "interaction")) {
    return("interaction")
  }
  if (type %in% c("modifiers", "subgroups")) {
    return("subgroups")
  }
  stop(
    "`type` must be `\"modifiers\"` or `\"interaction\"`.",
    call. = FALSE
  )
}

#' Keep only the requested forest type block(s) after stacking
#' @noRd
forestplot_filter_type <- function(resultsA,
                                       resultsB,
                                       row_labels,
                                       type = c("modifiers", "interaction", "subgroups"),
                                       include_combined = TRUE) {
  type <- forestplot_resolve_type(type)
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
    keep <- !forestplot_heading_is_interaction(heading)
    if (!comparison_forest_layer_on(include_combined)) {
      keep <- keep & !forestplot_heading_is_combined(heading)
    }
    keep
  } else {
    forestplot_heading_is_interaction(heading)
  }
  if (!any(idx)) {
    stop(
      sprintf(
        "`type = \"%s\"` produced an empty forest. ",
        type
      ),
      "Check `include_overall` / `include_prediction` / `include_combined` ",
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
                                                  type = c("modifiers", "interaction", "subgroups"),
                                                  fit_analysis = NULL) {
  type <- forestplot_resolve_type(type)
  annotations <- forestplot_normalize_annotations(annotations)
  # `subgroups` keeps the full footnote stack (overall IF + both heterogeneity
  # notes). Interaction keeps its slice; DA interaction uses subgroup-specific
  # tau from the separate subgroup RE-MAs.
  if (identical(type, "subgroups")) {
    return(annotations)
  }
  het <- annotations$heterogeneity
  het_key <- if (identical(fit_analysis, "da") && identical(type, "interaction")) {
    "combined"
  } else {
    "interaction"
  }
  het <- het[intersect(names(het), het_key)]
  annotations$heterogeneity <- het
  annotations$ifrac <- character()
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
