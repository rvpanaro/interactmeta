# ---- forestplot_comparison.R ----
#' Forest plot for a fitted CAMS model
#'
#' Layers observed study rows with optional grey **shrinkage** from
#' [shrinkage()], then **Overall** / **Prediction** summaries. The Combined
#' section (when enabled) shows only those population summaries for [combined()],
#' not trial-level rows. Drawn via the bayesmeta / forestplot
#' backend (see [forestplot.interactmeta]).
#'
#' @param fit A CAMS fit from [cams()] (`interactmeta`).
#' @param model_names Optional character vector of model row labels. When `NULL`
#'   (default), built from `include_shrinkage`, `include_overall`, and
#'   `include_prediction`.
#' @param include_shrinkage If `TRUE` (default), interleave grey shrinkage rows from
#'   [shrinkage()].
#' @param include_overall If `TRUE` (default), append overall-effect diamonds.
#' @param include_prediction If `TRUE` (default), append prediction marginal rows.
#' @param include_combined If `TRUE` (default), append combined-population rows from
#'   [combined()].
#' @details Panel rows carry `row_kind`: `"study"`, `"shrinkage"`, `"overall"`,
#'   `"prediction"`, and `"combined"`.
#' @param prob Coverage probability for intervals.
#' @param ifrac How the overall information fraction ω is pooled via [ifrac()]
#'   for annotations and for Overall / Prediction / Combined diamonds:
#'   `"interaction"` (default) or `"overall"`. Shrinkage greys use each trial's
#'   own \eqn{\pi_j}; see [shrinkage()].
#' @param contrast_sep Separator between subgroup labels on interaction rows.
#' @param max_ifrac_draws Optional cap for unweighted posterior draws.
#'   Leave `NULL` (default) for Bayesian DIRECT supports, which cannot be thinned.
#' @param exponentiate Logical; if `TRUE`, estimates are shown on the ratio scale.
#' @param n_var Column for study size labels; `FALSE` (default) omits them.
#' @param group `"subgroup"` (default) or `"trial"` nesting.
#' @param show_annotations If `TRUE` (default), annotate under the x-axis with
#'   overall-IF quantiles and heterogeneity credibility intervals.
#' @param annotation_digits Decimal places in annotation text (default `2`).
#' @param ... Passed to the forestplot renderer (see [forestplot.interactmeta]).
#'   Also accepts options such as `xlim`.
#'
#' @return A `interactmeta_forestplot` object (see [forestplot.interactmeta]).
#'
#' @noRd
comparison_forestplot_from_fit <- function(fit,
                                    model_names = NULL,
                                    include_shrinkage = TRUE,
                                    include_overall = TRUE,
                                    include_prediction = TRUE,
                                    include_combined = TRUE,
                                    exponentiate = FALSE,
                                    n_var = FALSE,
                                    group = "subgroup",
                                    show_annotations = TRUE,
                                    annotation_digits = 2,
                                    prob = 0.95,
                                    ifrac = c("interaction", "overall"),
                                    contrast_sep = "-",
                                    max_ifrac_draws = NULL,
                                    ...) {
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be a CAMS fit from `cams()`.", call. = FALSE)
  }

  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st) || is.null(st$ifrac_data)) {
    stop(
      "Fit has no `ifrac_data` in attribute `\"interactmeta\"`; refit with `cams()`.",
      call. = FALSE
    )
  }

  dots <- list(...)
  if ("exponentiate" %in% names(dots)) {
    exponentiate <- dots$exponentiate
    dots$exponentiate <- NULL
  }
  exponentiate <- isTRUE(exponentiate)
  type <- "modifiers"
  if ("type" %in% names(dots)) {
    type <- forestplot_resolve_type(as.character(dots$type)[[1L]])
    dots$type <- NULL
  } else {
    type <- forestplot_resolve_type(type)
  }
  if ("result" %in% names(dots)) {
    warning("`result` is deprecated; use `type` instead.", call. = FALSE)
    type <- forestplot_resolve_type(as.character(dots$result)[[1L]])
    dots$result <- NULL
  }
  # Legacy printplot / quiet ignored — draw via print() / forestplot_draw().
  dots$printplot <- NULL
  dots$quiet <- NULL
  if ("n_var" %in% names(dots)) {
    n_var <- dots$n_var
    dots$n_var <- NULL
  }
  if ("group" %in% names(dots)) {
    group <- dots$group
    dots$group <- NULL
  }
  if ("show_annotations" %in% names(dots)) {
    show_annotations <- dots$show_annotations
    dots$show_annotations <- NULL
  }
  if ("show_heterogeneity_footnote" %in% names(dots)) {
    # Back-compat alias
    show_annotations <- isTRUE(dots$show_heterogeneity_footnote)
    dots$show_heterogeneity_footnote <- NULL
  }
  if ("annotation_digits" %in% names(dots)) {
    annotation_digits <- dots$annotation_digits
    dots$annotation_digits <- NULL
  }
  if ("heterogeneity_digits" %in% names(dots)) {
    annotation_digits <- dots$heterogeneity_digits
    dots$heterogeneity_digits <- NULL
  }
  group <- match.arg(group, c("subgroup", "trial"))
  validate_forestplot_args(
    prob = prob,
    which_panel = dots$which_panel,
    xlim = dots$xlim,
    xticks = dots$xticks,
    fill = if (!is.null(dots$fill)) dots$fill else "white",
    n_var = n_var
  )
  ifrac <- match.arg(ifrac)
  model_names <- comparison_forest_resolve_model_names(
    model_names = model_names,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction
  )

  effect_bundle <- comparison_forest_effect_bundle(
    fit = fit,
    ifrac = ifrac,
    prob = prob,
    max_support = max_ifrac_draws,
    contrast_sep = contrast_sep,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
    show_annotations = show_annotations,
    annotation_digits = annotation_digits,
    model_names = model_names
  )

  render_args <- comparison_forest_build_fit_panels(
    fit = fit,
    prob = prob,
    contrast_sep = contrast_sep,
    max_ifrac_draws = max_ifrac_draws,
    model_names = model_names,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
    n_var = n_var,
    group = group,
    ifrac = ifrac,
    effect_bundle = effect_bundle
  )

  annotations <- effect_bundle$annotations %||% list(
    lines = character(),
    ifrac = character(),
    subgroup_ifrac = character(),
    heterogeneity = character()
  )

  fp <- do.call(
    forestplot_render,
    c(
      render_args,
      list(
        prob = prob,
        contrast_sep = contrast_sep,
        exponentiate = exponentiate,
        annotations = annotations,
        group = group,
        type = type
      ),
      dots
    )
  )
  attr(fp, "footnote") <- if (
    length(annotations$lines) ||
      length(annotations$ifrac) ||
      length(annotations$heterogeneity)
  ) {
    paste(
      c(
        annotations$lines,
        unname(annotations$ifrac),
        unname(annotations$heterogeneity)
      ),
      collapse = "\n"
    )
  } else {
    NULL
  }
  attr(fp, "annotations") <- annotations
  attr(fp, "forestplot_fit_render") <- render_args
  attr(fp, "forestplot_fit_options") <- list(
    prob = prob,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate,
    ifrac = ifrac,
    group = group,
    show_annotations = isTRUE(show_annotations),
    dots = dots
  )
  fp
}

#' Shared forest quantities computed once (ifrac, shrinkage, heterogeneity)
#' @noRd
comparison_forest_effect_bundle <- function(fit,
                                      ifrac = c("interaction", "overall"),
                                      prob = 0.95,
                                      max_support = NULL,
                                      contrast_sep = "-",
                                      include_shrinkage = TRUE,
                                      include_overall = TRUE,
                                      include_prediction = TRUE,
                                      include_combined = TRUE,
                                      show_annotations = TRUE,
                                      annotation_digits = 2,
                                      model_names = NULL) {
  ifrac <- match.arg(ifrac)
  oi <- interactmeta_eval_ifrac(fit, ifrac = ifrac)
  need_shrinkage <- comparison_forest_layer_on(include_shrinkage)
  rot <- if (need_shrinkage) {
    comparison_shrinkage_rotation_bundle(
      fit,
      max_support = max_support,
      contrast_sep = contrast_sep,
      ifrac = ifrac,
      overall_ifrac = oi
    )
  } else {
    NULL
  }
  shrinkage_effects <- if (!is.null(rot)) {
    list(
      subgroups = comparison_shrinkage_effect_from_bundle(
        rot, type = "modifiers", prob = prob, summary = TRUE, exponentiate = FALSE
      ),
      interaction = comparison_shrinkage_effect_from_bundle(
        rot, type = "interaction", prob = prob, summary = TRUE, exponentiate = FALSE
      ),
      overall = comparison_shrinkage_effect_from_bundle(
        rot, type = "combined", prob = prob, summary = TRUE, exponentiate = FALSE
      )
    )
  } else {
    NULL
  }

  summary_blocks <- comparison_forest_summary_blocks(
    fit,
    model_names = model_names,
    prob = prob,
    ifrac = ifrac,
    max_support = max_support,
    contrast_sep = contrast_sep,
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
    overall_ifrac = oi
  )

  het <- if (isTRUE(show_annotations)) {
    tryCatch(
      heterogeneity(fit, component = "all", prob = prob, summary = TRUE),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  annotations <- if (isTRUE(show_annotations)) {
    forestplot_fit_annotations(
      fit = fit,
      ifrac = ifrac,
      prob = prob,
      digits = as.integer(annotation_digits)[1L],
      overall_ifrac = oi,
      heterogeneity_obj = het
    )
  } else {
    list(
      lines = character(),
      ifrac = character(),
      subgroup_ifrac = character(),
      heterogeneity = character()
    )
  }

  list(
    overall_ifrac = oi,
    ifrac = ifrac,
    rotation = rot,
    shrinkage = shrinkage_effects,
    summary_blocks = summary_blocks,
    heterogeneity = het,
    annotations = annotations
  )
}

#' Build forest-plot annotation lines (heterogeneity + overall IF)
#'
#' @return A list with:
#'   - `lines`: unused footer lines (kept for compatibility)
#'   - `ifrac`: named character vector of overall-IF notes by subgroup
#'     (drawn under the x-axis)
#'   - `subgroup_ifrac`: short IF labels (legacy; not drawn on headers)
#'   - `heterogeneity`: named `overall` / `interaction` notes (drawn under
#'     the x-axis)
#' @noRd
forestplot_fit_annotations <- function(fit,
                                       ifrac = c("interaction", "overall"),
                                       prob = 0.95,
                                       digits = 2,
                                       show_heterogeneity = TRUE,
                                       show_ifrac = FALSE,
                                       overall_ifrac = NULL,
                                       heterogeneity_obj = NULL) {
  ifrac <- match.arg(ifrac)
  lines <- character()
  ifrac_notes <- character()
  subgroup_ifrac <- character()
  het_notes <- character()
  plugin_fit <- isTRUE(interactmeta_ifrac_is_plugin(fit))

  if (isTRUE(show_ifrac)) {
    if_parts <- forestplot_overall_ifrac_annotations(
      fit = fit,
      ifrac = ifrac,
      prob = prob,
      digits = digits,
      overall_ifrac = overall_ifrac
    )
    ifrac_notes <- if_parts$ifrac
    subgroup_ifrac <- if_parts$subgroup_ifrac
  }

  if (isTRUE(show_heterogeneity)) {
    h <- heterogeneity_obj %||% tryCatch(
      heterogeneity(fit, component = "all", prob = prob, summary = TRUE),
      error = function(e) NULL
    )
    if (!is.null(h) && !is.null(h$summary) && nrow(h$summary)) {
      summ <- h$summary
      if ("type" %in% names(summ)) {
        ord <- order(match(as.character(summ$type), c("intercept", "slope"), nomatch = 99L))
        summ <- summ[ord, , drop = FALSE]
      }
      for (i in seq_len(nrow(summ))) {
        piece <- heterogeneity_footnote(
          summ[i, , drop = FALSE],
          prob = prob,
          digits = digits,
          include_ci = !isTRUE(plugin_fit)
        )
        if (!nzchar(piece)) {
          next
        }
        if (!interactmeta_is_bayesian_fit(fit) &&
              identical(toupper(as.character(fit$method)[[1L]]), "ML")) {
          piece <- paste0(piece, " (MLE)")
        }
        ty <- as.character(summ$type[[i]])
        key <- if (identical(ty, "slope")) "interaction" else "combined"
        het_notes[[key]] <- piece
      }
    }
  }

  list(
    lines = lines,
    ifrac = ifrac_notes,
    subgroup_ifrac = subgroup_ifrac,
    heterogeneity = het_notes
  )
}

#' Combined-IF estimate (and interval, when it is a posterior) per non-reference subgroup
#' @noRd
forestplot_overall_ifrac_annotations <- function(fit,
                                                 ifrac = c("interaction", "overall"),
                                                 prob = 0.95,
                                                 digits = 2,
                                                 overall_ifrac = NULL) {
  ifrac <- match.arg(ifrac)
  st <- attr(fit, "interactmeta", exact = TRUE)
  vars <- character()
  if (!is.null(st$ifrac_vars) && length(st$ifrac_vars)) {
    vars <- as.character(st$ifrac_vars)
  } else if (!is.null(st$ifrac_vars_std) && length(st$ifrac_vars_std)) {
    vars <- as.character(st$ifrac_vars_std)
  } else if (!is.null(st$ifrac_data)) {
    vars <- detect_ifrac_vars(st$ifrac_data)
  } else if (!is.null(st$ifrac_var) && nzchar(st$ifrac_var)) {
    vars <- st$ifrac_var
  }
  vars <- unique(vars[nzchar(vars)])

  ifrac_notes <- character()
  subgroup_ifrac <- character()
  fmt <- function(v) {
    format(round(as.numeric(v), digits), nsmall = digits, trim = TRUE)
  }
  if_prefix <- forestplot_ifrac_annotation_prefix(ifrac)
  plugin_fit <- isTRUE(interactmeta_ifrac_is_plugin(fit))

  for (iv in vars) {
    oi <- if (!is.null(overall_ifrac) && length(vars) == 1L) {
      overall_ifrac
    } else {
      tryCatch(
        interactmeta_eval_ifrac(
          fit, ifrac_var = iv, ifrac = ifrac, summary = FALSE, prob = prob
        ),
        error = function(e) NULL
      )
    }
    if (is.null(oi)) {
      next
    }
    summ <- attr(oi, "summary", exact = TRUE)
    if (is.null(summ) || !nrow(summ)) {
      next
    }
    sg <- attr(oi, "subgroup_label", exact = TRUE) %||%
      overall_ifrac_subgroup_label(iv, interactmeta = st)
    if (is.null(sg) || !nzchar(as.character(sg))) {
      sg <- interactmeta_ifrac_subgroup_label(iv) %||% iv
    }
    sg <- as.character(sg)[[1L]]
    est <- fmt(summ$estimate[[1L]])
    sg_tag <- paste0(" for ", sg)
    if (isTRUE(plugin_fit)) {
      note <- if (!interactmeta_is_bayesian_fit(fit)) {
        if (identical(ifrac, "interaction")) {
          " (at fitted \u03c4_\u03b3)"
        } else {
          " (at fitted \u03c4)"
        }
      } else {
        " (\u03c4 MLE plug-in)"
      }
      ifrac_notes[[sg]] <- paste0(
        if_prefix,
        sg_tag,
        ": ",
        est,
        note
      )
      subgroup_ifrac[[sg]] <- est
    } else {
      ifrac_notes[[sg]] <- paste0(
        if_prefix,
        sg_tag,
        ": ",
        est,
        " (",
        fmt(summ$lci[[1L]]),
        ", ",
        fmt(summ$uci[[1L]]),
        ")"
      )
      subgroup_ifrac[[sg]] <- paste0(
        est,
        " [",
        fmt(summ$lci[[1L]]),
        ", ",
        fmt(summ$uci[[1L]]),
        "]"
      )
    }
  }

  list(
    lines = unname(ifrac_notes),
    ifrac = ifrac_notes,
    subgroup_ifrac = subgroup_ifrac
  )
}

#' Axis-footnote label for the reporting information fraction
#' @noRd
forestplot_ifrac_annotation_prefix <- function(ifrac = c("interaction", "overall")) {
  match.arg(ifrac)
  "Information fraction (\u03c0)"
}

#' Reporting IF for the Combined `Overall` diamond IF column
#'
#' Percent estimate with a credible interval when `ifrac()` is posterior
#' (Bayesian); point estimate only for plug-in metafor IF.
#' @noRd
forestplot_reporting_if_cell <- function(fit,
                                         prob = 0.95,
                                         digits = 1L,
                                         ifrac = c("interaction", "overall"),
                                         overall_ifrac = NULL,
                                         ifrac_var = NULL,
                                         reference_subgroup = FALSE) {
  ifrac <- match.arg(ifrac)
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(ifrac_var) || !nzchar(ifrac_var)) {
    vars <- character()
    if (!is.null(st$ifrac_vars) && length(st$ifrac_vars)) {
      vars <- as.character(st$ifrac_vars)
    } else if (!is.null(st$ifrac_vars_std) && length(st$ifrac_vars_std)) {
      vars <- as.character(st$ifrac_vars_std)
    } else if (!is.null(st$ifrac_data)) {
      vars <- detect_ifrac_vars(st$ifrac_data)
    } else if (!is.null(st$ifrac_var) && nzchar(st$ifrac_var)) {
      vars <- st$ifrac_var
    }
    vars <- unique(vars[nzchar(vars)])
    ifrac_var <- vars[[1L]]
  }
  if (is.null(ifrac_var) || !nzchar(ifrac_var)) {
    return("")
  }
  oi <- overall_ifrac
  if (is.null(oi)) {
    oi <- tryCatch(
      interactmeta_eval_ifrac(
        fit,
        ifrac_var = ifrac_var,
        ifrac = ifrac,
        summary = FALSE,
        prob = prob
      ),
      error = function(e) NULL
    )
  }
  if (is.null(oi)) {
    return("")
  }
  summ <- attr(oi, "summary", exact = TRUE)
  if (is.null(summ) || !nrow(summ)) {
    draws <- if (is.numeric(oi)) oi else NULL
    if (is.null(draws) || !length(draws)) {
      return("")
    }
    w <- interactmeta_ifrac_draw_weights(fit, draws)
    summ <- overall_ifrac_summary_df(draws, prob = prob, weight = w)
  }
  if (is.null(summ) || !nrow(summ)) {
    return("")
  }
  if (isTRUE(reference_subgroup)) {
    draws <- if (is.numeric(oi)) oi else NULL
    if (!is.null(draws) && length(draws)) {
      w <- interactmeta_ifrac_draw_weights(fit, draws)
      summ <- overall_ifrac_summary_df(
        1 - as.numeric(draws),
        prob = prob,
        weight = w
      )
    } else if (is.finite(summ$estimate[[1L]])) {
      summ$estimate[[1L]] <- 1 - as.numeric(summ$estimate[[1L]])
      if (all(c("lci", "uci") %in% names(summ)) &&
            is.finite(summ$lci[[1L]]) &&
            is.finite(summ$uci[[1L]])) {
        lci <- 1 - as.numeric(summ$uci[[1L]])
        uci <- 1 - as.numeric(summ$lci[[1L]])
        summ$lci[[1L]] <- lci
        summ$uci[[1L]] <- uci
      }
    }
  }
  forestplot_format_reporting_if_summary(summ, fit = fit, digits = digits)
}

#' @noRd
forestplot_format_reporting_if_summary <- function(summ, fit, digits = 1L) {
  if (is.null(summ) || !nrow(summ)) {
    return("")
  }
  fmt_pct <- function(v) {
    forestplot_format_ifrac(as.numeric(v), digits = digits)[[1L]]
  }
  est <- fmt_pct(summ$estimate[[1L]])
  if (!is.finite(as.numeric(summ$estimate[[1L]]))) {
    return("")
  }
  plugin_fit <- isTRUE(interactmeta_ifrac_is_plugin(fit))
  if (isTRUE(plugin_fit) ||
        !all(c("lci", "uci") %in% names(summ)) ||
        !is.finite(summ$lci[[1L]]) ||
        !is.finite(summ$uci[[1L]])) {
    return(est)
  }
  paste0(
    est,
    " [",
    fmt_pct(summ$lci[[1L]]),
    ", ",
    fmt_pct(summ$uci[[1L]]),
    "]"
  )
}

#' Treat include-layer flags as on unless explicitly `FALSE`
#' @noRd
comparison_forest_layer_on <- function(x) {
  !isFALSE(x)
}

#' Build `model_names` from include flags when `model_names` is not set
#' @noRd
comparison_forest_resolve_model_names <- function(model_names,
                                            include_shrinkage,
                                            include_overall,
                                            include_prediction) {
  layers <- c("Shrinkage", "Overall", "Prediction")
  flags <- c(
    comparison_forest_layer_on(include_shrinkage),
    comparison_forest_layer_on(include_overall),
    comparison_forest_layer_on(include_prediction)
  )
  enabled_by_flag <- layers[flags]
  if (is.null(model_names)) {
    return(enabled_by_flag)
  }

  requested <- as.character(model_names)
  unknown <- setdiff(requested, layers)
  if (length(unknown)) {
    warning(
      "Unrecognized `model_names` dropped: ",
      paste(unknown, collapse = ", "),
      ". Known layers: ",
      paste(layers, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  out <- intersect(requested, layers)
  if (comparison_forest_layer_on(include_shrinkage)) {
    out <- unique(c("Shrinkage", out))
  }
  if (comparison_forest_layer_on(include_overall)) {
    out <- unique(c(out, "Overall"))
  }
  if (comparison_forest_layer_on(include_prediction)) {
    out <- unique(c(out, "Prediction"))
  }
  if (isFALSE(include_shrinkage)) {
    out <- setdiff(out, "Shrinkage")
  }
  if (isFALSE(include_overall)) {
    out <- setdiff(out, "Overall")
  }
  if (isFALSE(include_prediction)) {
    out <- setdiff(out, "Prediction")
  }
  ord <- match(out, layers)
  out[order(ord)]
}

#' Point scale for shrinkage rows relative to observed rows in CAMS forest plots
#' @noRd
comparison_forest_shrinkage_point_scale <- function() 0.15

#' Hex colour for forest study markers
#' @noRd
comparison_forest_marker_colour <- function() {
  interactmeta_plot_constants()$forest_study
}

#' Hex colour for forest shrinkage markers (`gray45`)
#' @noRd
comparison_forest_shrinkage_colour <- function() {
  interactmeta_plot_constants()$forest_shrinkage
}

#' Apply gray45 styling to CAMS forest shrinkage rows
#' @noRd
comparison_forest_apply_shrinkage_colours <- function(resultsA,
                                                resultsB,
                                                row_labels,
                                                subgroup_shrinkage_row_index,
                                                interaction_shrinkage_row_index,
                                                ref) {
  force(row_labels)
  force(ref)
  if (!length(subgroup_shrinkage_row_index) && !length(interaction_shrinkage_row_index)) {
    return(list(resultsA = resultsA, resultsB = resultsB))
  }

  n <- nrow(resultsA)
  if (!"fill" %in% names(resultsA)) {
    resultsA$fill <- rep(NA_character_, n)
    resultsA$colour <- rep(NA_character_, n)
  }
  if (!"fill" %in% names(resultsB)) {
    resultsB$fill <- rep(NA_character_, n)
    resultsB$colour <- rep(NA_character_, n)
  }

  shr_hex <- comparison_forest_shrinkage_colour()

  if (length(subgroup_shrinkage_row_index)) {
    resultsA$fill[subgroup_shrinkage_row_index] <- shr_hex
    resultsA$colour[subgroup_shrinkage_row_index] <- shr_hex
    resultsB$fill[subgroup_shrinkage_row_index] <- NA_character_
    resultsB$colour[subgroup_shrinkage_row_index] <- NA_character_
  }

  if (length(interaction_shrinkage_row_index)) {
    int_yi <- resultsB$yi[interaction_shrinkage_row_index]
    resultsB$fill[interaction_shrinkage_row_index] <- shr_hex
    resultsB$colour[interaction_shrinkage_row_index] <- shr_hex
    resultsA$fill[interaction_shrinkage_row_index] <- NA_character_
    resultsA$colour[interaction_shrinkage_row_index] <- NA_character_
    na_b <- is.na(int_yi)
    resultsB$fill[interaction_shrinkage_row_index[na_b]] <- NA_character_
    resultsB$colour[interaction_shrinkage_row_index[na_b]] <- NA_character_
  }

  list(resultsA = resultsA, resultsB = resultsB)
}

#' @noRd
comparison_forest_diamond_shape <- function() {
  18L
}

#' Placeholder row label for forest panels (must not be `""`)
#'
#' Empty strings collapse distinct forest rows. Uses a zero-width space so the
#' trial-name column stays blank without showing as a visible glyph or
#' `<U+00A0>` escape in graphics devices.
#' @noRd
comparison_forest_blank_row_label <- function() {
  "\u200b"
}

#' @noRd
comparison_forest_is_blank_row_label <- function(x) {
  is.na(x) | !nzchar(x) |
    x %in% c(comparison_forest_blank_row_label(), "\u00a0", "\u00A0", "\u200b")
}

#' Build model-based forest panel data frames from a CAMS fit
#'
#' Overall / Prediction / Combined diamonds are evaluated at the overall IF
#' \eqn{\omega}. Shrinkage greys are trial-specific CAMS reconstructions
#' from [shrinkage()] at each trial's \eqn{\pi_j} (not the FE linpred).
#' @noRd
comparison_forest_build_fit_panels <- function(fit,
                                         prob,
                                         contrast_sep,
                                         max_ifrac_draws,
                                         model_names,
                                         include_shrinkage = TRUE,
                                         include_overall = TRUE,
                                         include_prediction = TRUE,
                                         include_combined = TRUE,
                                         n_var = FALSE,
                                         group = "subgroup",
                                         ifrac = c("interaction", "overall"),
                                         effect_bundle = NULL) {
  group <- match.arg(group, c("subgroup", "trial"))
  ifrac <- match.arg(ifrac)
  st <- attr(fit, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()

  model_names <- comparison_forest_resolve_model_names(
    model_names = model_names,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction
  )
  data <- validate_forestplot_data(st$ifrac_data)
  panels <- split_ifrac_for_forest(
    data,
    contrast_sep = contrast_sep,
    group = group
  )
  ref <- if (interactmeta_is_non_cams_analysis(fit) && !is.null(fit$reference)) {
    as.character(fit$reference)[[1L]]
  } else {
    interactmeta_reference(fit)
  }
  # Keep data$reference aligned with the fit attribute
  if ("reference" %in% names(data)) {
    interactmeta_assert_reference_agrees(
      interactmeta_reference_from_column(data$reference),
      ref,
      where = "forest ifrac_data"
    )
  }
  aligned <- align_interaction_to_subgroups(
    panels$subgroup,
    panels$interaction,
    reference = ref
  )

  ifrac_disp <- forestplot_row_information_fraction(
    aligned$subgroup,
    contrast_sep = contrast_sep,
    digits = 2
  )

  resultsA <- forestplot_panel_df(
    aligned$subgroup,
    prob = prob,
    ifrac_display = ifrac_disp
  )
  resultsA <- forestplot_attach_trial_size(
    resultsA,
    aligned$subgroup,
    trial_col = "trial",
    subgroup_col = "subgroup",
    n_var = n_var,
    contrast_sep = contrast_sep
  )
  resultsB <- forestplot_panel_df(aligned$interaction, prob = prob)
  n_obs <- nrow(resultsA)

  row_labels <- forestplot_study_row_labels(
    trial = aligned$subgroup$trial,
    subgroup = aligned$subgroup$subgroup,
    group = group
  )

  resultsA$fill <- "white"
  resultsA$colour <- comparison_forest_marker_colour()
  resultsA$bold <- FALSE
  resultsB$fill <- "white"
  resultsB$colour <- comparison_forest_marker_colour()
  resultsB$bold <- FALSE
  resultsA$row_kind <- "study"
  resultsB$row_kind <- "study"

  ix_obs <- seq_len(n_obs)
  ix_sg_shr <- integer(0)
  ix_int_shr <- integer(0)
  ix_shr <- integer(0)
  overall_shrinkage <- NULL

  if (comparison_forest_layer_on(include_shrinkage)) {
    shrinkage <- comparison_forest_shrinkage_data(
      fit,
      prob = prob,
      contrast_sep = contrast_sep,
      max_support = max_ifrac_draws,
      ifrac = ifrac,
      group = group,
      shrinkage_effects = effect_bundle$shrinkage
    )
    overall_shrinkage <- shrinkage$overall
    inter <- interleave_comparison_shrinkage_rows(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      shrinkage = shrinkage,
      prob = prob,
      ref = ref,
      contrast_sep = contrast_sep,
      ifrac_display = ifrac_disp
    )
    resultsA <- inter$resultsA
    resultsB <- inter$resultsB
    row_labels <- inter$row_labels
    ix_obs <- inter$obs_row_index
    ix_sg_shr <- inter$subgroup_shrinkage_row_index
    ix_int_shr <- inter$interaction_shrinkage_row_index
    ix_shr <- c(ix_sg_shr, ix_int_shr)
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
    resultsA$row_kind[ix_shr] <- "shrinkage"
    resultsB$row_kind[ix_shr] <- "shrinkage"
  }

  resultsA$shape <- NA_integer_
  resultsB$shape <- NA_integer_
  resultsA$shape[ix_obs] <- 22L
  resultsB$shape[ix_obs] <- 21L
  if (length(ix_sg_shr)) {
    resultsA$shape[ix_sg_shr] <- 15L
    resultsB$shape[ix_sg_shr] <- NA_integer_
  }
  if (length(ix_int_shr)) {
    resultsA$shape[ix_int_shr] <- NA_integer_
    resultsB$shape[ix_int_shr] <- 16L
  }
  resultsA$point_scale <- 1
  resultsB$point_scale <- 1
  if (length(ix_shr)) {
    shr_scale <- comparison_forest_shrinkage_point_scale()
    resultsA$point_scale[ix_shr] <- shr_scale
    resultsB$point_scale[ix_shr] <- shr_scale
  }

  marginal_names <- setdiff(model_names, "Shrinkage")
  blocks <- list()
  need_marginal <- length(marginal_names) &&
    (comparison_forest_layer_on(include_overall) ||
      comparison_forest_layer_on(include_prediction) ||
      comparison_forest_layer_on(include_combined))
  if (need_marginal) {
    blocks <- effect_bundle$summary_blocks %||% comparison_forest_summary_blocks(
      fit,
      model_names = marginal_names,
      prob = prob,
      ifrac = ifrac,
      max_support = max_ifrac_draws,
      contrast_sep = contrast_sep,
      include_overall = include_overall,
      include_prediction = include_prediction,
      include_combined = include_combined,
      overall_ifrac = effect_bundle$overall_ifrac
    )
  }
  if (length(blocks)) {
    if (!"colour_label" %in% names(row_labels)) {
      row_labels$colour_label <- row_labels$label
    }
    extra <- append_comparison_forest_blocks(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      blocks = blocks,
      prob = prob,
      ref = ref,
      contrast_sep = contrast_sep
    )
    resultsA <- extra$resultsA
    resultsB <- extra$resultsB
    row_labels <- extra$row_labels
  }

  styled <- comparison_forest_apply_effect_styles(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    ref = ref,
    contrast_sep = contrast_sep
  )
  resultsA <- styled$resultsA
  resultsB <- styled$resultsB
  if (!is.null(overall_shrinkage)) {
    attr(resultsA, "overall_shrinkage") <- overall_shrinkage
  }

  list(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    n_study = n_obs,
    ref = ref,
    obs_row_index = ix_obs,
    include_shrinkage = comparison_forest_layer_on(include_shrinkage),
    include_prediction = comparison_forest_layer_on(include_prediction),
    include_combined = comparison_forest_layer_on(include_combined),
    overall_shrinkage = overall_shrinkage
  )
}

#' Diamond styling for Overall / Prediction / Combined
#'
#' Summary diamonds use blue (`forestplot_default_summary_colour()`) for
#' every block — the same whether `group` is `"subgroup"` or `"trial"`.
#' Shrinkage rows are styled separately via
#' `comparison_forest_apply_shrinkage_colours()` (`gray45`). `marker` is only a
#' fallback when a heading is empty.
#' @noRd
comparison_forest_apply_effect_styles <- function(resultsA,
                                            resultsB,
                                            row_labels,
                                            ref,
                                            contrast_sep = "-",
                                            marker = NULL,
                                            bold_summary = FALSE) {
  force(row_labels)
  force(ref)
  force(contrast_sep)
  if (!"row_kind" %in% names(resultsA)) {
    return(list(resultsA = resultsA, resultsB = resultsB))
  }

  diamond <- comparison_forest_diamond_shape()
  fallback <- if (is.null(marker) || !nzchar(as.character(marker)[[1L]])) {
    forestplot_default_summary_colour()
  } else {
    as.character(marker)[[1L]]
  }
  overall_ix <- which(resultsA$row_kind == "overall")
  pred_ix <- which(resultsA$row_kind == "prediction")
  combined_ix <- which(resultsA$row_kind == "combined")
  if (!length(overall_ix) && !length(pred_ix) && !length(combined_ix)) {
    return(list(resultsA = resultsA, resultsB = resultsB))
  }

  headings <- if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    rep("", nrow(resultsA))
  }
  if (length(headings) < nrow(resultsA)) {
    headings <- rep_len(headings, nrow(resultsA))
  }
  # Before stacking, marginal rows often have blank heading and the subgroup
  # in colour_label — use that so colours are correct either side of stack.
  empty_h <- !nzchar(headings) | is.na(headings)
  if (any(empty_h) && "colour_label" %in% names(row_labels)) {
    cl <- as.character(row_labels$colour_label)
    if (length(cl) < nrow(resultsA)) {
      cl <- rep_len(cl, nrow(resultsA))
    }
    use_cl <- empty_h & nzchar(cl) & !is.na(cl)
    headings[use_cl] <- cl[use_cl]
  }
  sum_cols <- forestplot_summary_diamond_colours(
    headings = headings,
    reference = ref
  )
  empty_h <- !nzchar(headings) | is.na(headings)
  if (any(empty_h)) {
    sum_cols[empty_h] <- fallback
  }

  paint_summary_ix <- function(ix) {
    if (!length(ix)) {
      return()
    }
    cols <- sum_cols[ix]
    resultsA$shape[ix] <<- diamond
    resultsA$fill[ix] <<- cols
    resultsA$colour[ix] <<- cols
    resultsA$bold[ix] <<- isTRUE(bold_summary)
    resultsB$bold[ix] <<- isTRUE(bold_summary)
  }

  if (length(overall_ix)) {
    paint_summary_ix(overall_ix)
    resultsB$shape[overall_ix] <- NA_integer_
    resultsB$fill[overall_ix] <- NA_character_
    resultsB$colour[overall_ix] <- NA_character_
    ok_b <- overall_ix[is.finite(resultsB$yi[overall_ix])]
    if (length(ok_b)) {
      resultsB$shape[ok_b] <- diamond
      resultsB$fill[ok_b] <- sum_cols[ok_b]
      resultsB$colour[ok_b] <- sum_cols[ok_b]
    }
  }

  if (length(pred_ix)) {
    has_a <- is.finite(resultsA$yi[pred_ix])
    has_b <- is.finite(resultsB$yi[pred_ix])
    paint_summary_ix(pred_ix)
    miss_a <- pred_ix[!has_a]
    resultsA$fill[miss_a] <- NA_character_
    resultsA$colour[miss_a] <- NA_character_
    resultsA$shape[miss_a] <- NA_integer_
    resultsB$shape[pred_ix] <- NA_integer_
    resultsB$fill[pred_ix] <- NA_character_
    resultsB$colour[pred_ix] <- NA_character_
    ok_b <- pred_ix[has_b]
    if (length(ok_b)) {
      resultsB$shape[ok_b] <- diamond
      resultsB$fill[ok_b] <- sum_cols[ok_b]
      resultsB$colour[ok_b] <- sum_cols[ok_b]
    }
  }

  if (length(combined_ix)) {
    overall_col <- forestplot_overall_summary_colour()
    resultsA$shape[combined_ix] <- diamond
    resultsB$shape[combined_ix] <- NA_integer_
    resultsA$fill[combined_ix] <- overall_col
    resultsA$colour[combined_ix] <- overall_col
    resultsB$fill[combined_ix] <- NA_character_
    resultsB$colour[combined_ix] <- NA_character_
    resultsA$bold[combined_ix] <- isTRUE(bold_summary)
    resultsB$bold[combined_ix] <- isTRUE(bold_summary)
  }

  list(resultsA = resultsA, resultsB = resultsB)
}

#' Shrinkage rows from [shrinkage()]
#' @noRd
comparison_forest_shrinkage_data <- function(fit,
                                       prob,
                                       contrast_sep = "-",
                                       max_support = NULL,
                                       ifrac = c("interaction", "overall"),
                                       group = "subgroup",
                                       shrinkage_effects = NULL) {
  group <- match.arg(group, c("subgroup", "trial"))
  st <- attr(fit, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var
  subgroup_var <- st$subgroup_var
  ref <- st$reference

  sg_dat <- comparison_forest_subgroup_rows(st$ifrac_data, trial_var, subgroup_var, contrast_sep)
  if (is.null(shrinkage_effects)) {
    rot <- comparison_shrinkage_rotation_bundle(
      fit,
      max_support = max_support,
      contrast_sep = contrast_sep,
      ifrac = match.arg(ifrac)
    )
    shrinkage_effects <- list(
      subgroups = comparison_shrinkage_effect_from_bundle(
        rot, type = "modifiers", prob = prob, summary = TRUE, exponentiate = FALSE
      ),
      interaction = comparison_shrinkage_effect_from_bundle(
        rot, type = "interaction", prob = prob, summary = TRUE, exponentiate = FALSE
      ),
      overall = comparison_shrinkage_effect_from_bundle(
        rot, type = "combined", prob = prob, summary = TRUE, exponentiate = FALSE
      )
    )
  }
  me_sg <- shrinkage_effects$subgroups
  me_int <- shrinkage_effects$interaction
  me_ov <- shrinkage_effects$overall
  sg <- interactmeta_effects_merge_yivi(
    template = sg_dat,
    summary = me_sg$summary,
    trial_var = trial_var,
    label_var = subgroup_var,
    type = "Shrinkage",
    use_posterior_qi = TRUE
  )
  sg <- comparison_forest_order_to_template(
    sg,
    split_ifrac_for_forest(
      validate_forestplot_data(st$ifrac_data),
      contrast_sep = contrast_sep,
      group = group
    )$subgroup,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )

  int_template <- build_contrast_rows_from_subgroups(
    sg,
    contrast_sep = contrast_sep,
    reference = ref
  )
  int <- interactmeta_effects_merge_yivi(
    template = int_template,
    summary = me_int$summary,
    trial_var = trial_var,
    label_var = subgroup_var,
    summary_label_var = "interaction",
    type = "Shrinkage",
    use_posterior_qi = TRUE
  )

  ov <- comparison_forest_overall_shrinkage_yivi(
    summary = me_ov$summary,
    trial_var = trial_var
  )

  list(subgroup = sg, interaction = int, overall = ov)
}

#' Trial-level overall shrinkage `yi` / `vi` for the Combined forest block
#'
#' Uses the median and posterior QI of reconstructed overall (`m_j`) draws.
#' @noRd
comparison_forest_overall_shrinkage_yivi <- function(summary, trial_var = "study") {
  if (is.null(summary) || !nrow(summary)) {
    return(data.frame(
      trial = character(0),
      yi = numeric(0),
      vi = numeric(0),
      ci.lb = numeric(0),
      ci.ub = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  s <- summary
  if ("type" %in% names(s)) {
    s <- s[as.character(s$type) == "Shrinkage" | is.na(s$type), , drop = FALSE]
  }
  if (!nrow(s) || !("trial" %in% names(s) || trial_var %in% names(s))) {
    return(data.frame(
      trial = character(0),
      yi = numeric(0),
      vi = numeric(0),
      ci.lb = numeric(0),
      ci.ub = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  tr_col <- if ("trial" %in% names(s)) "trial" else trial_var
  yi <- if ("estimate" %in% names(s)) {
    as.numeric(s$estimate)
  } else if ("mean" %in% names(s)) {
    as.numeric(s$mean)
  } else {
    as.numeric(s$yi)
  }
  vi <- if ("sd" %in% names(s)) {
    pmax(as.numeric(s$sd)^2, 1e-12, na.rm = TRUE)
  } else {
    rep(1e-12, nrow(s))
  }
  ci.lb <- if ("lci" %in% names(s)) as.numeric(s$lci) else rep(NA_real_, nrow(s))
  ci.ub <- if ("uci" %in% names(s)) as.numeric(s$uci) else rep(NA_real_, nrow(s))
  data.frame(
    trial = as.character(s[[tr_col]]),
    yi = yi,
    vi = vi,
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    stringsAsFactors = FALSE
  )
}

#' Build one Overall / Prediction block (Overall rows, then Prediction rows)
#' @noRd
comparison_forest_summary_blocks <- function(fit,
                                       model_names,
                                       prob,
                                       ifrac = c("interaction", "overall"),
                                       max_support = NULL,
                                       contrast_sep = "-",
                                       include_overall = TRUE,
                                       include_prediction = TRUE,
                                       include_combined = TRUE,
                                       overall_ifrac = NULL) {
  ifrac <- match.arg(ifrac)
  st <- attr(fit, "interactmeta", exact = TRUE)
  ref <- st$reference
  model_names <- comparison_forest_resolve_model_names(
    model_names = setdiff(model_names, "Shrinkage"),
    include_shrinkage = FALSE,
    include_overall = include_overall,
    include_prediction = include_prediction
  )
  show_overall <- comparison_forest_layer_on(include_overall) && "Overall" %in% model_names
  show_prediction <- comparison_forest_layer_on(include_prediction) && "Prediction" %in% model_names
  if (!length(model_names) || (!show_overall && !show_prediction)) {
    return(list())
  }

  oi <- overall_ifrac %||% interactmeta_eval_ifrac(fit, ifrac = match.arg(ifrac))
  overall_label <- sprintf("%.3f", stats::median(oi))
  subgroup_order <- comparison_forest_subgroup_order(st$ifrac_data, contrast_sep)

  empty_sg <- data.frame(
    trial = character(0),
    subgroup = character(0),
    yi = numeric(0),
    vi = numeric(0),
    stringsAsFactors = FALSE
  )
  overall_sg <- empty_sg
  overall_int <- empty_sg
  pred_sg <- empty_sg
  pred_int <- empty_sg

  if (show_overall) {
    overall_sg <- comparison_forest_marginal_subgroup_block(
      fit, oi, "Overall", overall_label, prob, max_support, contrast_sep,
      subgroup_order = subgroup_order,
      ifrac = ifrac
    )
    overall_int <- comparison_forest_marginal_interaction_block(
      fit, oi, "Overall", overall_label, prob, max_support, contrast_sep,
      ifrac = ifrac
    )
  }
  if (show_prediction) {
    pred_sg <- comparison_forest_marginal_subgroup_block(
      fit, oi, "Prediction", overall_label, prob, max_support, contrast_sep,
      subgroup_order = subgroup_order,
      ifrac = ifrac
    )
    pred_int <- comparison_forest_marginal_interaction_block(
      fit, oi, "Prediction", overall_label, prob, max_support, contrast_sep,
      ifrac = ifrac
    )
  }

  combined <- comparison_forest_interleave_overall_prediction(
    overall_sg = overall_sg,
    pred_sg = pred_sg,
    overall_int = overall_int,
    pred_int = pred_int,
    reference = ref,
    include_overall = show_overall,
    include_prediction = show_prediction
  )
  ref_chr <- as.character(ref)[[1L]]
  reporting_if_for_row <- function(colour_lab) {
    forestplot_reporting_if_cell(
      fit = fit,
      prob = prob,
      digits = 1L,
      ifrac = ifrac,
      overall_ifrac = oi,
      reference_subgroup = identical(as.character(colour_lab)[[1L]], ref_chr)
    )
  }
  is_overall <- combined$subgroup$row_kind == "overall"
  combined$subgroup$ifrac <- ifelse(
    is_overall,
    vapply(
      seq_len(nrow(combined$subgroup)),
      function(i) {
        if (!is_overall[[i]]) {
          return("")
        }
        reporting_if_for_row(combined$subgroup$colour_label[[i]])
      },
      character(1)
    ),
    ""
  )

  if (comparison_forest_layer_on(include_combined)) {
    # Always reserve Overall + Prediction Combined slots when either summary
    # layer is shown, so switching one off does not shrink the figure.
    overall_types <- character(0)
    if (show_overall || show_prediction) {
      overall_types <- c("Overall", "Prediction")
    }
    if (length(overall_types)) {
      overall_extra <- comparison_forest_overall_rows(
        fit = fit,
        overall_ifrac = oi,
        types = overall_types,
        overall_label = overall_label,
        prob = prob,
        max_support = max_support,
        ifrac = ifrac,
        include_types = c(
          if (show_overall) "Overall",
          if (show_prediction) "Prediction"
        )
      )
      # Combined Overall / Prediction: no IF column label (use footnotes for ω).
      overall_extra$subgroup$ifrac <- ""
      if (nrow(overall_extra$subgroup)) {
        combined$subgroup <- dplyr::bind_rows(
          combined$subgroup,
          overall_extra$subgroup
        )
        combined$panel_b <- dplyr::bind_rows(
          combined$panel_b,
          overall_extra$panel_b
        )
      }
    }
  }

  if (show_prediction &&
        !any(
          combined$subgroup$row_kind == "prediction" &
            is.finite(combined$subgroup$yi),
          na.rm = TRUE
        )) {
    warning(
      "No prediction rows were added to the forest plot; check that the fit supports new-trial prediction (trial random effects).",
      call. = FALSE
    )
  }

  blocks <- list(list(
    name = "",
    subgroup = combined$subgroup,
    panel_b = combined$panel_b,
    interaction = NULL
  ))
  attr(blocks, "overall_ifrac") <- oi
  blocks
}

#' @noRd
comparison_forest_marginal_subgroup_block <- function(fit,
                                                overall_ifrac,
                                                type,
                                                overall_label,
                                                prob,
                                                max_support,
                                                contrast_sep,
                                                subgroup_order = NULL,
                                                ifrac = c("interaction", "overall")) {
  sg_me <- modifiers(
    fit,
    at = overall_ifrac,
    ifrac = ifrac,
    type = type,
    prob = prob,
    max_support = max_support,
    summary = TRUE,
    contrast_sep = contrast_sep,
    exponentiate = FALSE
  )
  comparison_forest_overall_subgroup_block(
    sg_me$summary,
    type = type,
    overall_label = overall_label,
    subgroup_order = subgroup_order
  )
}

#' @noRd
comparison_forest_marginal_interaction_block <- function(fit,
                                                 overall_ifrac,
                                                 type,
                                                 overall_label,
                                                 prob,
                                                 max_support,
                                                 contrast_sep,
                                                 ifrac = c("interaction", "overall")) {
  int_me <- interactions(
    fit,
    at = overall_ifrac,
    ifrac = ifrac,
    type = type,
    prob = prob,
    max_support = max_support,
    summary = TRUE,
    contrast_sep = contrast_sep,
    exponentiate = FALSE
  )
  comparison_forest_overall_interaction_block(
    int_me$summary,
    type = type,
    overall_label = overall_label
  )
}

#' Combined marginal rows from [combined()] (left panel only)
#'
#' `types` is the full slot order (usually `Overall` then `Prediction`).
#' `include_types` is which of those get real estimates; the rest are `NA`
#' placeholders that keep vertical space when a layer is switched off.
#' @noRd
comparison_forest_overall_rows <- function(fit,
                                            overall_ifrac,
                                            types,
                                            overall_label,
                                            prob,
                                            max_support,
                                            ifrac = c("interaction", "overall"),
                                            include_types = types) {
  empty <- data.frame(
    trial = character(0),
    subgroup = character(0),
    yi = numeric(0),
    vi = numeric(0),
    stringsAsFactors = FALSE
  )
  if (!length(types)) {
    return(list(subgroup = empty, panel_b = empty))
  }
  include_types <- as.character(include_types %||% character())

  me <- if (length(include_types)) {
    combined(
      fit,
      at = overall_ifrac,
      ifrac = ifrac,
      type = include_types,
      prob = prob,
      max_support = max_support,
      summary = TRUE,
      exponentiate = FALSE
    )
  } else {
    NULL
  }
  s <- if (!is.null(me)) me$summary else NULL
  rows <- lapply(types, function(ty) {
    row <- if (!is.null(s) && ty %in% include_types) {
      s[s$type == ty, , drop = FALSE]
    } else {
      s[0L, , drop = FALSE]
    }
    if (!nrow(row) || !ty %in% include_types) {
      return(data.frame(
        trial = overall_label,
        subgroup = ty,
        colour_label = forestplot_combined_block_label(),
        row_kind = "combined",
        yi = NA_real_,
        vi = NA_real_,
        ci.lb = NA_real_,
        ci.ub = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      trial = overall_label,
      subgroup = ty,
      colour_label = forestplot_combined_block_label(),
      row_kind = "combined",
      yi = comparison_forest_summary_point(row)[[1L]],
      vi = pmax(row$sd[[1L]]^2, 1e-12, na.rm = TRUE),
      ci.lb = comparison_forest_interval_col(row, "lci")[[1L]],
      ci.ub = comparison_forest_interval_col(row, "uci")[[1L]],
      stringsAsFactors = FALSE
    )
  })
  sg <- dplyr::bind_rows(rows)
  if (!nrow(sg)) {
    return(list(subgroup = empty, panel_b = empty))
  }
  panel_b <- data.frame(
    yi = rep(NA_real_, nrow(sg)),
    vi = rep(NA_real_, nrow(sg)),
    ci.lb = rep(NA_real_, nrow(sg)),
    ci.ub = rep(NA_real_, nrow(sg)),
    stringsAsFactors = FALSE
  )
  list(subgroup = sg, panel_b = panel_b)
}

#' Map one overall interaction contrast onto non-reference subgroup rows
#' @noRd
comparison_forest_align_marginal_interaction <- function(subgroup_df,
                                                   interaction_df,
                                                   reference) {
  if (!nrow(subgroup_df)) {
    return(data.frame(
      subgroup = character(0),
      yi = numeric(0),
      vi = numeric(0),
      ci.lb = numeric(0),
      ci.ub = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  out <- data.frame(
    subgroup = subgroup_df$subgroup,
    yi = NA_real_,
    vi = NA_real_,
    ci.lb = NA_real_,
    ci.ub = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!nrow(interaction_df)) {
    return(out)
  }
  int_yi <- interaction_df$yi[[1L]]
  int_vi <- interaction_df$vi[[1L]]
  int_lci <- comparison_forest_interval_col(interaction_df, "ci.lb")[[1L]]
  int_uci <- comparison_forest_interval_col(interaction_df, "ci.ub")[[1L]]
  ref_chr <- as.character(reference)
  non_ref <- nzchar(as.character(out$subgroup)) &
    as.character(out$subgroup) != ref_chr
  out$yi[non_ref] <- int_yi
  out$vi[non_ref] <- int_vi
  out$ci.lb[non_ref] <- int_lci
  out$ci.ub[non_ref] <- int_uci
  out
}

#' Read an optional precomputed interval column
#' @noRd
comparison_forest_interval_col <- function(x, name) {
  if (name %in% names(x)) {
    return(as.numeric(x[[name]]))
  }
  rep(NA_real_, nrow(x))
}

#' Interleave marginal rows by subgroup (pooled then prediction per level)
#'
#' When one of `include_overall` / `include_prediction` is `FALSE` but the other
#' is `TRUE`, the omitted layer still gets an all-`NA` placeholder row so the
#' forest keeps the same vertical rhythm as the full plot.
#' @noRd
comparison_forest_interleave_overall_prediction <- function(overall_sg,
                                                     pred_sg,
                                                     overall_int,
                                                     pred_int,
                                                     reference,
                                                     include_overall = TRUE,
                                                     include_prediction = TRUE) {
  want_overall <- isTRUE(include_overall)
  want_pred <- isTRUE(include_prediction)
  # Reserve both slots whenever either layer is shown.
  reserve_pair <- want_overall || want_pred
  template <- if (want_overall && nrow(overall_sg) > 0L) {
    overall_sg
  } else if (want_pred && nrow(pred_sg) > 0L) {
    pred_sg
  } else {
    overall_sg
  }
  if (!nrow(template) || !reserve_pair) {
    return(list(
      subgroup = template[0L, , drop = FALSE],
      panel_b = data.frame(yi = numeric(0), vi = numeric(0), stringsAsFactors = FALSE)
    ))
  }
  if (want_overall && want_pred &&
      nrow(overall_sg) > 0L && nrow(pred_sg) > 0L &&
      nrow(pred_sg) != nrow(overall_sg)) {
    stop("Overall and Prediction subgroup rows must match.", call. = FALSE)
  }

  overall_b <- if (want_overall && nrow(overall_sg) > 0L) {
    comparison_forest_align_marginal_interaction(
      overall_sg,
      overall_int,
      reference = reference
    )
  } else {
    NULL
  }
  pred_b <- if (want_pred && nrow(pred_sg) > 0L) {
    comparison_forest_align_marginal_interaction(
      pred_sg,
      pred_int,
      reference = reference
    )
  } else {
    NULL
  }
  rows <- vector("list", 0L)
  n <- nrow(template)

  blank_slot <- function(tr, lab, colour_lab, row_kind) {
    data.frame(
      trial = tr,
      subgroup = lab,
      colour_label = colour_lab,
      row_kind = row_kind,
      yi = NA_real_,
      vi = NA_real_,
      ci.lb = NA_real_,
      ci.ub = NA_real_,
      b_yi = NA_real_,
      b_vi = NA_real_,
      b_ci.lb = NA_real_,
      b_ci.ub = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  for (i in seq_len(n)) {
    lab <- as.character(template$subgroup[[i]])
    tr <- template$trial[[i]]
    colour_lab <- if ("colour_label" %in% names(template)) {
      as.character(template$colour_label[[i]])
    } else {
      lab
    }
    if (is.na(colour_lab) || !nzchar(colour_lab)) {
      colour_lab <- lab
    }

    if (want_overall && nrow(overall_sg) > 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        trial = tr,
        subgroup = lab,
        colour_label = colour_lab,
        row_kind = "overall",
        yi = overall_sg$yi[[i]],
        vi = overall_sg$vi[[i]],
        ci.lb = comparison_forest_interval_col(overall_sg, "ci.lb")[[i]],
        ci.ub = comparison_forest_interval_col(overall_sg, "ci.ub")[[i]],
        b_yi = overall_b$yi[[i]],
        b_vi = overall_b$vi[[i]],
        b_ci.lb = comparison_forest_interval_col(overall_b, "ci.lb")[[i]],
        b_ci.ub = comparison_forest_interval_col(overall_b, "ci.ub")[[i]],
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1L]] <- blank_slot(
        tr = tr,
        lab = lab,
        colour_lab = colour_lab,
        row_kind = "overall"
      )
    }

    if (want_pred && nrow(pred_sg) > 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        trial = tr,
        subgroup = "",
        colour_label = colour_lab,
        row_kind = "prediction",
        yi = pred_sg$yi[[i]],
        vi = pred_sg$vi[[i]],
        ci.lb = comparison_forest_interval_col(pred_sg, "ci.lb")[[i]],
        ci.ub = comparison_forest_interval_col(pred_sg, "ci.ub")[[i]],
        b_yi = pred_b$yi[[i]],
        b_vi = pred_b$vi[[i]],
        b_ci.lb = comparison_forest_interval_col(pred_b, "ci.lb")[[i]],
        b_ci.ub = comparison_forest_interval_col(pred_b, "ci.ub")[[i]],
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1L]] <- blank_slot(
        tr = tr,
        lab = "",
        colour_lab = colour_lab,
        row_kind = "prediction"
      )
    }
  }

  out <- dplyr::bind_rows(rows)
  panel_b <- data.frame(
    yi = out$b_yi,
    vi = out$b_vi,
    ci.lb = out$b_ci.lb,
    ci.ub = out$b_ci.ub,
    stringsAsFactors = FALSE
  )
  out[c("b_yi", "b_vi", "b_ci.lb", "b_ci.ub")] <- NULL
  list(subgroup = out, panel_b = panel_b)
}

#' Point estimate column on a marginal summary table
#' @noRd
comparison_forest_summary_point <- function(s) {
  if ("estimate" %in% names(s)) {
    return(as.numeric(s$estimate))
  }
  as.numeric(s$mean)
}

#' Level labels on a [modifiers()] summary table
#' @noRd
comparison_forest_summary_level <- function(s) {
  if ("level" %in% names(s)) {
    return(as.character(s$level))
  }
  if ("subgroup" %in% names(s)) {
    return(as.character(s$subgroup))
  }
  stop("Marginal summary is missing a level / subgroup column.", call. = FALSE)
}

#' Interaction labels on an [interactions()] summary table
#' @noRd
comparison_forest_summary_interaction <- function(s) {
  if ("interaction" %in% names(s)) {
    return(as.character(s$interaction))
  }
  if ("subgroup" %in% names(s)) {
    return(as.character(s$subgroup))
  }
  stop("Interaction summary is missing an interaction column.", call. = FALSE)
}

#' @noRd
comparison_forest_overall_subgroup_block <- function(summary,
                                               type,
                                               overall_label,
                                               subgroup_order = NULL) {
  s <- summary[summary$type == type, , drop = FALSE]
  out <- data.frame(
    trial = overall_label,
    subgroup = comparison_forest_summary_level(s),
    yi = comparison_forest_summary_point(s),
    vi = pmax(s$sd^2, 1e-12, na.rm = TRUE),
    ci.lb = comparison_forest_interval_col(s, "lci"),
    ci.ub = comparison_forest_interval_col(s, "uci"),
    row_kind = "subgroup_effect",
    stringsAsFactors = FALSE
  )
  if (is.null(subgroup_order)) {
    return(order_forest_rows(out))
  }
  ord <- match(as.character(subgroup_order), as.character(out$subgroup))
  if (any(is.na(ord))) {
    stop(
      "Overall subgroup rows do not match forest subgroup levels.",
      call. = FALSE
    )
  }
  out[ord, , drop = FALSE]
}

#' Subgroup order used within each trial in CAMS forest plots
#' @noRd
comparison_forest_subgroup_order <- function(ifrac_data, contrast_sep) {
  sg <- split_ifrac_for_forest(
    validate_forestplot_data(ifrac_data),
    contrast_sep = contrast_sep
  )$subgroup
  if (is.factor(sg$subgroup)) {
    return(levels(sg$subgroup))
  }
  unique(as.character(sg$subgroup))
}

#' @noRd
comparison_forest_overall_interaction_block <- function(summary, type, overall_label) {
  s <- summary[summary$type == type, , drop = FALSE]
  data.frame(
    trial = overall_label,
    subgroup = comparison_forest_summary_interaction(s),
    yi = comparison_forest_summary_point(s),
    vi = pmax(s$sd^2, 1e-12, na.rm = TRUE),
    ci.lb = comparison_forest_interval_col(s, "lci"),
    ci.ub = comparison_forest_interval_col(s, "uci"),
    row_kind = "interaction_effect",
    stringsAsFactors = FALSE
  )
}

#' @noRd
comparison_forest_subgroup_rows <- function(dat, trial_var, subgroup_var, contrast_sep) {
  out <- dat[
    is_ifrac_component_row(dat[[subgroup_var]], contrast_sep = contrast_sep),
    ,
    drop = FALSE
  ]
  if (!nrow(out)) {
    stop("No subgroup rows found in `ifrac_data` for the forest plot.", call. = FALSE)
  }
  out
}

#' Reorder forest rows to match a trial-by-subgroup template
#' @noRd
comparison_forest_order_to_template <- function(df, template, trial_var = "study", subgroup_var = "level") {
  if (!nrow(df) || !nrow(template)) {
    return(df)
  }
  row_key <- function(tr, sg) paste(tr, as.character(sg), sep = "\r")
  keys <- row_key(template$trial, template$subgroup)
  df_keys <- row_key(df[[trial_var]], df[[subgroup_var]])
  idx <- match(keys, df_keys)
  if (any(is.na(idx))) {
    stop(
      "Forest rows do not align with the trial-by-subgroup template.",
      call. = FALSE
    )
  }
  df[idx, , drop = FALSE]
}

#' Align columns before binding forest panel rows
#' @noRd
comparison_forest_align_row_cols <- function(row, template) {
  if (!nrow(template)) {
    return(row)
  }
  cols <- union(names(template), names(row))
  for (nm in setdiff(cols, names(row))) {
    if (nm %in% c("fill", "colour")) {
      row[[nm]] <- NA_character_
    } else if (nm == "bold") {
      row[[nm]] <- FALSE
    } else if (nm == "row_kind") {
      row[[nm]] <- NA_character_
    } else if (nm == "shape") {
      row[[nm]] <- NA_integer_
    } else if (nm %in% c("ifrac", "trial_size")) {
      row[[nm]] <- ""
    } else {
      row[[nm]] <- NA
    }
  }
  row[, cols, drop = FALSE]
}

#' Blank effect columns on a forest panel row
#' @noRd
comparison_forest_blank_panel_row <- function(row) {
  row$yi <- NA_real_
  row$vi <- NA_real_
  if ("ci.lb" %in% names(row)) {
    row$ci.lb <- NA_real_
    row$ci.ub <- NA_real_
  }
  row
}

#' Interleave observed and shrinkage rows
#'
#' Within each heading block, pairs each observed row with its subgroup
#' shrinkage (and interaction shrinkage when present): study, then its shrinkage
#' row (grey), …
#' @noRd
interleave_comparison_shrinkage_rows <- function(resultsA,
                                           resultsB,
                                           row_labels,
                                           shrinkage,
                                           prob,
                                           ref,
                                           contrast_sep,
                                           ifrac_display = NULL) {
  n <- nrow(row_labels)
  sg <- shrinkage$subgroup
  if (!nrow(sg) || nrow(sg) != n) {
    stop(
      "Shrinkage rows must match the number of observed subgroup rows.",
      call. = FALSE
    )
  }

  shr_ifrac <- rep("", nrow(sg))
  sg_sum <- forestplot_panel_df(sg, prob = prob, ifrac_display = shr_ifrac)
  if ("row_kind" %in% names(resultsA)) {
    sg_sum$row_kind <- "shrinkage"
  }
  sg_sum$ifrac <- ""

  int_df <- shrinkage$interaction
  if (is.null(int_df) || !nrow(int_df)) {
    int_df <- build_contrast_rows_from_subgroups(
      sg,
      contrast_sep = contrast_sep,
      reference = ref
    )
  }
  int_aligned <- align_interaction_to_subgroups(sg, int_df, reference = ref)
  int_sum <- forestplot_panel_df(int_aligned$interaction, prob = prob)
  if ("row_kind" %in% names(resultsB)) {
    int_sum$row_kind <- "shrinkage"
  }

  trials <- unique(row_labels$subheading)
  outA <- list()
  outB <- list()
  out_lab <- list()
  ifrac_out <- if (!is.null(ifrac_display)) list() else NULL
  obs_row_index <- integer(0)
  subgroup_shrinkage_row_index <- integer(0)
  interaction_shrinkage_row_index <- integer(0)
  row_n <- 0L
  ref_chr <- as.character(ref)

  push_row <- function(panel_a, panel_b, lab_row, ifrac_val = NULL) {
    row_n <<- row_n + 1L
    outA[[row_n]] <<- panel_a
    outB[[row_n]] <<- panel_b
    out_lab[[row_n]] <<- lab_row
    if (!is.null(ifrac_out)) {
      ifrac_out[[row_n]] <<- if (is.null(ifrac_val)) "" else ifrac_val
    }
    row_n
  }

  for (tr in trials) {
    ii <- which(row_labels$subheading == tr)

    for (j in seq_along(ii)) {
      i <- ii[[j]]
      obsA <- resultsA[i, , drop = FALSE]
      obsB <- resultsB[i, , drop = FALSE]
      row_template <- comparison_forest_align_row_cols(sg_sum[i, , drop = FALSE], obsA)
      obsA <- comparison_forest_align_row_cols(obsA, row_template)
      obsB <- comparison_forest_align_row_cols(obsB, row_template)

      obs_lab <- row_labels[i, , drop = FALSE]
      colour_lab <- forestplot_row_colour_label(obs_lab)
      if (!"colour_label" %in% names(obs_lab) || is.na(obs_lab$colour_label[[1]])) {
        obs_lab$colour_label <- colour_lab
      }
      lab <- as.character(obs_lab$label[[1]])

      obs_idx <- push_row(
        obsA,
        obsB,
        obs_lab,
        if (!is.null(ifrac_display)) ifrac_display[[i]] else NULL
      )
      obs_row_index <- c(obs_row_index, obs_idx)

      # Subgroup shrinkage immediately after the matching observed row
      shrA <- comparison_forest_align_row_cols(sg_sum[i, , drop = FALSE], obsA)
      shrB <- comparison_forest_align_row_cols(comparison_forest_blank_panel_row(obsB), obsB)
      shr_lab <- row_labels[i, , drop = FALSE]
      shr_lab$label <- lab
      shr_lab$colour_label <- colour_lab
      sg_idx <- push_row(shrA, shrB, shr_lab, "")
      subgroup_shrinkage_row_index <- c(subgroup_shrinkage_row_index, sg_idx)

      # Interaction shrinkage for non-reference rows, still paired with the study
      show_int_shr <- nzchar(colour_lab) &&
        colour_lab != ref_chr &&
        is.finite(int_sum$yi[[i]])
      if (show_int_shr) {
        intA <- comparison_forest_align_row_cols(comparison_forest_blank_panel_row(obsA), obsA)
        intB <- comparison_forest_align_row_cols(int_sum[i, , drop = FALSE], obsB)
        int_lab <- row_labels[i, , drop = FALSE]
        int_lab$label <- lab
        int_lab$colour_label <- colour_lab
        int_idx <- push_row(intA, intB, int_lab, "")
        interaction_shrinkage_row_index <- c(
          interaction_shrinkage_row_index,
          int_idx
        )
      }
    }
  }

  resultsA <- do.call(rbind, outA)
  resultsB <- do.call(rbind, outB)
  if (!is.null(ifrac_display)) {
    if (length(ifrac_display) != n) {
      stop("Information fraction display length does not match observed rows.", call. = FALSE)
    }
    resultsA$ifrac <- unlist(ifrac_out, use.names = FALSE)
  }

  row_labels_out <- do.call(rbind, out_lab)
  if (!"colour_label" %in% names(row_labels_out)) {
    row_labels_out$colour_label <- row_labels_out$label
  }

  list(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels_out,
    n_obs = n,
    obs_row_index = obs_row_index,
    subgroup_shrinkage_row_index = subgroup_shrinkage_row_index,
    interaction_shrinkage_row_index = interaction_shrinkage_row_index,
    shrinkage_row_index = c(subgroup_shrinkage_row_index, interaction_shrinkage_row_index)
  )
}

#' @noRd
append_comparison_forest_blocks <- function(resultsA,
                                      resultsB,
                                      row_labels,
                                      blocks,
                                      prob,
                                      ref,
                                      contrast_sep) {
  for (block in blocks) {
    extra <- append_comparison_forest_block(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      subgroup = block$subgroup,
      interaction = block$interaction,
      block_name = block$name,
      prob = prob,
      ref = ref,
      contrast_sep = contrast_sep,
      panel_b = block$panel_b
    )
    resultsA <- extra$resultsA
    resultsB <- extra$resultsB
    row_labels <- extra$row_labels
  }
  list(resultsA = resultsA, resultsB = resultsB, row_labels = row_labels)
}

#' @noRd
append_comparison_forest_block <- function(resultsA,
                                     resultsB,
                                     row_labels,
                                     subgroup,
                                     interaction,
                                     block_name,
                                     prob,
                                     ref,
                                     contrast_sep,
                                     panel_b = NULL) {
  if (!nrow(subgroup)) {
    return(list(resultsA = resultsA, resultsB = resultsB, row_labels = row_labels))
  }

  n <- nrow(row_labels)
  n_add <- nrow(subgroup)

  sg_sum <- forestplot_panel_df(subgroup, prob = prob)
  if ("row_kind" %in% names(subgroup)) {
    sg_sum$row_kind <- subgroup$row_kind
  }
  if ("ifrac" %in% names(subgroup)) {
    sg_sum$ifrac <- as.character(subgroup$ifrac)
  }
  group_lab <- if ("colour_label" %in% names(subgroup)) {
    as.character(subgroup$colour_label)
  } else {
    as.character(subgroup$subgroup)
  }
  group_lab[is.na(group_lab)] <- ""

  sg_lab <- as.character(subgroup$subgroup)
  sg_lab[is.na(sg_lab)] <- ""

  row_kind_chr <- if ("row_kind" %in% names(subgroup)) {
    rk <- as.character(subgroup$row_kind)
    rk[is.na(rk)] <- ""
    rk
  } else {
    rep("", n_add)
  }

  block_lab <- if (!is.null(block_name) && nzchar(as.character(block_name)[[1]])) {
    as.character(block_name)[[1]]
  } else {
    ""
  }

  # Axis detail: block name, else Overall/Prediction from subgroup for overall
  # rows, else title-cased row_kind. Avoids NBSP placeholders and duplicate
  # "Younger" / "All" axis ticks.
  detail_lab <- if (nzchar(block_lab)) {
    rep(block_lab, n_add)
  } else {
    ifelse(
      row_kind_chr == "combined" & nzchar(sg_lab),
      sg_lab,
      ifelse(
        nzchar(row_kind_chr),
        paste0(toupper(substr(row_kind_chr, 1L, 1L)), substring(row_kind_chr, 2L)),
        ""
      )
    )
  }
  detail_lab[detail_lab == group_lab] <- ""

  extra_labels <- data.frame(
    subheading = group_lab,
    heading = "",
    label = detail_lab,
    colour_label = subgroup$colour_label,
    stringsAsFactors = FALSE
  )

  pad_forest_rows <- function(df, n_add) {
    if (!n_add) {
      return(df)
    }
    blank <- df[rep(1L, n_add), , drop = FALSE]
    for (col in names(blank)) {
      if (is.character(blank[[col]])) {
        blank[[col]] <- ""
      } else {
        blank[[col]] <- NA
      }
    }
    if ("row_kind" %in% names(blank)) {
      blank$row_kind <- NA_character_
    }
    if ("point_scale" %in% names(blank)) {
      blank$point_scale <- 1
    }
    rbind(df, blank)
  }

  fill_forest_block <- function(target, row_idx, source) {
    shared <- intersect(names(target), names(source))
    target[row_idx, shared] <- source[, shared, drop = FALSE]
    target
  }

  resultsA <- pad_forest_rows(resultsA, n_add)
  resultsB <- pad_forest_rows(resultsB, n_add)
  row_labels <- rbind(row_labels, extra_labels)

  row_idx <- (n + 1L):(n + n_add)
  resultsA <- fill_forest_block(resultsA, row_idx, sg_sum)

  if (!is.null(panel_b) && nrow(panel_b)) {
    b_sum <- forestplot_panel_df(panel_b, prob = prob)
    if ("row_kind" %in% names(subgroup)) {
      b_sum$row_kind <- subgroup$row_kind
    }
    resultsB <- fill_forest_block(resultsB, row_idx, b_sum)
  } else if (nrow(interaction)) {
    int_panel <- comparison_forest_interaction_panel_for_block(
      subgroup = subgroup,
      interaction = interaction,
      prob = prob,
      ref = ref
    )
    resultsB <- fill_forest_block(resultsB, row_idx, int_panel)
  } else {
    resultsB <- fill_forest_block(resultsB, row_idx, resultsA[row_idx, , drop = FALSE])
  }

  list(resultsA = resultsA, resultsB = resultsB, row_labels = row_labels)
}

#' Match interaction effect rows to a subgroup summary block
#' @noRd
comparison_forest_interaction_panel_for_block <- function(subgroup,
                                                    interaction,
                                                    prob,
                                                    ref) {
  aligned <- align_interaction_to_subgroups(subgroup, interaction, reference = ref)
  int_rows <- forestplot_panel_df(aligned$interaction, prob = prob)
  int_rows
}
