#' Shared IF / contribution-plot helpers
#' @noRd
detect_ifrac_vars <- function(dat) {
  cols <- grep("^ifrac_", names(dat), value = TRUE)
  cols <- cols[!grepl("^ifrac2_", cols)]
  sort(cols)
}

resolve_bubbleplot_settings <- function(object) {
  st <- attr(object, "interactmeta", exact = TRUE)
  raw_dat <- if (!is.null(st) && !is.null(st$ifrac_data)) {
    st$ifrac_data
  } else {
    object$data
  }
  ifrac_vars <- detect_ifrac_vars(raw_dat)

  if (!is.null(st) && !is.null(st$ifrac_var)) {
    return(list(
      ifrac_var = st$ifrac_var,
      ifrac_vars = ifrac_vars,
      ifrac2_var = st$ifrac2_var,
      trial_var = st$trial_var,
      subgroup_var = st$subgroup_var,
      effects = st$effects,
      ifrac_data = st$ifrac_data,
      reference = st$reference
    ))
  }

  st <- interactmeta_settings(object)
  st$ifrac_vars <- ifrac_vars
  st$ifrac_data <- NULL
  st$reference <- NULL
  if ("reference" %in% names(raw_dat)) {
    ref <- unique(raw_dat$reference)
    ref <- ref[!is.na(ref) & nzchar(ref)]
    if (length(ref) == 1L) {
      st$reference <- ref[[1]]
    }
  }
  st
}

resolve_bubbleplot_raw <- function(object, settings) {
  if (!is.null(settings$ifrac_data)) {
    return(as.data.frame(settings$ifrac_data))
  }
  as.data.frame(object$data)
}

resolve_bubbleplot_limits <- function(ifrac_vals, y_vals, xlim, ylim, x_breaks, y_breaks) {
  ifrac_vals <- ifrac_vals[is.finite(ifrac_vals)]
  # y_vals kept for call-site compatibility; treatment-effect limits are not
  # auto-derived so ggplot2 can use its defaults when ylim is NULL.
  invisible(y_vals)

  # Information fraction is typically in [0, 1], but do not force the axis to
  # span that full range — pad from the observed values only.
  if (is.null(xlim) && length(ifrac_vals)) {
    rng <- range(ifrac_vals)
    pad <- max(0.02, 0.05 * diff(rng))
    xlim <- c(max(0, rng[1] - pad), rng[2] + pad)
  }
  # Treatment-effect ylim stays NULL unless the user sets it, so ggplot2 can
  # choose default limits from the geoms.
  if (is.null(x_breaks) && !is.null(xlim)) {
    x_breaks <- pretty(xlim, n = 5)
  }
  if (is.null(y_breaks) && !is.null(ylim)) {
    y_breaks <- pretty(ylim, n = 5)
  }

  list(xlim = xlim, ylim = ylim, x_breaks = x_breaks, y_breaks = y_breaks)
}

#' ggplot2 breaks: keep auto ticks when the user did not set them
#'
#' `scale_*_continuous(breaks = NULL)` draws no tick labels. Use a waiver
#' instead so ggplot2 chooses a default.
#' @noRd
interactmeta_continuous_breaks <- function(breaks) {
  if (is.null(breaks)) {
    ggplot2::waiver()
  } else {
    breaks
  }
}

#' IF evaluation grid for conditional-effect ribbons
#'
#' Uses the observed information-fraction range (padded), not a fixed 0–1 span.
#' @noRd
bubbleplot_ifrac_grid <- function(ifrac_vals, resolution = 200) {
  n <- max(2L, as.integer(resolution)[[1L]])
  ifrac_vals <- ifrac_vals[is.finite(ifrac_vals)]
  if (!length(ifrac_vals)) {
    return(seq(0, 1, length.out = n))
  }
  rng <- range(ifrac_vals)
  pad <- max(0.02, 0.05 * diff(rng))
  lo <- max(0, rng[[1L]] - pad)
  hi <- min(1, rng[[2L]] + pad)
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    return(seq(0, 1, length.out = n))
  }
  seq(lo, hi, length.out = n)
}

#' Require a CAMS model fit (`interactmeta`)
#' @noRd
check_comparison_fit <- function(object = NULL) {
  if (!is.null(object)) {
    interactmeta_abort_non_cams(object, "This method")
    if (is_comparison_fit(object)) {
      return(invisible(NULL))
    }
    stop("`object` must be an `interactmeta` from `cam()`.", call. = FALSE)
  }
  invisible(NULL)
}

bubbleplot_has_ce <- function(ce_df) {
  is.data.frame(ce_df) &&
    nrow(ce_df) > 0L &&
    all(c("estimate__", "lower__", "upper__", "ifrac_x") %in% names(ce_df))
}

bubbleplot_ce_estimates <- function(ce_df) {
  if (!bubbleplot_has_ce(ce_df)) {
    return(numeric())
  }
  c(ce_df$estimate__, ce_df$lower__, ce_df$upper__)
}

prepare_conditional_effects_df <- function(ce_df, subgroup_var, ifrac_var, subgroup_labels) {
  sub_col <- setdiff(
    names(ce_df),
    c(ifrac_var, "estimate__", "se__", "lower__", "upper__", "cond__")
  )
  sub_col <- sub_col[grepl(subgroup_var, sub_col, fixed = TRUE)][1L]
  if (is.na(sub_col)) {
    sub_col <- "effect2__"
    if (!(sub_col %in% names(ce_df))) {
      sub_col <- names(ce_df)[2L]
    }
  }

  lvls <- unique(as.character(ce_df[[sub_col]]))
  display_labels <- resolve_subgroup_display_labels(lvls, subgroup_labels)

  ce_df$display <- factor(as.character(ce_df[[sub_col]]), levels = lvls, labels = display_labels)
  ce_df
}

prepare_plot_raw_df <- function(raw, subgroup_var, ifrac_var, study_var, subgroup_labels, level_order) {
  lvls <- unique(as.character(raw[[subgroup_var]]))
  display_labels <- resolve_subgroup_display_labels(lvls, subgroup_labels)

  raw$display <- factor(as.character(raw[[subgroup_var]]), levels = lvls, labels = display_labels)
  if (!is.null(level_order)) {
    raw$display <- factor(raw$display, levels = level_order)
  }

  if ("se" %in% names(raw)) {
    raw$precision <- 1 / pmax(raw$se, 1e-6)
  } else if ("precision" %in% names(raw)) {
    raw$precision <- raw$precision
  } else {
    raw$precision <- 1 / pmax(sqrt(raw$vi), 1e-6)
  }

  raw
}

resolve_subgroup_display_labels <- function(levels, subgroup_labels) {
  if (is.null(subgroup_labels)) {
    return(levels)
  }
  if (length(subgroup_labels) != length(levels)) {
    stop(
      "`subgroup_labels` must have the same length as subgroup levels in the data.",
      call. = FALSE
    )
  }
  subgroup_labels
}

#' Wald interval for observed study effects; optionally exponentiate
#'
#' Intervals are formed on the log scale (`yi ± z * se`) and then
#' exponentiated with the point estimate when `exponentiate = TRUE`.
#' @noRd
attach_observed_effect_interval <- function(raw_df,
                                            prob = 0.95,
                                            exponentiate = FALSE) {
  z <- stats::qnorm((1 - as.numeric(prob)[[1L]]) / 2, lower.tail = FALSE)
  if ("se" %in% names(raw_df)) {
    se <- as.numeric(raw_df$se)
  } else if ("vi" %in% names(raw_df)) {
    se <- sqrt(pmax(as.numeric(raw_df$vi), 0))
  } else {
    se <- NULL
  }
  if (!is.null(se)) {
    raw_df$yi_lci <- raw_df$yi - z * se
    raw_df$yi_uci <- raw_df$yi + z * se
  }
  if (isTRUE(exponentiate)) {
    raw_df$yi <- exp(raw_df$yi)
    if ("yi_lci" %in% names(raw_df)) {
      raw_df$yi_lci <- exp(raw_df$yi_lci)
      raw_df$yi_uci <- exp(raw_df$yi_uci)
    }
  }
  raw_df
}

#' Null-effect reference line (0 on the log scale, 1 on the ratio scale)
#' @noRd
add_effects_null_line <- function(p_out, y_vals, exponentiate, const) {
  ref <- if (isTRUE(exponentiate)) 1 else 0
  y_rng <- range(y_vals[is.finite(y_vals)], na.rm = TRUE)
  if (!(length(y_rng) == 2L && is.finite(y_rng[[1L]]) && is.finite(y_rng[[2L]]) &&
        y_rng[[1L]] <= ref && y_rng[[2L]] >= ref)) {
    return(p_out)
  }
  p_out +
    ggplot2::geom_hline(
      yintercept = ref,
      colour = const$zero_line,
      linewidth = const$line_hair,
      inherit.aes = FALSE
    )
}

trial_yi_segments <- function(raw, study_var, ifrac_var) {
  raw |>
    dplyr::group_by(.data[[study_var]]) |>
    dplyr::filter(dplyr::n() == 2L) |>
    dplyr::summarise(
      ifrac_x = .data[[ifrac_var]][1],
      yi_1 = .data$yi[1],
      yi_2 = .data$yi[2],
      .groups = "drop"
    )
}

resolve_bubbleplot_level_keys <- function(raw_df, display_levels, subgroup_var) {
  if (!subgroup_var %in% names(raw_df)) {
    return(display_levels)
  }
  key_tab <- dplyr::distinct(
    raw_df,
    dplyr::pick(dplyr::all_of(c(subgroup_var, "display")))
  )
  key_tab[[subgroup_var]][match(display_levels, as.character(key_tab$display))]
}

resolve_ifrac_plot_style <- function(display_levels,
                                    level_keys = NULL,
                                    reference = NULL,
                                    subgroup_labels,
                                    legend_math,
                                    colors,
                                    fills) {
  lv <- display_levels
  if (is.null(level_keys)) {
    level_keys <- lv
  }

  sg <- subgroup_colour_map(
    level_keys,
    reference = reference,
    palette = interactmeta_palette_name()
  )
  default_colors <- stats::setNames(sg$colour, lv)
  default_fills <- stats::setNames(sg$fill, lv)

  clrs <- if (is.null(colors)) default_colors else colors
  fill_vals <- if (is.null(fills)) default_fills else fills
  fill_clrs <- stats::setNames(fill_vals, names(fill_vals))

  legend_labels <- if (is.null(legend_math)) {
    lv
  } else {
    labels <- legend_math[lv]
    labels[is.na(labels)] <- lv[is.na(labels)]
    interactmeta_parse_legend_labels(labels, fallback = lv)
  }

  list(
    colors = clrs,
    fill_colors = fill_clrs,
    legend_labels = legend_labels,
    legend_title = "Subgroup"
  )
}
