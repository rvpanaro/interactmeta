# ---- contribution.R ----
#' CAMS subgroup and overall effects on an information-fraction grid
#'
#' Evaluates the **fitted** CAMS posterior at a grid of information fractions
#' without refitting. Each posterior row already stores the joint pairing of
#' \eqn{(\alpha,\theta)} (overall block) with the independent interaction
#' draw \eqn{\gamma}. For draw \eqn{s} and grid value \eqn{\pi_k}:
#' \deqn{M^{(s)}(\pi_k)=\alpha^{(s)}+\theta^{(s)}\pi_k,}
#' \deqn{\mu_A^{(s)}(\pi_k)=M^{(s)}(\pi_k)-\pi_k\gamma^{(s)},}
#' \deqn{\mu_B^{(s)}(\pi_k)=M^{(s)}(\pi_k)+(1-\pi_k)\gamma^{(s)}.}
#' Stored `b_Intercept` is \eqn{\alpha}, stored `b_ifrac` is
#' \eqn{\beta=\theta-\gamma}, and stored `b_subgroup*` is \eqn{\gamma}, so
#' \eqn{\theta=\beta+\gamma}. Subgroup lines are parallel with vertical
#' gap \eqn{\gamma}.
#'
#' Pointwise `estimate` / `lower` / `upper` use mixture intervals from the
#' DIRECT-node Gaussian mixture (weights `.weight`, within-node SDs
#' `.cond_*`). Alpha and theta stay on the joint overall-block support
#' (never convolved with each other).
#'
#' @param fit A [`interactmeta`] from [cams()].
#' @param if_grid Numeric information-fraction grid. If `NULL`, spans the
#'   observed study IF range with a small margin (see [plot_contribution()]).
#' @param level Coverage probability for `lower` / `upper` (default `0.95`).
#' @param ifrac_var `ifrac_*` column used to build the default grid.
#' @param resolution Grid length when `if_grid` is `NULL`.
#' @param max_support Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports; omit it.
#' @param exponentiate If `TRUE`, summarize on the exponential scale.
#' @param interval.type `"central"` or `"shortest"` mixture interval on
#'   Bayesian DIRECT nodes (default: `fit$interval.type`, usually `"central"`).
#'   Pass `"shortest"` for the narrowest equal-probability interval.
#'
#' @return A data frame with one row per grid point and subgroup:
#'   `ifrac`, `subgroup`, `estimate`, `lower`, `upper`, plus
#'   `overall_estimate`, `overall_lower`, `overall_upper` for
#'   \eqn{M(\pi)=\alpha+\theta\pi} at that grid point.
#'
#' @seealso [modifiers()], [combined()], [plot_contribution()]
#' @name predict_cams_by_if
#' @export
predict_cams_by_if <- function(fit,
                               if_grid = NULL,
                               level = 0.95,
                               ifrac_var = NULL,
                               resolution = 200,
                               max_support = NULL,
                               exponentiate = FALSE,
                               interval.type = NULL) {
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an `interactmeta` from `cams()`.", call. = FALSE)
  }
  level <- as.numeric(level)[[1L]]
  if (!is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a coverage probability in (0, 1).", call. = FALSE)
  }
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% fit$interval.type
  )

  settings <- resolve_bubbleplot_settings(fit)
  if (is.null(ifrac_var)) {
    ifrac_var <- settings$ifrac_vars[[1L]] %||% settings$ifrac_var
  }
  if (is.null(if_grid)) {
    raw <- resolve_bubbleplot_raw(fit, settings)
    if (is.null(ifrac_var) || !ifrac_var %in% names(raw)) {
      stop("Could not resolve an `ifrac_*` column for the default IF grid.", call. = FALSE)
    }
    if_grid <- bubbleplot_ifrac_grid(raw[[ifrac_var]], resolution = resolution)
  }
  if_grid <- as.numeric(if_grid)
  if (!length(if_grid) || any(!is.finite(if_grid))) {
    stop("`if_grid` must be a finite numeric vector.", call. = FALSE)
  }

  design <- comparison_marginal_design(fit)
  post <- comparison_support_df(fit)
  idx <- comparison_marginal_draw_index(
    nrow(post),
    max_support = max_support,
    weight = if (".weight" %in% names(post)) post$.weight else NULL
  )
  post <- post[idx, , drop = FALSE]
  sg_coef <- comparison_subgroup_coef_terms(fit, design)
  ref <- as.character(design$reference)[[1L]]

  rows <- lapply(if_grid, function(pi) {
    m_draws <- comparison_marginal_overall_effect_draws(
      post = post,
      omega = pi,
      design = design,
      fit = fit,
      type = "Overall"
    )
    m_qi <- predict_cam_if_qi(
      x = m_draws,
      post = post,
      a_alpha = 1,
      a_theta = pi,
      a_gamma = 0,
      level = level,
      exponentiate = exponentiate,
      fit = fit,
      interval.type = interval.type
    )
    sg_rows <- lapply(design$subgroup_levels, function(lbl) {
      mu <- comparison_marginal_fe_draws(
        post = post,
        evaluation_if = pi,
        design = design,
        subgroup_coef = sg_coef[[lbl]]
      )
      x_nonref <- as.numeric(!identical(as.character(lbl), ref))
      a_gamma <- ifelse(x_nonref == 0, -pi, 1 - pi)
      qi <- predict_cam_if_qi(
        x = mu,
        post = post,
        a_alpha = 1,
        a_theta = pi,
        a_gamma = a_gamma,
        level = level,
        exponentiate = exponentiate,
        fit = fit,
        interval.type = interval.type
      )
      data.frame(
        ifrac = pi,
        level = lbl,
        estimate = qi$estimate,
        lower = qi$lci,
        upper = qi$uci,
        overall_estimate = m_qi$estimate,
        overall_lower = m_qi$lci,
        overall_upper = m_qi$uci,
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(sg_rows)
  })
  out <- dplyr::bind_rows(rows)
  attr(out, "level") <- level
  attr(out, "reference") <- ref
  attr(out, "ifrac_var") <- ifrac_var
  if (interactmeta_is_bayesian_fit(fit)) {
    attr(out, "interval.rule") <- if (identical(interval.type, "central")) {
      "mixture_equal_tail"
    } else {
      "mixture_shortest"
    }
  }
  attr(out, "interval.type") <- interval.type
  class(out) <- c("interactmeta_cam_if_pred", "data.frame")
  out
}

#' Pointwise CAMS interval at one information fraction.
#'
#' Mixture interval at one IF (shortest or central on DIRECT nodes).
#' Alpha and theta stay jointly paired.
#' @noRd
predict_cam_if_qi <- function(x, post, a_alpha, a_theta, a_gamma, level, exponentiate,
                               fit = NULL, interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% if (is.list(fit)) fit$interval.type else NULL
  )
  if (!is.null(fit) && !interactmeta_is_bayesian_fit(fit)) {
    lin <- interactmeta_metafor_lincomb(
      fit,
      L_ov = c(a_alpha, a_theta),
      L_int = a_gamma,
      exponentiate = exponentiate,
      prob = level
    )
    if (!is.null(lin)) {
      return(lin)
    }
  }
  cond_sd <- if (interactmeta_has_cond_moments(post)) {
    interactmeta_lincomb_cond_sd(post, a_alpha, a_theta, a_gamma)
  } else {
    NULL
  }
  w <- if (".weight" %in% names(post)) post$.weight else NULL
  comparison_draw_moments(
    x,
    prob = level,
    weight = w,
    cond_sd = cond_sd,
    exponentiate = exponentiate,
    interval.type = interval.type
  )
}

#' Contribution-adjusted CAMS meta-regression over information fraction
#'
#' Core contribution-adjusted display: population subgroup curves
#' \eqn{\mu_A(\pi)} / \eqn{\mu_B(\pi)} with pointwise credible ribbons from
#' [predict_cams_by_if()], the overall CAMS line
#' \eqn{M(\pi)=\alpha+\theta\pi}, and observed study effects at each trial's
#' own \eqn{\pi_j}. This is the plot of the CAMS identities
#' \eqn{\mu_B(\pi)-\mu_A(\pi)=\gamma} and
#' \eqn{M(\pi)=(1-\pi)\mu_A(\pi)+\pi\mu_B(\pi)}.
#'
#' @param fit A [`interactmeta`] / `interactmeta` from [cams()].
#' @param ifrac_var,prob,type,max_support,contrast_sep Passed to
#'   [predict_cams_by_if()] and [ifrac()].
#' @param exponentiate If `TRUE`, draw the CAMS curves, study points, and
#'   intervals on the exponential (ratio) scale. Observed `yi` are exponentiated
#'   after the log-scale Wald interval is formed.
#' @param subgroup_var,trial_var Optional subgroup and trial column names.
#'   Defaults to the names stored on `fit`.
#' @param at_ifrac Optional dotted reference IF; default `NULL` (no marker).
#' @param oi_support Optional overall-IF draws (see [ifrac()]).
#' @param effects,re_formula,resolution Passed to the conditional-effects helper
#'   (`re_formula = NA` by default, population effects only).
#' @param subgroup_labels,legend_math,colors,fills,ylab Plot styling.
#' @param xlab Axis label; default `"Information fraction (IF) (%)"`.
#' @param xlim,ylim,x_breaks,y_breaks Axis limits and breaks (`xlim` / `x_breaks`
#'   are on the **percent** scale, 0–100).
#' @param show_overall_line Logical; dashed \eqn{M(\pi)=\alpha+\theta\pi} line
#'   (legend \dQuote{All} when drawn).
#' @param ribbon_interval.type `"central"` or `"shortest"` credible interval
#'   for subgroup ribbon limits on the DIRECT mixture (default `"central"`).
#' @param point_size_range Range for precision-scaled study points.
#' @param show_footnote If `TRUE` (default), add a caption below the plot with
#'   the reporting information fraction (percent, interval), and study count.
#'   Use `FALSE` to omit it.
#' @param digits Decimal places in the footnote.
#' @param base_size Base **ggplot2** theme size.
#' @param ... Ignored.
#'
#' @return A **ggplot2** object, invisibly.
#'
#' @seealso [predict_cams_by_if()], [ifrac()]
#' @examples
#' \dontrun{
#' plot(fit)
#' }
#' @keywords internal
plot_contribution <- function(fit,
                                    ifrac_var = NULL,
                                    prob = 0.95,
                                    type = "Overall",
                                    max_support = NULL,
                                    contrast_sep = "-",
                                    exponentiate = FALSE,
                                    subgroup_var = NULL,
                                    at_ifrac = NULL,
                                    oi_support = NULL,
                                    trial_var = NULL,
                                    effects = NULL,
                                    re_formula = NA,
                                    resolution = 200,
                                    subgroup_labels = NULL,
                                    legend_math = NULL,
                                    colors = NULL,
                                    fills = NULL,
                                    xlab = NULL,
                                    ylab = "Treatment effect",
                                    xlim = NULL,
                                    ylim = NULL,
                                    x_breaks = NULL,
                                    y_breaks = NULL,
                                    show_overall_line = TRUE,
                                    ribbon_interval.type = "central",
                                    point_size_range = c(2.5, 6.5),
                                    show_footnote = TRUE,
                                    digits = 2,
                                    base_size = interactmeta_base_size(default = 11L),
                                    ...) {
  check_comparison_fit(fit)
  p <- build_contribution_plot(
    fit = fit,
    ifrac_var = ifrac_var,
    prob = prob,
    type = type,
    max_support = max_support,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate,
    subgroup_var = subgroup_var,
    at_ifrac = at_ifrac,
    oi_support = oi_support,
    trial_var = trial_var,
    effects = effects,
    re_formula = re_formula,
    resolution = resolution,
    subgroup_labels = subgroup_labels,
    legend_math = legend_math,
    colors = colors,
    fills = fills,
    xlab = xlab,
    ylab = ylab,
    xlim = xlim,
    ylim = ylim,
    x_breaks = x_breaks,
    y_breaks = y_breaks,
    show_overall_line = show_overall_line,
    ribbon_interval.type = ribbon_interval.type,
    point_size_range = point_size_range,
    show_footnote = show_footnote,
    digits = digits,
    base_size = base_size
  )
  print(p)
  invisible(p)
}

#' @noRd
build_contribution_plot <- function(fit,
                                          ifrac_var = NULL,
                                          prob = 0.95,
                                          type = "Overall",
                                          max_support = NULL,
                                          contrast_sep = "-",
                                          exponentiate = FALSE,
                                          subgroup_var = NULL,
                                          at_ifrac = NULL,
                                          oi_support = NULL,
                                          trial_var = NULL,
                                          effects = NULL,
                                          re_formula = NA,
                                          resolution = 200,
                                          subgroup_labels = NULL,
                                          legend_math = NULL,
                                          colors = NULL,
                                          fills = NULL,
                                          xlab = NULL,
                                          ylab = "Treatment effect",
                                          xlim = NULL,
                                          ylim = NULL,
                                          x_breaks = NULL,
                                          y_breaks = NULL,
                                          show_overall_line = TRUE,
                                          ribbon_interval.type = "central",
                                          point_size_range = c(2.5, 6.5),
                                          show_footnote = TRUE,
                                          digits = 2,
                                          base_size = interactmeta_base_size(default = 11L)) {
  ribbon_interval.type <- interactmeta_match_interval_type(ribbon_interval.type)
  settings <- resolve_bubbleplot_settings(fit)
  if (is.null(ifrac_var)) {
    ifrac_var <- settings$ifrac_vars[[1]]
    if (is.null(ifrac_var)) {
      ifrac_var <- settings$ifrac_var
    }
  }
  if (is.null(subgroup_var)) {
    subgroup_var <- settings$subgroup_var
  }
  study_var <- if (is.null(trial_var)) settings$trial_var else trial_var
  raw <- resolve_bubbleplot_raw(fit, settings)

  pred <- predict_cams_by_if(
    fit,
    if_grid = NULL,
    level = prob,
    ifrac_var = ifrac_var,
    resolution = resolution,
    max_support = max_support,
    exponentiate = exponentiate,
    interval.type = ribbon_interval.type
  )
  ce_df <- as.data.frame(pred)
  ce_df$ifrac_x <- ce_df$ifrac
  ce_df$estimate__ <- ce_df$estimate
  ce_df$lower__ <- ce_df$lower
  ce_df$upper__ <- ce_df$upper
  ce_df[[ifrac_var]] <- ce_df$ifrac
  src_col <- interactmeta_modifier_summary_col(ce_df)
  ce_df[[subgroup_var]] <- ce_df[[src_col]]
  ce_df <- prepare_conditional_effects_df(
    ce_df,
    subgroup_var = subgroup_var,
    ifrac_var = ifrac_var,
    subgroup_labels = subgroup_labels
  )
  if (!nrow(ce_df)) {
    stop(
      "Could not compute CAMS IF predictions; check the fit and `ifrac_var`.",
      call. = FALSE
    )
  }

  display_levels <- levels(ce_df$display)
  if (is.null(display_levels)) {
    display_levels <- unique(as.character(ce_df$display))
  }

  raw_df <- prepare_plot_raw_df(
    raw,
    subgroup_var = subgroup_var,
    ifrac_var = ifrac_var,
    study_var = study_var,
    subgroup_labels = subgroup_labels,
    level_order = display_levels
  )
  raw_df$ifrac_x <- raw_df[[ifrac_var]]
  raw_df <- attach_observed_effect_interval(
    raw_df,
    prob = prob,
    exponentiate = exponentiate
  )
  seg_df <- trial_yi_segments(raw_df, study_var = study_var, ifrac_var = "ifrac_x")

  overall_df <- unique(ce_df[, c("ifrac_x", "overall_estimate", "overall_lower", "overall_upper")])
  names(overall_df)[names(overall_df) == "overall_estimate"] <- "estimate__"

  level_keys <- resolve_bubbleplot_level_keys(
    raw_df = raw_df,
    display_levels = display_levels,
    subgroup_var = subgroup_var
  )
  style <- resolve_ifrac_plot_style(
    display_levels = display_levels,
    level_keys = level_keys,
    reference = settings$reference,
    subgroup_labels = subgroup_labels,
    legend_math = legend_math,
    colors = colors,
    fills = fills
  )

  xlab_user <- xlab
  if (is.null(xlab_user)) {
    xlab <- bubbleplot_default_xlab(ifrac_var)
  } else {
    xlab <- interactmeta_as_ggplot_label(xlab_user)
  }

  y_vals <- c(
    bubbleplot_ce_estimates(ce_df),
    raw_df$yi,
    raw_df$yi_lci,
    raw_df$yi_uci
  )
  limits <- resolve_bubbleplot_limits(
    ifrac_vals = c(ce_df$ifrac_x, raw_df$ifrac_x, at_ifrac),
    y_vals = y_vals,
    xlim = xlim,
    ylim = ylim,
    x_breaks = x_breaks,
    y_breaks = y_breaks
  )

  ce_df$ifrac_x <- interactmeta_ifrac_fraction_to_axis_pct(ce_df$ifrac_x)
  raw_df$ifrac_x <- interactmeta_ifrac_fraction_to_axis_pct(raw_df$ifrac_x)
  if (nrow(seg_df)) {
    seg_df$ifrac_x <- interactmeta_ifrac_fraction_to_axis_pct(seg_df$ifrac_x)
    seg_df$ifrac_xend <- seg_df$ifrac_x
  }
  if (nrow(overall_df)) {
    overall_df$ifrac_x <- interactmeta_ifrac_fraction_to_axis_pct(overall_df$ifrac_x)
  }
  at_ifrac_pct <- NULL
  if (!is.null(at_ifrac) && length(at_ifrac) == 1L && is.finite(at_ifrac)) {
    at_ifrac_pct <- interactmeta_ifrac_fraction_to_axis_pct(at_ifrac)
  }

  const <- interactmeta_plot_constants()
  p_out <- ggplot2::ggplot(ce_df, ggplot2::aes(
    x = .data$ifrac_x,
    y = .data$estimate__,
    colour = .data$display,
    fill = .data$display
  ))

  p_out <- add_effects_null_line(
    p_out,
    y_vals = y_vals,
    exponentiate = exponentiate,
    const = const
  )

  if (!is.null(at_ifrac_pct)) {
    p_out <- p_out +
      ggplot2::geom_vline(
        xintercept = at_ifrac_pct,
        colour = scales::alpha(const$overall_line_colour, 0.55),
        linewidth = 0.4,
        linetype = "dotted",
        inherit.aes = FALSE
      )
  }

  p_out <- p_out +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower__, ymax = .data$upper__),
      alpha = const$ribbon_alpha,
      colour = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = const$line_strong, show.legend = FALSE)

  overall_lty_key <- "__overall_mp__"
  show_overall_legend <- isTRUE(show_overall_line) && nrow(overall_df) > 0L
  if (show_overall_legend) {
    overall_df$overall_lty <- overall_lty_key
    p_out <- p_out +
      ggplot2::geom_line(
        data = overall_df,
        ggplot2::aes(
          x = .data$ifrac_x,
          y = .data$estimate__,
          linetype = .data$overall_lty
        ),
        colour = const$overall_line_colour,
        linewidth = const$line_mid,
        inherit.aes = FALSE
      )
  }

  if (nrow(seg_df)) {
    p_out <- p_out +
      ggplot2::geom_segment(
        data = seg_df,
        ggplot2::aes(
          x = .data$ifrac_x,
          xend = .data$ifrac_x,
          y = .data$yi_1,
          yend = .data$yi_2
        ),
        colour = const$segment_colour,
        linewidth = const$line_mid,
        inherit.aes = FALSE
      )
  }

  if (all(c("yi_lci", "yi_uci") %in% names(raw_df))) {
    p_out <- p_out +
      ggplot2::geom_linerange(
        data = raw_df,
        ggplot2::aes(
          x = .data$ifrac_x,
          ymin = .data$yi_lci,
          ymax = .data$yi_uci,
          colour = .data$display
        ),
        linewidth = const$line_hair,
        alpha = 0.55,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  p_out <- p_out +
    ggplot2::geom_point(
      data = raw_df,
      ggplot2::aes(
        x = .data$ifrac_x,
        y = .data$yi,
        colour = .data$display,
        size = .data$precision
      ),
      shape = 21,
      fill = NA,
      stroke = 0.55,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = style$colors,
      name = style$legend_title,
      na.value = const$na_colour,
      labels = style$legend_labels
    ) +
    ggplot2::scale_fill_manual(
      values = style$fill_colors,
      name = style$legend_title,
      na.value = const$na_colour,
      labels = style$legend_labels
    ) +
    ggplot2::scale_size_continuous(range = point_size_range, guide = "none") +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::scale_x_continuous(
      breaks = interactmeta_continuous_breaks(limits$x_breaks),
      labels = interactmeta_ifrac_axis_pct_labels,
      expand = interactmeta_effects_axis_expansion()
    ) +
    ggplot2::scale_y_continuous(
      breaks = interactmeta_continuous_breaks(limits$y_breaks),
      expand = interactmeta_effects_axis_expansion()
    ) +
    ifrac_effects_theme(base_size = base_size)

  if (isTRUE(show_footnote)) {
    oi_obj <- oi_support
    if (is.null(oi_obj)) {
      oi_obj <- ifrac(fit, ifrac_var = ifrac_var, prob = prob, summary = FALSE)
    }
    bundle <- overall_ifrac_bundle(oi_obj, prob = prob)
    cap_txt <- overall_ifrac_footnote_bundle(bundle, digits = digits)
    p_out <- p_out +
      ggplot2::labs(caption = forestplot_footnote_as_plotmath(cap_txt))
  }

  coord_lim <- list()
  if (!is.null(limits$ylim)) {
    coord_lim$ylim <- limits$ylim
  }
  if (!is.null(limits$xlim)) {
    coord_lim$xlim <- limits$xlim
  }
  if (length(coord_lim)) {
    p_out <- p_out + do.call(ggplot2::coord_cartesian, coord_lim)
  }

  p_out <- p_out + interactmeta_effects_legend_guide(
    base_size = base_size,
    colour_order = 1L
  )
  if (show_overall_legend) {
    p_out <- p_out +
      ggplot2::scale_linetype_manual(
        name = NULL,
        values = stats::setNames("22", overall_lty_key),
        labels = stats::setNames("All", overall_lty_key),
        guide = ggplot2::guide_legend(
          order = 2L,
          override.aes = list(
            colour = const$overall_line_colour,
            linewidth = const$line_mid
          )
        )
      )
  }
  p_out
}

# ---- bubbleplot.R ----
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

#' Map stored IF (0–1) to contribution-plot x-axis (percent)
#' @noRd
interactmeta_ifrac_fraction_to_axis_pct <- function(x) {
  100 * x
}

#' Percent tick labels for contribution-plot x-axis
#' @noRd
interactmeta_ifrac_axis_pct_labels <- function(x) {
  paste0(format(as.numeric(x), trim = TRUE, digits = 4), "%")
}

resolve_bubbleplot_limits <- function(ifrac_vals, y_vals, xlim, ylim, x_breaks, y_breaks) {
  ifrac_vals <- ifrac_vals[is.finite(ifrac_vals)]
  # y_vals kept for call-site compatibility; treatment-effect limits are not
  # auto-derived so ggplot2 can use its defaults when ylim is NULL.
  invisible(y_vals)

  # User `xlim` / `x_breaks` are on the percent (0–100) axis; convert for padding.
  xlim_frac <- if (!is.null(xlim)) xlim / 100 else NULL
  x_breaks_frac <- if (!is.null(x_breaks)) x_breaks / 100 else NULL

  # Information fraction is typically in [0, 1], but do not force the axis to
  # span that full range — pad from the observed values only.
  if (is.null(xlim_frac) && length(ifrac_vals)) {
    rng <- range(ifrac_vals)
    pad <- max(0.02, 0.05 * diff(rng))
    xlim_frac <- c(max(0, rng[1] - pad), rng[2] + pad)
  }
  xlim <- if (!is.null(xlim_frac)) interactmeta_ifrac_fraction_to_axis_pct(xlim_frac) else NULL
  # Treatment-effect ylim stays NULL unless the user sets it, so ggplot2 can
  # choose default limits from the geoms.
  if (is.null(x_breaks_frac) && !is.null(xlim_frac)) {
    x_breaks_frac <- pretty(xlim_frac, n = 5)
  }
  x_breaks <- if (!is.null(x_breaks_frac)) {
    interactmeta_ifrac_fraction_to_axis_pct(x_breaks_frac)
  } else {
    NULL
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
    stop("`object` must be an `interactmeta` from `cams()`.", call. = FALSE)
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

# ---- plot_labels.R ----
#' Convert a label to a ggplot plotmath object (expression or character)
#'
#' Character strings with plotmath markup (`~`, `[`, `^`, `==`) are parsed;
#' existing expressions and plain text without markup are passed through.
#'
#' @param label Character scalar, expression, or language object.
#' @noRd
interactmeta_as_ggplot_label <- function(label) {
  if (is.null(label)) {
    return(NULL)
  }
  if (inherits(label, "expression")) {
    return(label)
  }
  if (is.language(label)) {
    return(as.expression(label))
  }
  if (!is.character(label) || length(label) != 1L) {
    return(label)
  }
  if (!nzchar(label)) {
    return(label)
  }
  if (grepl("~|\\[|\\]|\\^|==|\\*", label)) {
    parsed <- tryCatch(parse(text = label), error = function(e) NULL)
    if (!is.null(parsed)) {
      return(parsed)
    }
    # Operators such as `>=` in subgroup names are not valid plotmath;
    # never leak the source tildes as literal characters.
    return(gsub("~", " ", label, fixed = TRUE))
  }
  label
}

#' Resolve display subgroup label for overall-IF plots
#' @noRd
overall_ifrac_subgroup_label <- function(ifrac_var, interactmeta = NULL) {
  if (!is.null(interactmeta) && !is.null(interactmeta$ifrac_var) && nzchar(interactmeta$ifrac_var)) {
    lab <- interactmeta_ifrac_subgroup_label(interactmeta$ifrac_var)
    if (!is.null(lab) && nzchar(lab)) {
      return(lab)
    }
  }
  interactmeta_ifrac_subgroup_label(ifrac_var)
}

#' Subgroup label from an ifrac column name (e.g. `ifrac_Older` -> `Older`)
#' @noRd
interactmeta_ifrac_subgroup_label <- function(ifrac_var) {
  if (is.null(ifrac_var) || !nzchar(ifrac_var)) {
    return(NULL)
  }
  sub("^ifrac2?_", "", ifrac_var)
}

#' Default x-axis label for overall IF plots
#'
#' Plain character so subgroup names with operators (`>=65`) are not
#' parsed as plotmath (which would leave `Information~fraction~...` on the
#' axis).
#' @noRd
overall_ifrac_default_xlab <- function(subgroup_label = NULL) {
  base <- if (is.null(subgroup_label) || !nzchar(subgroup_label)) {
    "Information fraction"
  } else {
    sprintf("Information fraction (%s)", subgroup_label)
  }
  paste0(base, " (%)")
}

#' Default x-axis label for [plot_contribution()] (percent-scale axis)
#' @noRd
bubbleplot_default_xlab <- function(ifrac_var = NULL) {
  paste0("Information fraction (IF)", " (%)")
}

#' Parse legend labels that may use plotmath strings
#' @noRd
interactmeta_parse_legend_labels <- function(labels, fallback) {
  if (is.null(labels)) {
    return(fallback)
  }
  out <- vapply(
    seq_along(labels),
    function(i) {
      lab <- labels[[i]]
      if (is.na(lab)) {
        return(fallback[[i]])
      }
      parsed <- interactmeta_as_ggplot_label(lab)
      if (inherits(parsed, "expression")) {
        return(parsed[[1]])
      }
      if (is.language(parsed)) {
        return(parsed)
      }
      parsed
    },
    FUN.VALUE = quote(x),
    USE.NAMES = FALSE
  )
  stats::setNames(as.list(out), names(labels))
}
