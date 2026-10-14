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
#' Pointwise `estimate` / `lower` / `upper` use shortest intervals from the
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
#' @param max_draws Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports; omit it.
#' @param exponentiate If `TRUE`, summarize on the exponential scale.
#'
#' @return A data frame with one row per grid point and subgroup:
#'   `ifrac`, `subgroup`, `estimate`, `lower`, `upper`, plus
#'   `overall_estimate`, `overall_lower`, `overall_upper` for
#'   \eqn{M(\pi)=\alpha+\theta\pi} at that grid point.
#'
#' @seealso [subgroups()], [combined()], [plot_contribution()]
#' @name predict_cams_by_if
#' @export
predict_cams_by_if <- function(fit,
                               if_grid = NULL,
                               level = 0.95,
                               ifrac_var = NULL,
                               resolution = 200,
                               max_draws = NULL,
                               exponentiate = FALSE) {
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an `interactmeta` from `cams()`.", call. = FALSE)
  }
  level <- as.numeric(level)[[1L]]
  if (!is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a coverage probability in (0, 1).", call. = FALSE)
  }

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
  post <- comparison_posterior_df(fit)
  idx <- comparison_marginal_draw_index(
    nrow(post),
    max_draws = max_draws,
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
      fit = fit
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
        fit = fit
      )
      data.frame(
        ifrac = pi,
        subgroup = lbl,
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
  if (identical(fit$engine, "bayesmeta")) {
    attr(out, "interval.rule") <- "mixture_shortest"
  }
  class(out) <- c("interactmeta_cam_if_pred", "interactmeta_cams_if_pred", "data.frame")
  out
}

#' @rdname predict_cams_by_if
#' @export
predict_cam_by_if <- predict_cams_by_if

#' Pointwise CAMS interval at one information fraction.
#'
#' Public rule: shortest intervals of the DIRECT-node mixture (weights and
#' within-node SDs). Alpha and theta stay jointly paired.
#' @noRd
predict_cam_if_qi <- function(x, post, a_alpha, a_theta, a_gamma, level, exponentiate,
                               fit = NULL) {
  if (!is.null(fit) && identical(fit$engine, "metafor")) {
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
    exponentiate = exponentiate
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
#' @param ifrac_var,prob,type,max_draws,contrast_sep Passed to
#'   [predict_cams_by_if()] and [ifrac()].
#' @param exponentiate If `TRUE`, draw the CAMS curves, study points, and
#'   intervals on the exponential (ratio) scale. Observed `yi` are exponentiated
#'   after the log-scale Wald interval is formed.
#' @param at_ifrac Optional dotted reference IF; default `NULL` (no marker).
#' @param oi_draws Optional overall-IF draws (see [ifrac()]).
#' @param effects,re_formula,resolution Passed to the conditional-effects helper
#'   (`re_formula = NA` by default, population effects only).
#' @param subgroup_labels,legend_math,colors,fills,ylab Plot styling.
#' @param xlab Axis label; default is the information-fraction label.
#' @param xlim,ylim,x_breaks,y_breaks Axis limits and breaks.
#' @param show_overall_line Logical; dashed \eqn{M(\pi)=\alpha+\theta\pi} line.
#' @param point_size_range Range for precision-scaled study points.
#' @param show_footnote If `TRUE` (default), add a caption below the plot with
#'   the overall information fraction, its interval, study count, and reference.
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
#' plot_contribution(fit)
#' }
#' @export
plot_contribution <- function(fit,
                                    ifrac_var = NULL,
                                    prob = 0.95,
                                    type = "Overall",
                                    max_draws = NULL,
                                    contrast_sep = "-",
                                    exponentiate = FALSE,
                                    subgroup_var = NULL,
                                    at_ifrac = NULL,
                                    oi_draws = NULL,
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
    max_draws = max_draws,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate,
    subgroup_var = subgroup_var,
    at_ifrac = at_ifrac,
    oi_draws = oi_draws,
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
                                          max_draws = NULL,
                                          contrast_sep = "-",
                                          exponentiate = FALSE,
                                          subgroup_var = NULL,
                                          at_ifrac = NULL,
                                          oi_draws = NULL,
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
                                          point_size_range = c(2.5, 6.5),
                                          show_footnote = TRUE,
                                          digits = 2,
                                          base_size = interactmeta_base_size(default = 11L)) {
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

  pred <- predict_cam_by_if(
    fit,
    if_grid = NULL,
    level = prob,
    ifrac_var = ifrac_var,
    resolution = resolution,
    max_draws = max_draws,
    exponentiate = exponentiate
  )
  ce_df <- as.data.frame(pred)
  ce_df$ifrac_x <- ce_df$ifrac
  ce_df$estimate__ <- ce_df$estimate
  ce_df$lower__ <- ce_df$lower
  ce_df$upper__ <- ce_df$upper
  ce_df[[ifrac_var]] <- ce_df$ifrac
  ce_df[[subgroup_var]] <- ce_df$subgroup
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

  if (!is.null(at_ifrac) && length(at_ifrac) == 1L && is.finite(at_ifrac)) {
    p_out <- p_out +
      ggplot2::geom_vline(
        xintercept = at_ifrac,
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

  if (isTRUE(show_overall_line) && nrow(overall_df)) {
    p_out <- p_out +
      ggplot2::geom_line(
        data = overall_df,
        ggplot2::aes(x = .data$ifrac_x, y = .data$estimate__),
        colour = const$overall_line_colour,
        linewidth = const$line_mid,
        linetype = "22",
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
      expand = interactmeta_effects_axis_expansion()
    ) +
    ggplot2::scale_y_continuous(
      breaks = interactmeta_continuous_breaks(limits$y_breaks),
      expand = interactmeta_effects_axis_expansion()
    ) +
    ifrac_effects_theme(base_size = base_size)

  if (isTRUE(show_footnote)) {
    oi_obj <- oi_draws
    if (is.null(oi_obj)) {
      oi_obj <- ifrac(fit, ifrac_var = ifrac_var, prob = prob, summary = FALSE)
    }
    bundle <- overall_ifrac_bundle(oi_obj, prob = prob)
    p_out <- p_out +
      ggplot2::labs(caption = overall_ifrac_footnote_bundle(bundle, digits = digits))
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

  p_out + interactmeta_effects_legend_guide(base_size = base_size)
}
