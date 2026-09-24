#' Default base font size for interactmeta ggplot figures
#'
#' @param override User-supplied `base_size`; when non-`NULL`, returned as-is.
#' @param default Fallback when `override` is `NULL` and no option is set.
#' @details
#' ggplot2 figures (contribution plots, overall-IF density) use
#' `base_size` via [interactmeta_theme_report()] / [interactmeta_ggplot_theme()]. Set a session
#' default with `options(interactmeta.base_size = 11)`.
#'
#' Forest plots use a separate default from [interactmeta_forest_base_size()];
#' see the *Plot typography* section in [interactmeta-package].
#' @noRd
interactmeta_base_size <- function(override = NULL, default = 11L) {
  if (!is.null(override)) {
    return(override)
  }
  getOption("interactmeta.base_size", default)
}

#' Default sans-serif family for interactmeta ggplot figures
#' @noRd
interactmeta_base_family <- function(override = NULL) {
  if (!is.null(override) && nzchar(as.character(override)[[1L]])) {
    return(as.character(override)[[1L]])
  }
  getOption("interactmeta.base_family", "Helvetica")
}

#' Default base font size for forest plots
#'
#' @param override User-supplied `base_size`; when non-`NULL`, returned as-is.
#' @param default Fallback when `override` is `NULL` and no option is set.
#' @details
#' Default `base_size` for forest plots. Also drives CI / diamond
#' line width via [forestplot_base_line_size()]. Set a session default with
#' `options(interactmeta.forest_base_size = 9)`.
#' @noRd
interactmeta_forest_base_size <- function(override = NULL, default = 9L) {
  if (!is.null(override)) {
    return(override)
  }
  getOption("interactmeta.forest_base_size", default)
}

#' Default study-row point size for forest plots
#' @noRd
interactmeta_forest_pointsize <- function(default = 7) {
  getOption("interactmeta.forest_pointsize", default)
}

#' Default point stroke for forest plots
#' @noRd
interactmeta_forest_stroke <- function(default = 4L) {
  getOption("interactmeta.forest_stroke", default)
}

#' ggplot2 point size for effect-plot markers at a given `base_size`
#' @noRd
interactmeta_effects_point_size <- function(base_size = interactmeta_base_size()) {
  base_size * (3 / 12)
}

#' ggplot2 legend key size for effect plots at a given `base_size`
#' @noRd
interactmeta_effects_legend_point_size <- function(base_size = interactmeta_base_size()) {
  base_size * (3.5 / 12)
}

#' Shared visual constants for interactmeta figures
#'
#' Subgroup identity (reference / non-reference) stays in
#' [interactmeta_palette()] / [subgroup_colour_map()]. Everything else —
#' chrome, density, forests, method-comparison accents — should read from
#' here so plots do not invent one-off hex or line weights.
#' @noRd
interactmeta_plot_constants <- function() {
  pal <- interactmeta_palette()
  shrinkage_rgb <- grDevices::col2rgb("gray45")
  list(
    ribbon_alpha = 0.18,
    density_alpha = 0.55,
    line_hair = 0.35,
    line_mid = 0.7,
    line_strong = 1,
    segment_colour = "#D8D8D8",
    overall_line_colour = "#8A8A8A",
    na_colour = "#B0B0B0",
    text_muted = "#4A4A4A",
    text_caption = "#7A7A7A",
    strip_fill = "#F7F7F7",
    strip_border = "#E2E2E2",
    axis_line = "#2A2A2A",
    grid_major = "#E8E8E8",
    zero_line = "#D0D0D0",
    null_line = "#9A9A9A",
    density_fill = unname(pal[["slate"]]),
    density_point = "#3D4A54",
    forest_study = "#000000",
    forest_shrinkage = grDevices::rgb(
      shrinkage_rgb[1L, 1L],
      shrinkage_rgb[2L, 1L],
      shrinkage_rgb[3L, 1L],
      maxColorValue = 255
    ),
    forest_summary = unname(pal[["navy"]]),
    forest_text = "#000000",
    forest_rule = "#BEBEBE"
  )
}

#' Discrete colours for DA / SWADA / CAMS comparison plots
#'
#' Not subgroup identity: methods use palette accents; published rows use
#' muted chrome. Unknown approach labels (e.g. `CAMS — bayesmeta`) recycle
#' the same accents rather than a hardcoded approach name.
#' @noRd
interactmeta_approach_colours <- function(palette = NULL, approaches = NULL) {
  pal <- interactmeta_palette(palette = palette)
  const <- interactmeta_plot_constants()
  named <- c(
    Published = const$text_muted,
    DA = unname(pal[["navy"]]),
    AD = unname(pal[["teal"]]),
    SWADA = unname(pal[["teal"]]),
    CAMS = unname(pal[["navy_mid"]])
  )
  if (is.null(approaches)) {
    return(named)
  }
  approaches <- unique(as.character(approaches))
  approaches <- approaches[nzchar(approaches) & !is.na(approaches)]
  if (!length(approaches)) {
    return(named)
  }
  accents <- c(
    unname(pal[["navy"]]),
    unname(pal[["teal"]]),
    unname(pal[["navy_mid"]]),
    unname(pal[["slate"]])
  )
  out <- character(length(approaches))
  names(out) <- approaches
  extra <- 0L
  for (i in seq_along(approaches)) {
    nm <- approaches[[i]]
    if (nm %in% names(named)) {
      out[[i]] <- named[[nm]]
    } else {
      extra <- extra + 1L
      out[[i]] <- accents[[((extra - 1L) %% length(accents)) + 1L]]
    }
  }
  out
}

#' Report / dashboard theme for interactmeta ggplot figures
#'
#' Consistent sans-serif typography, muted greys, light horizontal gridlines
#' only, and a right-hand legend with comfortable padding. Use this (or
#' [interactmeta_ggplot_theme()]) so contribution plots share
#' one look across a report.
#'
#' Set session defaults with `options(interactmeta.base_size = 11)` and
#' `options(interactmeta.base_family = "Helvetica")`.
#'
#' @param base_size Base font size in points (axis titles ~ same size; title
#'   ~1.25×; caption ~0.95×). Default `11`.
#' @param base_family Font family (default Helvetica via
#'   `options(interactmeta.base_family)`).
#' @param variant `"effects"` (contribution plot: horizontal grid only),
#'   `"density"` (heterogeneity / IF), or `"minimal"`.
#'
#' @return A **ggplot2** theme object.
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
#'   geom_point() +
#'   interactmeta_theme_report()
#' }
interactmeta_theme_report <- function(base_size = interactmeta_base_size(default = 11L),
                              base_family = NULL,
                              variant = c("effects", "density", "minimal")) {
  variant <- match.arg(variant)
  base_family <- interactmeta_base_family(base_family)
  const <- interactmeta_plot_constants()

  title_size <- max(14, min(16, base_size * 1.25))
  axis_title_size <- max(11, min(12, base_size))
  # Footnotes / captions: closer to axis text than the old ~0.8× shrink.
  caption_size <- max(10, min(12, base_size * 0.95))

  theme <- ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, colour = const$axis_line),
      axis.text = ggplot2::element_text(
        size = base_size * 0.9,
        colour = const$axis_line,
        family = base_family
      ),
      axis.title = ggplot2::element_text(
        size = axis_title_size,
        face = "plain",
        colour = const$axis_line,
        family = base_family
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 10, unit = "pt")
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 10, unit = "pt")
      ),
      axis.line = ggplot2::element_line(colour = const$axis_line, linewidth = const$line_hair),
      axis.ticks = ggplot2::element_line(colour = const$axis_line, linewidth = 0.3),
      axis.ticks.length = grid::unit(0.12, "cm"),
      plot.title = ggplot2::element_text(
        face = "bold",
        colour = const$axis_line,
        family = base_family,
        hjust = 0,
        size = title_size,
        margin = ggplot2::margin(b = 6, unit = "pt")
      ),
      plot.subtitle = ggplot2::element_text(
        colour = const$text_muted,
        family = base_family,
        size = base_size * 0.9,
        hjust = 0,
        margin = ggplot2::margin(b = 8, unit = "pt")
      ),
      plot.caption = ggplot2::element_text(
        size = caption_size,
        colour = const$forest_text,
        family = base_family,
        hjust = 0,
        margin = ggplot2::margin(t = 10, unit = "pt")
      ),
      plot.margin = ggplot2::margin(12, 14, 10, 12, unit = "pt"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "top",
      legend.title = ggplot2::element_text(
        size = base_size * 0.85,
        face = "bold",
        colour = const$axis_line,
        family = base_family
      ),
      legend.text = ggplot2::element_text(
        size = base_size * 0.9,
        colour = const$axis_line,
        family = base_family
      ),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(t = 0, r = 4, b = 4, l = 8, unit = "pt"),
      legend.box.margin = ggplot2::margin(l = 8, unit = "pt"),
      legend.spacing.x = grid::unit(6, "pt"),
      legend.spacing.y = grid::unit(4, "pt")
    )

  if (variant == "effects") {
    theme <- theme +
      ggplot2::theme(
        # Horizontal guides only — no vertical grid competing with IF markers.
        panel.grid.major.y = ggplot2::element_line(
          colour = const$grid_major,
          linewidth = 0.25
        ),
        panel.grid.major.x = ggplot2::element_blank(),
        legend.key.size = grid::unit(base_size * 0.55, "pt"),
        legend.key.height = grid::unit(base_size * 0.55, "pt"),
        legend.key.width = grid::unit(base_size * 1.0, "pt"),
        strip.background = ggplot2::element_rect(
          fill = const$strip_fill,
          colour = const$strip_border,
          linewidth = 0.3
        ),
        strip.text = ggplot2::element_text(
          face = "bold",
          size = base_size * 0.9,
          colour = const$axis_line,
          family = base_family,
          margin = ggplot2::margin(t = 4, b = 4, unit = "pt")
        ),
        panel.spacing = grid::unit(1.0, "lines")
      )
  }

  if (variant == "density") {
    theme <- theme +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          margin = ggplot2::margin(t = 3, b = 0, unit = "pt")
        ),
        axis.title.x = ggplot2::element_text(
          margin = ggplot2::margin(t = 6, unit = "pt")
        ),
        legend.position = "none"
      )
  }

  theme
}

#' @rdname interactmeta_theme_report
#' @export
subgroupmeta_theme_report <- interactmeta_theme_report

#' Unified **ggplot2** theme for interactmeta figures
#'
#' Alias of [interactmeta_theme_report()] kept for internal call sites.
#' @noRd
interactmeta_ggplot_theme <- function(base_size = interactmeta_base_size(),
                              base_family = NULL,
                              variant = c("effects", "density", "minimal")) {
  interactmeta_theme_report(
    base_size = base_size,
    base_family = base_family,
    variant = variant
  )
}

#' Theme for conditional-effects / bubble plots
#' @noRd
ifrac_effects_theme <- function(base_size = interactmeta_base_size()) {
  interactmeta_theme_report(base_size = base_size, variant = "effects")
}

#' Legend guide defaults for subgroup-coloured effect plots
#' @noRd
interactmeta_effects_legend_guide <- function(base_size = interactmeta_base_size()) {
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      title.position = "top",
      title.hjust = 0,
      ncol = 1L,
      override.aes = list(
        size = interactmeta_effects_legend_point_size(base_size),
        alpha = 1,
        stroke = 0.55,
        fill = NA,
        linewidth = 0,
        linetype = 0
      )
    ),
    fill = "none"
  )
}

#' Axis expansion for effect plots
#' @noRd
interactmeta_effects_axis_expansion <- function() {
  ggplot2::expansion(mult = c(0.04, 0.04))
}

#' Add a null-effect reference line when the y-axis spans zero
#' @noRd
interactmeta_add_zero_line <- function(p, ylim = NULL, const = NULL) {
  if (is.null(const)) {
    const <- interactmeta_plot_constants()
  }
  if (is.null(ylim)) {
    built <- ggplot2::ggplot_build(p)
    ylim <- built$layout$panel_params[[1]]$ylim
  }
  if (length(ylim) == 2L && is.finite(ylim[[1L]]) && is.finite(ylim[[2L]]) &&
      ylim[[1L]] < 0 && ylim[[2L]] > 0) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        colour = const$zero_line,
        linewidth = const$line_hair,
        linetype = "solid",
        inherit.aes = FALSE
      )
  }
  p
}

#' Final polish layer for contribution effect plots
#' @noRd
interactmeta_finish_effects_plot <- function(p, limits = NULL) {
  p <- p + interactmeta_effects_legend_guide()
  ylim <- if (!is.null(limits)) limits$ylim else NULL
  interactmeta_add_zero_line(p, ylim = ylim)
}

#' Default plot margins for forest plots
#' @noRd
interactmeta_forestplot_margin <- function() {
  ggplot2::margin(t = 4, r = 10, b = 6, l = 8, unit = "mm")
}
