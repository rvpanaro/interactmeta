# ---- interactmeta_palette.R ----
#' Default discrete colours for subgroup plots
#'
#' First two are ColorBrewer Set1 red / blue for the reference /
#' non-reference pair. Remaining hues follow the colourblind-safe Okabe–Ito
#' set. This is subgroup identity only — not for heterogeneity, forests, or
#' method-comparison series.
#' @noRd
interactmeta_okabe_ito_colours <- function(n = 7L) {
  cols <- c(
    "#E41A1C", # red (reference)
    "#377EB8", # blue (non-reference)
    "#009E73", # bluish green
    "#E69F00", # orange
    "#0072B2", # blue (darker)
    "#CC79A7", # reddish purple
    "#F0E442"  # yellow
  )
  n <- max(1L, as.integer(n)[[1L]])
  if (n <= length(cols)) {
    cols[seq_len(n)]
  } else {
    grDevices::colorRampPalette(cols)(n)
  }
}

#' Resolve a paletteer discrete palette name
#'
#' @param palette Optional `"package::palette"` string, or `"okabeito"` for the
#'   built-in Okabe–Ito colours. When `NULL`, uses the `"interactmeta.palette"` option
#'   (default `"okabeito"`).
#' @return A single palette name.
#' @keywords internal
#' @noRd
interactmeta_palette_name <- function(palette = NULL) {
  if (is.null(palette)) {
    getOption("interactmeta.palette", "okabeito")
  } else {
    as.character(palette)[[1L]]
  }
}

#' @keywords internal
#' @noRd
interactmeta_palette_grey_name <- function(palette = NULL) {
  if (is.null(palette)) {
    getOption("interactmeta.palette_grey", "RColorBrewer::Greys")
  } else {
    palette
  }
}

#' @keywords internal
#' @noRd
interactmeta_normalize_hex <- function(x) {
  x <- as.character(x)
  sub("(#[0-9A-Fa-f]{6}).*", "\\1", x)
}

#' @keywords internal
#' @noRd
interactmeta_paletteer_d <- function(palette, n) {
  pal <- as.character(palette)[[1L]]
  n <- max(1L, as.integer(n)[[1L]])
  if (identical(pal, "okabeito") || identical(pal, "OkabeIto") ||
        grepl("okabeito", pal, ignore.case = TRUE)) {
    return(interactmeta_okabe_ito_colours(n))
  }
  interactmeta_normalize_hex(paletteer::paletteer_d(pal, n = n))
}

#' Colour palette for interact figures
#'
#' Named colours for figures. The default discrete scale (`"okabeito"`) uses
#' **red** for the reference subgroup and **blue** for the non-reference, then
#' colourblind-safe Okabe–Ito hues for additional levels
#' (`options(interactmeta.palette = "okabeito")`). That pair is subgroup
#' identity only. Any
#' [paletteer](https://emilhvitfeldt.github.io/paletteer/) discrete palette name
#' in `"package::palette"` form is also accepted.
#'
#' @param palette Optional palette name (`"okabeito"` or `"package::palette"`).
#'   See [paletteer::palettes_d].
#' @param n Number of colours to request from the palette (minimum 6).
#'
#' @return A named character vector of hex colours.
#' @export
#'
#' @examples
#' interactmeta_palette()
interactmeta_palette <- function(palette = NULL, n = 6L) {
  pal <- interactmeta_palette_name(palette)
  cols <- interactmeta_paletteer_d(pal, max(as.integer(n), 6L))
  # Prefer blue accents for lines (avoid green reading as a third subgroup).
  use_okabe <- identical(pal, "okabeito") || grepl("okabeito", pal, ignore.case = TRUE)
  accent <- if (use_okabe) {
    list(
      line = cols[[min(2L, length(cols))]], # blue
      deep = cols[[min(5L, length(cols))]], # darker blue
      mid = cols[[min(4L, length(cols))]]   # orange
    )
  } else {
    list(
      line = cols[[min(3L, length(cols))]],
      deep = cols[[min(4L, length(cols))]],
      mid = cols[[min(5L, length(cols))]]
    )
  }
  c(
    teal = accent$line,
    navy = accent$deep,
    teal_light = scales::alpha(accent$line, 0.35),
    mint = scales::alpha(accent$deep, 0.2),
    navy_mid = accent$mid,
    slate = "#6B7C8A",
    white = "#FFFFFF"
  )
}

#' Map subgroup labels to fill and line colours (internal)
#'
#' @param labels Character vector of subgroup labels.
#' @param reference Reference subgroup label.
#' @param palette Optional palette name; defaults to the `"interactmeta.palette"` option
#'   (`"okabeito"`).
#' @return A list with `fill` and `colour` character vectors.
#' @keywords internal
#' @name subgroup_colour_map
subgroup_colour_map <- function(labels, reference = NULL, palette = NULL) {
  levels <- unique(as.character(labels))
  levels <- levels[nzchar(levels)]
  if (length(levels) == 0L) {
    return(list(fill = character(0), colour = character(0)))
  }

  if (!is.null(reference) && length(reference) == 1L && reference %in% levels) {
    ord <- c(reference, setdiff(levels, reference))
  } else {
    ord <- levels
  }

  pal <- if (is.null(palette)) {
    interactmeta_palette_name(NULL)
  } else {
    as.character(palette)[[1L]]
  }
  n_ord <- length(ord)
  line_cols <- subgroup_palette_colours(pal, n_ord)
  # Opaque fills; ribbon / area geoms apply `interactmeta_plot_constants()$ribbon_alpha`.
  fill_cols <- line_cols

  fill <- colour <- character(n_ord)
  names(fill) <- names(colour) <- ord
  fill[] <- fill_cols
  colour[] <- line_cols

  list(
    fill = unname(fill[labels]),
    colour = unname(colour[labels])
  )
}

#' Discrete colours from a paletteer / built-in palette
#' @noRd
subgroup_palette_colours <- function(palette, n) {
  n <- max(1L, as.integer(n)[[1L]])
  is_grey <- grepl("Greys|Grays|grey|gray", palette, ignore.case = TRUE)
  if (!is_grey) {
    return(interactmeta_paletteer_d(palette, n))
  }
  # Brewer Greys starts near-white; skip light shades so lines stay readable.
  n_pool <- max(9L, n + 5L)
  pool <- interactmeta_paletteer_d(palette, n_pool)
  start <- max(2L, ceiling(n_pool * 0.4))
  usable <- pool[start:n_pool]
  idx <- unique(pmax(1L, pmin(
    length(usable),
    round(seq(1, length(usable), length.out = n))
  )))
  if (length(idx) < n) {
    idx <- seq_len(min(n, length(usable)))
  }
  usable[idx]
}

interaction_row_colours <- function(labels, yi, reference = NULL, contrast_sep = "-", palette = NULL) {
  n <- length(labels)
  if (!n) {
    return(list(fill = character(0), colour = character(0)))
  }

  labels <- as.character(labels)
  levels <- unique(labels)
  sg <- subgroup_colour_map(levels, reference = reference, palette = palette)
  fill_by_level <- stats::setNames(sg$fill, levels)
  colour_by_level <- stats::setNames(sg$colour, levels)

  out_fill <- out_colour <- rep(NA_character_, n)
  has_yi <- !is.na(yi)

  for (i in which(has_yi)) {
    lab <- labels[i]
    if (!nzchar(lab)) {
      next
    }

    non_ref <- NULL
    if (grepl(contrast_sep, lab, fixed = TRUE)) {
      parts <- strsplit(lab, contrast_sep, fixed = TRUE)[[1]]
      if (length(parts) >= 2L) {
        non_ref <- parts[[1]]
      }
    } else if (!is.null(reference) && lab != reference) {
      non_ref <- lab
    } else if (is.null(reference)) {
      non_ref <- lab
    }

    if (!length(non_ref) || !nzchar(non_ref)) {
      next
    }

    j <- match(non_ref, names(fill_by_level), nomatch = NA_integer_)
    if (!is.na(j)) {
      out_fill[i] <- fill_by_level[j]
      out_colour[i] <- colour_by_level[j]
    }
  }

  list(fill = out_fill, colour = out_colour)
}

#' Grey-scale fills matching [subgroup_colour_map()] roles (reference vs non-reference)
#' @noRd
subgroup_shrinkage_colour_map <- function(labels, reference = NULL, palette = NULL) {
  labels <- as.character(labels)
  if (!length(labels)) {
    return(list(fill = character(0), colour = character(0)))
  }
  levels <- unique(labels[nzchar(labels)])
  if (length(levels) == 0L) {
    return(list(fill = rep(NA_character_, length(labels)), colour = rep(NA_character_, length(labels))))
  }

  if (!is.null(reference) && length(reference) == 1L && reference %in% levels) {
    ord <- c(reference, setdiff(levels, reference))
  } else {
    ord <- levels
  }

  grey_pal <- interactmeta_palette_grey_name(palette)
  n_ord <- length(ord)
  n_grey <- max(9L, n_ord * 3L)
  grey <- interactmeta_paletteer_d(grey_pal, n_grey)

  fill <- colour <- character(n_ord)
  names(fill) <- names(colour) <- ord
  fill[] <- grey[seq_len(n_ord) + 1L]
  colour[] <- grey[rev(seq_len(n_ord)) + (n_grey - n_ord)]

  list(
    fill = unname(fill[labels]),
    colour = unname(colour[labels])
  )
}

#' @rdname interactmeta_palette
#' @export
subgroupmeta_palette <- interactmeta_palette

# ---- interactmeta_theme.R ----
#' Default base font size for interact ggplot figures
#'
#' @param override User-supplied `base_size`; when non-`NULL`, returned as-is.
#' @param default Fallback when `override` is `NULL` and no option is set.
#' @details
#' ggplot2 figures (contribution plots, overall-IF density) use
#' `base_size` via [interactmeta_theme_report()] / `interactmeta_ggplot_theme()`. Set a session
#' default with `options(interactmeta.base_size = 11)`.
#'
#' Forest plots use a separate default from `interactmeta_forest_base_size()`;
#' see the *Plot typography* section in [interact-package].
#' @noRd
interactmeta_base_size <- function(override = NULL, default = 11L) {
  if (!is.null(override)) {
    return(override)
  }
  getOption("interactmeta.base_size", default)
}

#' Default sans-serif family for interact ggplot figures
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
#' line width via `forestplot_base_line_size()`. Set a session default with
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

#' ggplot2 legend key size for effect plots at a given `base_size`
#' @noRd
interactmeta_effects_legend_point_size <- function(base_size = interactmeta_base_size()) {
  base_size * (3.5 / 12)
}

#' Shared visual constants for interact figures
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
#' muted chrome. Unknown approach labels (e.g. `CAMS — Bayesian`) recycle
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

#' Report / dashboard theme for interact ggplot figures
#'
#' Consistent sans-serif typography, muted greys, light horizontal gridlines
#' only, and a right-hand legend with comfortable padding. Use this (or
#' `interactmeta_ggplot_theme()`) so contribution plots share
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

#' Theme for conditional-effects / bubble plots
#' @noRd
ifrac_effects_theme <- function(base_size = interactmeta_base_size()) {
  interactmeta_theme_report(base_size = base_size, variant = "effects")
}

#' Legend guide defaults for subgroup-coloured effect plots
#' @noRd
interactmeta_effects_legend_guide <- function(
    base_size = interactmeta_base_size(),
    colour_order = NULL) {
  colour_guide <- list(
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
  )
  if (!is.null(colour_order)) {
    colour_guide$order <- colour_order
  }
  ggplot2::guides(
    colour = do.call(ggplot2::guide_legend, colour_guide),
    fill = "none"
  )
}

#' Axis expansion for effect plots
#' @noRd
interactmeta_effects_axis_expansion <- function() {
  ggplot2::expansion(mult = c(0.04, 0.04))
}

#' Default plot margins for forest plots
#' @noRd
interactmeta_forestplot_margin <- function() {
  ggplot2::margin(t = 4, r = 10, b = 6, l = 8, unit = "mm")
}
