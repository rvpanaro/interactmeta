#' Recommended figure dimensions for a forest plot
#'
#' Returns width and height in inches scaled to the number of layout rows
#' (studies, pooled / prediction diamonds, subgroup headers), plus footnote
#' lines and axis / favours chrome. Full forests use about **0.28 inches per
#' content row**, **0.22 inches per footnote line**, and ~2.2 inches of chrome,
#' capped at 24 inches.
#'
#' Use the returned sizes with `png()` / `pdf()` or **knitr** chunk options
#' (`fig.width`, `fig.height`). Compute dimensions *before* opening the device
#' or in a prior chunk, then pass them into the export / chunk options — the
#' active device height is what prevents row overlap / cut-off.
#'
#' @param x A `forestplot()` result.
#' @param width Target figure width in inches. When `NULL`, uses a size derived
#'   from the number of rows and panels.
#'
#' @return A list with elements `width`, `height`, `aspect` (height / width),
#'   `n_rows` (content / layout rows), `n_footnotes`, `n_chrome_in` (axis +
#'   favours band in inches), and `formula` (character description of the
#'   height calculation).
#'
#' @examples
#' \dontrun{
#' fp <- forestplot(fit)
#' dims <- forestplot_height(fp, width = 10)
#' png("forestplot.png", width = dims$width, height = dims$height,
#'     units = "in", res = 300)
#' forestplot_draw(fp)
#' dev.off()
#' }
#'
#' @seealso [forestplot_draw()], [forestplot()]
#' @export
forestplot_height <- function(x, width = 9.5) {
  dims <- forestplot_figure_inches(x)
  if (is.null(width)) {
    width <- dims$width
  }
  if (!is.finite(width) || width <= 0) {
    stop("`width` must be a positive finite number.", call. = FALSE)
  }
  # Height is driven by layout rows + footnotes, not by aspect with width.
  height <- dims$height
  height <- max(3, min(24, height))
  list(
    width = width,
    height = height,
    aspect = height / width,
    n_rows = dims$n_rows,
    n_footnotes = dims$n_footnotes %||% 0L,
    n_chrome_in = dims$n_chrome_in %||% NA_real_,
    formula = dims$formula %||% NA_character_
  )
}

#' Draw a forest plot on the active graphics device
#'
#' Forest plots are drawn with [metafor::forest()] / [forestplot::forestplot()].
#' When `typography_auto = TRUE` was requested, sizes may be re-fit if the
#' device has changed.
#'
#' @param x A `forestplot()` result.
#'
#' @return `x`, invisibly.
#'
#' @seealso [forestplot_height()], [forestplot()]
#' @export
forestplot_draw <- function(x) {
  if (forestplot_is_metafor_engine(x)) {
    x <- forestplot_refit_typography(x)
    if (forestplot_has_graphics_device()) {
      grDevices::dev.hold()
      on.exit(grDevices::dev.flush(), add = TRUE)
    }
    forestplot_metafor_draw(x)
    return(invisible(x))
  }

  # Legacy ggplot / grid objects (if any)
  x <- forestplot_refit_typography(x)
  grob <- forestplot_as_grob(x)
  dims <- forestplot_grob_inches(grob)
  if (!all(is.finite(unlist(dims))) || dims$width <= 0 || dims$height <= 0) {
    stop("Forest plot has non-finite grob dimensions.", call. = FALSE)
  }

  grid::grid.newpage()
  dev_size <- grDevices::dev.size("in")
  scale <- min(dev_size[1] / dims$width, dev_size[2] / dims$height)
  w_npc <- (dims$width * scale) / dev_size[1]
  h_npc <- (dims$height * scale) / dev_size[2]

  x_left <- (1 - w_npc) / 2
  y_bottom <- (1 - h_npc) / 2

  vp <- grid::viewport(
    x = x_left,
    y = y_bottom,
    width = w_npc,
    height = h_npc,
    just = c("left", "bottom")
  )
  grid::pushViewport(vp)
  on.exit(grid::popViewport(), add = TRUE)
  grid::grid.draw(grob)
  invisible(x)
}

#' @export
print.interactmeta_forestplot <- function(x, ...) {
  forestplot_maybe_print(x)
  invisible(x)
}

#' @noRd
forestplot_is_fitted <- function(fp) {
  !is.null(attr(fp, "forestplot_fit_render", exact = TRUE))
}

#' @noRd
forestplot_as_grob <- function(x) {
  if (forestplot_is_metafor_engine(x)) {
    stop(
      "metafor-backed forest plots are drawn with `forestplot_draw()` / ",
      "`print()` (base graphics), not as grid grobs.",
      call. = FALSE
    )
  }
  if (inherits(x, "grob")) {
    return(x)
  }
  if (inherits(x, "ggplot")) {
    return(ggplot2::ggplotGrob(x))
  }
  if (inherits(x, "patchwork")) {
    return(patchwork::patchworkGrob(x))
  }
  if (is.list(x) && !is.null(x$figure) && inherits(x$figure, "grob")) {
    return(x$figure)
  }
  if (is.list(x) && !is.null(x$plot) && inherits(x$plot, "ggplot")) {
    gt <- ggplot2::ggplotGrob(x$plot)
    return(forestplot_set_panel_height(gt, attr(x, "panel.height", exact = TRUE)))
  }
  stop(
    "`x` must be a forest plot object, ggplot, patchwork, or grid grob.",
    call. = FALSE
  )
}

#' Figure size in inches for forestplot_height()
#' @noRd
forestplot_figure_inches <- function(x) {
  if (forestplot_is_metafor_engine(x)) {
    n_rows <- forestplot_n_layout_rows(x)
    n_panels <- forestplot_n_panels(x)
    n_footnotes <- forestplot_n_footnote_lines(x)
    has_favours <- isTRUE(attr(x, "has_favours", exact = TRUE)) ||
      (!is.null(x$specs$favours) && length(x$specs$favours) >= 2L)
    # Wide enough for study labels + forest + Effect text; leave room so
    # wrapped footnotes are not clipped on the right.
    width <- if (n_panels <= 1L) {
      12.5
    } else {
      5.8 * n_panels + 1.6
    }
    # Content rows (studies + diamonds + subgroup headers + table header):
    #   0.28 in each. Footnote lines sit in a tight band below the axis: 0.14 in
    #   each. Chrome covers axis ticks, xlab, and favours (~2.0–2.3 in).
    row_pitch <- 0.28
    foot_pitch <- 0.14
    # Reserve favours chrome even when favours are off so recommended height
    # stays stable across `favours` / `type` toggles.
    force(has_favours)
    chrome <- 2.3
    # +1 row of slack so the last Overall / Prediction diamond is not clipped.
    height <- row_pitch * (n_rows + 1L) + foot_pitch * n_footnotes + chrome
    height <- max(3.8, min(28, height))
    formula <- sprintf(
      "height = 0.28*(%d+1) + 0.14*%d + %.1f = %.2f",
      n_rows, n_footnotes, chrome, height
    )
    # Forests are naturally wider than tall. Never inflate height to meet an
    # aspect floor: that stretches row pitch when the caller then uses a
    # narrower `fig.width`. Only widen a skyscraper (height / width > 2.8).
    aspect <- height / width
    if (aspect > 2.8) {
      width <- height / 2.8
    }
    return(list(
      width = width,
      height = height,
      n_rows = n_rows,
      n_footnotes = n_footnotes,
      n_chrome_in = chrome,
      formula = formula
    ))
  }
  out <- forestplot_grob_inches(forestplot_as_grob(x))
  out$n_rows <- NA_integer_
  out$n_footnotes <- 0L
  out$n_chrome_in <- NA_real_
  out$formula <- NA_character_
  out
}

#' Number of wrapped footnote lines for height sizing
#' @noRd
forestplot_n_footnote_lines <- function(x, wrap_at = 96L) {
  specs <- x$specs
  if (is.null(specs)) {
    return(0L)
  }
  length(forestplot_bayesmeta_footnote_lines(specs, wrap_at = wrap_at))
}

#' Number of forest panels on a forestplot object
#' @noRd
forestplot_n_panels <- function(x) {
  panels <- attr(x, "panels", exact = TRUE)
  if (is.list(panels) && length(panels)) {
    return(length(panels))
  }
  if (!is.null(x$specs$panels)) {
    return(length(x$specs$panels))
  }
  2L
}

#' Layout row count for height sizing (panel rows + subgroup headers / gaps)
#'
#' Counts panel data rows plus unique non-empty headings (and a table header
#' row). Prefer this over the legacy metafor ylim span so bayesmeta forests
#' size correctly when Overall includes both Pooled and Prediction.
#'
#' @noRd
forestplot_n_layout_rows <- function(x) {
  specs <- x$specs
  # Prefer draw panels from specs when multiple axis panels are present.
  draw_panels <- specs$panels
  if (is.list(draw_panels) && length(draw_panels) && is.data.frame(draw_panels[[1L]]) &&
        !is.null(specs$panel_roles) &&
        length(specs$panel_roles) == length(draw_panels) &&
        length(draw_panels) > 1L) {
    rl <- specs$row_labels
    total <- 0L
    for (i in seq_along(draw_panels)) {
      rl_i <- if (is.list(rl) && !is.data.frame(rl)) rl[[i]] else rl
      total <- total + forestplot_n_layout_rows_one(draw_panels[[i]], rl_i)
    }
    return(max(1L, total))
  }

  panels <- attr(x, "panels", exact = TRUE)
  if (is.null(panels) || !length(panels)) {
    panels <- draw_panels
  }
  if (!is.list(panels) || !length(panels) || !is.data.frame(panels[[1L]])) {
    return(10L)
  }
  panel <- if (!is.null(draw_panels) && length(draw_panels) && is.data.frame(draw_panels[[1L]])) {
    draw_panels[[1L]]
  } else {
    panels[[1L]]
  }
  row_labels <- specs$row_labels
  if (is.list(row_labels) && !is.data.frame(row_labels)) {
    row_labels <- row_labels[[1L]]
  }
  if (is.null(row_labels) && !is.null(attr(x, "row_labels", exact = TRUE))) {
    row_labels <- attr(x, "row_labels", exact = TRUE)
  }
  forestplot_n_layout_rows_one(panel, row_labels)
}

#' @noRd
forestplot_n_layout_rows_one <- function(panel, row_labels) {
  n <- max(1L, nrow(panel))
  n_headers <- 0L
  if (is.data.frame(row_labels) && "heading" %in% names(row_labels)) {
    h <- as.character(row_labels$heading)
    h <- unique(h[nzchar(h) & !is.na(h)])
    n_headers <- length(h)
  }
  # +1 table header row
  as.integer(n + n_headers + 1L)
}

#' Force absolute height of ggplot panel cell(s) in a gtable
#' @noRd
forestplot_set_panel_height <- function(gt, panel.height) {
  if (is.null(panel.height) || !inherits(gt, "gtable")) {
    return(gt)
  }
  h <- if (inherits(panel.height, "unit")) {
    panel.height
  } else if (is.numeric(panel.height) && length(panel.height) == 1L && is.finite(panel.height)) {
    grid::unit(panel.height, "mm")
  } else {
    return(gt)
  }
  panel_rows <- unique(gt$layout$t[grepl("^panel", gt$layout$name)])
  for (r in panel_rows) {
    gt$heights[r] <- h
  }
  gt
}

#' Force absolute width of ggplot panel cell(s) in a gtable
#' @noRd
forestplot_set_panel_width <- function(gt, panel.width) {
  if (is.null(panel.width) || !inherits(gt, "gtable")) {
    return(gt)
  }
  w <- if (inherits(panel.width, "unit")) {
    panel.width
  } else if (is.numeric(panel.width) && length(panel.width) == 1L && is.finite(panel.width)) {
    grid::unit(panel.width, "mm")
  } else {
    return(gt)
  }
  panel_cols <- unique(gt$layout$l[grepl("^panel", gt$layout$name)])
  for (c in panel_cols) {
    gt$widths[c] <- w
  }
  gt
}

#' @noRd
forestplot_grob_inches <- function(x) {
  grob <- if (inherits(x, "grob")) x else forestplot_as_grob(x)
  list(
    width = grid::convertWidth(grid::grobWidth(grob), "in", valueOnly = TRUE),
    height = grid::convertHeight(grid::grobHeight(grob), "in", valueOnly = TRUE)
  )
}

#' Active (or fallback) device size in inches
#' @noRd
forestplot_device_inches <- function() {
  if (forestplot_has_graphics_device()) {
    return(as.numeric(grDevices::dev.size("in")))
  }
  fallback <- getOption("interactmeta.forest_device_inches", c(12, 8.5))
  fallback <- as.numeric(fallback)
  if (length(fallback) != 2L || any(!is.finite(fallback)) || any(fallback <= 0)) {
    return(c(12, 8.5))
  }
  fallback
}

#' Typography and panel height tuned to display size
#'
#' Chooses `base_size`, `pointsize`, and `panel.height` so row labels stay
#' readable on the target device (about 7–12 pt after draw scaling).
#' @noRd
forestplot_display_metrics <- function(n_rows,
                                       device_width_in = NULL,
                                       device_height_in = NULL,
                                       base_size = NULL,
                                       pointsize = NULL,
                                       panel.height = NULL,
                                       has_favours = FALSE,
                                       show_ci_text = FALSE,
                                       n_forest_panels = 2L) {
  n_rows <- max(1L, as.integer(n_rows)[1])
  dev <- forestplot_device_inches()
  if (is.null(device_width_in)) {
    device_width_in <- dev[[1]]
  }
  if (is.null(device_height_in)) {
    device_height_in <- dev[[2]]
  }
  device_width_in <- as.numeric(device_width_in)[1]
  device_height_in <- as.numeric(device_height_in)[1]
  if (!is.finite(device_width_in) || device_width_in <= 0) {
    device_width_in <- 8
  }
  if (!is.finite(device_height_in) || device_height_in <= 0) {
    device_height_in <- 6
  }

  # Chrome: title, x-axis, margins, optional favours+arrows strip
  chrome_in <- 0.95 + if (isTRUE(has_favours)) 0.75 else 0
  usable_h_mm <- max(28, (device_height_in - chrome_in) * 25.4)

  row_mm <- usable_h_mm / n_rows
  # Keep rows tall enough for readable type
  row_mm <- max(3.8, min(8.5, row_mm))

  if (is.null(panel.height)) {
    panel_mm <- min(360, max(36, row_mm * n_rows))
    panel.height <- grid::unit(panel_mm, "mm")
  } else if (is.numeric(panel.height) && length(panel.height) == 1L) {
    panel_mm <- as.numeric(panel.height)
    row_mm <- panel_mm / n_rows
    panel.height <- grid::unit(panel_mm, "mm")
  } else if (inherits(panel.height, "unit")) {
    panel_mm <- grid::convertHeight(panel.height, "mm", valueOnly = TRUE)
    if (is.finite(panel_mm) && panel_mm > 0) {
      row_mm <- panel_mm / n_rows
    }
  }

  # ~55% of row pitch in points
  row_pt <- row_mm * (72 / 25.4)
  auto_base <- max(7, min(12, round(row_pt * 0.55, 1)))

  # Narrow multi-column layouts (CI text) need slightly smaller type
  n_forest_panels <- max(1L, as.integer(n_forest_panels)[1])
  col_units <- 1 + n_forest_panels * (1 + isTRUE(show_ci_text))
  min_width <- 1.55 * col_units
  if (device_width_in < min_width) {
    auto_base <- max(6.5, round(auto_base * device_width_in / min_width, 1))
  }

  if (is.null(base_size)) {
    base_size <- auto_base
  }
  base_size <- as.numeric(base_size)[1]
  if (!is.finite(base_size) || base_size <= 0) {
    base_size <- 9
  }

  if (is.null(pointsize)) {
    pointsize <- max(4.5, min(8, base_size * 0.75))
  }
  pointsize <- as.numeric(pointsize)[1]
  if (!is.finite(pointsize) || pointsize <= 0) {
    pointsize <- max(4.5, base_size * 0.75)
  }

  list(
    base_size = base_size,
    pointsize = pointsize,
    panel.height = panel.height,
    row_mm = row_mm,
    device_inches = c(device_width_in, device_height_in)
  )
}

#' @noRd
forestplot_default_panel_height <- function(n_rows) {
  metrics <- forestplot_display_metrics(n_rows = n_rows)
  metrics$panel.height
}

#' @noRd
forestplot_resolve_panel_height <- function(panels, panel.height) {
  if (!is.null(panel.height)) {
    return(panel.height)
  }
  n_rows <- max(vapply(panels, nrow, integer(1L)))
  forestplot_default_panel_height(n_rows)
}

#' Apply display metrics inside forestplot_render()
#' @noRd
forestplot_apply_display_metrics <- function(panels,
                                             base_size,
                                             pointsize,
                                             panel.height,
                                             typography_auto = FALSE,
                                             has_favours = FALSE,
                                             show_ci_text = FALSE) {
  n_rows <- max(1L, max(vapply(panels, nrow, integer(1L))))
  metrics <- forestplot_display_metrics(
    n_rows = n_rows,
    device_width_in = NULL,
    device_height_in = NULL,
    base_size = if (isTRUE(typography_auto)) NULL else base_size,
    pointsize = if (isTRUE(typography_auto)) NULL else pointsize,
    panel.height = panel.height,
    has_favours = has_favours,
    show_ci_text = show_ci_text,
    n_forest_panels = length(panels)
  )
  list(
    base_size = metrics$base_size,
    pointsize = metrics$pointsize,
    panel.height = metrics$panel.height,
    device_inches = metrics$device_inches
  )
}

#' Re-render forest typography when the active device size changes
#' @noRd
forestplot_refit_typography <- function(x) {
  if (!inherits(x, "interactmeta_forestplot")) {
    return(x)
  }
  if (!isTRUE(attr(x, "typography_auto", exact = TRUE))) {
    return(x)
  }
  if (!forestplot_has_graphics_device()) {
    return(x)
  }

  panels <- attr(x, "panels", exact = TRUE)
  n_rows <- forestplot_n_layout_rows(x)
  n_panels <- forestplot_n_panels(x)
  dev <- as.numeric(grDevices::dev.size("in"))
  prev <- attr(x, "display_inches", exact = TRUE)
  if (
    is.numeric(prev) && length(prev) == 2L &&
      abs(prev[[1]] - dev[[1]]) < 0.2 &&
      abs(prev[[2]] - dev[[2]]) < 0.2
  ) {
    return(x)
  }

  has_favours <- isTRUE(attr(x, "has_favours", exact = TRUE))
  show_ci_text <- isTRUE(attr(x, "show_ci_text", exact = TRUE))

  metrics <- forestplot_display_metrics(
    n_rows = n_rows,
    device_width_in = dev[[1]],
    device_height_in = dev[[2]],
    has_favours = has_favours,
    show_ci_text = show_ci_text,
    n_forest_panels = n_panels
  )

  if (forestplot_is_fitted(x)) {
    render_args <- attr(x, "forestplot_fit_render", exact = TRUE)
    opts <- attr(x, "forestplot_fit_options", exact = TRUE)
    if (is.null(render_args) || is.null(opts)) {
      return(x)
    }
    dots <- opts$dots %||% list()
    override <- list(
      prob = opts$prob,
      contrast_sep = opts$contrast_sep,
      exponentiate = opts$exponentiate %||% FALSE,
      base_size = metrics$base_size,
      pointsize = metrics$pointsize,
      panel.height = metrics$panel.height,
      typography_auto = TRUE,
      annotations = attr(x, "annotations", exact = TRUE),
      group = opts$group %||% "subgroup"
    )
    # Drop names we set explicitly so do.call does not see duplicates.
    if (length(dots)) {
      dots <- dots[setdiff(names(dots), c(names(override), "show_annotations", "annotation_digits", "layout", "printplot", "quiet"))]
    }
    fp <- do.call(
      forestplot_render,
      c(render_args, override, dots)
    )
    if (!is.null(attr(x, "footnote", exact = TRUE))) {
      attr(fp, "footnote") <- attr(x, "footnote", exact = TRUE)
    }
    if (!is.null(attr(x, "annotations", exact = TRUE))) {
      attr(fp, "annotations") <- attr(x, "annotations", exact = TRUE)
    }
  } else {
    snap <- attr(x, "forestplot_render", exact = TRUE)
    if (!is.list(snap)) {
      return(x)
    }
    snap$base_size <- metrics$base_size
    snap$pointsize <- metrics$pointsize
    snap$panel.height <- metrics$panel.height
    snap$typography_auto <- TRUE
    snap$printplot <- NULL
    snap$quiet <- NULL
    fp <- do.call(forestplot_render, snap)
  }

  attr(fp, "typography_auto") <- TRUE
  attr(fp, "display_inches") <- metrics$device_inches
  attr(fp, "has_favours") <- has_favours
  attr(fp, "show_ci_text") <- show_ci_text
  fp
}

#' @noRd
forestplot_add_class <- function(x) {
  class(x) <- unique(c("interactmeta_forestplot", class(x)))
  x
}

#' @noRd
forestplot_has_graphics_device <- function() {
  grDevices::dev.cur() > 1L
}

#' @noRd
forestplot_maybe_print <- function(x) {
  # Never call dev.new() / open a window — draw on the active device only.
  # In RStudio the Plots pane is used; knitr already has a device open.
  forestplot_draw(x)
  invisible(x)
}
