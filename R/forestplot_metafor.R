#' Forest-plot renderer via bayesmeta / forestplot
#'
#' Specs packaging for `forestplot_bayesmeta_draw()`. Drawing uses
#' [bayesmeta::forestplot.bayesmeta] (single block) or [forestplot::forestplot]
#' with the same conventions (multi-block CAMS forests).
#'
#' @noRd
forestplot_render_interactmeta_panels <- function(panels,
                                          row_labels,
                                          xlim,
                                          xticks,
                                          xlab,
                                          col_right_heading,
                                          estimate_col = "yi",
                                          base_size = interactmeta_forest_base_size(),
                                          pointsize = interactmeta_forest_pointsize(),
                                          digits = 2,
                                          exponentiate = FALSE,
                                          nullval = 0,
                                          plot.margin = interactmeta_forestplot_margin(),
                                          panel.height = NULL,
                                          show_left = TRUE,
                                          left_names = character(),
                                          left_headings = character(),
                                          show_ci_text = FALSE,
                                          favours = NULL,
                                          forest_args = list(),
                                          annotations = NULL,
                                          group = "subgroup",
                                          reference = NULL,
                                          panel_roles = NULL) {
  n <- length(panels)
  group <- match.arg(group, c("subgroup", "trial"))
  xlab <- as.list(xlab)
  if (length(xlab) < n) {
    xlab <- rep_len(xlab, n)
  }
  if (!is.list(col_right_heading)) {
    col_right_heading <- as.list(rep_len(list(col_right_heading), n))
  } else if (length(col_right_heading) < n) {
    col_right_heading <- rep_len(col_right_heading, n)
  }
  if (!is.list(xlim)) {
    xlim <- rep(list(xlim), n)
  } else if (length(xlim) < n) {
    xlim <- rep_len(xlim, n)
  }
  if (!is.list(xticks)) {
    xticks <- rep(list(xticks), n)
  } else if (length(xticks) < n) {
    xticks <- rep_len(xticks, n)
  }
  if (!is.list(forest_args)) {
    forest_args <- list()
  }
  annotations <- forestplot_normalize_annotations(annotations)
  if (is.character(reference) && length(reference)) {
    reference <- as.character(reference)[[1L]]
    if (!nzchar(reference) || is.na(reference)) {
      reference <- NULL
    }
  } else if (!is.null(reference)) {
    reference <- as.character(reference)[[1L]]
  }
  if (is.null(panel_roles)) {
    panel_roles <- rep("subgroups", n)
  } else {
    panel_roles <- as.character(panel_roles)
    if (length(panel_roles) < n) {
      panel_roles <- rep_len(panel_roles, n)
    }
  }
  # One panel: keep a data.frame so `fp$specs$row_labels$heading` works.
  # Several panels: one label table per panel.
  if (is.data.frame(row_labels)) {
    if (n > 1L) {
      row_labels <- rep(list(row_labels), n)
    }
  } else if (is.list(row_labels) && length(row_labels) >= n) {
    if (n == 1L && is.data.frame(row_labels[[1L]])) {
      row_labels <- row_labels[[1L]]
    }
  } else {
    stop("`row_labels` must be a data frame or a list of panel label tables.", call. = FALSE)
  }

  specs <- list(
    panels = panels,
    row_labels = row_labels,
    xlim = xlim,
    xticks = xticks,
    xlab = xlab,
    col_right_heading = col_right_heading,
    estimate_col = estimate_col,
    base_size = base_size,
    pointsize = pointsize,
    digits = digits,
    exponentiate = isTRUE(exponentiate),
    nullval = nullval,
    show_left = isTRUE(show_left),
    left_names = left_names,
    left_headings = left_headings,
    show_ci_text = isTRUE(show_ci_text),
    favours = favours,
    forest_args = forest_args,
    annotations = annotations,
    group = group,
    reference = reference,
    panel_roles = panel_roles
  )

  out <- list(
    engine = "bayesmeta",
    specs = specs,
    figure = structure(list(engine = "bayesmeta"), class = "interactmeta_forestplot_bayesmeta"),
    plots = lapply(seq_len(n), function(i) {
      list(panel = i, engine = "bayesmeta", role = panel_roles[[i]])
    })
  )
  if (n == 1L) {
    out$plot <- out$plots[[1L]]
  }
  out
}

#' Draw forest plot — bayesmeta / forestplot backend
#' @noRd
forestplot_metafor_draw <- function(x) {
  forestplot_bayesmeta_draw(x)
}

#' Whether a forestplot result draws via bayesmeta (or legacy metafor tag)
#' @noRd
forestplot_is_metafor_engine <- function(x) {
  inherits(x, "interactmeta_forestplot") &&
    (identical(x$engine, "bayesmeta") ||
       identical(x$engine, "metafor") ||
       inherits(x$figure, "interactmeta_forestplot_bayesmeta") ||
       inherits(x$figure, "interactmeta_forestplot_metafor"))
}

#' @noRd
forestplot_is_bayesmeta_engine <- function(x) {
  inherits(x, "interactmeta_forestplot") &&
    (identical(x$engine, "bayesmeta") ||
       inherits(x$figure, "interactmeta_forestplot_bayesmeta"))
}

#' Soft white pad at the far right device edge only
#'
#' Used for vignette / PNG crop safety. Never pad the left edge — that clipped
#' the first character of study and subgroup labels.
#' @noRd
forestplot_metafor_pad_device_edges <- function(pad_in = 0.04) {
  pad_in <- suppressWarnings(as.numeric(pad_in)[1L])
  if (!is.finite(pad_in) || pad_in <= 0) {
    return(invisible(NULL))
  }
  din <- graphics::par("din")
  if (length(din) < 2L || any(!is.finite(din)) || any(din <= 0)) {
    return(invisible(NULL))
  }
  y0 <- graphics::grconvertY(0, from = "inches", to = "user")
  y1 <- graphics::grconvertY(din[[2L]], from = "inches", to = "user")
  x_right0 <- graphics::grconvertX(din[[1L]] - pad_in, from = "inches", to = "user")
  x_right1 <- graphics::grconvertX(din[[1L]], from = "inches", to = "user")
  graphics::rect(
    x_right0, y0, x_right1, y1,
    col = "white", border = NA, xpd = NA
  )
  invisible(NULL)
}

#' One panel following metafor's subgroup-forest recipe
#'
#' Study / shrinkage rows → [metafor::forest()] with `rows=`;
#' block headers → [graphics::text()];
#' pooled / prediction / overall → [metafor::addpoly()].
#'
#' @noRd
forestplot_metafor_one_panel <- function(panel,
                                         row_labels,
                                         xlim,
                                         xticks,
                                         xlab = "",
                                         title = NULL,
                                         estimate_col = "yi",
                                         base_size = interactmeta_forest_base_size(),
                                         pointsize = interactmeta_forest_pointsize(),
                                         digits = 2,
                                         exponentiate = FALSE,
                                         nullval = 0,
                                         show_ilab = FALSE,
                                         left_names = character(),
                                         left_headings = character(),
                                         show_ci_text = FALSE,
                                         favours = NULL,
                                         forest_args = list(),
                                         annotations = NULL,
                                         group = "subgroup",
                                         device_height_in = NULL) {
  n <- nrow(panel)
  if (!n) {
    stop("Forest panel data is empty.", call. = FALSE)
  }
  if (!is.list(forest_args)) {
    forest_args <- list()
  }
  group <- match.arg(group, c("subgroup", "trial"))
  # Layout no longer reserves in-plot note rows; IF / het go under the axis.
  annotations <- forestplot_normalize_annotations(annotations)

  # Favours go under the x-axis (via xlab), not on reserved in-plot rows.
  layout <- forestplot_metafor_subgroup_layout(
    panel,
    row_labels,
    estimate_col,
    heterogeneity = NULL,
    ifrac = NULL,
    favours = NULL,
    group = group
  )
  n_layout_rows <- max(1, ceiling(layout$ylim_hi - layout$ylim_lo))
  # Prefer the active figure-band height over the full device so cex shrinks
  # when rows are packed into a short band.
  if (is.null(device_height_in) || !is.finite(device_height_in) || device_height_in <= 0) {
    pin_h <- tryCatch(graphics::par("pin")[[2L]], error = function(e) NA_real_)
    mai_h <- tryCatch(sum(graphics::par("mai")[c(1L, 3L)]), error = function(e) NA_real_)
    if (is.finite(pin_h) && is.finite(mai_h)) {
      device_height_in <- pin_h + mai_h
    }
  }
  cex <- forestplot_metafor_cex(
    base_size,
    n_layout_rows = n_layout_rows,
    device_height_in = device_height_in
  )

  # Keep estimates on the log / original scale; transform the axis like the
  # metafor vignette (atransf = exp) when exponentiate = TRUE.
  yi <- layout$study$yi
  ci.lb <- layout$study$ci.lb
  ci.ub <- layout$study$ci.ub
  slab <- layout$study$slab
  rows <- layout$study$rows
  pch <- layout$study$pch
  psize <- layout$study$psize
  col <- layout$study$col
  pt <- suppressWarnings(as.numeric(pointsize)[1L])
  if (!is.finite(pt) || pt <= 0) {
    pt <- as.numeric(interactmeta_forest_pointsize())
  }
  # metafor draws symbols at cex * psize; scale by user/default pointsize.
  psize <- psize * max(0.5, min(2.0, pt / 10))
  # Summary diamonds must drive the axis when study yi is empty.
  poly_x <- unlist(lapply(layout$polys, function(p) {
    se <- sqrt(max(as.numeric(p$vi), 0))
    z <- stats::qnorm(0.975)
    c(p$yi, p$yi - z * se, p$yi + z * se)
  }), use.names = FALSE)
  if (isTRUE(exponentiate)) {
    atransf <- exp
    refline <- 0
    if (length(xlim) == 2L && all(is.finite(xlim)) && all(xlim > 0)) {
      alim <- log(pmax(xlim, 1e-8))
      at <- if (!is.null(xticks) && length(xticks)) {
        log(pmax(as.numeric(xticks), 1e-8))
      } else {
        log(pmax(
          forestplot_nice_axis(xlim, exponentiate = TRUE, nullval = 1)$xticks,
          1e-8
        ))
      }
    } else {
      finite_x <- c(yi, ci.lb, ci.ub, poly_x)
      finite_x <- finite_x[is.finite(finite_x)]
      alim <- if (length(finite_x)) {
        range(c(finite_x, 0), na.rm = TRUE)
      } else {
        c(-1, 1)
      }
      nice <- forestplot_nice_axis(alim, exponentiate = FALSE, nullval = 0)
      alim <- nice$xlim
      at <- nice$xticks
    }
  } else {
    atransf <- NULL
    refline <- if (is.null(nullval)) 0 else nullval
    if (length(xlim) != 2L || any(!is.finite(xlim))) {
      finite_x <- c(yi, ci.lb, ci.ub, poly_x)
      finite_x <- finite_x[is.finite(finite_x)]
      alim <- if (length(finite_x)) {
        range(c(finite_x, refline), na.rm = TRUE)
      } else {
        c(-1, 1)
      }
      nice <- forestplot_nice_axis(alim, exponentiate = FALSE, nullval = refline)
      alim <- nice$xlim
      at <- if (is.null(xticks) || !length(xticks)) nice$xticks else xticks
    } else {
      alim <- xlim
      at <- if (is.null(xticks) || !length(xticks)) {
        forestplot_nice_axis(alim, exponentiate = FALSE, nullval = refline)$xticks
      } else {
        xticks
      }
    }
  }

  # Keep every tick inside alim so metafor does not clamp/dedupe endpoints
  # (which drops labels). Prefer the tick range when we have an explicit grid.
  at <- as.numeric(at)
  at <- at[is.finite(at)]
  if (length(at)) {
    alim <- range(c(as.numeric(alim), at), na.rm = TRUE)
  }

  # Axis title + Favours are drawn after forest() so they sit on separate
  # margin lines (metafor length-3 xlab puts all three on one overlapping line).
  if (is.null(xlab) || !length(xlab) || all(!nzchar(as.character(xlab)))) {
    xlab <- forestplot_default_xlab(exponentiate)
  }
  xlab_title <- forestplot_metafor_xlab(xlab, favours = NULL)
  xlab_arg <- ""
  # Right header uses the classic IV / Random / CI label when title is default.
  header <- forestplot_metafor_header(
    if (is.null(title) || !length(title) ||
          identical(as.character(title)[[1L]], "Estimate") ||
          grepl("^Effect", as.character(title)[[1L]])) {
      "IV, Random, 95% CI"
    } else {
      title
    },
    slab_label = "Study"
  )

  is_placeholder <- isTRUE(layout$study$placeholder)

  ilab <- NULL
  ilab.lab <- NULL
  if (
    isTRUE(show_ilab) &&
      length(left_names) &&
      length(layout$study$index) &&
      !is_placeholder
  ) {
    ilab_full <- forestplot_metafor_ilab(panel, left_names)
    if (!is.null(ilab_full)) {
      ilab <- ilab_full[layout$study$index, , drop = FALSE]
      ilab.lab <- if (length(left_headings) == ncol(ilab)) {
        as.character(left_headings)
      } else {
        as.character(left_names)
      }
    }
  }

  # Drop missing outcomes ourselves. metafor::forest(na.action = "na.omit")
  # remaps remaining `rows` downward, which then collides with addpoly()
  # coordinates that still use the original layout.
  keep <- !is.na(yi) & !is.na(ci.lb) & !is.na(ci.ub)
  if (is_placeholder || !any(keep)) {
    # Invisible stub so forest() can set up axes without painting over diamonds.
    yi <- 0
    ci.lb <- 0
    ci.ub <- 0
    pch <- NA_integer_
    col <- NA_character_
  } else if (!all(keep)) {
    yi <- yi[keep]
    ci.lb <- ci.lb[keep]
    ci.ub <- ci.ub[keep]
    slab <- slab[keep]
    rows <- rows[keep]
    pch <- pch[keep]
    psize <- psize[keep]
    col <- col[keep]
    if (!is.null(ilab)) {
      ilab <- ilab[keep, , drop = FALSE]
    }
  }

  # Placeholder rows must not print Effect-column text (addpoly still can).
  forest_annotate <- isTRUE(show_ci_text) && !is_placeholder
  # addpoly() still draws CI text when show_ci_text is on — reserve right
  # pad for that column even if forest() itself does not annotate.
  annotate_pad <- isTRUE(show_ci_text)

  plot_xlim <- forestplot_metafor_plot_xlim(
    alim = alim,
    has_ilab = !is.null(ilab),
    annotate = annotate_pad
  )
  # Modest inset for Effect; nudge study labels inward so the first glyph is
  # never clipped at usr[1].
  textpos <- forestplot_metafor_textpos(
    plot_xlim,
    right_shift = 0.06,
    left_shift = 0.12
  )

  ilab.xpos <- NULL
  if (!is.null(ilab)) {
    ilab.xpos <- forestplot_metafor_ilab_xpos(
      xlim = plot_xlim,
      alim = alim,
      yi = yi,
      ci.lb = ci.lb,
      n_cols = ncol(ilab),
      left_names = left_names
    )
  }

  # Annotations (IF / heterogeneity) are drawn under the x-axis, not in-plot.
  ylim_lo <- layout$ylim_lo
  # metafor centres column headers one row above the rule — too tight for
  # tall left labels. When headers wrap or
  # group = "trial" (first content is a block heading at ylim_hi), reserve
  # extra top space and redraw the header chrome after forest().
  n_header_lines <- max(
    1L,
    forestplot_metafor_header_lines(ilab.lab),
    forestplot_metafor_header_lines(header)
  )
  # Gap from first content / trial heading up to the header rule.
  content_gap <- if (identical(group, "trial")) {
    0.55
  } else {
    0.40
  }
  # Always reserve room for column headers above the rule (we redraw them).
  header_body <- as.numeric(n_header_lines)
  pad_top <- max(0.45, 0.45 * n_header_lines + 0.15)
  # metafor `top` is empty rows above the highest content — keep tight.
  top <- header_body + pad_top + 0.25
  ylim <- c(
    ylim_lo,
    layout$ylim_hi + content_gap + header_body + pad_top
  )

  # Study labels are re-drawn after forest() with xpd = NA (passing blank
  # slabs here avoids clipped-then-doubled text from metafor's clip region).
  args <- list(
    x = yi,
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    slab = if (is_placeholder) slab else rep(" ", length(slab)),
    xlim = plot_xlim,
    alim = alim,
    at = at,
    ylim = ylim,
    top = top,
    rows = rows,
    refline = refline,
    xlab = xlab_arg,
    # Suppress metafor header / ilab.lab chrome — we redraw cleanly below
    # (avoids orphaned dashes next to the Effect column header).
    header = FALSE,
    annotate = forest_annotate,
    # digits[1] = annotations; digits[2] = axis tick labels (need enough
    # decimals so every tick — e.g. 0.125 on a ratio axis — prints a number).
    digits = {
      ann_d <- as.integer(digits)[1L]
      if (!is.finite(ann_d) || ann_d < 0L) {
        ann_d <- 2L
      }
      axis_d <- if (isTRUE(exponentiate)) {
        max(ann_d, 3L)
      } else {
        max(ann_d, 1L)
      }
      c(ann_d, axis_d)
    },
    pch = pch,
    psize = psize,
    col = col,
    ilab = ilab,
    ilab.lab = NULL,
    ilab.xpos = ilab.xpos,
    textpos = textpos,
    # Study CI end-caps off; modest diamond height via length-2 efac.
    efac = forestplot_default_efac(),
    cex = cex,
    cex.lab = max(cex, 0.85),
    # Tick numbers smaller than body text so every label fits.
    cex.axis = max(0.45, min(0.75, cex * 0.65))
  )
  if (!is.null(atransf)) {
    args$atransf <- atransf
  }
  # Suppress metafor's header rule (lty[2]) — we draw our own inset rule.
  args$lty <- c("solid", "blank")

  # User `...` from forestplot() / forestplot.interactmeta() overrides defaults,
  # but layout-critical arguments stay under interactmeta control.
  args <- forestplot_merge_metafor_args(args, forest_args, which = "forest")

  # Keep header rule suppressed even if the user passes lty.
  if (!is.null(args$lty)) {
    lty_v <- args$lty
    if (is.character(lty_v) || is.numeric(lty_v)) {
      lty_v <- rep_len(lty_v, 2L)
      lty_v[[2L]] <- if (is.character(lty_v)) "blank" else 0
      args$lty <- lty_v
    }
  }

  suppressWarnings(do.call(metafor::forest, args))

  # forest() restores the pre-compensated `mar` on exit. Leave it — changing
  # mar after plotting breaks the user coordinate mapping used for Favours /
  # axis-title annotations below.

  # Re-paint study labels with xpd = NA so the first glyph cannot be clipped
  # at the plot-region edge (metafor draws them with the default clip).
  if (!is_placeholder && length(layout$study$slab)) {
    forestplot_metafor_repaint_slabs(
      slab = layout$study$slab,
      rows = layout$study$rows,
      textpos = textpos,
      cex = cex,
      plot_xlim = plot_xlim,
      ilab.xpos = ilab.xpos
    )
  }

  # Draw column headers ourselves (metafor header/ilab.lab were suppressed).
  forestplot_metafor_redraw_column_headers(
    ylim = ylim,
    top = top,
    ylim_hi = layout$ylim_hi,
    content_gap = content_gap,
    header = header,
    ilab.lab = ilab.lab,
    ilab.xpos = ilab.xpos,
    plot_xlim = plot_xlim,
    textpos = textpos,
    cex = cex,
    draw_rule = TRUE
  )

  # Full-height vertical null line through the content band.
  if (is.finite(refline)) {
    y_line <- layout$ylim_hi + content_gap
    graphics::segments(
      refline,
      ylim_lo,
      refline,
      y_line,
      lwd = 0.8,
      xpd = TRUE
    )
  }

  # Subgroup / block headers (metafor: text(..., pos = 4, font = 4))
  op <- graphics::par(cex = cex, font = 4, xpd = NA)
  on.exit(graphics::par(op), add = TRUE)
  if (length(layout$headers$label)) {
    graphics::text(
      forestplot_metafor_safe_label_x(textpos),
      layout$headers$row,
      pos = 4,
      offset = 0.35,
      layout$headers$label,
      xpd = NA
    )
  }
  graphics::par(op)

  # Summary polygons under each block / at the bottom
  if (length(layout$polys)) {
    for (poly in layout$polys) {
      poly_args <- list(
        x = poly$yi,
        vi = poly$vi,
        rows = poly$row,
        # Labels re-drawn below with xpd = NA (avoids first-glyph clip).
        mlab = " ",
        cex = cex,
        efac = forestplot_default_poly_efac(),
        col = poly$col,
        border = poly$col,
        # Bold left labels (Overall / Prediction / …) and right CI text.
        fonts = forestplot_default_poly_fonts()
      )
      if (!is.null(atransf)) {
        poly_args$atransf <- atransf
      }
      if (isTRUE(show_ci_text)) {
        poly_args$annotate <- TRUE
        poly_args$digits <- as.integer(digits)[1L]
      } else {
        poly_args$annotate <- FALSE
      }
      poly_args <- forestplot_merge_metafor_args(
        poly_args,
        forest_args,
        which = "addpoly"
      )
      suppressWarnings(do.call(metafor::addpoly, poly_args))
    }
    # Bold summary labels (Overall / Prediction / …) without clip.
    op_lab <- graphics::par(cex = cex, font = 2, xpd = NA)
    on.exit(graphics::par(op_lab), add = TRUE)
    for (poly in layout$polys) {
      if (!nzchar(poly$mlab)) {
        next
      }
      graphics::text(
        forestplot_metafor_safe_label_x(textpos),
        poly$row,
        labels = poly$mlab,
        pos = 4,
        offset = 0.35,
        xpd = NA
      )
    }
    graphics::par(op_lab)
  }

  # Heterogeneity / IF notes are drawn under the x-axis (not in-plot).

  # In-plot favours rows (legacy layout path) — axis favours drawn separately.
  if (length(layout$favours_rows) && !is.null(favours) && length(favours) >= 2L) {
    forestplot_metafor_draw_favours(
      favours = favours,
      rows = layout$favours_rows,
      alim = alim,
      refline = refline,
      cex = cex,
      exponentiate = exponentiate
    )
  }

  # Axis title on its own margin line; Favours+arrows on the line below;
  # annotation footnotes under that.
  forestplot_metafor_draw_axis_title(
    xlab = xlab_title,
    alim = alim,
    cex = max(cex, 0.85)
  )
  if (!is.null(favours) && length(favours) >= 2L) {
    forestplot_metafor_draw_axis_favours(
      favours = favours,
      alim = alim,
      refline = refline,
      cex = cex
    )
  }
  forestplot_metafor_draw_axis_annotations(
    lines = forestplot_annotation_lines(annotations),
    plot_xlim = plot_xlim,
    has_favours = !is.null(favours) && length(favours) >= 2L,
    cex = cex
  )

  invisible(NULL)
}

#' Normalize forest annotation payload
#' @noRd
forestplot_normalize_annotations <- function(annotations) {
  empty <- list(
    lines = character(),
    ifrac = character(),
    subgroup_ifrac = character(),
    heterogeneity = character()
  )
  if (is.null(annotations)) {
    return(empty)
  }
  if (is.character(annotations)) {
    empty$lines <- annotations[nzchar(annotations)]
    return(empty)
  }
  if (!is.list(annotations)) {
    return(empty)
  }
  lines <- annotations$lines %||% character()
  lines <- as.character(lines)
  lines <- lines[!is.na(lines) & nzchar(lines)]

  clean_named <- function(x) {
    x <- x %||% character()
    if (!length(x)) {
      return(character())
    }
    if (!is.null(names(x))) {
      x <- x[nzchar(names(x)) & !is.na(x) & nzchar(as.character(x))]
    } else {
      x <- x[!is.na(x) & nzchar(as.character(x))]
    }
    x
  }

  list(
    lines = lines,
    ifrac = clean_named(annotations$ifrac),
    subgroup_ifrac = clean_named(annotations$subgroup_ifrac),
    heterogeneity = clean_named(annotations$heterogeneity)
  )
}

#' Group forest footnotes into rows
#'
#' All notes (IF, heterogeneity, interaction heterogeneity, and free-form
#' `lines`) share a single row, joined with `" \u00b7 "` when drawn.
#' @noRd
forestplot_annotation_footnote_rows <- function(annotations) {
  annotations <- forestplot_normalize_annotations(annotations)
  clean <- function(x) {
    x <- as.character(x)
    x[!is.na(x) & nzchar(x)]
  }
  if_notes <- character()
  if (length(annotations$ifrac)) {
    if_notes <- clean(unname(annotations$ifrac))
  }
  het_cells <- character()
  het <- annotations$heterogeneity
  if (length(het)) {
    if (!is.null(names(het))) {
      for (nm in c("combined", "interaction")) {
        if (nm %in% names(het)) {
          het_cells <- c(het_cells, clean(het[[nm]][[1L]]))
        }
      }
      extra <- setdiff(names(het), c("combined", "interaction"))
      if (length(extra)) {
        het_cells <- c(het_cells, clean(unname(het[extra])))
      }
    } else {
      het_cells <- clean(het)
    }
  }
  extra_lines <- character()
  if (length(annotations$lines)) {
    extra_lines <- clean(annotations$lines)
  }
  cells <- c(if_notes, het_cells, extra_lines)
  if (!length(cells)) {
    return(list())
  }
  list(cells)
}

#' Flatten annotation payload to footnote lines (axis footer order)
#'
#' All notes are joined on one line (`" \u00b7 "`).
#' @noRd
forestplot_annotation_lines <- function(annotations) {
  rows <- forestplot_annotation_footnote_rows(annotations)
  if (!length(rows)) {
    return(character())
  }
  vapply(
    rows,
    function(cells) paste(as.character(cells), collapse = " \u00b7 "),
    character(1L),
    USE.NAMES = FALSE
  )
}

#' Draw annotation footnotes under the x-axis (below title / Favours)
#' @noRd
forestplot_metafor_draw_axis_annotations <- function(lines,
                                                     plot_xlim,
                                                     has_favours = FALSE,
                                                     cex = 1) {
  lines <- as.character(lines)
  lines <- lines[!is.na(lines) & nzchar(lines)]
  if (!length(lines)) {
    return(invisible(NULL))
  }
  # Below Favours band when present; otherwise just under the axis title.
  start_line <- if (isTRUE(has_favours)) 5.9 else 3.2
  cex_lab <- max(0.55, cex * 0.78)
  op <- graphics::par(xpd = NA)
  on.exit(graphics::par(op), add = TRUE)
  x_mid <- mean(plot_xlim)
  if (!is.finite(x_mid)) {
    x_mid <- plot_xlim[[1L]]
  }
  for (i in seq_along(lines)) {
    y <- graphics::grconvertY(-(start_line + (i - 1L) * 0.78), from = "lines", to = "user")
    if (!is.finite(y)) {
      next
    }
    graphics::text(
      x_mid,
      y,
      labels = forestplot_footnote_as_plotmath(lines[[i]]),
      adj = c(0.5, 1),
      cex = cex_lab,
      xpd = NA
    )
  }
  invisible(NULL)
}

#' Draw multi-line annotation footer under a metafor forest panel
#'
#' Legacy in-plot footer (kept for tests / callers). New forests use
#' `forestplot_metafor_draw_axis_annotations()`.
#' @noRd
forestplot_metafor_draw_annotations <- function(lines, xlim, ylim_lo, cex = 1) {
  lines <- as.character(lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) {
    return(invisible(NULL))
  }
  x_left <- xlim[[1L]]
  y0 <- ylim_lo + 0.35
  ys <- y0 + 0.85 * (rev(seq_along(lines)) - 1L)
  op <- graphics::par(cex = cex * 1.05, font = 1)
  on.exit(graphics::par(op), add = TRUE)
  graphics::text(
    x_left,
    ys,
    labels = as.expression(lapply(lines, forestplot_footnote_as_plotmath)),
    pos = 4,
    xpd = TRUE
  )
  invisible(NULL)
}

#' Merge user metafor args into interactmeta defaults
#'
#' User `...` from [forestplot()] / `forest.interactmeta()` overrides presentation
#' defaults. Layout-critical data vectors stay under interactmeta control.
#'
#' @param which `"forest"` or `"addpoly"` — reserved names that must not be
#'   overridden by the user (layout / data vectors).
#' @noRd
forestplot_merge_metafor_args <- function(args, user = list(), which = c("forest", "addpoly")) {
  which <- match.arg(which)
  if (!length(user)) {
    return(args)
  }
  # Data / geometry owned by interactmeta subgroup layout (metafor compliance otherwise).
  reserved <- if (identical(which, "forest")) {
    c("x", "vi", "sei", "ci.lb", "ci.ub", "slab", "rows", "ylim")
  } else {
    c("x", "vi", "sei", "ci.lb", "ci.ub", "row", "rows", "mlab")
  }
  user <- user[setdiff(names(user), reserved)]
  user <- user[nzchar(names(user))]
  if (!length(user)) {
    return(args)
  }
  utils::modifyList(args, user, keep.null = TRUE)
}

#' Build rows / headers / addpoly layout for one forest panel
#'
#' Mirrors the metafor subgroup example: study rows occupy `rows=`, headers sit
#' above each block, and summary diamonds use half-row positions for `addpoly`.
#' Summary rows (`pooled` / `prediction`) are attached under the matching
#' subgroup header even when they appear later in the panel data.
#' Heterogeneity and overall-IF notes are no longer reserved in-plot (they are
#' drawn under the x-axis). When `favours` is length 2, left/right labels may
#' still be reserved on in-plot rows for legacy callers.
#'
#' @noRd
forestplot_metafor_subgroup_layout <- function(panel,
                                               row_labels,
                                               estimate_col = "yi",
                                               heterogeneity = NULL,
                                               ifrac = NULL,
                                               favours = NULL,
                                               group = "subgroup") {
  group <- match.arg(group, c("subgroup", "trial"))
  n <- nrow(panel)
  est <- panel[[estimate_col]]
  if (is.null(est)) {
    est <- panel$yi
  }
  ci.lb <- panel$ci.lb
  ci.ub <- panel$ci.ub
  vi <- if ("vi" %in% names(panel)) {
    as.numeric(panel$vi)
  } else {
    ((ci.ub - ci.lb) / (2 * stats::qnorm(0.975)))^2
  }

  heading <- if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    rep("", n)
  }
  label <- if ("label" %in% names(row_labels)) {
    as.character(row_labels$label)
  } else {
    rep("", n)
  }
  heading[is.na(heading)] <- ""
  label[is.na(label)] <- ""
  label[comparison_forest_is_blank_row_label(label)] <- ""
  heading[comparison_forest_is_blank_row_label(heading)] <- ""

  if (any(nzchar(heading))) {
    last <- ""
    for (i in seq_len(n)) {
      if (nzchar(heading[[i]])) {
        last <- heading[[i]]
      } else if (nzchar(last)) {
        heading[[i]] <- last
      }
    }
  }

  row_kind <- if ("row_kind" %in% names(panel)) {
    as.character(panel$row_kind)
  } else {
    rep("study", n)
  }
  row_kind[is.na(row_kind) | !nzchar(row_kind)] <- "study"

  pch <- if ("shape" %in% names(panel)) {
    as.integer(panel$shape)
  } else {
    rep(15L, n)
  }
  pch[is.na(pch)] <- 15L

  col <- if ("colour" %in% names(panel)) {
    as.character(panel$colour)
  } else {
    rep(comparison_forest_marker_colour(), n)
  }
  col[is.na(col) | !nzchar(col)] <- comparison_forest_marker_colour()

  pscale <- if ("point_scale" %in% names(panel)) {
    ps <- as.numeric(panel$point_scale)
    ps[is.na(ps)] <- 1
    ps
  } else {
    rep(1, n)
  }
  psize <- pmax(0.4, pmin(1.5, pscale))

  is_point <- row_kind %in% c("study", "shrinkage")
  is_poly <- row_kind %in% c("overall", "prediction")
  is_combined <- row_kind == "combined"

  # Study / shrinkage blocks in first-seen heading order
  point_idx <- which(is_point)
  block_order <- unique(heading[point_idx])
  if (!length(block_order)) {
    # summaries only — Combined is drawn separately via combined_idx
    block_order <- unique(heading[is_poly])
  }
  block_order <- block_order[nzchar(block_order) & !is.na(block_order)]
  # Combined study / shrinkage are drawn with the Combined diamonds below.
  combined_point_idx <- point_idx[heading[point_idx] == "Combined"]
  block_order <- setdiff(block_order, "Combined")

  # When group = "trial", study headings are trials but subgroup pooled /
  # prediction rows keep subgroup headings. Append those summary-only blocks
  # (with both pooled and prediction) after trial blocks and before Interaction.
  poly_headings <- unique(heading[is_poly])
  poly_headings <- poly_headings[nzchar(poly_headings) & !is.na(poly_headings)]
  extra_sg <- setdiff(poly_headings, c(block_order, "Interaction", "Combined"))
  if (length(extra_sg)) {
    # Preserve subgroup order from the poly rows as they appear in the panel.
    extra_sg <- unique(heading[is_poly][heading[is_poly] %in% extra_sg])
    if ("Interaction" %in% block_order) {
      ix <- match("Interaction", block_order)
      block_order <- c(
        block_order[seq_len(ix - 1L)],
        extra_sg,
        block_order[ix:length(block_order)]
      )
    } else {
      block_order <- c(block_order, extra_sg)
    }
  }
  # Interaction may have only pooled/prediction rows (e.g. estimates_only).
  if ("Interaction" %in% poly_headings && !"Interaction" %in% block_order) {
    block_order <- c(block_order, "Interaction")
  }

  poly_by_hdr <- split(which(is_poly), heading[is_poly])
  combined_idx <- which(is_combined)

  headers <- list(label = character(), row = numeric())
  notes <- list(label = character(), row = numeric())
  favours_rows <- numeric()
  study_idx <- integer()
  study_rows <- numeric()
  polys <- list()

  show_favours <- !is.null(favours) && length(favours) >= 2L
  has_interaction_pred <- "Interaction" %in% names(poly_by_hdr) &&
    any(row_kind[poly_by_hdr[["Interaction"]]] == "prediction", na.rm = TRUE)
  has_combined_pred <- length(combined_idx) > 0L &&
    any(tolower(label[combined_idx]) == "prediction", na.rm = TRUE)

  # Last subgroup block that has a prediction diamond (for favours when Int/Combined absent).
  sg_blocks <- setdiff(block_order, c("Interaction", "Combined"))
  sg_pred_blocks <- sg_blocks[
    vapply(sg_blocks, function(h) {
      ix <- poly_by_hdr[[h]]
      length(ix) > 0L && any(row_kind[ix] == "prediction", na.rm = TRUE)
    }, logical(1L))
  ]
  last_subgroup_pred_hdr <- if (length(sg_pred_blocks)) {
    sg_pred_blocks[[length(sg_pred_blocks)]]
  } else {
    NA_character_
  }

  het <- heterogeneity %||% character()
  het_interaction <- if ("interaction" %in% names(het)) {
    as.character(het[["interaction"]])
  } else {
    NA_character_
  }
  het_overall <- if ("combined" %in% names(het)) {
    as.character(het[["combined"]])
  } else {
    NA_character_
  }
  has_het_interaction <- !is.na(het_interaction) && nzchar(het_interaction)
  has_het_overall <- !is.na(het_overall) && nzchar(het_overall)

  ifrac_notes <- ifrac %||% character()
  if (!is.null(names(ifrac_notes))) {
    ifrac_notes <- ifrac_notes[nzchar(names(ifrac_notes)) & nzchar(as.character(ifrac_notes))]
  } else {
    ifrac_notes <- character()
  }
  ifrac_remaining <- ifrac_notes

  # Within each block, draw pooled then prediction so notes follow prediction.
  order_polys <- function(idx) {
    if (!length(idx)) {
      return(idx)
    }
    rk <- row_kind[idx]
    idx[order(match(rk, c("overall", "prediction"), nomatch = 99L))]
  }

  add_poly <- function(j, hdr) {
    mlab <- label[[j]]
    colour_lab <- if ("colour_label" %in% names(row_labels)) {
      as.character(row_labels$colour_label[[j]])
    } else {
      NA_character_
    }
    if (identical(row_kind[[j]], "combined")) {
      sg <- if (!is.na(colour_lab) && nzchar(colour_lab)) {
        colour_lab
      } else if (nzchar(mlab)) {
        mlab
      } else {
        ""
      }
      mlab <- if (nzchar(sg)) paste("Combined:", sg) else "Combined"
    } else {
      if (!nzchar(mlab)) {
        mlab <- paste0(
          toupper(substr(row_kind[[j]], 1L, 1L)),
          substring(row_kind[[j]], 2L)
        )
      }
      if (nzchar(hdr) && nzchar(mlab)) {
        mlab <- paste0(mlab, " (", hdr, ")")
      }
    }
    yi_j <- est[[j]]
    vi_j <- vi[[j]]
    if (!is.finite(vi_j) || vi_j <= 0) {
      if (is.finite(ci.lb[[j]]) && is.finite(ci.ub[[j]])) {
        vi_j <- ((ci.ub[[j]] - ci.lb[[j]]) / (2 * stats::qnorm(0.975)))^2
      } else {
        vi_j <- 1e-6
      }
    }
    if (!is.finite(yi_j)) {
      return(NULL)
    }
    list(
      yi = yi_j,
      vi = max(vi_j, 1e-12),
      mlab = mlab,
      col = col[[j]]
    )
  }

  # Estimate height (top padding + blocks + overall + in-plot notes + bottom)
  # Extra full rows:
  # - after the column header / before first content (all group modes)
  # - before Interaction / Combined headers when they follow other blocks
  # - when group = "subgroup": between consecutive subgroup blocks
  # Summary-only subgroup blocks omit the alone heading (name stays on mlab).
  subgroup_block_gap <- if (identical(group, "subgroup")) 1 else 0
  height <- 1.25 + 1L
  for (i_hdr in seq_along(block_order)) {
    hdr <- block_order[[i_hdr]]
    n_pt <- sum(heading[point_idx] == hdr)
    n_py <- length(poly_by_hdr[[hdr]])
    is_int_ov <- forestplot_heading_is_interaction(hdr) ||
      forestplot_heading_is_combined(hdr)
    skip_heading <- !is_int_ov && identical(n_pt, 0L)
    if (
      i_hdr > 1L &&
        is_int_ov
    ) {
      height <- height + 1L
    }
    # Header row only when studies (or Interaction / Combined) need it.
    height <- height + if (skip_heading) 0L else 1L
    height <- height + n_pt + (if (n_py > 0) n_py + 0.25 else 0) + 0.40
    if (hdr %in% names(ifrac_notes)) {
      height <- height + 1L
    }
    if (identical(hdr, "Interaction") && has_het_interaction) {
      height <- height + 1L
    }
    if (identical(hdr, "Interaction") && show_favours && has_interaction_pred) {
      height <- height + 1L
    }
    if (
      show_favours &&
        !has_interaction_pred &&
        !has_combined_pred &&
        identical(hdr, last_subgroup_pred_hdr)
    ) {
      height <- height + 1L
    }
    if (
      subgroup_block_gap > 0 &&
        i_hdr < length(block_order) &&
        !forestplot_heading_is_interaction(hdr) &&
        !forestplot_heading_is_combined(hdr)
    ) {
      nxt <- block_order[[i_hdr + 1L]]
      if (
        !forestplot_heading_is_interaction(nxt) &&
          !forestplot_heading_is_combined(nxt)
      ) {
        height <- height + subgroup_block_gap
      }
    }
  }
  if (length(combined_idx) || length(combined_point_idx)) {
    # Spacer before Combined header + header + points + polys + trailing.
    # If Combined is the only content (no block_order), the post-header gap
    # already covers the lead-in — still add the mid-plot Combined spacer when
    # subgroup / interaction blocks precede it.
    ov_spacer <- if (length(block_order)) 1L else 0L
    height <- height + ov_spacer + 1L + length(combined_point_idx) + length(combined_idx) + 0.25
  }
  if (show_favours && has_combined_pred) {
    height <- height + 1L
  }
  if (has_het_overall) {
    height <- height + 1L
  }
  # Notes with no matching block: reserve trailing rows
  unmatched_ifrac <- setdiff(names(ifrac_notes), block_order)
  height <- height + length(unmatched_ifrac)
  if (has_het_interaction && !"Interaction" %in% block_order) {
    height <- height + 1L
  }

  y <- as.numeric(height)
  ylim_hi <- y
  # Gap under the column header before the first content block.
  y <- y - 1

  for (i_hdr in seq_along(block_order)) {
    hdr <- block_order[[i_hdr]]
    n_pt <- sum(heading[point_idx] == hdr)
    is_int_ov <- forestplot_heading_is_interaction(hdr) ||
      forestplot_heading_is_combined(hdr)
    skip_heading <- !is_int_ov && identical(n_pt, 0L)
    if (i_hdr > 1L && is_int_ov) {
      y <- y - 1
    }
    if (!skip_heading) {
      headers$label <- c(headers$label, if (nzchar(hdr)) hdr else " ")
      headers$row <- c(headers$row, y)
      y <- y - 1
    }

    pt_ii <- point_idx[heading[point_idx] == hdr]
    if (length(pt_ii)) {
      rows_b <- y - seq_along(pt_ii) + 1
      study_idx <- c(study_idx, pt_ii)
      study_rows <- c(study_rows, rows_b)
      y <- y - length(pt_ii)
    }

    py_ii <- order_polys(poly_by_hdr[[hdr]])
    if (length(py_ii)) {
      for (j in py_ii) {
        p <- add_poly(j, hdr)
        if (!is.null(p)) {
          p$row <- y
          polys[[length(polys) + 1L]] <- p
        }
        y <- y - 1
      }
      y <- y - 0.25
    }

    # Favours left/right under Interaction prediction.
    if (
      identical(hdr, "Interaction") &&
        show_favours &&
        has_interaction_pred
    ) {
      favours_rows <- c(favours_rows, y)
      y <- y - 1
    }

    # Subgroups-only (or any plot without Interaction/Combined): favours under
    # the last subgroup prediction block so the feature still appears.
    if (
      show_favours &&
        !has_interaction_pred &&
        !has_combined_pred &&
        identical(hdr, last_subgroup_pred_hdr)
    ) {
      favours_rows <- c(favours_rows, y)
      y <- y - 1
    }

    # Combined IF for this subgroup sits after its prediction (IF) diamond.
    if (hdr %in% names(ifrac_remaining)) {
      notes$label <- c(notes$label, ifrac_remaining[[hdr]])
      notes$row <- c(notes$row, y)
      y <- y - 1
      ifrac_remaining <- ifrac_remaining[setdiff(names(ifrac_remaining), hdr)]
    }

    # Heterogeneity for the interaction block sits after its prediction diamond.
    if (identical(hdr, "Interaction") && has_het_interaction) {
      notes$label <- c(notes$label, het_interaction)
      notes$row <- c(notes$row, y)
      y <- y - 1
      has_het_interaction <- FALSE
    }

    y <- y - 0.40
    # Full empty row between consecutive subgroup blocks.
    if (
      subgroup_block_gap > 0 &&
        i_hdr < length(block_order) &&
        !forestplot_heading_is_interaction(hdr) &&
        !forestplot_heading_is_combined(hdr)
    ) {
      nxt <- block_order[[i_hdr + 1L]]
      if (
        !forestplot_heading_is_interaction(nxt) &&
          !forestplot_heading_is_combined(nxt)
      ) {
        y <- y - subgroup_block_gap
      }
    }
  }

  if (length(combined_idx) || length(combined_point_idx)) {
    if (length(block_order)) {
      y <- y - 1
    }
    headers$label <- c(headers$label, "Combined")
    headers$row <- c(headers$row, y)
    y <- y - 1

    if (length(combined_point_idx)) {
      # Preserve input order (study then its shrinkage per trial).
      for (j in combined_point_idx) {
        study_idx <- c(study_idx, j)
        study_rows <- c(study_rows, y)
        y <- y - 1
      }
    }

    if (length(combined_idx)) {
      ov_ordered <- combined_idx[order(match(
        label[combined_idx],
        c("Overall", "Prediction"),
        nomatch = 99L
      ))]
      for (j in ov_ordered) {
        p <- add_poly(j, "Combined")
        if (!is.null(p)) {
          p$row <- y
          polys[[length(polys) + 1L]] <- p
        }
        y <- y - 1
      }
    }
    y <- y - 0.25
  }

  # Favours left/right under Combined prediction.
  if (show_favours && has_combined_pred) {
    favours_rows <- c(favours_rows, y)
    y <- y - 1
  }

  # Combined treatment-effect heterogeneity after Combined prediction (or at end).
  if (has_het_overall) {
    notes$label <- c(notes$label, het_overall)
    notes$row <- c(notes$row, y)
    y <- y - 1
  }

  # Fallback notes when no matching heading block existed.
  if (length(ifrac_remaining)) {
    for (nm in names(ifrac_remaining)) {
      notes$label <- c(notes$label, ifrac_remaining[[nm]])
      notes$row <- c(notes$row, y)
      y <- y - 1
    }
  }
  if (has_het_interaction) {
    notes$label <- c(notes$label, het_interaction)
    notes$row <- c(notes$row, y)
    y <- y - 1
  }

  # Leave a small gap above the x-axis (avoid large empty band under the plot).
  ylim_lo <- y - 0.5

  slab <- label[study_idx]
  # Shrinkage sits under its study — omit the repeated study name.
  slab[row_kind[study_idx] == "shrinkage"] <- ""
  slab[!nzchar(slab)] <- " "

  # Edge case: no study rows (e.g. estimates_only) — forest() still needs one
  # input row. Park it below the content so its CI / annotation cannot land on
  # a pooled or prediction diamond (ylim_hi - 2 collides with the first block's
  # prediction row when there are no study points).
  if (!length(study_idx)) {
    return(list(
      study = list(
        index = 1L,
        rows = ylim_lo - 1,
        yi = NA_real_,
        ci.lb = NA_real_,
        ci.ub = NA_real_,
        slab = " ",
        pch = NA_integer_,
        psize = 1,
        col = NA_character_,
        placeholder = TRUE
      ),
      headers = headers,
      notes = notes,
      favours_rows = favours_rows,
      polys = polys,
      ylim_lo = ylim_lo,
      ylim_hi = ylim_hi
    ))
  }

  list(
    study = list(
      index = study_idx,
      rows = study_rows,
      yi = est[study_idx],
      ci.lb = ci.lb[study_idx],
      ci.ub = ci.ub[study_idx],
      slab = slab,
      pch = pch[study_idx],
      psize = psize[study_idx],
      col = col[study_idx]
    ),
    headers = headers,
    notes = notes,
    favours_rows = favours_rows,
    polys = polys,
    ylim_lo = ylim_lo,
    ylim_hi = ylim_hi
  )
}

#' Draw “Favours …” labels left/right of the null line
#'
#' Labels sit in the plot region on reserved rows (Interaction / Combined
#' prediction), centered between the null line and each axis limit.
#'
#' @noRd
forestplot_metafor_draw_favours <- function(favours,
                                            rows,
                                            alim,
                                            refline,
                                            cex = 1,
                                            exponentiate = FALSE) {
  force(exponentiate)
  favours <- as.character(favours)
  rows <- as.numeric(rows)
  rows <- rows[is.finite(rows)]
  if (length(favours) < 2L || !length(rows)) {
    return(invisible(NULL))
  }
  alim <- as.numeric(alim)
  refline <- as.numeric(refline)[[1L]]
  if (length(alim) != 2L || any(!is.finite(alim)) || !is.finite(refline)) {
    return(invisible(NULL))
  }

  # Forest coordinates stay on the model scale (log when atransf = exp).
  x_left <- mean(c(alim[[1L]], refline))
  x_right <- mean(c(refline, alim[[2L]]))
  if (!is.finite(x_left) || !is.finite(x_right)) {
    return(invisible(NULL))
  }

  op <- graphics::par(cex = cex * 0.7, font = 1)
  on.exit(graphics::par(op), add = TRUE)
  for (y in rows) {
    graphics::text(x_left, y, labels = favours[[1L]], xpd = TRUE)
    graphics::text(x_right, y, labels = favours[[2L]], xpd = TRUE)
  }
  invisible(NULL)
}

#' Re-draw study labels without clipping the first character
#' @noRd
forestplot_metafor_repaint_slabs <- function(slab,
                                             rows,
                                             textpos,
                                             cex,
                                             plot_xlim,
                                             ilab.xpos = NULL) {
  force(plot_xlim)
  force(ilab.xpos)
  slab <- as.character(slab)
  rows <- as.numeric(rows)
  ok <- is.finite(rows) & nzchar(slab) & slab != " "
  if (!any(ok)) {
    return(invisible(NULL))
  }
  x_left <- forestplot_metafor_safe_label_x(textpos)
  if (!is.finite(x_left)) {
    return(invisible(NULL))
  }
  op <- graphics::par(xpd = NA, cex = cex, font = 1)
  on.exit(graphics::par(op), add = TRUE)
  graphics::text(
    x_left,
    rows[ok],
    labels = slab[ok],
    pos = 4,
    offset = 0.35,
    xpd = NA
  )
  invisible(NULL)
}

#' User-x for left labels, inset ≥0.2 in from the plot clip edge
#' @noRd
forestplot_metafor_safe_label_x <- function(textpos) {
  x <- as.numeric(textpos)[[1L]]
  if (!is.finite(x)) {
    return(NA_real_)
  }
  # Guarantee clearance from the plot-region clip edge in inches.
  x_edge_in <- graphics::grconvertX(graphics::par("usr")[[1L]], "user", "inches")
  x_safe <- graphics::grconvertX(x_edge_in + 0.22, "inches", "user")
  if (is.finite(x_safe)) {
    max(x, x_safe)
  } else {
    x
  }
}

#' Draw the x-axis title below the tick marks
#' @noRd
forestplot_metafor_draw_axis_title <- function(xlab, alim, cex = 1) {
  xlab <- as.character(xlab)[[1L]]
  if (!nzchar(xlab)) {
    return(invisible(NULL))
  }
  alim <- as.numeric(alim)
  at <- if (length(alim) == 2L && all(is.finite(alim))) {
    mean(alim)
  } else {
    mean(graphics::par("usr")[1:2])
  }
  op <- graphics::par(xpd = NA)
  on.exit(graphics::par(op), add = TRUE)
  graphics::mtext(
    xlab,
    side = 1,
    line = 1.7,
    at = at,
    cex = cex,
    font = 1
  )
  invisible(NULL)
}

#' Draw Favours labels + outward arrows under the axis title
#'
#' Axis title is on `line = 2.0`. Arrows and labels sit further down so they
#' never share a line with the title. X positions use the forest `alim` /
#' `refline` (model / log scale, same as `par("usr")`).
#'
#' @noRd
forestplot_metafor_draw_axis_favours <- function(favours,
                                                 alim,
                                                 refline,
                                                 cex = 1) {
  favours <- as.character(favours)
  if (length(favours) < 2L) {
    return(invisible(NULL))
  }
  alim <- as.numeric(alim)
  refline <- as.numeric(refline)[[1L]]
  if (length(alim) != 2L || any(!is.finite(alim)) || !is.finite(refline)) {
    return(invisible(NULL))
  }

  cex_lab <- max(0.5, cex * 0.65)
  op <- graphics::par(xpd = NA, cex = cex_lab, font = 1)
  on.exit(graphics::par(op), add = TRUE)

  y_arrow <- graphics::grconvertY(-3.6, from = "lines", to = "user")
  y_lab <- graphics::grconvertY(-4.7, from = "lines", to = "user")
  if (!is.finite(y_arrow) || !is.finite(y_lab)) {
    return(invisible(NULL))
  }

  arr_len <- min(0.07, graphics::par("cin")[[1L]] * 0.28)
  # Leave a clear gap around the null so the two arrows do not meet.
  span_left <- refline - alim[[1L]]
  span_right <- alim[[2L]] - refline
  gap_frac <- 0.50
  x0_left <- refline - gap_frac * span_left
  x0_right <- refline + gap_frac * span_right
  graphics::arrows(
    x0_left, y_arrow, alim[[1L]], y_arrow,
    length = arr_len, angle = 25, code = 2, lwd = 0.7, xpd = NA
  )
  graphics::arrows(
    x0_right, y_arrow, alim[[2L]], y_arrow,
    length = arr_len, angle = 25, code = 2, lwd = 0.7, xpd = NA
  )

  # Center labels in each half; shrink so they stay clear of the null gap.
  x_left <- mean(c(alim[[1L]], x0_left))
  x_right <- mean(c(x0_right, alim[[2L]]))
  max_w_left <- 0.88 * abs(x0_left - alim[[1L]])
  max_w_right <- 0.88 * abs(alim[[2L]] - x0_right)
  cex_left <- cex_lab
  cex_right <- cex_lab
  while (
    graphics::strwidth(favours[[1L]], cex = cex_left) > max_w_left &&
      cex_left > 0.4
  ) {
    cex_left <- cex_left * 0.9
  }
  while (
    graphics::strwidth(favours[[2L]], cex = cex_right) > max_w_right &&
      cex_right > 0.4
  ) {
    cex_right <- cex_right * 0.9
  }
  graphics::text(
    x_left,
    y_lab,
    labels = favours[[1L]],
    adj = 0.5,
    cex = cex_left,
    xpd = NA
  )
  graphics::text(
    x_right,
    y_lab,
    labels = favours[[2L]],
    adj = 0.5,
    cex = cex_right,
    xpd = NA
  )
  invisible(NULL)
}

#' Axis title for metafor::forest (Favours are drawn separately)
#' @noRd
forestplot_metafor_xlab <- function(xlab, favours = NULL) {
  force(favours) # accepted for backward-compatible callers; ignored
  if (is.null(xlab) || !length(xlab) || all(!nzchar(as.character(xlab)))) {
    return("")
  }
  as.character(xlab)[[1L]]
}

#' @noRd
forestplot_metafor_header <- function(title, slab_label = "Study") {
  right <- if (is.null(title) || !length(title)) {
    "Estimate"
  } else {
    paste(as.character(unlist(title)), collapse = " ")
  }
  right <- gsub("\n", " ", right, fixed = TRUE)
  c(slab_label, right)
}

#' Max number of lines in character labels that may contain `\\n`
#' @noRd
forestplot_metafor_header_lines <- function(labs) {
  if (is.null(labs) || !length(labs)) {
    return(1L)
  }
  n <- vapply(
    as.character(labs),
    function(s) {
      if (is.na(s) || !nzchar(s)) {
        1L
      } else {
        length(strsplit(s, "\n", fixed = TRUE)[[1L]])
      }
    },
    integer(1L)
  )
  max(1L, max(n))
}

#' Redraw column headers above the rule (bottom-aligned, no clipping)
#'
#' metafor centres `ilab.lab` / `header` one row above the rule, which clips
#' multiline labels. Wipe that chrome and redraw with `adj` y = 0 so text sits
#' fully above a rule that clears the first content / trial heading.
#' @noRd
forestplot_metafor_redraw_column_headers <- function(ylim,
                                                    top,
                                                    ylim_hi,
                                                    content_gap = 1,
                                                    header,
                                                    ilab.lab = NULL,
                                                    ilab.xpos = NULL,
                                                    plot_xlim,
                                                    textpos = NULL,
                                                    cex = 1,
                                                    draw_rule = TRUE) {
  usr <- graphics::par("usr")
  y_line_default <- ylim[[2L]] - (top - 1)
  y_line <- ylim_hi + content_gap
  # Keep the new rule below metafor's header text so the wipe is clean.
  y_line <- min(y_line, y_line_default - 0.05)
  # Bottom-align headers with a clear gap above the rule (do not use pos=,
  # which vertically centres text and lets it touch the hline).
  y_gap <- if (isTRUE(draw_rule)) {
    max(0.40, 0.35 * cex)
  } else {
    max(0.25, 0.25 * cex)
  }
  y_text <- y_line + y_gap

  # Wipe any residual metafor header ink (full width, generous vertical band).
  graphics::rect(
    usr[[1L]] - diff(usr[1:2]) * 0.02,
    y_line_default - 0.35,
    usr[[2L]] + diff(usr[1:2]) * 0.02,
    ylim[[2L]] + 1.5,
    col = "white",
    border = NA,
    xpd = NA
  )
  # Inset rule — avoid endpoint ticks at the clip edge that look like
  # orphaned dashes beside the Study / Effect headers.
  if (isTRUE(draw_rule)) {
    x0 <- forestplot_metafor_safe_label_x(textpos)
    x1 <- if (
      !is.null(textpos) && length(textpos) >= 2L && is.finite(textpos[[2L]])
    ) {
      textpos[[2L]]
    } else {
      plot_xlim[[2L]] - 0.015 * diff(plot_xlim)
    }
    graphics::segments(
      x0,
      y_line,
      x1,
      y_line,
      lty = 1,
      xpd = FALSE
    )
  }

  op <- graphics::par(cex = cex, font = 2)
  on.exit(graphics::par(op), add = TRUE)

  # Small inset so labels are not flush against the plot edge.
  x_pad <- 0.02 * diff(plot_xlim)
  x_left <- if (
    !is.null(textpos) && length(textpos) >= 1L && is.finite(textpos[[1L]])
  ) {
    textpos[[1L]]
  } else {
    plot_xlim[[1L]] + x_pad
  }
  x_right <- if (
    !is.null(textpos) && length(textpos) >= 2L && is.finite(textpos[[2L]])
  ) {
    textpos[[2L]] - x_pad
  } else {
    plot_xlim[[2L]] - x_pad
  }

  # Study (left) — bottom-left; inset so the first glyph is not clipped.
  if (length(header) >= 1L && nzchar(header[[1L]])) {
    graphics::text(
      forestplot_metafor_safe_label_x(textpos),
      y_text,
      labels = header[[1L]],
      adj = c(0, 0),
      xpd = NA
    )
  }
  # Estimate / Effect (right) — bottom-right at textpos[2]
  if (length(header) >= 2L && nzchar(header[[2L]])) {
    graphics::text(
      x_right,
      y_text,
      labels = header[[2L]],
      adj = c(1, 0),
      xpd = TRUE
    )
  }
  if (!is.null(ilab.lab) && length(ilab.lab) && !is.null(ilab.xpos)) {
    n <- min(length(ilab.lab), length(ilab.xpos))
    for (l in seq_len(n)) {
      if (!nzchar(ilab.lab[[l]])) {
        next
      }
      graphics::text(
        ilab.xpos[[l]],
        y_text,
        labels = ilab.lab[[l]],
        adj = c(0.5, 0),
        xpd = TRUE
      )
    }
  }
  invisible(NULL)
}

#' @noRd
forestplot_metafor_ilab <- function(panel, left_names) {
  cols <- intersect(as.character(left_names), names(panel))
  if (!length(cols)) {
    return(NULL)
  }
  mat <- as.data.frame(panel[cols], stringsAsFactors = FALSE)
  for (nm in names(mat)) {
    mat[[nm]] <- as.character(mat[[nm]])
    mat[[nm]][is.na(mat[[nm]])] <- ""
  }
  as.matrix(mat)
}

#' Character expansion for metafor forest text / markers
#'
#' Scales down when the active figure height is too short for the number of
#' layout rows so study and pooled labels do not paint on top of each other.
#' @noRd
forestplot_metafor_cex <- function(base_size,
                                   n_layout_rows = NULL,
                                   device_height_in = NULL) {
  bs <- suppressWarnings(as.numeric(base_size)[1L])
  if (!is.finite(bs) || bs <= 0) {
    bs <- 9
  }
  cex <- max(0.55, min(1.15, bs / 11))
  n <- suppressWarnings(as.integer(n_layout_rows)[1L])
  if (!is.finite(n) || n < 1L) {
    return(cex)
  }
  dh <- if (is.null(device_height_in)) {
    forestplot_device_inches()[[2L]]
  } else {
    suppressWarnings(as.numeric(device_height_in)[1L])
  }
  if (!is.finite(dh) || dh <= 0) {
    dh <- 12
  }
  # Usable vertical inches for row content (minus approximate axis / header).
  usable <- max(0.75, dh - 1.35)
  inches_per_row <- usable / n
  # Target ~0.22 in per row at cex = 1; shrink when rows are packed tighter.
  cex <- min(cex, max(0.45, inches_per_row / 0.22))
  cex
}

#' Horizontal positions for left annotation columns (`ilab`)
#'
#' Shifts the information-fraction column slightly left of metafor's default
#' (`0.65` of the slab-to-forest gap) so it sits closer to the study labels.
#'
#' @noRd
forestplot_metafor_ilab_xpos <- function(xlim,
                                         alim,
                                         yi,
                                         ci.lb,
                                         n_cols,
                                         left_names = character()) {
  xlim <- as.numeric(xlim)
  forest_left <- min(c(yi, ci.lb, alim), na.rm = TRUE)
  if (!is.finite(forest_left)) {
    forest_left <- alim[[1L]]
  }
  dist <- forest_left - xlim[[1L]]
  if (!is.finite(dist) || dist <= 0) {
    dist <- diff(xlim) * 0.35
  }

  # metafor defaults: 1 col 0.75; 2 cols c(0.65, 0.85); 3 cols c(0.6, 0.75, 0.9)
  fracs <- if (n_cols <= 1L) {
    0.75
  } else if (n_cols == 2L) {
    c(0.65, 0.85)
  } else if (n_cols == 3L) {
    c(0.6, 0.75, 0.9)
  } else {
    seq(0.5, 0.9, length.out = n_cols)
  }

  # Nudge the ifrac column left (toward study labels).
  if ("ifrac" %in% left_names) {
    ix <- match("ifrac", left_names)
    if (!is.na(ix) && ix <= length(fracs)) {
      fracs[[ix]] <- fracs[[ix]] - 0.08
    }
  }
  fracs <- pmax(0.35, pmin(0.95, fracs))
  xlim[[1L]] + dist * fracs
}

#' Wider plot-region xlim so slab / annotations fit (metafor convention)
#' @noRd
forestplot_metafor_plot_xlim <- function(alim,
                                         has_ilab = FALSE,
                                         annotate = TRUE) {
  alim <- as.numeric(alim)
  span <- diff(range(alim, na.rm = TRUE))
  if (!is.finite(span) || span <= 0) {
    span <- 1
  }
  # Generous left pad for study labels / ilab; keep enough forest width so
  # the plot does not look left-heavy. Right pad must fit the Effect \[CI\]
  # column (header + values) without clipping — but stay close to `alim` so
  # there is no large empty band between the forest and Effect text.
  left_pad <- if (isTRUE(has_ilab)) 4.2 * span else 2.8 * span
  right_pad <- if (isTRUE(annotate)) 2.5 * span else 0.45 * span
  left_pad <- max(left_pad, if (isTRUE(has_ilab)) 3.2 else 2.0)
  right_pad <- max(right_pad, if (isTRUE(annotate)) 1.9 else 0.35)
  c(alim[[1L]] - left_pad, alim[[2L]] + right_pad)
}

#' Left / right text columns for slabs and Effect annotations
#'
#' Nudges the right (Effect) column left of `xlim[2]` so the header and CI
#' text clear the plot / device edge. Optionally nudges the left (study)
#' column inward so the first character is not clipped at the device edge.
#' @noRd
forestplot_metafor_textpos <- function(xlim,
                                       right_shift = 0.06,
                                       left_shift = 0) {
  xlim <- as.numeric(xlim)
  span <- diff(xlim)
  if (!is.finite(span) || span <= 0) {
    span <- 1
  }
  rs <- suppressWarnings(as.numeric(right_shift)[1L])
  if (!is.finite(rs) || rs <= 0) {
    rs <- 0
  }
  ls <- suppressWarnings(as.numeric(left_shift)[1L])
  if (!is.finite(ls) || ls <= 0) {
    ls <- 0
  }
  shift_r <- max(0, min(0.22, rs)) * span
  shift_l <- max(0, min(0.20, ls)) * span
  c(xlim[[1L]] + shift_l, xlim[[2L]] - shift_r)
}
