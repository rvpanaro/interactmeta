
#' Blank repeated group headings so trials (or marginal groups) appear once
#'
#' Consecutive identical non-empty headings keep only the first occurrence,
#' which groups forest rows visually by study / trial (or by Younger / Older
#' for pooled–prediction blocks).
#' @noRd
forestplot_group_headings <- function(heading) {
  heading <- as.character(heading)
  heading[is.na(heading)] <- ""
  if (!length(heading)) {
    return(heading)
  }
  out <- heading
  if (length(heading) > 1L) {
    same_as_prev <- heading[-1L] == heading[-length(heading)] & nzchar(heading[-1L])
    out[c(FALSE, same_as_prev)] <- ""
  }
  out
}

#' Compose forest row labels from heading + subgroup
#'
#' When `group_by_heading = TRUE` (default for display), the trial / group name
#' is shown once at the start of each consecutive block and subgroup labels are
#' left-aligned under a fixed-width trial column.
#' @noRd
forestplot_compose_row_label <- function(heading, label, group_by_heading = FALSE) {
  heading <- as.character(heading)
  label <- as.character(label)
  heading[is.na(heading)] <- ""
  label[is.na(label)] <- ""

  if (!isTRUE(group_by_heading)) {
    return(ifelse(
      nzchar(heading) & nzchar(label),
      paste(heading, label, sep = " / "),
      ifelse(nzchar(heading), heading, label)
    ))
  }

  heading_show <- forestplot_group_headings(heading)
  # Pad to the full heading width so subgroup-only rows stay aligned
  width <- max(1L, max(nchar(heading, type = "width"), na.rm = TRUE))
  trial_col <- ifelse(
    nzchar(heading_show),
    format(heading_show, width = width, justify = "left"),
    strrep(" ", width)
  )
  ifelse(
    nzchar(label),
    paste0(trial_col, "  ", label),
    trimws(trial_col, which = "right")
  )
}

#' Marker / null-line colours (default forest chrome)
#'
#' Summary diamonds use blue via
#' `forestplot_summary_diamond_colour()` (same for every block role and for
#' `group = "subgroup"` / `group = "trial"`). `summary` here is only a fallback.
#' @noRd
forestplot_style_colours <- function(...) {
  const <- interactmeta_plot_constants()
  list(
    marker = comparison_forest_marker_colour(),
    summary = forestplot_default_summary_colour(),
    null_line = const$axis_line,
    axis = const$axis_line
  )
}

#' Default summary diamond colour (palette navy)
#'
#' Used for all Overall / Prediction / Combined diamonds. Not subgroup identity.
#' @noRd
forestplot_default_summary_colour <- function() {
  interactmeta_plot_constants()$forest_summary
}

#' Reference-subgroup summary diamond colour
#' @noRd
forestplot_reference_summary_colour <- forestplot_default_summary_colour

#' Non-reference subgroup summary diamond colour
#' @noRd
forestplot_nonref_summary_colour <- forestplot_default_summary_colour

#' Interaction-block summary diamond colour
#' @noRd
forestplot_interaction_summary_colour <- forestplot_default_summary_colour

#' Combined-block summary diamond colour
#' @noRd
forestplot_overall_summary_colour <- forestplot_default_summary_colour

#' Summary diamond / prediction-bar colour for a forest block heading
#'
#' All summary diamonds use palette navy, independent of `group`
#' and block role (subgroups, Interaction, Combined).
#'
#' @param heading Block heading (unused for colour; kept for call compatibility).
#' @param reference Optional reference subgroup label (unused for colour).
#' @param subgroup_levels Unused; kept for call compatibility.
#' @param palette Unused; kept for call compatibility.
#' @noRd
forestplot_summary_diamond_colour <- function(heading,
                                              reference = NULL,
                                              subgroup_levels = NULL,
                                              palette = NULL) {
  force(heading)
  force(reference)
  force(subgroup_levels)
  force(palette)
  forestplot_default_summary_colour()
}

#' Vectorised `forestplot_summary_diamond_colour()` over headings
#' @noRd
forestplot_summary_diamond_colours <- function(headings,
                                               reference = NULL,
                                               subgroup_levels = NULL,
                                               palette = NULL) {
  headings <- as.character(headings %||% character())
  vapply(
    headings,
    function(h) {
      forestplot_summary_diamond_colour(
        heading = h,
        reference = reference,
        subgroup_levels = subgroup_levels,
        palette = palette
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

#' forestplot drawer that forces a summary fill / line colour
#' @noRd
forestplot_summary_draw_fn <- function(col, bar = FALSE) {
  fill <- as.character(col %||% forestplot_nonref_summary_colour())[[1L]]
  force(fill)
  if (isTRUE(bar)) {
    function(..., col, y.offset = 0.5) {
      forestplot::fpDrawBarCI(..., y.offset = 0.5, col = fill)
    }
  } else {
    function(..., col, y.offset = 0.5) {
      forestplot::fpDrawSummaryCI(..., y.offset = 0.5, col = fill)
    }
  }
}

#' Study CI end-caps off; modest diamond height for [metafor::forest()]
#' @noRd
forestplot_default_efac <- function() {
  c(0, 1.5)
}

#' Modest diamond height for [metafor::addpoly()]
#' @noRd
forestplot_default_poly_efac <- function() {
  1.5
}

#' Bold fonts for addpoly labels and annotations
#'
#' metafor stores `fonts` as a named length-2 vector: names are font families,
#' values are font faces (`1` plain, `2` bold). Unspecified `fonts` inherit
#' plain faces from the forest plot, so bold must be set explicitly.
#' @noRd
forestplot_default_poly_fonts <- function(family = NULL) {
  if (is.null(family)) {
    family <- graphics::par("family")
  }
  family <- as.character(family)[[1L]]
  if (!length(family) || is.na(family)) {
    family <- ""
  }
  stats::setNames(c(2L, 2L), c(family, family))
}

#' Default / validate axis favour labels
#'
#' For `type = "interaction"`, defaults name the non-reference and reference
#' subgroups (same left/right roles as Control / Experimental). Other types
#' keep `c("Control", "Experimental")`. Pass length-2 labels or `FALSE` to
#' override.
#'
#' @noRd
resolve_forestplot_favours <- function(favours,
                                       ...,
                                       exponentiate = FALSE,
                                       type = NULL,
                                       reference = NULL,
                                       nonref_subgroup = NULL) {
  force(exponentiate)
  if (isFALSE(favours)) {
    return(NULL)
  }
  if (!is.null(favours)) {
    favours <- as.character(favours)
    if (length(favours) != 2L) {
      stop("`favours` must be length 2, FALSE, or NULL.", call. = FALSE)
    }
    return(favours)
  }
  if (identical(as.character(type %||% "")[[1L]], "interaction")) {
    int_fav <- forestplot_interaction_favours(
      reference = reference,
      nonref_subgroup = nonref_subgroup
    )
    if (!is.null(int_fav)) {
      # Leading spaces nudge the right-hand label (same as Experimental).
      return(c(int_fav[[1L]], paste0("   ", int_fav[[2L]])))
    }
  }
  # Leading spaces nudge "Experimental" right under the forest column
  # (bayesmeta / forestplot favours band alignment).
  c("Control", "   Experimental")
}

#' Interaction-plot favour labels (subgroup comparison)
#'
#' Left / right match Control / Experimental: reference on the left (negative
#' contrast), non-reference on the right (positive contrast).
#' @noRd
forestplot_interaction_favours <- function(reference = NULL,
                                           nonref_subgroup = NULL) {
  ref <- as.character(reference %||% "")[[1L]]
  nonref <- as.character(nonref_subgroup %||% "")[[1L]]
  if (!nzchar(ref) || !nzchar(nonref)) {
    return(NULL)
  }
  c(ref, nonref)
}

#' Reference / first non-reference subgroup labels from a CAMS fit
#' @noRd
forestplot_fit_reference_nonref <- function(fit) {
  if (interactmeta_is_non_cams_analysis(fit) && !is.null(fit$reference)) {
    return(list(
      reference = as.character(fit$reference)[[1L]],
      nonref = as.character(fit$nonref %||% "")[[1L]]
    ))
  }
  st <- attr(fit, "interactmeta", exact = TRUE)
  empty <- list(reference = "", nonref = "")
  if (is.null(st)) {
    return(empty)
  }
  ref <- as.character(st$reference %||% "")[[1L]]
  if (is.na(ref)) {
    ref <- ""
  }
  nonref <- ""
  ifrac_orig <- st$coef_map$ifrac_orig %||% st$ifrac_vars
  if (length(ifrac_orig)) {
    nonref <- sub("^ifrac_", "", as.character(ifrac_orig)[[1L]])
  } else if (length(st$subgroup_levels)) {
    others <- setdiff(as.character(st$subgroup_levels), ref)
    if (length(others)) {
      nonref <- others[[1L]]
    }
  }
  list(reference = ref, nonref = nonref)
}

#' Default markers: black study squares, navy summary diamonds
#'
#' Summary diamonds use `interactmeta_plot_constants()$forest_summary`.
#' Shrinkage stays `forest_shrinkage` (gray45).
#' @noRd
forestplot_apply_default_colours <- function(resultsA, resultsB) {
  styled <- forestplot_apply_journal_colours(
    resultsA,
    resultsB,
    marker = comparison_forest_marker_colour(),
    diamond_studies = FALSE,
    square_studies = TRUE
  )
  paint_summary <- function(df) {
    if (!nrow(df) || !"row_kind" %in% names(df)) {
      return(df)
    }
    ix <- which(df$row_kind %in% c("overall", "prediction", "combined"))
    if (!length(ix)) {
      return(df)
    }
    col <- forestplot_default_summary_colour()
    if (!"colour" %in% names(df)) {
      df$colour <- rep(NA_character_, nrow(df))
    }
    if (!"fill" %in% names(df)) {
      df$fill <- rep(NA_character_, nrow(df))
    }
    ok <- ix[!is.na(df$shape[ix])]
    if (length(ok)) {
      df$colour[ok] <- col
      df$fill[ok] <- col
    }
    df
  }
  list(
    resultsA = paint_summary(styled$resultsA),
    resultsB = paint_summary(styled$resultsB)
  )
}

#' @noRd
identical_row_kind_shrinkage <- function(row_kind) {
  rk <- as.character(row_kind)
  !is.na(rk) & rk == "shrinkage"
}

#' Shared journal-style colour / shape paint
#' @noRd
forestplot_apply_journal_colours <- function(resultsA,
                                             resultsB,
                                             marker = comparison_forest_marker_colour(),
                                             diamond_studies = FALSE,
                                             square_studies = FALSE) {
  paint <- function(df) {
    n <- nrow(df)
    if (!n) {
      return(df)
    }
    if (!"colour" %in% names(df)) {
      df$colour <- rep(NA_character_, n)
    }
    if (!"fill" %in% names(df)) {
      df$fill <- rep(NA_character_, n)
    }
    paint_ix <- if ("shape" %in% names(df)) {
      !is.na(df$shape)
    } else if ("plot_yi" %in% names(df)) {
      is.finite(df$yi) | is.finite(df$plot_yi)
    } else {
      is.finite(df$yi)
    }
    # Shrinkage keeps forest_shrinkage from comparison_forest_apply_shrinkage_colours();
    # pooled / prediction / overall diamonds stay black.
    if ("row_kind" %in% names(df)) {
      rk <- as.character(df$row_kind)
      paint_ix <- paint_ix &
        !identical_row_kind_shrinkage(rk) &
        !(rk %in% c("overall", "prediction", "combined") & !is.na(rk))
    }
    if ("shape" %in% names(df)) {
      shapes <- as.integer(df$shape)
      study_like <- paint_ix & shapes %in% c(15L, 16L, 19L, 20L, 21L, 22L, 23L, 24L, 25L)
      if (isTRUE(diamond_studies)) {
        # BMJ uses filled diamonds for study and summary markers
        df$shape[study_like] <- 18L
      } else if (isTRUE(square_studies)) {
        # metafor: filled squares for studies; leave diamonds (18) for summaries
        df$shape[study_like & shapes != 18L] <- 15L
      }
    }
    df$colour[paint_ix] <- marker
    shapes <- if ("shape" %in% names(df)) {
      as.integer(df$shape)
    } else {
      rep(22L, n)
    }
    # Solid pch 15–20 use colour only; 21–25 use fill + border colour
    solid <- paint_ix & shapes %in% c(15L, 16L, 17L, 18L, 19L, 20L)
    bordered <- paint_ix & !solid
    df$fill[solid] <- marker
    df$fill[bordered] <- if (
      identical(marker, "black") ||
        identical(marker, "#000000") ||
        identical(marker, comparison_forest_marker_colour())
    ) {
      "white"
    } else {
      marker
    }
    df
  }
  list(resultsA = paint(resultsA), resultsB = paint(resultsB))
}

# ---- forestplot_layout.R ----
#' Recommended figure dimensions for a forest plot
#'
#' Returns width and height in inches scaled to the number of layout rows
#' (studies, pooled / prediction diamonds, subgroup headers), plus footnote
#' lines and axis / favours chrome. Full forests use about **0.28 inches per
#' content row**, **0.14 inches per footnote line**, and ~2.3 inches of chrome,
#' capped at 24 inches.
#'
#' Height follows row count; the **auto-recommended width** (when `width =
#' NULL`) is at least **2:1 width-to-height** so exports stay landscape.
#' Passing an explicit `width` keeps that width unchanged (height still follows
#' content). Set `width_height = FALSE` to skip widening the auto width.
#'
#' Use the returned sizes with `png()` / `pdf()` or **knitr** chunk options
#' (`fig.width`, `fig.height`). Compute dimensions *before* opening the device
#' or in a prior chunk, then pass them into the export / chunk options — the
#' active device height is what prevents row overlap / cut-off.
#'
#' @param x A `forestplot()` result.
#' @param width Target figure width in inches. When `NULL`, uses a size derived
#'   from the number of rows and panels (including the width-to-height rule).
#' @param width_height Target **width / height** ratio for the auto-recommended
#'   width. Default `NULL` uses [forestplot_default_width_height_ratio()] (`2`).
#'   `FALSE` disables widening.
#'
#' @return A list with elements `width`, `height`, `aspect` (width / height;
#'   `2` is a 2:1 landscape figure),
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
forestplot_height <- function(x, width = 9.5, width_height = NULL) {
  dims <- forestplot_figure_inches(x, width_height = width_height)
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
#' @keywords internal
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

#' @param ... Ignored.
#' @rdname forestplot_draw
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

#' Default width-to-height ratio for recommended forest export size
#' @export
forestplot_default_width_height_ratio <- function() {
  2
}

#' Widen `width` so width / height meets `width_height` (never shrinks height)
#' @noRd
forestplot_fit_width_to_height <- function(width, height, width_height = NULL) {
  if (identical(width_height, FALSE)) {
    return(width)
  }
  wh <- width_height %||% forestplot_default_width_height_ratio()
  wh <- suppressWarnings(as.numeric(wh)[1L])
  if (!is.finite(wh) || wh <= 0) {
    return(width)
  }
  max(as.numeric(width)[1L], as.numeric(height)[1L] * wh)
}

#' Figure size in inches for forestplot_height()
#' @noRd
forestplot_figure_inches <- function(x, width_height = NULL) {
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
    # +1 row of slack so the last Combined / Prediction diamond is not clipped.
    height <- row_pitch * (n_rows + 1L) + foot_pitch * n_footnotes + chrome
    height <- max(3.8, min(28, height))
    formula <- sprintf(
      "height = 0.28*(%d+1) + 0.14*%d + %.1f = %.2f",
      n_rows, n_footnotes, chrome, height
    )
    width0 <- width
    width <- forestplot_fit_width_to_height(
      width = width,
      height = height,
      width_height = width_height
    )
    if (is.finite(width) && is.finite(width0) && width > width0 + 1e-9) {
      wh <- width_height %||% forestplot_default_width_height_ratio()
      formula <- paste0(
        formula,
        sprintf(
          "; width = max(%.2f, %.2g*%.2f) = %.2f",
          width0, wh, height, width
        )
      )
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
#' size correctly when Combined includes both Overall and Prediction.
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
