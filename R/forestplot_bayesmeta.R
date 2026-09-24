#' Draw a bayesmeta / forestplot-backed forest on the active device
#'
#' CAMS forests always use [forestplot::forestplot()] so every `type`
#' (`"subgroups"`, `"interaction"`, `"overall"`) shares the same aesthetics:
#' black study boxes; blue summary diamonds for every block;
#' prediction bars; optional grey shrinkage overlays; and IF columns.
#' (Standalone [as_bayesmeta()] objects can still be drawn with
#' [bayesmeta::forestplot.bayesmeta].)
#'
#' @noRd
forestplot_bayesmeta_draw <- function(x) {
  check_forestplot_deps()
  specs <- x$specs
  if (is.null(specs) || !is.list(specs$panels) || !length(specs$panels)) {
    stop("Forest plot is missing bayesmeta draw specifications.", call. = FALSE)
  }

  n <- length(specs$panels)
  roles <- as.character(specs$panel_roles %||% rep("subgroups", n))
  if (length(roles) < n) {
    roles <- rep_len(roles, n)
  }
  foot_lines <- forestplot_bayesmeta_footnote_rows(specs)
  fps <- vector("list", n)
  row_weights <- numeric(n)
  for (i in seq_len(n)) {
    panel_i <- specs$panels[[i]]
    rl_i <- forestplot_bayesmeta_panel_row_labels(specs, i)
    blocks <- forestplot_bayesmeta_blocks(
      panel_i,
      rl_i,
      estimate_col = specs$estimate_col
    )
    specs_i <- forestplot_bayesmeta_panel_specs(specs, i)
    fp_args <- forestplot_bayesmeta_forestplot_args(blocks, specs_i)
    fp <- forestplot_bayesmeta_build_fp(fp_args)
    fav_i <- forestplot_bayesmeta_panel_favours(specs, i)
    fp <- forestplot_bayesmeta_apply_favours(fp, fav_i)
    fps[[i]] <- fp
    row_weights[[i]] <- forestplot_bayesmeta_panel_row_weight(panel_i, rl_i)
  }

  if (identical(n, 1L)) {
    forestplot_bayesmeta_print_with_footnotes(fps[[1L]], foot_lines)
  } else {
    forestplot_bayesmeta_print_panels(
      fps,
      row_weights,
      foot_lines,
      roles = roles
    )
  }
  invisible(x)
}

#' @noRd
forestplot_bayesmeta_panel_row_labels <- function(specs, i) {
  rl <- specs$row_labels
  if (is.data.frame(rl)) {
    return(rl)
  }
  if (is.list(rl) && length(rl) >= i) {
    return(rl[[i]])
  }
  stop("Forest plot is missing row labels for panel ", i, ".", call. = FALSE)
}

#' @noRd
forestplot_bayesmeta_panel_specs <- function(specs, i) {
  out <- specs
  pick <- function(x) {
    if (is.list(x) && !is.data.frame(x) && length(x) >= i) {
      x[[i]]
    } else {
      x
    }
  }
  out$xlim <- pick(specs$xlim)
  out$xticks <- pick(specs$xticks)
  out$xlab <- pick(specs$xlab)
  out$col_right_heading <- pick(specs$col_right_heading)
  out$row_labels <- forestplot_bayesmeta_panel_row_labels(specs, i)
  out$panels <- list(specs$panels[[i]])
  # Secondary panels keep IF/N off via drop_left_on_secondary already.
  if (i > 1L) {
    out$show_left <- FALSE
    out$left_names <- character()
    out$left_headings <- character()
  }
  out$panel_role <- {
    roles <- specs$panel_roles %||% "subgroups"
    as.character(roles)[[min(i, length(roles))]]
  }
  out
}

#' @noRd
forestplot_bayesmeta_panel_favours <- function(specs, i) {
  fav <- specs$favours
  if (is.null(fav) || isFALSE(fav)) {
    return(NULL)
  }
  if (is.character(fav) && length(fav) == 2L) {
    return(fav)
  }
  if (is.list(fav) && length(fav) >= i) {
    return(fav[[i]])
  }
  NULL
}

#' @noRd
forestplot_bayesmeta_panel_row_weight <- function(panel, row_labels) {
  n <- max(1L, nrow(panel))
  n_headers <- 0L
  if (is.data.frame(row_labels) && "heading" %in% names(row_labels)) {
    h <- unique(as.character(row_labels$heading))
    h <- h[nzchar(h) & !is.na(h)]
    n_headers <- length(h)
  }
  # Header row + content + chrome for axis / favours.
  as.numeric(n + n_headers + 3L)
}

#' Print several forest panels with optional footnotes
#'
#' Stacks panels vertically in one column. Footnotes sit under the stack.
#' `roles` is accepted for caller compatibility but ignored.
#' @noRd
forestplot_bayesmeta_print_panels <- function(fps,
                                              row_weights,
                                              foot_lines,
                                              roles = NULL) {
  n <- length(fps)
  if (!n) {
    return(invisible(NULL))
  }
  row_weights <- as.numeric(row_weights)
  if (length(row_weights) != n) {
    row_weights <- rep(1, n)
  }
  row_weights[!is.finite(row_weights) | row_weights <= 0] <- 1

  forestplot_bayesmeta_with_interactmeta_favors({
    grid::grid.newpage()
    forestplot_bayesmeta_print_stack_layout(fps, row_weights, foot_lines)
  })
  invisible(fps)
}

#' Spacing for annotation footnotes under a forestplot panel
#'
#' Returns layout units: gap above the footnote band, band height, and the
#' per-line pitch used when placing text inside the band. Pitch is kept tight
#' (near the footnote `cex`) so consecutive annotation lines have little gap.
#' @noRd
forestplot_bayesmeta_footnote_layout <- function(n_foot) {
  n_foot <- max(0L, as.integer(n_foot)[[1L]])
  line_pitch <- 0.82
  list(
    # Negative gap tucks the band up toward the favours / axis title.
    gap = grid::unit(-0.45, "lines"),
    band = grid::unit(n_foot * line_pitch + 0.02, "lines"),
    line_pitch = line_pitch
  )
}

#' Coerce footnote payload to a list of character rows
#' @noRd
forestplot_footnote_as_rows <- function(foot_lines) {
  if (is.null(foot_lines) || !length(foot_lines)) {
    return(list())
  }
  if (is.list(foot_lines) && !is.character(foot_lines)) {
    lapply(foot_lines, function(x) {
      x <- as.character(x)
      x[!is.na(x) & nzchar(x)]
    })
  } else {
    as.list(as.character(foot_lines))
  }
}

#' Draw footnote rows in the current viewport
#'
#' Notes on a row (IF, heterogeneity, interaction heterogeneity) are joined
#' with `" \u00b7 "` and drawn as one centred line.
#' @noRd
forestplot_bayesmeta_draw_footnote_texts <- function(foot_lines, line_pitch = 0.82) {
  rows <- forestplot_footnote_as_rows(foot_lines)
  gp <- grid::gpar(
    cex = 0.72,
    fontface = "plain",
    col = interactmeta_plot_constants()$forest_text
  )
  for (j in seq_along(rows)) {
    cells <- rows[[j]]
    cells <- as.character(cells)
    cells <- cells[!is.na(cells) & nzchar(cells)]
    if (!length(cells)) {
      next
    }
    y <- grid::unit(1, "npc") +
      grid::unit(0.28, "lines") -
      grid::unit((j - 1) * line_pitch, "lines")
    grid::grid.text(
      forestplot_footnote_as_plotmath(paste(cells, collapse = " \u00b7 ")),
      x = grid::unit(0.5, "npc"),
      y = y,
      just = c("center", "bottom"),
      gp = gp
    )
  }
  invisible(NULL)
}

#' Vertical stack of forest panels (multi-panel layout)
#' @noRd
forestplot_bayesmeta_print_stack_layout <- function(fps, row_weights, foot_lines) {
  n <- length(fps)
  n_foot <- length(foot_lines)
  units <- list()
  for (i in seq_len(n)) {
    units[[length(units) + 1L]] <- grid::unit(row_weights[[i]], "null")
    if (i < n) {
      units[[length(units) + 1L]] <- grid::unit(0.45, "lines")
    }
  }
  if (n_foot) {
    foot_lay <- forestplot_bayesmeta_footnote_layout(n_foot)
    units[[length(units) + 1L]] <- foot_lay$gap
    units[[length(units) + 1L]] <- foot_lay$band
  }
  heights <- do.call(grid::unit.c, units)
  lay <- grid::grid.layout(nrow = length(heights), ncol = 1L, heights = heights)
  grid::pushViewport(grid::viewport(layout = lay))
  row_i <- 1L
  for (i in seq_len(n)) {
    grid::pushViewport(grid::viewport(layout.pos.row = row_i, layout.pos.col = 1L))
    print(fps[[i]], new_page = FALSE)
    grid::popViewport()
    row_i <- row_i + 1L
    if (i < n) {
      row_i <- row_i + 1L
    }
  }
  if (n_foot) {
    row_i <- row_i + 1L
    grid::pushViewport(grid::viewport(
      layout.pos.row = row_i,
      layout.pos.col = 1L,
      name = "interactmetaFootnotes"
    ))
    forestplot_bayesmeta_draw_footnote_texts(foot_lines, foot_lay$line_pitch)
    grid::popViewport()
  }
  grid::popViewport()
  invisible(NULL)
}

#' Build a forestplot object without flashing a blank graphics device
#'
#' [forestplot::forestplot()] calls [grid::grid.newpage()] when no device is
#' open (`dev.cur() == 1`), which opens a visible empty window before we print.
#' Construct on an invisible null PDF device in that case, then draw once.
#' @noRd
forestplot_bayesmeta_build_fp <- function(fp_args) {
  if (!is.list(fp_args)) {
    fp_args <- as.list(fp_args)
  }
  # We control paging in forestplot_bayesmeta_print_with_footnotes().
  fp_args$new_page <- FALSE
  opened_null <- FALSE
  # dev.cur() is a named integer; identical(..., 1L) is FALSE because of names.
  if (!forestplot_has_graphics_device()) {
    grDevices::pdf(NULL)
    opened_null <- TRUE
  }
  on.exit(
    {
      if (opened_null && forestplot_has_graphics_device()) {
        try(grDevices::dev.off(), silent = TRUE)
      }
    },
    add = TRUE
  )
  fp <- do.call(forestplot::forestplot, fp_args)
  if (opened_null) {
    grDevices::dev.off()
    opened_null <- FALSE
  }
  forestplot_bayesmeta_add_graph_separator(
    fp,
    graph.pos = fp_args$graph.pos
  )
}

#' Vertical rule at the left of the forest graph (label | forest)
#' @noRd
forestplot_bayesmeta_add_graph_separator <- function(fp, graph.pos) {
  if (is.null(fp) || is.null(graph.pos) || is.character(graph.pos)) {
    return(fp)
  }
  graph.pos <- as.integer(graph.pos)[[1L]]
  if (!is.finite(graph.pos) || graph.pos < 1L) {
    return(fp)
  }
  if (!is.function(get0("fp_add_lines", envir = asNamespace("forestplot"), inherits = FALSE))) {
    return(fp)
  }
  v_name <- paste0("v_", graph.pos)
  args <- list(x = fp)
  args[[v_name]] <- grid::gpar(col = interactmeta_plot_constants()$forest_rule)
  do.call(forestplot::fp_add_lines, args)
}

#' `Pooled` / `Prediction` labels with a subgroup / block heading in parentheses
#'
#' Used when the heading is not drawn alone (summary-only subgroup blocks or
#' `group = "trial"`), so the block remains identifiable on the summary row.
#'
#' @param allow_block_headings If `TRUE`, also attach Interaction / Overall
#'   headings (e.g. `"Pooled (Overall)"`). Default skips those block names.
#' @noRd
forestplot_summary_group_label <- function(label,
                                           heading,
                                           attach_heading = FALSE,
                                           allow_block_headings = FALSE) {
  lab <- as.character(label %||% "")[[1L]]
  if (!isTRUE(attach_heading)) {
    return(lab)
  }
  if (
    !identical(lab, "Pooled") &&
      !identical(lab, "Prediction")
  ) {
    return(lab)
  }
  hdr <- as.character(heading %||% "")[[1L]]
  if (!nzchar(hdr)) {
    return(lab)
  }
  if (
    !isTRUE(allow_block_headings) &&
      (forestplot_heading_is_interaction(hdr) ||
         forestplot_heading_is_overall(hdr))
  ) {
    return(lab)
  }
  paste0(lab, " (", hdr, ")")
}

#' @rdname forestplot_summary_group_label
#' @noRd
forestplot_trial_group_summary_label <- function(label,
                                                 heading,
                                                 group = "subgroup") {
  forestplot_summary_group_label(
    label = label,
    heading = heading,
    attach_heading = identical(group, "trial")
  )
}

#' Attach left / right Favours labels under the forest graph
#'
#' Uses [forestplot::fp_set_favors()] when available so forestplot reserves
#' axis space for the favours band. Drawing itself is handled by
#' [forestplot_interactmeta_plot_graph_favors()] (wider arrow gap, non-overlapping
#' labels). Falls back to embedding in `xlab` only if `fp_set_favors` is missing.
#' @noRd
forestplot_bayesmeta_apply_favours <- function(fp, favours) {
  if (is.null(fp) || is.null(favours) || length(favours) < 2L) {
    return(fp)
  }
  low <- as.character(favours[[1L]])
  high <- as.character(favours[[2L]])
  if (!nzchar(low) && !nzchar(high)) {
    return(fp)
  }
  if (is.function(get0("fp_set_favors", envir = asNamespace("forestplot"), inherits = FALSE))) {
    return(
      forestplot::fp_set_favors(
        fp,
        low = low,
        high = high,
        arrows = TRUE,
        position = "outside",
        txt_gp = grid::gpar(cex = 0.65, col = interactmeta_plot_constants()$forest_text)
      )
    )
  }
  # Older forestplot: append to xlab (may look shifted for long labels).
  xlab <- fp$xlab %||% ""
  xlab <- as.character(xlab)[[1L]]
  fp$xlab <- paste0(xlab, "\n", low, "                              ", high)
  fp
}

#' Favours band with a center gap and equal-sized labels
#'
#' Drop-in replacement for [forestplot:::plotGraphFavors()]: arrows leave a
#' clear gap at the null. Left / right labels share one cex (the larger of the
#' two per-side fits) so Control / Experimental stay the same size.
#' @noRd
forestplot_interactmeta_plot_graph_favors <- function(favors, txt_gp, col) {
  gp_list <- as.list(txt_gp)
  gp_list$col <- interactmeta_plot_constants()$forest_text
  if (!is.null(favors$txt_gp)) {
    for (n in names(favors$txt_gp)) {
      gp_list[[n]] <- favors$txt_gp[[n]]
    }
  }

  cex <- gp_list$cex
  if (is.null(cex) || !is.numeric(cex) || is.na(cex[[1L]]) || cex[[1L]] <= 0) {
    cex <- 0.65
  } else {
    cex <- cex[[1L]]
  }

  low <- as.character(favors$low %||% "")[[1L]]
  high <- as.character(favors$high %||% "")[[1L]]

  # Keep a clear band around the null (npc fraction of the forest column).
  center_gap <- 0.26
  left_limit <- 0.5 - center_gap / 2
  right_limit <- 0.5 + center_gap / 2
  left_tip_npc <- 0.02
  right_tip_npc <- 0.98
  left_base_npc <- left_limit
  right_base_npc <- right_limit

  # Fit each label inside its half, then use the larger cex for both.
  fit_cex <- function(label, max_npc, start_cex) {
    c <- start_cex
    while (c > 0.4) {
      gp_list$cex <- c
      w <- tryCatch(
        as.numeric(grid::convertWidth(
          grid::grobWidth(grid::textGrob(label, gp = do.call(grid::gpar, gp_list))),
          "npc",
          valueOnly = TRUE
        )),
        error = function(e) NA_real_
      )
      if (!is.finite(w) || w <= max_npc) {
        return(c)
      }
      c <- c * 0.9
    }
    c
  }
  max_w_left <- max(0.12, left_limit - left_tip_npc - 0.01)
  max_w_right <- max(0.12, right_tip_npc - right_limit - 0.01)
  cex_left <- fit_cex(low, max_w_left, cex)
  cex_right <- fit_cex(high, max_w_right, cex)
  shared_cex <- max(cex_left, cex_right)

  arrow_y <- if (identical(favors$position, "outside")) {
    # Sit a bit below mid-band so arrows clear the axis ticks / xlab.
    outside_arrow_y <- 0.52 - 0.18 * (shared_cex - 0.65)
    grid::unit(min(max(outside_arrow_y, 0.38), 0.62), "npc")
  } else {
    grid::unit(0.2, "npc")
  }

  label_y_nudge <- favors$label_y_nudge
  if (is.null(label_y_nudge)) {
    label_y_nudge <- grid::unit(0.32, "line")
  } else if (!grid::is.unit(label_y_nudge)) {
    label_y_nudge <- grid::unit(as.numeric(label_y_nudge)[[1L]], "mm")
  }
  label_y <- if (isTRUE(favors$arrows)) {
    arrow_y + label_y_nudge
  } else {
    grid::unit(0.35, "npc") + label_y_nudge
  }

  # Same size on both sides.
  gp_left <- gp_list
  gp_left$cex <- shared_cex
  gp_right <- gp_list
  gp_right$cex <- shared_cex
  grid::grid.text(
    low,
    x = grid::unit(left_tip_npc, "npc"),
    y = label_y,
    hjust = 0,
    vjust = 0,
    gp = do.call(grid::gpar, gp_left)
  )
  grid::grid.text(
    high,
    x = grid::unit(right_tip_npc, "npc"),
    y = label_y,
    hjust = 1,
    vjust = 0,
    gp = do.call(grid::gpar, gp_right)
  )

  if (isTRUE(favors$arrows)) {
    arrow_gp <- list(col = gp_list$col)
    if (!is.null(gp_list$lwd)) {
      arrow_gp$lwd <- gp_list$lwd
    }
    if (!is.null(favors$arrow_gp)) {
      for (n in names(favors$arrow_gp)) {
        arrow_gp[[n]] <- favors$arrow_gp[[n]]
      }
    }
    if (is.null(arrow_gp$lwd)) {
      arrow_gp$lwd <- 0.65
    }
    segments_gp <- do.call(grid::gpar, arrow_gp)
    arrow_spec <- grid::arrow(
      type = "open",
      ends = "last",
      length = grid::unit(1.5, "mm")
    )
    grid::grid.segments(
      x0 = grid::unit(left_base_npc, "npc"),
      y0 = arrow_y,
      x1 = grid::unit(left_tip_npc, "npc"),
      y1 = arrow_y,
      gp = segments_gp,
      arrow = arrow_spec
    )
    grid::grid.segments(
      x0 = grid::unit(right_base_npc, "npc"),
      y0 = arrow_y,
      x1 = grid::unit(right_tip_npc, "npc"),
      y1 = arrow_y,
      gp = segments_gp,
      arrow = arrow_spec
    )
  }
}

#' Temporarily use interactmeta favours drawing inside forestplot
#' @noRd
forestplot_bayesmeta_with_interactmeta_favors <- function(expr) {
  ns <- asNamespace("forestplot")
  if (!exists("plotGraphFavors", envir = ns, inherits = FALSE)) {
    return(force(expr))
  }
  orig <- get("plotGraphFavors", envir = ns, inherits = FALSE)
  was_locked <- bindingIsLocked("plotGraphFavors", ns)
  if (was_locked) {
    unlockBinding("plotGraphFavors", ns)
  }
  assign("plotGraphFavors", forestplot_interactmeta_plot_graph_favors, envir = ns)
  on.exit(
    {
      assign("plotGraphFavors", orig, envir = ns)
      if (was_locked) {
        lockBinding("plotGraphFavors", ns)
      }
    },
    add = TRUE
  )
  force(expr)
}

#' Print a forestplot object with footnotes below the axis / favours band
#'
#' Always starts with [grid::grid.newpage()]. The no-footnote path used to
#' call `print(fp)` alone; without clearing the device, a follow-up plot with
#' `show_annotations = FALSE` overlaid the previous forest (ghosted / double
#' drawing) while leftover footnotes from the prior draw stayed visible.
#' @noRd
forestplot_bayesmeta_print_with_footnotes <- function(fp, foot_lines) {
  if (is.null(fp)) {
    return(invisible(NULL))
  }
  n_foot <- length(foot_lines)

  forestplot_bayesmeta_with_interactmeta_favors({
    grid::grid.newpage()
    if (!n_foot) {
      # Pass new_page through forestplot's print() so it does not clear again.
      print(fp, new_page = FALSE)
    } else {
      # Row 1: forest; row 2: overlap; row 3: footnotes (cdot-separated).
      foot_lay <- forestplot_bayesmeta_footnote_layout(n_foot)
      heights <- grid::unit.c(
        grid::unit(1, "null"),
        foot_lay$gap,
        foot_lay$band
      )
      lay <- grid::grid.layout(nrow = 3L, ncol = 1L, heights = heights)
      grid::pushViewport(grid::viewport(layout = lay))
      grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
      print(fp, new_page = FALSE)
      grid::popViewport()
      grid::pushViewport(grid::viewport(
        layout.pos.row = 3L,
        layout.pos.col = 1L,
        name = "interactmetaFootnotes"
      ))
      forestplot_bayesmeta_draw_footnote_texts(foot_lines, foot_lay$line_pitch)
      grid::popViewport()
      grid::popViewport()
    }
  })
  invisible(fp)
}

#' @noRd
forestplot_bayesmeta_draw_footnotes_only <- function(foot_lines) {
  rows <- forestplot_footnote_as_rows(foot_lines)
  if (!length(rows)) {
    return(invisible(NULL))
  }
  # Fallback when forestplot.bayesmeta already drew on the device: try to add
  # text in the bottom margin without wiping the figure.
  tryCatch(
    {
      grid::seekViewport("axis_margin")
      forestplot_bayesmeta_draw_footnote_texts(rows)
    },
    error = function(e) invisible(NULL)
  )
  invisible(NULL)
}

# --- bayesmeta object construction ------------------------------------------

#' Build a minimal bayesmeta object from study / summary estimates
#'
#' Enough fields for [bayesmeta::forestplot.bayesmeta] (not a full fitted MA).
#'
#' @param y Numeric study estimates (log scale).
#' @param sigma Numeric study standard errors.
#' @param labels Character study labels.
#' @param mu Named numeric with `median`, `lower`, `upper`, and optional `sd`
#'   for the pooled mean.
#' @param prediction Named numeric with `median`, `lower`, `upper`, and optional
#'   `sd` for the prediction interval (`theta` in bayesmeta).
#' @param shrinkage Optional data frame / matrix with columns (or rows)
#'   `median`, `lower`, `upper` (and optional `sd`) — one row per study.
#' @param tau Optional named numeric with `median`, `lower`, `upper` for the
#'   heterogeneity footnote.
#'
#' @return An object of class `"bayesmeta"`.
#' @export
as_bayesmeta <- function(y,
                         sigma,
                         labels,
                         mu,
                         prediction = NULL,
                         shrinkage = NULL,
                         tau = NULL) {
  y <- as.numeric(y)
  sigma <- as.numeric(sigma)
  labels <- as.character(labels)
  k <- length(y)
  if (length(sigma) != k || length(labels) != k) {
    stop("`y`, `sigma`, and `labels` must have the same length.", call. = FALSE)
  }
  if (k < 1L) {
    stop("Need at least one study estimate.", call. = FALSE)
  }

  pick <- function(x, what) {
    if (is.null(x)) {
      return(NA_real_)
    }
    nm <- names(x)
    if (!is.null(nm) && what %in% nm) {
      return(as.numeric(x[[what]]))
    }
    # aliases
    alt <- switch(
      what,
      lower = c("lower", "lci", "ci.lb", "95% lower"),
      upper = c("upper", "uci", "ci.ub", "95% upper"),
      median = c("median", "estimate", "yi", "mean"),
      sd = c("sd", "se", "sigma"),
      what
    )
    hit <- intersect(alt, nm)
    if (length(hit)) {
      return(as.numeric(x[[hit[[1L]]]]))
    }
    NA_real_
  }

  summary <- matrix(
    NA_real_,
    nrow = 6L,
    ncol = 3L,
    dimnames = list(
      c("mode", "median", "mean", "sd", "95% lower", "95% upper"),
      c("tau", "mu", "theta")
    )
  )
  for (par in c("tau", "mu", "theta")) {
    src <- switch(par, tau = tau, mu = mu, theta = prediction)
    summary["median", par] <- pick(src, "median")
    summary["mean", par] <- pick(src, "median")
    summary["mode", par] <- pick(src, "median")
    summary["sd", par] <- pick(src, "sd")
    summary["95% lower", par] <- pick(src, "lower")
    summary["95% upper", par] <- pick(src, "upper")
  }
  # Prediction defaults to mu if omitted (bayesmeta always has theta column).
  if (is.null(prediction) || !is.finite(summary["median", "theta"])) {
    summary[, "theta"] <- summary[, "mu"]
  }
  if (is.null(tau) || !is.finite(summary["median", "tau"])) {
    summary[, "tau"] <- c(0, 0, 0, 0, 0, 0)
  }

  theta <- matrix(
    NA_real_,
    nrow = 8L,
    ncol = k,
    dimnames = list(
      c("y", "sigma", "mode", "median", "mean", "sd", "95% lower", "95% upper"),
      labels
    )
  )
  theta["y", ] <- y
  theta["sigma", ] <- sigma
  if (!is.null(shrinkage)) {
    sh <- as.data.frame(shrinkage)
    # Accept bayesmeta-style column names or ci.lb/ci.ub.
    med <- if ("median" %in% names(sh)) {
      sh$median
    } else if ("yi" %in% names(sh)) {
      sh$yi
    } else if ("estimate" %in% names(sh)) {
      sh$estimate
    } else {
      sh[[1L]]
    }
    lo <- if ("95% lower" %in% names(sh)) {
      sh[["95% lower"]]
    } else if ("lower" %in% names(sh)) {
      sh$lower
    } else if ("ci.lb" %in% names(sh)) {
      sh$ci.lb
    } else if ("lci" %in% names(sh)) {
      sh$lci
    } else {
      rep(NA_real_, k)
    }
    hi <- if ("95% upper" %in% names(sh)) {
      sh[["95% upper"]]
    } else if ("upper" %in% names(sh)) {
      sh$upper
    } else if ("ci.ub" %in% names(sh)) {
      sh$ci.ub
    } else if ("uci" %in% names(sh)) {
      sh$uci
    } else {
      rep(NA_real_, k)
    }
    sd_sh <- if ("sd" %in% names(sh)) sh$sd else rep(NA_real_, k)
    theta["median", ] <- as.numeric(med)
    theta["mean", ] <- as.numeric(med)
    theta["mode", ] <- as.numeric(med)
    theta["sd", ] <- as.numeric(sd_sh)
    theta["95% lower", ] <- as.numeric(lo)
    theta["95% upper", ] <- as.numeric(hi)
  } else {
    # No shrinkage: mirror quoted estimates so forestplot.bayesmeta can still
    # run with shrinkage = FALSE.
    theta["median", ] <- y
    theta["mean", ] <- y
    theta["mode", ] <- y
    theta["sd", ] <- sigma
    q <- stats::qnorm(0.975)
    theta["95% lower", ] <- y - q * sigma
    theta["95% upper", ] <- y + q * sigma
  }

  out <- list(
    y = y,
    sigma = sigma,
    labels = labels,
    k = k,
    summary = summary,
    theta = theta,
    call = match.call()
  )
  class(out) <- "bayesmeta"
  out
}

#' @noRd
forestplot_as_bayesmeta_block <- function(block) {
  studies <- block$studies
  if (!nrow(studies)) {
    stop("Block has no study rows for as_bayesmeta().", call. = FALSE)
  }
  se <- forestplot_bayesmeta_se(studies)
  mu <- block$pooled
  pred <- block$prediction
  shr <- block$shrinkage
  as_bayesmeta(
    y = studies$yi,
    sigma = se,
    labels = studies$label,
    mu = mu,
    prediction = pred,
    shrinkage = shr,
    tau = block$tau
  )
}

# --- block parsing ----------------------------------------------------------

#' @noRd
forestplot_bayesmeta_blocks <- function(panel, row_labels, estimate_col = "yi") {
  n <- nrow(panel)
  if (!n) {
    return(list())
  }
  est <- panel[[estimate_col]]
  if (is.null(est)) {
    est <- panel$yi
  }
  panel$.yi <- as.numeric(est)
  panel$.ci.lb <- as.numeric(panel$ci.lb)
  panel$.ci.ub <- as.numeric(panel$ci.ub)
  if ("vi" %in% names(panel)) {
    panel$.vi <- as.numeric(panel$vi)
  } else {
    q <- stats::qnorm(0.975)
    panel$.vi <- ((panel$.ci.ub - panel$.ci.lb) / (2 * q))^2
  }

  heading <- if (is.data.frame(row_labels) && "heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    rep("", n)
  }
  label <- if (is.data.frame(row_labels) && "label" %in% names(row_labels)) {
    as.character(row_labels$label)
  } else {
    rep("", n)
  }
  heading[is.na(heading)] <- ""
  label[is.na(label)] <- ""
  # Carry headings forward within a block.
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

  # Unique headings in order of first appearance.
  heads <- heading
  heads[!nzchar(heads)] <- "(Overall)"
  uniq <- unique(heads)

  lapply(uniq, function(h) {
    ii <- which(heads == h)
    kind_i <- row_kind[ii]
    lab_i <- label[ii]
    pan_i <- panel[ii, , drop = FALSE]

    study_ii <- which(kind_i %in% c("study", "shrinkage"))
    # Pair study + following shrinkage into one bayesmeta study row.
    studies <- list()
    shrinkage_rows <- list()
    j <- 1L
    while (j <= length(study_ii)) {
      idx <- study_ii[[j]]
      if (identical(kind_i[[idx]], "shrinkage")) {
        j <- j + 1L
        next
      }
      st <- data.frame(
        yi = pan_i$.yi[[idx]],
        ci.lb = pan_i$.ci.lb[[idx]],
        ci.ub = pan_i$.ci.ub[[idx]],
        vi = pan_i$.vi[[idx]],
        label = lab_i[[idx]],
        stringsAsFactors = FALSE
      )
      # IF / N for labeltext extras
      if ("ifrac" %in% names(pan_i)) {
        st$ifrac <- pan_i$ifrac[[idx]]
      }
      if ("trial_size" %in% names(pan_i)) {
        st$trial_size <- pan_i$trial_size[[idx]]
      }
      sh <- NULL
      if (j < length(study_ii)) {
        nxt <- study_ii[[j + 1L]]
        if (identical(kind_i[[nxt]], "shrinkage")) {
          sh <- data.frame(
            median = pan_i$.yi[[nxt]],
            lower = pan_i$.ci.lb[[nxt]],
            upper = pan_i$.ci.ub[[nxt]],
            stringsAsFactors = FALSE
          )
          j <- j + 1L
        }
      }
      studies[[length(studies) + 1L]] <- st
      shrinkage_rows[[length(shrinkage_rows) + 1L]] <- sh
      j <- j + 1L
    }
    studies_df <- if (length(studies)) {
      dplyr::bind_rows(studies)
    } else {
      data.frame(
        yi = numeric(), ci.lb = numeric(), ci.ub = numeric(), vi = numeric(),
        label = character(), stringsAsFactors = FALSE
      )
    }
    shr_df <- if (length(shrinkage_rows) && !all(vapply(shrinkage_rows, is.null, logical(1)))) {
      dplyr::bind_rows(lapply(shrinkage_rows, function(z) {
        if (is.null(z)) {
          data.frame(median = NA_real_, lower = NA_real_, upper = NA_real_)
        } else {
          z
        }
      }))
    } else {
      NULL
    }

    as_sum_at <- function(i0) {
      list(
        median = pan_i$.yi[[i0]],
        lower = pan_i$.ci.lb[[i0]],
        upper = pan_i$.ci.ub[[i0]],
        sd = sqrt(max(pan_i$.vi[[i0]], 0)),
        label = {
          lab0 <- lab_i[[i0]]
          if (!nzchar(lab0) && "subgroup" %in% names(pan_i)) {
            lab0 <- as.character(pan_i$subgroup[[i0]])
            if (is.na(lab0)) lab0 <- ""
          }
          lab0
        },
        bar = {
          k0 <- kind_i[[i0]]
          lab0 <- lab_i[[i0]]
          if (!nzchar(lab0) && "subgroup" %in% names(pan_i)) {
            lab0 <- as.character(pan_i$subgroup[[i0]])
            if (is.na(lab0)) lab0 <- ""
          }
          identical(k0, "prediction") ||
            grepl("^prediction$", lab0, ignore.case = TRUE)
        }
      )
    }

    # All summary diamonds/bars in panel order (keeps Overall Pooled + Prediction).
    summaries <- list()
    for (idx in seq_along(kind_i)) {
      if (kind_i[[idx]] %in% c("pooled", "prediction", "overall")) {
        summaries[[length(summaries) + 1L]] <- as_sum_at(idx)
      }
    }

    # Convenience slots (first of each kind) for simple bayesmeta conversion.
    pooled <- NULL
    prediction <- NULL
    overall <- NULL
    for (s in summaries) {
      lab_l <- tolower(s$label %||% "")
      if (isTRUE(s$bar) || identical(lab_l, "prediction")) {
        if (is.null(prediction)) prediction <- s
      } else if (identical(lab_l, "pooled") || identical(lab_l, "mean")) {
        if (is.null(pooled)) pooled <- s
      } else {
        if (is.null(overall)) overall <- s
        if (is.null(pooled) && !grepl("overall", lab_l, fixed = TRUE)) {
          pooled <- s
        }
      }
    }
    # Subgroup blocks: kinds are pooled/prediction without relying on labels.
    if (is.null(pooled)) {
      hit <- which(kind_i == "pooled")
      if (length(hit)) pooled <- as_sum_at(hit[[1L]])
    }
    if (is.null(prediction)) {
      hit <- which(kind_i == "prediction")
      if (length(hit)) prediction <- as_sum_at(hit[[1L]])
    }
    if (is.null(overall)) {
      hit <- which(kind_i == "overall")
      if (length(hit)) overall <- as_sum_at(hit[[1L]])
    }

    list(
      heading = if (identical(h, "(Overall)")) "" else h,
      studies = studies_df,
      shrinkage = shr_df,
      summaries = summaries,
      pooled = pooled,
      prediction = prediction,
      overall = overall,
      tau = NULL
    )
  })
}

#' Whether a block looks like a plain bayesmeta forest (studies + mean [+ pred])
#'
#' Used by [as_bayesmeta()] helpers. CAMS drawing always goes through
#' [forestplot::forestplot()] so aesthetics match across `type` values.
#' @noRd
forestplot_bayesmeta_block_is_simple <- function(block) {
  nrow(block$studies) >= 1L &&
    !is.null(block$pooled) &&
    is.null(block$overall) &&
    length(block$summaries) <= 2L
}

#' @noRd
forestplot_bayesmeta_se <- function(studies) {
  se <- sqrt(pmax(as.numeric(studies$vi), 0))
  bad <- !is.finite(se) | se <= 0
  if (any(bad)) {
    q <- stats::qnorm(0.975)
    se[bad] <- pmax(
      (as.numeric(studies$ci.ub[bad]) - as.numeric(studies$ci.lb[bad])) / (2 * q),
      .Machine$double.eps
    )
  }
  se
}

#' Keep at least one non-summary row for forestplot vertical lines
#'
#' [forestplot::forestplot()] builds a default vertical rule with
#' `rows = which(!is.summary)[1]:n`. When every row is a summary (header,
#' spacers, diamonds), that start is `NA` and print fails with
#' `NA/NaN argument`. Flip the first empty spacer (or the first row)
#' so Pooled / Prediction-only forests still draw.
#'
#' @noRd
forestplot_ensure_nonsummary_row <- function(is_summary, mean_m) {
  is_summary <- as.logical(is_summary)
  if (!length(is_summary) || any(!is_summary, na.rm = TRUE)) {
    return(is_summary)
  }
  mean_vec <- if (is.matrix(mean_m)) {
    as.numeric(mean_m[, 1L])
  } else {
    as.numeric(mean_m)
  }
  n <- length(is_summary)
  if (length(mean_vec) < n) {
    mean_vec <- c(mean_vec, rep(NA_real_, n - length(mean_vec)))
  }
  # Prefer a spacer under the header so the vertical rule does not
  # start in the title row.
  idx <- if (n > 1L) seq.int(2L, n) else 1L
  spacer_i <- idx[is_summary[idx] & !is.finite(mean_vec[idx])]
  i <- if (length(spacer_i)) spacer_i[[1L]] else idx[[1L]]
  is_summary[[i]] <- FALSE
  is_summary
}

# --- forestplot::forestplot argument builders -------------------------------

#' @noRd
forestplot_bayesmeta_call_args <- function(specs, block) {
  has_shr <- !is.null(block$shrinkage) && nrow(block$shrinkage) &&
    any(is.finite(block$shrinkage$median))
  args <- list(
    exponentiate = isTRUE(specs$exponentiate),
    prediction = !is.null(block$prediction),
    shrinkage = has_shr,
    heterogeneity = FALSE,
    digits = specs$digits %||% 2,
    plot = TRUE,
    zero = specs$nullval %||% if (isTRUE(specs$exponentiate)) 1 else 0
  )
  xlab <- specs$xlab[[1L]]
  if (!is.null(xlab) && !identical(xlab, FALSE)) {
    args$xlab <- as.character(xlab)[[1L]]
  }
  # Prefer compact column header over a large centered title.
  args$title <- NULL
  # clip / xticks only when the user supplied them (bayesmeta default: omit).
  xlim <- if (is.list(specs$xlim)) specs$xlim[[1L]] else specs$xlim
  if (is.numeric(xlim) && length(xlim) == 2L && all(is.finite(xlim))) {
    args$clip <- xlim
    xt <- if (is.list(specs$xticks)) specs$xticks[[1L]] else specs$xticks
    if (is.numeric(xt) && length(xt)) {
      args$xticks <- xt
    }
  } else {
    xt <- if (is.list(specs$xticks)) specs$xticks[[1L]] else specs$xticks
    if (is.numeric(xt) && length(xt)) {
      args$xticks <- xt
    }
  }
  # User overrides from forest_args / ...
  fa <- specs$forest_args %||% list()
  # Drop metafor-only names.
  drop <- c(
    "atransf", "efac", "psize", "pch", "annotate", "header", "ilab",
    "ilab.xpos", "ilab.lab", "top", "annosym", "fonts", "shade", "lty",
    "rows", "ylim", "alim", "textpos", "slab", "vi", "sei", "ci.lb", "ci.ub"
  )
  fa <- fa[setdiff(names(fa), drop)]
  utils::modifyList(args, fa)
}

#' Build forestplot::forestplot args for one or more CAMS blocks
#' @noRd
forestplot_bayesmeta_forestplot_args <- function(blocks, specs) {
  exponentiate <- isTRUE(specs$exponentiate)
  digits <- as.integer(specs$digits %||% 2L)[[1L]]
  show_ci <- isTRUE(specs$show_ci_text)
  show_left <- isTRUE(specs$show_left)
  left_names <- specs$left_names %||% character()
  nullval <- specs$nullval %||% if (exponentiate) 1 else 0
  group <- match.arg(
    as.character(specs$group %||% "subgroup")[[1L]],
    c("subgroup", "trial")
  )
  reference <- specs$reference %||% NULL
  if (is.character(reference) && length(reference)) {
    reference <- as.character(reference)[[1L]]
    if (!nzchar(reference) || is.na(reference)) {
      reference <- NULL
    }
  } else {
    reference <- NULL
  }
  subgroup_levels <- unique(vapply(
    blocks,
    function(b) as.character(b$heading %||% "")[[1L]],
    character(1L)
  ))
  subgroup_levels <- subgroup_levels[
    nzchar(subgroup_levels) &
      !forestplot_heading_is_interaction(subgroup_levels) &
      !forestplot_heading_is_overall(subgroup_levels)
  ]

  label_cols <- list()
  mean_list <- list()
  lower_list <- list()
  upper_list <- list()
  is_summary <- logical()
  boxsize <- numeric()
  # Use a sentinel for "default summary drawer" — list[[i]] <- NULL drops the element.
  fn_sum <- list()
  fn_sum_missing <- list()
  has_any_shrinkage <- FALSE

  push_row <- function(lab_vec, mean_v, lower_v, upper_v, summary, box, sum_fn = NULL) {
    label_cols[[length(label_cols) + 1L]] <<- lab_vec
    mean_list[[length(mean_list) + 1L]] <<- mean_v
    lower_list[[length(lower_list) + 1L]] <<- lower_v
    upper_list[[length(upper_list) + 1L]] <<- upper_v
    is_summary <<- c(is_summary, summary)
    boxsize <<- c(boxsize, box)
    i <- length(fn_sum) + 1L
    if (is.null(sum_fn)) {
      fn_sum[[i]] <<- NA
      fn_sum_missing[[i]] <<- TRUE
    } else {
      fn_sum[[i]] <<- sum_fn
      fn_sum_missing[[i]] <<- FALSE
    }
  }

  push_spacer_row <- function() {
    push_row(
      lab_vec = rep("", n_lab),
      mean_v = c(NA_real_, NA_real_),
      lower_v = c(NA_real_, NA_real_),
      upper_v = c(NA_real_, NA_real_),
      summary = TRUE,
      # Finite boxsize: NA triggers forestplot warnings when shrinkage
      # adds a second estimate layer (max(info[i,], na.rm=TRUE) on all-NA).
      box = 0,
      sum_fn = NULL
    )
  }

  is_subgroup_heading <- function(h) {
    h <- as.character(h %||% "")[[1L]]
    nzchar(h) &&
      !forestplot_heading_is_interaction(h) &&
      !forestplot_heading_is_overall(h)
  }

  # Header row — Study (+ IF / N) and a single Effect column (estimate [CI]).
  # No large centered `title`, and no separate Estimate / 95% CI headers
  # (those duplicated Effect).
  effect_hdr <- {
    crh <- specs$col_right_heading[[1L]]
    if (is.null(crh) || identical(crh, FALSE) || !nzchar(as.character(crh)[[1L]])) {
      "Effect"
    } else {
      # Flatten newlines for the compact column header cell.
      gsub("\n", " ", as.character(crh)[[1L]], fixed = TRUE)
    }
  }
  hdr <- c("Study", effect_hdr)
  extra_hdr <- character()
  if (show_left && length(left_names)) {
    lh <- specs$left_headings %||% left_names
    extra_hdr <- as.character(lh)[seq_along(left_names)]
    hdr <- c("Study", extra_hdr, effect_hdr)
  }
  n_lab <- length(hdr)
  # Index of the effect / CI text column (last column).
  effect_col <- n_lab
  push_row(
    lab_vec = hdr,
    mean_v = c(NA_real_, NA_real_),
    lower_v = c(NA_real_, NA_real_),
    upper_v = c(NA_real_, NA_real_),
    summary = TRUE,
    box = 0,
    sum_fn = NULL
  )
  # Empty row under the column header (all forest types / group modes).
  push_spacer_row()

  fmt <- function(x) {
    if (!is.finite(x)) {
      return("")
    }
    sprintf(paste0("%.", digits, "f"), x)
  }
  fmt_ci <- function(lo, hi) {
    if (!is.finite(lo) || !is.finite(hi)) {
      return("")
    }
    paste0("[", fmt(lo), ", ", fmt(hi), "]")
  }
  maybe_exp <- function(x) {
    if (exponentiate) exp(x) else x
  }

  for (bi in seq_along(blocks)) {
    bl <- blocks[[bi]]
    n_st <- nrow(bl$studies)
    is_int_ov <- forestplot_heading_is_interaction(bl$heading) ||
      forestplot_heading_is_overall(bl$heading)
    # Summary-only subgroup blocks: drop the alone heading and put the name
    # on Pooled / Prediction instead, e.g. "Pooled (Younger)".
    skip_heading <- nzchar(bl$heading) && identical(n_st, 0L) && !is_int_ov
    attach_heading <- identical(group, "trial") || skip_heading

    # Empty row above Interaction / Overall when they follow other blocks
    # (the post-header spacer already covers the case when they lead).
    if (bi > 1L && is_int_ov) {
      push_spacer_row()
    }
    if (nzchar(bl$heading) && !skip_heading) {
      lab <- rep("", n_lab)
      lab[[1L]] <- bl$heading
      push_row(
        lab_vec = lab,
        mean_v = c(NA_real_, NA_real_),
        lower_v = c(NA_real_, NA_real_),
        upper_v = c(NA_real_, NA_real_),
        summary = TRUE,
        box = 0,
        sum_fn = NULL
      )
    }

    has_shr <- !is.null(bl$shrinkage) && nrow(bl$shrinkage) == n_st &&
      any(is.finite(bl$shrinkage$median))
    if (has_shr) {
      has_any_shrinkage <- TRUE
    }

    if (n_st) {
      for (i in seq_len(n_st)) {
        yi <- bl$studies$yi[[i]]
        lo <- bl$studies$ci.lb[[i]]
        hi <- bl$studies$ci.ub[[i]]
        yi_d <- maybe_exp(yi)
        lo_d <- maybe_exp(lo)
        hi_d <- maybe_exp(hi)

        lab <- rep("", n_lab)
        lab[[1L]] <- bl$studies$label[[i]]
        # Optional left columns (IF, N)
        col_i <- 2L
        if (show_left && length(left_names)) {
          for (nm in left_names) {
            val <- if (nm %in% names(bl$studies)) bl$studies[[nm]][[i]] else ""
            lab[[col_i]] <- if (is.null(val) || is.na(val)) "" else as.character(val)
            col_i <- col_i + 1L
          }
        }
        if (show_ci) {
          lab[[effect_col]] <- paste0(fmt(yi_d), " ", fmt_ci(lo_d, hi_d))
        }

        if (has_shr) {
          sy <- maybe_exp(bl$shrinkage$median[[i]])
          slo <- maybe_exp(bl$shrinkage$lower[[i]])
          shi <- maybe_exp(bl$shrinkage$upper[[i]])
          push_row(
            lab_vec = lab,
            mean_v = c(yi_d, sy),
            lower_v = c(lo_d, slo),
            upper_v = c(hi_d, shi),
            summary = FALSE,
            box = 0.25,
            sum_fn = NULL
          )
        } else {
          push_row(
            lab_vec = lab,
            mean_v = c(yi_d, NA_real_),
            lower_v = c(lo_d, NA_real_),
            upper_v = c(hi_d, NA_real_),
            summary = FALSE,
            box = 0.25,
            sum_fn = NULL
          )
        }
      }
    }

    add_summary <- function(sum_obj, default_lab = NULL, bar = NULL) {
      if (is.null(sum_obj)) {
        return()
      }
      yi_d <- maybe_exp(as.numeric(sum_obj[["median"]]))
      lo_d <- maybe_exp(as.numeric(sum_obj[["lower"]]))
      hi_d <- maybe_exp(as.numeric(sum_obj[["upper"]]))
      # Switched-off layers keep an empty row so plot height / rhythm match
      # the full (pooled + prediction) layout.
      if (!is.finite(yi_d) || !is.finite(lo_d) || !is.finite(hi_d)) {
        push_spacer_row()
        return()
      }
      lab <- rep("", n_lab)
      lab[[1L]] <- if (!is.null(sum_obj[["label"]]) && nzchar(sum_obj[["label"]])) {
        as.character(sum_obj[["label"]])
      } else if (!is.null(default_lab)) {
        default_lab
      } else {
        "mean"
      }
      # Put the block heading on Pooled / Prediction when it is not drawn alone.
      lab[[1L]] <- forestplot_summary_group_label(
        lab[[1L]],
        heading = bl$heading,
        attach_heading = attach_heading
      )
      if (show_ci) {
        lab[[effect_col]] <- paste0(fmt(yi_d), " ", fmt_ci(lo_d, hi_d))
      }
      use_bar <- if (!is.null(bar)) {
        isTRUE(bar)
      } else {
        isTRUE(sum_obj[["bar"]])
      }
      sum_col <- forestplot_summary_diamond_colour(
        heading = bl$heading,
        reference = reference,
        subgroup_levels = subgroup_levels
      )
      fn <- forestplot_summary_draw_fn(sum_col, bar = use_bar)
      push_row(
        lab_vec = lab,
        mean_v = c(yi_d, NA_real_),
        lower_v = c(lo_d, NA_real_),
        upper_v = c(hi_d, NA_real_),
        summary = TRUE,
        box = if (use_bar) 0.2 else 0.4,
        sum_fn = fn
      )
    }

    if (length(bl$summaries)) {
      for (sm in bl$summaries) {
        add_summary(sm)
      }
    } else {
      add_summary(bl$pooled, "Pooled", bar = FALSE)
      add_summary(bl$prediction, "Prediction", bar = TRUE)
      add_summary(bl$overall, "Overall", bar = FALSE)
    }

    # Empty row between consecutive subgroup blocks (group = "subgroup"),
    # including summary-only blocks that no longer draw a heading row.
    if (
      identical(group, "subgroup") &&
        bi < length(blocks) &&
        is_subgroup_heading(bl$heading) &&
        is_subgroup_heading(blocks[[bi + 1L]]$heading)
    ) {
      push_spacer_row()
    }
  }

  # Footnotes are drawn below the axis in forestplot_bayesmeta_draw(), not as
  # table rows (table rows sit above the x-axis and collide with ticks).

  labeltext <- do.call(rbind, label_cols)
  mean_m <- do.call(rbind, mean_list)
  lower_m <- do.call(rbind, lower_list)
  upper_m <- do.call(rbind, upper_list)

  if (!has_any_shrinkage) {
    mean_m <- mean_m[, 1L, drop = TRUE]
    lower_m <- lower_m[, 1L, drop = TRUE]
    upper_m <- upper_m[, 1L, drop = TRUE]
  }

  # forestplot::getVerticaLines() uses
  # `which(is.summary == FALSE)[1]:number_of_rows`. Summary-only
  # forests (Pooled / Prediction with estimates_only) have no study
  # row, so that start is NA unless a spacer is marked non-summary.
  is_summary <- forestplot_ensure_nonsummary_row(is_summary, mean_m)

  # Drawing functions: study squares; summary colours are baked into
  # each fn.ci_sum closure (role colours). col$summary is only a fallback.
  fp_const <- interactmeta_plot_constants()
  if (has_any_shrinkage) {
    fn.ci_norm <- list(
      function(...) forestplot::fpDrawPointCI(pch = 15, ...),
      function(...) forestplot::fpDrawPointCI(pch = 18, ...)
    )
    col <- forestplot::fpColors(
      box = c(fp_const$forest_study, fp_const$forest_shrinkage),
      lines = c(fp_const$forest_study, fp_const$forest_shrinkage),
      summary = forestplot_overall_summary_colour(),
      text = fp_const$forest_text,
      axes = fp_const$forest_text
    )
    legend <- c("quoted estimate", "shrinkage estimate")
  } else {
    fn.ci_norm <- function(...) forestplot::fpDrawPointCI(pch = 15, ...)
    col <- forestplot::fpColors(
      box = fp_const$forest_study,
      lines = fp_const$forest_study,
      summary = forestplot_overall_summary_colour(),
      text = fp_const$forest_text,
      axes = fp_const$forest_text
    )
    legend <- NULL
  }

  fn.ci_sum <- lapply(seq_along(fn_sum), function(i) {
    if (isTRUE(fn_sum_missing[[i]])) {
      function(y.offset, ...) forestplot::fpDrawSummaryCI(y.offset = 0.5, ...)
    } else {
      fn_sum[[i]]
    }
  })

  # Horizontal rules after header and before first summary in first block —
  # approximate bayesmeta's grey lines.
  hrzl_lines <- list(
    "2" = grid::gpar(col = fp_const$forest_rule)
  )

  xlab <- specs$xlab[[1L]]
  if (is.null(xlab) || identical(xlab, FALSE)) {
    role <- specs$panel_role %||% "subgroups"
    xlab <- forestplot_panel_xlab(role, exponentiate = exponentiate, xlab = NULL)
  } else {
    xlab <- as.character(xlab)[[1L]]
  }
  # Favours are applied after forestplot() via fp_set_favors() so they sit
  # under the left / right of the forest (not glued into xlab with spaces).

  args <- list(
    labeltext = labeltext,
    mean = mean_m,
    lower = lower_m,
    upper = upper_m,
    is.summary = is_summary,
    boxsize = boxsize,
    fn.ci_norm = fn.ci_norm,
    fn.ci_sum = fn.ci_sum,
    col = col,
    legend = legend,
    hrzl_lines = hrzl_lines,
    zero = nullval,
    xlab = xlab,
    # Forest sits before the Effect column (labels | forest | effect text).
    graph.pos = effect_col,
    # Fixed forest-band width so IF / N / Effect toggles do not reflow the
    # graph column when options change.
    graphwidth = grid::unit(55, "mm"),
    txt_gp = forestplot::fpTxtGp(
      label = grid::gpar(
        cex = max(0.55, (specs$base_size %||% 9) / 12),
        col = fp_const$forest_text
      ),
      ticks = grid::gpar(cex = 0.7, col = fp_const$forest_text),
      xlab = grid::gpar(cex = 0.8, col = fp_const$forest_text),
      title = grid::gpar(cex = 0.95, col = fp_const$forest_text)
    )
  )
  # No large centered title — Effect text lives in the header row only
  # (avoids duplicating Estimate / 95% CI + "Effect (95% interval)").
  # clip / xticks: only when supplied (bayesmeta leaves these to forestplot
  # defaults: clip = c(-Inf, Inf), xticks = NULL).
  xlim <- if (is.list(specs$xlim)) specs$xlim[[1L]] else specs$xlim
  if (is.numeric(xlim) && length(xlim) == 2L && all(is.finite(xlim))) {
    args$clip <- as.numeric(xlim)
    z <- args$zero
    if (length(z) == 1L && is.finite(z) && (z < min(xlim) || z > max(xlim))) {
      args$zero <- NA
    }
  }
  xticks <- if (is.list(specs$xticks)) specs$xticks[[1L]] else specs$xticks
  if (is.numeric(xticks) && length(xticks)) {
    args$xticks <- as.numeric(xticks)
  }
  if (exponentiate) {
    args$xlog <- TRUE
    # forestplot logs mean/lower/upper when xlog = TRUE; 0 or negative
    # values become -Inf and getTicks() errors.
    sanitize_pos <- function(x) {
      x[is.finite(x) & x <= 0] <- NA_real_
      x
    }
    args$mean <- sanitize_pos(args$mean)
    args$lower <- sanitize_pos(args$lower)
    args$upper <- sanitize_pos(args$upper)
  }

  fa <- specs$forest_args %||% list()
  drop <- c(
    "atransf", "efac", "psize", "pch", "annotate", "header", "ilab",
    "ilab.xpos", "ilab.lab", "top", "annosym", "fonts", "shade", "lty",
    "rows", "ylim", "alim", "textpos", "slab", "vi", "sei", "ci.lb", "ci.ub",
    "col", "cex", "cex.lab", "cex.axis", "title"
  )
  fa <- fa[setdiff(names(fa), c(drop, ""))]
  out <- utils::modifyList(args, fa)
  # Construction must not open / clear a device; we print once below.
  out$new_page <- FALSE
  out
}

#' Footnote rows for the forest table footer
#'
#' Multi-cell rows (IF + heterogeneity notes) are kept intact so they can be
#' joined with `" \u00b7 "` when drawn. Single-cell rows may wrap.
#' @noRd
forestplot_bayesmeta_footnote_rows <- function(specs, wrap_at = 96L) {
  rows <- forestplot_annotation_footnote_rows(specs$annotations)
  if (!length(rows)) {
    return(list())
  }
  wrap_at <- max(40L, as.integer(wrap_at)[[1L]])
  out <- list()
  for (cells in rows) {
    cells <- as.character(cells)
    cells <- cells[!is.na(cells) & nzchar(cells)]
    if (!length(cells)) {
      next
    }
    if (length(cells) > 1L) {
      out[[length(out) + 1L]] <- cells
    } else {
      wrapped <- strwrap(cells[[1L]], width = wrap_at, simplify = TRUE)
      for (w in wrapped) {
        if (nzchar(w)) {
          out[[length(out) + 1L]] <- w
        }
      }
    }
  }
  out
}

#' Flatten wrapped footnote rows to character lines
#' @noRd
forestplot_bayesmeta_footnote_lines <- function(specs, wrap_at = 96L) {
  rows <- forestplot_bayesmeta_footnote_rows(specs, wrap_at = wrap_at)
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

#' @noRd
forestplot_bayesmeta_draw_footer <- function(specs) {
  invisible(NULL)
}
