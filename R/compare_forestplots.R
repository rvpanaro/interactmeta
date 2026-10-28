# ---- compare_forestplots.R ----

#' Lay out several full forest plots side by side (or in a grid)
#'
#' Each argument must be a **`forestplot()`** result (assign to a variable;
#' do not rely on auto-print). Panels are captured with [grid::grid.grabExpr()]
#' and combined with **patchwork** for knitr / HTML / PDF output.
#'
#' When `align = TRUE` (default), plots share a common effect-axis range and
#' are captured at the same width and height so rows and x-axes line up.
#'
#' For pooled estimates across methods (DA, SWADA, CAMS), use
#' [compare_modifiers()] or [compare_interactions()] instead.
#'
#' @param ... Two or more objects of class `interactmeta_forestplot` from
#'   [forestplot()].
#' @param ncol Number of columns in the layout. Default: all panels in one row.
#' @param nrow Optional number of rows. If `NULL`, inferred from `ncol`.
#' @param widths Optional relative column widths passed to
#'   [patchwork::plot_layout()].
#' @param heights Optional relative row heights passed to
#'   [patchwork::plot_layout()].
#' @param align If `TRUE`, use a shared `xlim` / `xticks` and identical capture
#'   size for each panel (recommended for AD vs DA comparisons).
#' @param labels Optional character vector of panel captions (one label per
#'   forest, in layout order: left-to-right, then top-to-bottom). Rendered as
#'   **patchwork** panel tags (top-left of each panel). `NULL` (default) omits
#'   tags.
#' @param draw If `TRUE`, print the layout on the active graphics device.
#'
#' @return A **patchwork** object (invisibly if `draw = TRUE`).
#'
#' @seealso [forestplot()], [forestplot_draw()], [forestplot_height()],
#'   [compare_modifiers()]
#'
#' @examples
#' \dontrun{
#' fit <- da(my_data, reference = "ref")
#' fp_mod <- forestplot(fit, exponentiate = TRUE)
#' fp_int <- forestplot(fit, type = "interaction", exponentiate = TRUE)
#' compare_forestplots(fp_mod, fp_int, labels = c("Modifiers", "Interaction"))
#' }
#'
#' @export
compare_forestplots <- function(...,
                                ncol = NULL,
                                nrow = NULL,
                                widths = NULL,
                                heights = NULL,
                                align = TRUE,
                                labels = NULL,
                                draw = FALSE) {
  check_forestplot_deps()
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package `patchwork` is required for `compare_forestplots()`. ",
      "Install it with install.packages(\"patchwork\").",
      call. = FALSE
    )
  }

  fps <- list(...)
  if (length(fps) < 2L) {
    stop(
      "`compare_forestplots()` needs at least two `forestplot()` results.",
      call. = FALSE
    )
  }

  for (i in seq_along(fps)) {
    if (!inherits(fps[[i]], "interactmeta_forestplot")) {
      stop(
        sprintf(
          "Argument %d must be an `interactmeta_forestplot` from `forestplot()`.",
          i
        ),
        call. = FALSE
      )
    }
  }

  capture_w <- NULL
  capture_h <- NULL
  if (isTRUE(align)) {
    axis <- forestplot_fps_shared_axis(fps)
    if (!is.null(axis)) {
      fps <- lapply(fps, forestplot_fp_set_axis, xlim = axis$xlim, xticks = axis$xticks)
    }
    dims <- lapply(fps, forestplot_figure_inches)
    capture_w <- max(vapply(dims, `[[`, numeric(1), "width"))
    capture_h <- max(vapply(dims, `[[`, numeric(1), "height"))
    ncols <- if (is.null(ncol) && is.null(nrow)) length(fps) else as.integer(ncol %||% ceiling(length(fps) / as.integer(nrow)))
    if (!is.finite(ncols) || ncols < 1L) {
      ncols <- length(fps)
    }
    capture_w <- capture_w / ncols
  }

  grobs <- lapply(fps, forestplot_panel_grob, width_in = capture_w, height_in = capture_h)

  p <- Reduce(
    `+`,
    lapply(grobs, function(g) patchwork::wrap_elements(full = g))
  )

  n <- length(grobs)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- n
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / as.integer(nrow))
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / as.integer(ncol))
  }

  layout_args <- list(ncol = as.integer(ncol))
  if (!is.null(nrow)) {
    layout_args$nrow <- as.integer(nrow)
  }
  if (!is.null(widths)) {
    layout_args$widths <- widths
  }
  if (!is.null(heights)) {
    layout_args$heights <- heights
  }
  p <- p + do.call(patchwork::plot_layout, layout_args)

  if (!is.null(labels)) {
    labels <- as.character(labels)
    if (length(labels) != n) {
      stop(
        "`labels` must have one caption per forest (length ",
        n,
        "); you passed ",
        length(labels),
        ".",
        call. = FALSE
      )
    }
    p <- p + patchwork::plot_annotation(tag_levels = list(labels))
  }

  if (isTRUE(draw)) {
    print(p)
    return(invisible(p))
  }
  p
}

#' Shared effect-axis limits across several forestplot objects
#' @noRd
forestplot_fps_shared_axis <- function(fps) {
  panels <- lapply(fps, function(fp) {
    if (is.null(fp$specs$panels) || !length(fp$specs$panels)) {
      return(NULL)
    }
    fp$specs$panels[[1L]]
  })
  panels <- panels[!vapply(panels, is.null, logical(1L))]
  if (!length(panels)) {
    return(NULL)
  }
  exponentiate <- forestplot_fp_exponentiate(fps[[1L]])
  combined <- do.call(rbind, panels)
  resolve_forestplot_limits(
    panels = list(combined),
    exponentiate = exponentiate
  )
}

#' Whether a forestplot object uses ratio-scale display
#' @noRd
forestplot_fp_exponentiate <- function(fp) {
  snap <- attr(fp, "forestplot_fit_render", exact = TRUE)
  if (is.list(snap) && "exponentiate" %in% names(snap)) {
    return(isTRUE(snap$exponentiate))
  }
  isTRUE(fp$specs$exponentiate)
}

#' Assign axis limits on a forestplot object before drawing
#' @noRd
forestplot_fp_set_axis <- function(fp, xlim, xticks) {
  if (is.null(fp$specs) || !length(fp$specs$panels)) {
    return(fp)
  }
  n <- length(fp$specs$panels)
  fp$specs$xlim <- rep(list(as.numeric(xlim)), n)
  fp$specs$xticks <- rep(list(as.numeric(xticks)), n)
  fp
}

#' Capture one forest plot as a grid grob for layout
#' @noRd
forestplot_panel_grob <- function(fp, width_in = NULL, height_in = NULL) {
  if (is.finite(width_in) && is.finite(height_in) && width_in > 0 && height_in > 0) {
    grDevices::pdf(NULL, width = width_in, height = height_in)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  grid::grid.grabExpr(forestplot_draw(fp))
}
