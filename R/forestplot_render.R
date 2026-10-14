# Panel assembly / drawing lives in R/forestplot_metafor.R (metafor::forest).
# ggplot helpers below are retained only for label / colour utilities used elsewhere.

#' @noRd
forestplot_ggplot_panel <- function(panel,
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
                                    show_left = TRUE,
                                    left_names = character(),
                                    left_headings = character(),
                                    plot.margin = interactmeta_forestplot_margin(),
                                    favours = NULL,
                                    ...) {
  n <- nrow(panel)
  if (!n) {
    stop("Forest panel data is empty.", call. = FALSE)
  }

  est_raw <- panel[[estimate_col]]
  if (is.null(est_raw)) {
    est_raw <- panel$yi
  }
  lci_raw <- panel$ci.lb
  uci_raw <- panel$ci.ub

  if (isTRUE(exponentiate)) {
    est <- exp(est_raw)
    lci <- exp(lci_raw)
    uci <- exp(uci_raw)
  } else {
    est <- est_raw
    lci <- lci_raw
    uci <- uci_raw
  }

  row_id <- seq_len(n)
  lab_heading <- if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    rep("", n)
  }
  lab_label <- if ("label" %in% names(row_labels)) {
    as.character(row_labels$label)
  } else {
    rep("", n)
  }
  lab_heading[is.na(lab_heading)] <- ""
  lab_label[is.na(lab_label)] <- ""
  # Placeholder blanks (ZWSP / legacy NBSP) must not appear in axis text
  lab_label[comparison_forest_is_blank_row_label(lab_label)] <- ""
  lab_heading[comparison_forest_is_blank_row_label(lab_heading)] <- ""

  y_lab <- forestplot_compose_row_label(lab_heading, lab_label, group_by_heading = TRUE)
  y_lab[!nzchar(y_lab)] <- ""

  if (isTRUE(show_left) && length(left_names)) {
    y_lab <- forestplot_prefix_left_columns(
      y_lab = y_lab,
      panel = panel,
      left_names = left_names
    )
  }

  # Numeric y lets identical labels repeat (study vs shrinkage rows).
  # Row 1 at the top (y = n), last row at the bottom (y = 1).
  y_pos <- (n + 1) - row_id
  y_breaks <- y_pos
  y_labels <- y_lab

  shape <- if ("shape" %in% names(panel) && length(panel$shape) == n) {
    as.integer(panel$shape)
  } else {
    rep(22L, n)
  }
  colour <- if ("colour" %in% names(panel) && length(panel$colour) == n) {
    as.character(panel$colour)
  } else {
    rep(comparison_forest_marker_colour(), n)
  }
  fill <- if ("fill" %in% names(panel) && length(panel$fill) == n) {
    as.character(panel$fill)
  } else {
    rep("white", n)
  }
  stroke_v <- if ("stroke" %in% names(panel) && length(panel$stroke) == n) {
    as.numeric(panel$stroke)
  } else {
    rep(as.numeric(interactmeta_forest_stroke()), n)
  }
  pscale <- if ("point_scale" %in% names(panel) && length(panel$point_scale) == n) {
    ps <- as.numeric(panel$point_scale)
    ps[is.na(ps)] <- 1
    ps
  } else {
    rep(1, n)
  }
  shape[is.na(shape)] <- 22L
  colour[is.na(colour) | !nzchar(colour)] <- comparison_forest_marker_colour()
  fill[is.na(fill) | !nzchar(fill)] <- "white"
  stroke_v[!is.finite(stroke_v)] <- as.numeric(interactmeta_forest_stroke())

  pt_size <- suppressWarnings(as.numeric(pointsize)[1])
  if (!is.finite(pt_size) || pt_size <= 0) {
    pt_size <- as.numeric(interactmeta_forest_pointsize())
  }

  df <- data.frame(
    row_id = row_id,
    y = y_pos,
    estimate = est,
    lci = lci,
    uci = uci,
    shape = shape,
    colour = colour,
    fill = fill,
    stroke = stroke_v,
    pointsize = pscale * pt_size,
    stringsAsFactors = FALSE
  )

  df$ci_text <- forestplot_format_ci_text(
    estimate = df$estimate,
    lci = df$lci,
    uci = df$uci,
    digits = digits
  )

  if (length(xlim) != 2L || any(!is.finite(xlim))) {
    stop("`xlim` must be a length-2 finite numeric vector.", call. = FALSE)
  }
  if (is.null(xticks) || !length(xticks)) {
    xticks <- forestplot_nice_axis(
      xlim,
      exponentiate = FALSE,
      nullval = nullval
    )$xticks
  }

  line_size <- forestplot_base_line_size(base_size)
  const <- interactmeta_plot_constants()
  style_cols <- forestplot_style_colours()
  null_col <- style_cols$null_line
  null_lw <- line_size * 0.7
  p <- ggplot2::ggplot(df, ggplot2::aes(y = .data$y)) +
    ggplot2::geom_vline(
      xintercept = nullval,
      colour = null_col,
      linewidth = null_lw
    ) +
    ggplot2::geom_errorbar(
      data = df[is.finite(df$lci) & is.finite(df$uci), , drop = FALSE],
      ggplot2::aes(
        xmin = .data$lci,
        xmax = .data$uci,
        colour = .data$colour
      ),
      orientation = "y",
      width = 0.35,
      linewidth = line_size,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = df[is.finite(df$estimate), , drop = FALSE],
      ggplot2::aes(
        x = .data$estimate,
        shape = .data$shape,
        colour = .data$colour,
        fill = .data$fill,
        size = .data$pointsize,
        stroke = .data$stroke / 10
      ),
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(
      limits = xlim,
      breaks = xticks,
      expand = c(0, 0),
      trans = if (isTRUE(exponentiate)) "log" else "identity"
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = if (isTRUE(show_left)) y_labels else NULL,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::labs(
      x = if (is.null(xlab) || !nzchar(as.character(xlab)[[1]])) NULL else xlab,
      y = NULL,
      title = title
    ) +
    ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(colour = const$axis_line, linewidth = const$line_hair),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(colour = const$axis_line, linewidth = 0.3),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = base_size * 0.7,
        colour = const$axis_line
      ),
      axis.text.y = if (isTRUE(show_left)) {
        ggplot2::element_text(
          hjust = 1,
          size = base_size * 0.9,
          colour = const$axis_line,
          family = "sans",
          margin = ggplot2::margin(r = 4, unit = "pt")
        )
      } else {
        ggplot2::element_blank()
      },
      axis.title.x = ggplot2::element_text(
        size = base_size * 0.95,
        colour = const$axis_line,
        margin = ggplot2::margin(t = 4, unit = "pt")
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size,
        colour = const$axis_line,
        margin = ggplot2::margin(b = 4, unit = "pt")
      ),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = if (inherits(plot.margin, "margin")) {
        plot.margin
      } else {
        interactmeta_forestplot_margin()
      }
    )

  if (!is.null(favours) && length(favours) >= 2L) {
    p <- forestplot_add_favours(
      p,
      favours = favours,
      xlim = xlim,
      nullval = nullval,
      base_size = base_size,
      exponentiate = exponentiate
    )
  }

  p
}

#' Left meta column(s) for forest plots (IF, N, row labels)
#' @noRd
forestplot_ggplot_meta_panel <- function(panel,
                                         row_labels,
                                         left_names,
                                         left_headings = character(),
                                         base_size = interactmeta_forest_base_size(),
                                         plot.margin = interactmeta_forestplot_margin(),
                                         title_spacer = NULL,
                                         ...) {
  n <- nrow(panel)
  y_pos <- (n + 1) - seq_len(n)
  label_txt <- forestplot_meta_label_strings(panel, row_labels, left_names)

  df <- data.frame(
    y = y_pos,
    label = label_txt,
    stringsAsFactors = FALSE
  )
  const <- interactmeta_plot_constants()
  spacer <- if (!is.null(title_spacer) && nzchar(as.character(title_spacer)[[1]])) {
    as.character(title_spacer)[[1]]
  } else {
    " "
  }

  ggplot2::ggplot(df, ggplot2::aes(x = 0, y = .data$y, label = .data$label)) +
    ggplot2::geom_text(
      hjust = 0,
      size = base_size * 0.9 / ggplot2::.pt,
      colour = const$axis_line,
      family = "mono",
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = y_pos,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = spacer) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size,
        colour = "white",
        margin = ggplot2::margin(b = 4, unit = "pt")
      ),
      plot.margin = if (inherits(plot.margin, "margin")) {
        plot.margin
      } else {
        interactmeta_forestplot_margin()
      }
    )
}

#' Shared left-column strings for forest meta panel width / labels
#' @noRd
forestplot_meta_label_strings <- function(panel, row_labels, left_names) {
  n <- nrow(panel)
  lab_heading <- if ("heading" %in% names(row_labels)) {
    as.character(row_labels$heading)
  } else {
    rep("", n)
  }
  lab_label <- if ("label" %in% names(row_labels)) {
    as.character(row_labels$label)
  } else {
    rep("", n)
  }
  lab_heading[is.na(lab_heading)] <- ""
  lab_label[is.na(lab_label)] <- ""
  lab_label[comparison_forest_is_blank_row_label(lab_label)] <- ""
  lab_heading[comparison_forest_is_blank_row_label(lab_heading)] <- ""
  row_lab <- forestplot_compose_row_label(lab_heading, lab_label, group_by_heading = TRUE)

  cols <- left_names[left_names %in% names(panel)]
  meta_bits <- lapply(cols, function(nm) {
    x <- as.character(panel[[nm]])
    x[is.na(x)] <- ""
    width <- max(1L, max(nchar(x, type = "width"), na.rm = TRUE))
    format(x, width = width, justify = "right")
  })
  if (length(meta_bits)) {
    left_txt <- do.call(paste, c(meta_bits, list(sep = "  ")))
    paste0(left_txt, "   ", row_lab)
  } else {
    row_lab
  }
}

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

#' Prefix y labels with left-column metadata (IF, N) without pipe clutter
#' @noRd
forestplot_prefix_left_columns <- function(y_lab, panel, left_names) {
  out <- y_lab
  for (nm in left_names) {
    if (!nm %in% names(panel)) {
      next
    }
    prefix <- as.character(panel[[nm]])
    prefix[is.na(prefix)] <- ""
    keep <- nzchar(prefix)
    if (!any(keep)) {
      next
    }
    width <- max(nchar(prefix[keep], type = "width"), na.rm = TRUE)
    prefix_fmt <- ifelse(
      keep,
      format(prefix, width = width, justify = "right"),
      strrep(" ", width)
    )
    out <- paste0(prefix_fmt, "  ", out)
    break
  }
  out
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
forestplot_reference_summary_colour <- function() {
  forestplot_default_summary_colour()
}

#' Non-reference subgroup summary diamond colour
#' @noRd
forestplot_nonref_summary_colour <- function() {
  forestplot_default_summary_colour()
}

#' Interaction-block summary diamond colour
#' @noRd
forestplot_interaction_summary_colour <- function() {
  forestplot_default_summary_colour()
}

#' Combined-block summary diamond colour
#' @noRd
forestplot_overall_summary_colour <- function() {
  forestplot_default_summary_colour()
}

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

#' Format effect \[CI\] strings for the right-hand text column
#' @noRd
forestplot_format_ci_text <- function(estimate,
                                      lci,
                                      uci,
                                      digits = 2,
                                      ...) {
  fmt <- function(x) {
    ifelse(
      is.finite(x),
      format(round(x, digits), nsmall = digits, trim = TRUE),
      ""
    )
  }
  ok <- is.finite(estimate) & is.finite(lci) & is.finite(uci)
  out <- rep("", length(estimate))
  # Effect column format: "est [lower, upper]"
  out[ok] <- paste0(fmt(estimate[ok]), " [", fmt(lci[ok]), ", ", fmt(uci[ok]), "]")
  out
}

#' Right-hand effect \[CI\] text column
#' @noRd
forestplot_ggplot_ci_panel <- function(panel,
                                       estimate_col = "yi",
                                       digits = 2,
                                       exponentiate = FALSE,
                                       base_size = interactmeta_forest_base_size(),
                                       plot.margin = interactmeta_forestplot_margin(),
                                       title_spacer = NULL,
                                       ...) {
  n <- nrow(panel)
  y_pos <- (n + 1) - seq_len(n)

  est_raw <- panel[[estimate_col]]
  if (is.null(est_raw)) {
    est_raw <- panel$yi
  }
  lci_raw <- panel$ci.lb
  uci_raw <- panel$ci.ub
  if (isTRUE(exponentiate)) {
    est <- exp(est_raw)
    lci <- exp(lci_raw)
    uci <- exp(uci_raw)
  } else {
    est <- est_raw
    lci <- lci_raw
    uci <- uci_raw
  }

  label_txt <- forestplot_format_ci_text(
    estimate = est,
    lci = lci,
    uci = uci,
    digits = digits
  )

  df <- data.frame(
    y = y_pos,
    label = label_txt,
    stringsAsFactors = FALSE
  )
  const <- interactmeta_plot_constants()
  spacer <- if (!is.null(title_spacer) && nzchar(as.character(title_spacer)[[1]])) {
    as.character(title_spacer)[[1]]
  } else {
    " "
  }

  ggplot2::ggplot(df, ggplot2::aes(x = 1, y = .data$y, label = .data$label)) +
    ggplot2::geom_text(
      hjust = 1,
      size = base_size * 0.88 / ggplot2::.pt,
      colour = const$axis_line,
      family = "mono",
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = y_pos,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = spacer) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size,
        colour = "white",
        margin = ggplot2::margin(b = 4, unit = "pt")
      ),
      plot.margin = if (inherits(plot.margin, "margin")) {
        plot.margin
      } else {
        interactmeta_forestplot_margin()
      }
    )
}

#' Add “Favours …” labels under the forest axis
#' @noRd
forestplot_add_favours <- function(p,
                                   favours,
                                   xlim,
                                   nullval,
                                   base_size,
                                   exponentiate = FALSE) {
  favours <- as.character(favours)
  if (length(favours) < 2L) {
    return(p)
  }
  const <- interactmeta_plot_constants()
  if (isTRUE(exponentiate)) {
    # Geometric midpoints on the ratio scale
    x_left <- exp((log(xlim[[1]]) + log(nullval)) / 2)
    x_right <- exp((log(nullval) + log(xlim[[2]])) / 2)
  } else {
    x_left <- mean(c(xlim[[1]], nullval))
    x_right <- mean(c(nullval, xlim[[2]]))
  }
  if (!is.finite(x_left) || !is.finite(x_right)) {
    return(p)
  }

  # Extra bottom space so labels sit under tick labels
  cur_margin <- p$theme$plot.margin
  if (inherits(cur_margin, "margin")) {
    m <- as.numeric(cur_margin)
    # ggplot margin order: t, r, b, l
    bottom <- max(m[[3]], 8)
    p <- p + ggplot2::theme(
      plot.margin = ggplot2::margin(m[[1]], m[[2]], bottom, m[[4]], unit = "mm")
    )
  }

  p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::annotate(
      "text",
      x = x_left,
      y = -Inf,
      label = favours[[1]],
      vjust = 3.2,
      hjust = 0.5,
      size = base_size * 0.9 / ggplot2::.pt,
      colour = const$forest_text,
      fontface = "plain"
    ) +
    ggplot2::annotate(
      "text",
      x = x_right,
      y = -Inf,
      label = favours[[2]],
      vjust = 3.2,
      hjust = 0.5,
      size = base_size * 0.9 / ggplot2::.pt,
      colour = const$forest_text,
      fontface = "plain"
    )
}
