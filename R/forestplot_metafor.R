# ---- forestplot_metafor.R ----
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
    (interactmeta_is_bayesian_fit(x) ||
       !interactmeta_is_bayesian_fit(x) ||
       inherits(x$figure, "interactmeta_forestplot_bayesmeta") ||
       inherits(x$figure, "interactmeta_forestplot_metafor"))
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
  combined_lab <- forestplot_combined_block_label()
  combined_point_idx <- point_idx[forestplot_heading_is_combined(heading[point_idx])]
  block_order <- block_order[!forestplot_heading_is_combined(block_order)]

  # When group = "trial", study headings are trials but subgroup pooled /
  # prediction rows keep subgroup headings. Append those summary-only blocks
  # (with both pooled and prediction) after trial blocks and before Interaction.
  poly_headings <- unique(heading[is_poly])
  poly_headings <- poly_headings[nzchar(poly_headings) & !is.na(poly_headings)]
  extra_sg <- setdiff(
    poly_headings,
    c(block_order, "Interaction", "Combined", combined_lab)
  )
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
  sg_blocks <- block_order[
    !forestplot_heading_is_interaction(block_order) &
      !forestplot_heading_is_combined(block_order)
  ]
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
      mlab <- if (nzchar(sg)) {
        paste0(combined_lab, ": ", sg)
      } else {
        combined_lab
      }
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
    headers$label <- c(headers$label, combined_lab)
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
        p <- add_poly(j, combined_lab)
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

#' Axis title for metafor::forest (Favours are drawn separately)
#' @noRd
forestplot_metafor_xlab <- function(xlab, favours = NULL) {
  force(favours) # accepted for backward-compatible callers; ignored
  if (is.null(xlab) || !length(xlab) || all(!nzchar(as.character(xlab)))) {
    return("")
  }
  as.character(xlab)[[1L]]
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
