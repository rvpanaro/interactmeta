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
    panel_roles <- rep("modifiers", n)
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

#' Whether `x` is an interactmeta forestplot drawn on base graphics
#' @noRd
forestplot_is_metafor_engine <- function(x) {
  inherits(x, "interactmeta_forestplot")
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

#' Merge user metafor args into interactmeta defaults
#'
#' User `...` from [forestplot()] / `forest.interactmeta()` overrides presentation
