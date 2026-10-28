# ---- ifrac_design.R ----
#' Add information-fraction design columns
#'
#' Builds trial-level information fractions for non-reference subgroups. For
#' each trial, information is `1 / vi` on subgroup rows, and information
#' fractions are `IF_g = I_g / sum_h I_h`. For each non-reference subgroup `k`,
#' `ifrac_k = IF_k` (repeated on all rows in the trial).
#'
#' @param dat A data frame with at least `trial`, `subgroup`, `vi`, and `reference`.
#' @param trial_col Name of the trial column.
#' @param subgroup_col Name of the subgroup column.
#' @param vi_col Name of the variance column.
#' @param reference_col Name of the reference subgroup column (constant within trial).
#' @param include_contrast_rows If `FALSE` (default), rows whose `subgroup` label
#'   contains `contrast_sep` are excluded (e.g. `"g2-g1"`).
#' @param contrast_sep Separator used to identify derived contrast rows when
#'   `include_contrast_rows = FALSE`.
#' @param subgroup_levels Optional character vector of subgroup labels to keep.
#'   When provided, only these labels are used (overrides `contrast_sep` filtering).
#'
#' @return The input data (subgroup rows only, unless `include_contrast_rows = TRUE`)
#'   with columns `precision` and `ifrac_*` for each non-reference subgroup.
#'
#' @keywords internal
#' @export
add_ifrac_design <- function(dat,
                             trial_col = NULL,
                             subgroup_col = NULL,
                             vi_col = "vi",
                             reference_col = "reference",
                             include_contrast_rows = FALSE,
                             contrast_sep = "-",
                             subgroup_levels = NULL) {
  trial_col <- trial_col %||% interactmeta_resolve_study_col(dat)
  subgroup_col <- subgroup_col %||% interactmeta_resolve_level_col(dat)
  check_ifrac_design_cols(dat, trial_col, subgroup_col, vi_col, reference_col)

  out <- dat

  if (!include_contrast_rows) {
    if (!is.null(subgroup_levels)) {
      out <- out[out[[subgroup_col]] %in% subgroup_levels, , drop = FALSE]
    } else {
      keep <- is_ifrac_component_row(out[[subgroup_col]], contrast_sep = contrast_sep)
      out <- out[keep, , drop = FALSE]
    }
  }

  if (nrow(out) == 0L) {
    return(out)
  }

  vi <- out[[vi_col]]
  if (any(!is.finite(vi)) || any(vi <= 0)) {
    stop("`vi` must be positive and finite on all subgroup rows used for the design.", call. = FALSE)
  }

  ref_per_trial <- out |>
    dplyr::distinct(
      dplyr::across(dplyr::all_of(c(trial_col, reference_col)))
    )

  if (any(is.na(ref_per_trial[[reference_col]])) || any(!nzchar(ref_per_trial[[reference_col]]))) {
    stop("Each trial must have a non-missing `reference` value.", call. = FALSE)
  }

  n_trials <- dplyr::n_distinct(out[[trial_col]])
  if (nrow(ref_per_trial) != n_trials) {
    stop("Each trial must have exactly one reference subgroup label.", call. = FALSE)
  }

  out <- out |>
    dplyr::mutate(
      # Information: I_g = 1 / v_g. Trial IF of subgroup k is I_k / sum_h I_h.
      precision = 1 / .data[[vi_col]]
    )

  if_tbl <- out |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(trial_col, subgroup_col)))) |>
    dplyr::summarise(
      precision = sum(.data$precision),
      .groups = "drop"
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(trial_col))) |>
    dplyr::mutate(
      ifrac = .data$precision / sum(.data$precision)
    ) |>
    dplyr::ungroup()

  sum_if <- if_tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(trial_col))) |>
    dplyr::summarise(sum_ifrac = sum(.data$ifrac), .groups = "drop")

  if (any(abs(sum_if$sum_ifrac - 1) > 1e-8)) {
    stop("Within each trial, information fractions must sum to 1.", call. = FALSE)
  }

  non_ref <- if_tbl |>
    dplyr::left_join(
      ref_per_trial,
      by = trial_col
    ) |>
    dplyr::filter(.data[[subgroup_col]] != .data$reference) |>
    dplyr::select(
      dplyr::all_of(c(trial_col, subgroup_col)),
      "ifrac"
    )

  if (nrow(non_ref) > 0L) {
    ifrac_wide <- non_ref |>
      tidyr::pivot_wider(
        id_cols = dplyr::all_of(trial_col),
        names_from = dplyr::all_of(subgroup_col),
        values_from = "ifrac",
        names_prefix = "ifrac_"
      )

    out <- out |>
      dplyr::left_join(ifrac_wide, by = trial_col)

    non_ref_labels <- sort(unique(non_ref[[subgroup_col]]))
    for (lab in non_ref_labels) {
      icol <- paste0("ifrac_", lab)
      if (!icol %in% names(out)) {
        next
      }
      out[[icol]] <- as.numeric(out[[icol]])
    }
  }

  ifrac_k_cols <- sort(grep("^ifrac_", names(out), value = TRUE))
  ifrac_k_cols <- ifrac_k_cols[!grepl("^ifrac2_", ifrac_k_cols)]
  design_cols <- c("precision", ifrac_k_cols)
  design_cols <- design_cols[design_cols %in% names(out)]

  front <- setdiff(names(dat), c(design_cols, reference_col))
  tail <- if (reference_col %in% names(out)) reference_col else character(0)
  keep <- unique(c(front, design_cols, tail))
  keep <- keep[keep %in% names(out)]
  keep <- keep[!grepl("^ifrac2_", keep)]
  out <- out[, keep, drop = FALSE]
  coerce_ifrac_numeric_cols(out)
}

check_ifrac_design_cols <- function(dat, trial_col, subgroup_col, vi_col, reference_col) {
  need <- c(trial_col, subgroup_col, vi_col, reference_col)
  missing <- setdiff(need, names(dat))
  if (length(missing) > 0L) {
    stop(
      sprintf("Missing required column(s): %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Centered IF design columns used only for CAMS fitting
#'
#' `ifrac2_k = 1(subgroup == k) - ifrac_k`. Not returned by [ifrac()] or
#' [add_ifrac_design()].
#' @noRd
add_ifrac2_design <- function(dat, subgroup_col = NULL) {
  subgroup_col <- subgroup_col %||% interactmeta_resolve_level_col(dat)
  ifrac_cols <- grep("^ifrac_", names(dat), value = TRUE)
  ifrac_cols <- ifrac_cols[!grepl("^ifrac2_", ifrac_cols)]
  sg <- as.character(dat[[subgroup_col]])
  for (icol in ifrac_cols) {
    lab <- sub("^ifrac_", "", icol)
    i2col <- paste0("ifrac2_", lab)
    if (!i2col %in% names(dat)) {
      dat[[i2col]] <- as.numeric(sg == lab) - as.numeric(dat[[icol]])
    }
  }
  coerce_ifrac_numeric_cols(dat)
}

# ---- ifrac.R ----
#' Information fraction within each trial and pooled across trials
#'
#' **Design** (default): computes information fractions within each trial from
#' subgroup variances.
#'
#' Each component precision is defined as `1 / vi`.
#'
#' By default (`contrast = "all"`), all non-reference-vs-reference contrasts are
#' added per trial. Set `contrast = NULL` to return only subgroup rows with wide
#' `ifrac_<subgroup>` columns (one fraction per subgroup label).
#'
#' When `contrast` is provided, the information fraction for each comparison is
#' computed as the non-reference subgroup precision divided by the total subgroup
#' precision within the trial. By
#' convention, the reference subgroup is the label *after* `contrast_sep`, i.e.
#' in a label like `"g2-g1"`, `"g1"` is the reference subgroup.
#'
#' By default, components are the observed subgroup labels. When
#' `include_combined = TRUE`, one precision-weighted combined row per trial is
#' appended with `subgroup = "Combined"` (contribution-weighted mean $m_j$ on the
#' chosen scale). When `include_interaction = TRUE`, contrast rows are stacked
#' below subgroup rows with `subgroup` labels like `"g2-g1"` (reference after
#' `contrast_sep`), using `yi` and `vi` for the contrast estimate and variance.
#' Wide columns `ifrac_<non-reference>` are filled
#' on every row (no `NA`), plus `precision`, and a `reference` column (last
#' column). For each non-reference subgroup `k`, `ifrac_k` is the trial-level
#' information fraction. On combined
#' and contrast rows, `precision` is the sum of subgroup precisions in the trial
#' (or in that comparison for a single contrast).
#'
#' @param x Effect estimates. A numeric vector, or a data frame or
#'   [`escalc`](https://wviechtb.github.io/metafor/reference/escalc.html) object
#'   with columns `yi` and `vi`. If omitted, `yi` is used.
#' @param yi Either a numeric vector of subgroup effect estimates, or a data frame /
#'   `escalc` object containing columns `yi` and `vi`. When a data frame is
#'   supplied as `x` or `yi`, `vi` must not be passed separately; `trial` and
#'   `subgroup` name columns (defaults `"trial"` and `"subgroup"`).
#' @param vi Numeric vector. Variance of the subgroup estimate. Not used when the
#'   estimates are a data frame or `escalc` object.
#' @param trial Trial identifier. A vector (same length as the estimates), or when
#'   the estimates are a data frame the name of the trial column (default `"trial"`).
#' @param modifier Column name for the subgroup levels when the estimates are a
#'   data frame. A single string, supplied by name. Overrides `subgroup`. This is
#'   not the vector of subgroup labels; pass that vector in `subgroup`.
#' @param subgroup Subgroup identifier. A vector (same length as the estimates), or
#'   when the estimates are a data frame the name of the subgroup column (default
#'   `"subgroup"`).
#' @param big_var Positive number. Variance used to impute a missing
#'   `trial`-by-`subgroup` cell (`yi = 0`, `vi = big_var`). For an `escalc`
#'   object, a row that already has the placeholder `yi = 0` and `vi = 8` keeps
#'   that estimate, and `vi` is replaced by `big_var`. Defaults to `100^2`.
#' @param ... Further arguments passed to methods.
#' @param contrast Either `NULL` (subgroup rows only, wide `ifrac_<subgroup>`
#'   columns), a character vector of length 2 giving the two subgroup labels to
#'   contrast within each trial (computed as `label1 - label2`, where `label2` is
#'   the reference), or `"all"` (default) to add each non-reference subgroup vs
#'   the reference.
#' @param reference Optional. Used when `contrast = "all"`. Comparisons are
#'   generated only between each non-reference subgroup and the reference subgroup
#'   (e.g. `g2-g1`, `g3-g1`, ...). If `reference` is `NULL`, the reference defaults
#'   to the first level of a factor `subgroup`, otherwise the first label in
#'   sorted unique order. Can be either an integer index or a subgroup label
#'   string.
#'
#'   When supplied, extra rows per trial are added with `subgroup` labels
#'   `paste0(label1, contrast_sep, label2)`, and `yi`/`vi` computed as
#'   `yi1 - yi2` and `vi1 + vi2`. Missing `trial`-by-`subgroup` combinations use
#'   the imputed `yi=0` and `vi=big_var`.
#' @param contrast_sep String used between subgroup labels in the contrast row
#'   label. Defaults to `"-"`.
#' @param include_combined Logical. If `FALSE` (default), return subgroup rows only.
#'   If `TRUE`, append one precision-weighted combined row per trial with
#'   `subgroup = "Combined"`.
#' @param include_interaction Logical. If `FALSE` (default), return subgroup rows
#'   only while still using `contrast` to define the reference-centered
#'   `ifrac_*` design columns. If `TRUE` and `contrast` is not
#'   `NULL`, append interaction/contrast rows below subgroup rows.
#'
#' **Pooling** ( [`interactmeta`], or [ifrac()] data with `ifrac_*` columns):
#' pass `ifrac` to pool trial-level fractions across studies.
#' \deqn{\bar{\omega} = \sum_i \pi_i w_i / \sum_i w_i.}
#' With `ifrac = "interaction"` (default),
#' \eqn{w_i = 1/(\tau_\gamma^2 + \sigma_{Ai}^2 + \sigma_{Bi}^2)} using the
#' interaction-block SD \eqn{\tau_\gamma}.
#' With `ifrac = "overall"`, \eqn{w_i = 1/(\tau^2 + v_i)} using the
#' overall-block intercept SD \eqn{\tau} and
#' \eqn{v_i=\mathrm{Var}(m_i)=(1-\pi_i)^2\sigma_{Ai}^2+\pi_i^2\sigma_{Bi}^2}.
#' On Bayesian CAMS, \eqn{\tau} / \eqn{\tau_\gamma} are DIRECT posterior supports, so
#' the pooled IF is a distribution. On frequentist CAMS (`engine = "metafor"`),
#' pooling uses the fitted \eqn{\tau} or \eqn{\tau_\gamma}.
#' See argument details below and examples.
#'
#' @param ifrac For pooling only: `"interaction"` (default) or `"overall"`.
#' @param fit,tau,tau_gamma,ifrac_var,summary,prob Pooling options when `x` is a
#'   [`interactmeta`] or design data with `ifrac_*` columns (`fit` or `tau` /
#'   `tau_gamma` required for data-only pooling).
#'
#' @return For the **design** method, a tibble with one row per `trial`-by-`subgroup` combination (subgroup
#'   levels and contrast rows when requested), including `trial`, `subgroup`
#'   (a factor with the reference subgroup as the first level on subgroup rows),
#'   `yi`, `vi`, and information-fraction columns. Columns `yi`, `vi`, `precision`
#'   (when present), and all `ifrac_<...>` columns are numeric.
#'
#'   For **pooling**, a `interactmeta_overall_ifrac` object (draws with metadata), or a
#'   summary data frame when `summary = TRUE`.
#'
#' @examples
#' yi <- c(0.2, 0.1, 0.3, 0.15)
#' vi <- c(0.01, 0.04, 0.0225, 0.09)
#' trial <- c("A", "A", "B", "B")
#' subgroup <- c("g1", "g2", "g1", "g2")
#'
#' # Default: all non-reference vs reference contrasts (reference = first
#' # factor level, or first sorted label for character subgroups).
#' # Name `subgroup`: the fourth positional argument is `modifier`.
#' ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup)
#'
#' # Subgroup rows only, wide ifrac_<subgroup> for every label
#' ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, contrast = NULL)
#'
#' # Single contrast (g1 - g2, reference g2)
#' ifrac(
#'   yi = yi, vi = vi, trial = trial, subgroup = subgroup,
#'   contrast = c("g1", "g2")
#' )
#'
#' # Specify reference subgroup explicitly
#' ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, reference = "g1")
#'
#' # Append combined rows (one per trial)
#' ifrac(
#'   yi = yi, vi = vi, trial = trial, subgroup = subgroup,
#'   include_combined = TRUE
#' )
#'
#' # Append interaction rows at the bottom
#' ifrac(
#'   yi = yi, vi = vi, trial = trial, subgroup = subgroup,
#'   include_interaction = TRUE
#' )
#'
#' \dontrun{
#' library(interact)
#' e <- escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
#' d <- ifrac(e, reference = "Lower", include_interaction = TRUE)
#' fit <- cams(d, reference = "Lower")
#' forestplot(fit)
#' plot_contribution(fit)
#'
#' # Pool across trials (interaction-IF, the default)
#' ifrac(fit)
#' ifrac(fit, ifrac = "overall")
#' }
#'
#' @export
ifrac <- function(x, ...) {
  UseMethod("ifrac")
}

#' Call [ifrac()] when a local `ifrac` argument would mask the generic
#' @noRd
interactmeta_eval_ifrac <- function(...) {
  ifrac(...)
}

#' @rdname ifrac
#' @export
ifrac.default <- function(x,
                          vi = NULL,
                          trial = NULL,
                          modifier = NULL,
                          subgroup = NULL,
                          big_var = 100^2,
                          contrast = "all",
                          contrast_sep = "-",
                          reference = NULL,
                          include_combined = FALSE,
                          include_interaction = FALSE,
                          yi = NULL,
                          ...) {
  if (missing(x)) {
    if (!is.null(yi)) {
      x <- yi
    } else {
      stop("Supply `x` (effect estimates) or `yi`.", call. = FALSE)
    }
  }
  if (is.character(modifier) && length(modifier) == 1L && nzchar(modifier)) {
    subgroup <- modifier
  }
  ifrac_design(
    yi = x,
    vi = vi,
    trial = trial,
    subgroup = subgroup,
    big_var = big_var,
    contrast = contrast,
    contrast_sep = contrast_sep,
    reference = reference,
    include_combined = include_combined,
    include_interaction = include_interaction
  )
}

#' @rdname ifrac
#' @export
ifrac.data.frame <- function(x,
                             ifrac = NULL,
                             fit = NULL,
                             tau = NULL,
                             tau_gamma = NULL,
                             reference = NULL,
                             trial = "study",
                             modifier = NULL,
                             subgroup = "level",
                             ifrac_var = NULL,
                             contrast_sep = "-",
                             summary = FALSE,
                             prob = 0.95,
                             big_var = 100^2,
                             contrast = "all",
                             include_combined = FALSE,
                             include_interaction = FALSE,
                             yi = NULL,
                             ...) {
  if (is.character(modifier) && length(modifier) == 1L && nzchar(modifier)) {
    subgroup <- modifier
  }
  if (is.data.frame(x) || inherits(x, "escalc")) {
    subgroup <- interactmeta_resolve_modifier_column(
      modifier,
      subgroup,
      dat = as.data.frame(x)
    )
  }
  if (!is.null(yi)) {
    x <- yi
  }
  if (ifrac_requests_pool(
    x,
    ifrac = ifrac,
    fit = fit,
    tau = tau,
    tau_gamma = tau_gamma,
    summary = summary
  )) {
    return(overall_ifrac_from_data(
      x,
      fit = fit,
      tau = tau,
      tau_gamma = tau_gamma,
      ifrac = ifrac %||% "interaction",
      reference = reference,
      trial = trial,
      subgroup = subgroup,
      ifrac_var = ifrac_var,
      contrast_sep = contrast_sep,
      summary = summary,
      prob = prob
    ))
  }
  ifrac_design(
    yi = x,
    vi = NULL,
    trial = trial,
    subgroup = subgroup,
    big_var = big_var,
    contrast = contrast,
    contrast_sep = contrast_sep,
    reference = reference,
    include_combined = include_combined,
    include_interaction = include_interaction
  )
}

#' @noRd
ifrac_requests_pool <- function(x,
                                ifrac = NULL,
                                fit = NULL,
                                tau = NULL,
                                tau_gamma = NULL,
                                summary = FALSE) {
  if (!is.null(ifrac)) {
    return(TRUE)
  }
  if (!is.null(fit) || !is.null(tau) || !is.null(tau_gamma)) {
    return(TRUE)
  }
  isTRUE(summary) && length(detect_ifrac_vars(x)) > 0L
}

#' @noRd
ifrac_design <- function(yi, vi = NULL, trial = NULL, subgroup = NULL, big_var = 100^2, contrast = "all", contrast_sep = "-", reference = NULL, include_combined = FALSE, include_interaction = FALSE) {
  check_big_var(big_var)
  inputs <- coerce_ifrac_inputs(yi, vi, trial, subgroup, big_var = big_var)
  yi <- inputs$yi
  vi <- inputs$vi
  trial <- inputs$trial
  subgroup <- inputs$subgroup

  check_yi(yi)
  check_vi(vi)
  check_lengths(yi = yi, vi = vi, trial = trial, subgroup = subgroup)
  check_contrast(contrast)
  check_contrast_sep(contrast_sep)
  check_reference(reference)
  check_include_combined(include_combined)
  check_include_interaction(include_interaction)

  # Prefer factor level order (first level = default reference); else sort.
  subgroup_level_order <- inputs$subgroup_levels

  dat <- dplyr::tibble(
    trial = trial,
    subgroup = as.character(subgroup),
    yi = yi,
    vi = vi
  )
  if (!is.null(inputs$n)) {
    dat$n <- inputs$n
  }

  dat <- dat |>
    dplyr::mutate(.prec = 1 / .data$vi) |>
    dplyr::group_by(.data$trial, .data$subgroup)
  if ("n" %in% names(dat)) {
    dat <- dat |>
      dplyr::summarise(
        yi = sum(.data$yi * .data$.prec) / sum(.data$.prec),
        vi = 1 / sum(.data$.prec),
        n = sum(.data$n, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    dat <- dat |>
      dplyr::summarise(
        yi = sum(.data$yi * .data$.prec) / sum(.data$.prec),
        vi = 1 / sum(.data$.prec),
        .groups = "drop"
      )
  }

  observed <- unique(dat$subgroup)
  all_subgroups <- if (!is.null(subgroup_level_order) && length(subgroup_level_order)) {
    intersect(as.character(subgroup_level_order), observed)
  } else {
    sort(observed)
  }
  # Keep any observed labels missing from the preferred order (defensive).
  all_subgroups <- c(all_subgroups, setdiff(sort(observed), all_subgroups))

  complete_fill <- list(yi = 0, vi = big_var)
  if ("n" %in% names(dat)) {
    complete_fill$n <- 0
  }
  dat_complete <- dat |>
    tidyr::complete(
      trial,
      subgroup = all_subgroups,
      fill = complete_fill
    )

  contrast_rows <- NULL
  ref_label <- NULL
  if (!is.null(contrast)) {
    # Build contrast pairs
    if (identical(contrast, "all")) {
      if (length(all_subgroups) < 2L) {
        pairs <- tibble::tibble(sg1 = character(0), sg2 = character(0))
      } else {
        # Always use a single reference (sg2, after the separator). If not
        # provided, default to the first subgroup label (factor level order
        # when available, otherwise sorted).
        ref_label <- if (is.null(reference)) all_subgroups[[1]] else resolve_reference(reference, all_subgroups)
        others <- setdiff(all_subgroups, ref_label)
        pairs <- tibble::tibble(sg1 = others, sg2 = ref_label)
      }
    } else {
      ref_label <- contrast[[2]]
      pairs <- tibble::tibble(sg1 = contrast[[1]], sg2 = ref_label)
    }

    pairs <- pairs |>
      dplyr::mutate(contrast_label = paste0(.data$sg1, contrast_sep, .data$sg2))

    # Wide per trial for fast lookups
    wide_tv <- dat_complete |>
      dplyr::mutate(.prec = 1 / .data$vi) |>
      dplyr::select(trial, subgroup, yi, vi, .prec) |>
      tidyr::pivot_wider(names_from = "subgroup", values_from = c("yi", "vi", ".prec"))

    # Cross join (trial x pair) then compute contrasts
    contrast_rows <- tidyr::crossing(
      wide_tv |> dplyr::select(trial),
      pairs
    ) |>
      dplyr::left_join(wide_tv, by = "trial")

    if (nrow(contrast_rows) > 0) {
      yi1_cols <- paste0("yi_", contrast_rows$sg1)
      yi2_cols <- paste0("yi_", contrast_rows$sg2)
      vi1_cols <- paste0("vi_", contrast_rows$sg1)
      vi2_cols <- paste0("vi_", contrast_rows$sg2)
      prec1_cols <- paste0(".prec_", contrast_rows$sg1)
      prec2_cols <- paste0(".prec_", contrast_rows$sg2)

      idx <- seq_len(nrow(contrast_rows))
      contrast_rows$yi <- vapply(idx, function(i) contrast_rows[[yi1_cols[[i]]]][[i]] - contrast_rows[[yi2_cols[[i]]]][[i]], numeric(1L))
      contrast_rows$vi <- vapply(idx, function(i) contrast_rows[[vi1_cols[[i]]]][[i]] + contrast_rows[[vi2_cols[[i]]]][[i]], numeric(1L))
      contrast_rows$precision <- vapply(
        idx,
        function(i) contrast_rows[[prec1_cols[[i]]]][[i]] + contrast_rows[[prec2_cols[[i]]]][[i]],
        numeric(1L)
      )
    } else {
      contrast_rows$sg1 <- character(0)
      contrast_rows$yi <- numeric(0)
      contrast_rows$vi <- numeric(0)
      contrast_rows$precision <- numeric(0)
    }

    contrast_rows <- contrast_rows |>
      dplyr::transmute(
        trial = .data$trial,
        subgroup = .data$contrast_label,
        yi = .data$yi,
        vi = .data$vi,
        precision = .data$precision,
        sg1 = .data$sg1
      )
  }

  ifrac_wide <- NULL

  if (is.null(contrast)) {
    ifrac_wide <- dat_complete |>
      dplyr::mutate(.prec = 1 / .data$vi) |>
      dplyr::group_by(.data$trial, .data$subgroup) |>
      dplyr::summarise(.prec = sum(.data$.prec), .groups = "drop") |>
      dplyr::group_by(.data$trial) |>
      dplyr::mutate(ifrac = .data$.prec / sum(.data$.prec)) |>
      dplyr::select(-.prec) |>
      tidyr::pivot_wider(
        names_from = "subgroup",
        values_from = "ifrac",
        names_glue = "ifrac_{subgroup}"
      )
  }

  subgroup_cols <- c("trial", "subgroup", "yi", "vi")
  if ("n" %in% names(dat_complete)) {
    subgroup_cols <- c(subgroup_cols, "n")
  }
  subgroup_out <- dat_complete |>
    dplyr::select(dplyr::all_of(subgroup_cols))

  if (!is.null(contrast)) {
    subgroup_out$reference <- ref_label
    subgroup_out <- add_ifrac_design(
      subgroup_out,
      contrast_sep = contrast_sep,
      subgroup_levels = all_subgroups
    )
  }

  out <- subgroup_out
  derived_out <- NULL

  if (isTRUE(include_combined)) {
    derived_out <- build_ifrac_combined_rows(dat_complete)
    if (!is.null(contrast) && nrow(derived_out) > 0L) {
      ifrac_k_cols <- grep("^ifrac_", names(subgroup_out), value = TRUE)
      ifrac_k_cols <- ifrac_k_cols[!grepl("^ifrac2_", ifrac_k_cols)]
      design_by_trial <- subgroup_out |>
        dplyr::distinct(.data$trial, dplyr::across(dplyr::all_of(ifrac_k_cols)))
      derived_out <- derived_out |>
        dplyr::left_join(design_by_trial, by = "trial") |>
        dplyr::mutate(reference = ref_label)
    } else if (is.null(contrast) && !is.null(ifrac_wide) && nrow(derived_out) > 0L) {
      derived_out <- derived_out |>
        dplyr::left_join(ifrac_wide, by = "trial")
    }
  }

  if (isTRUE(include_interaction) && !is.null(contrast) && nrow(contrast_rows) > 0L) {
    ifrac_k_cols <- grep("^ifrac_", names(subgroup_out), value = TRUE)
    ifrac_k_cols <- ifrac_k_cols[!grepl("^ifrac2_", ifrac_k_cols)]
    design_by_trial <- subgroup_out |>
      dplyr::distinct(.data$trial, dplyr::across(dplyr::all_of(ifrac_k_cols)))

    contrast_out <- contrast_rows |>
      dplyr::select("trial", "subgroup", "yi", "vi", "precision") |>
      dplyr::left_join(design_by_trial, by = "trial") |>
      dplyr::mutate(reference = ref_label)
    if ("n" %in% names(subgroup_out)) {
      contrast_out$n <- NA_real_
    }

    derived_out <- if (is.null(derived_out)) {
      contrast_out
    } else {
      dplyr::bind_rows(derived_out, contrast_out)
    }
  }

  if (!is.null(derived_out) && nrow(derived_out) > 0L) {
    out <- dplyr::bind_rows(subgroup_out, derived_out)

    base_cols <- ifrac_output_base_cols(out)
    other_cols <- setdiff(names(out), base_cols)
    ref_pos <- match("reference", other_cols)
    if (!is.na(ref_pos)) {
      other_cols <- c(other_cols[-ref_pos], "reference")
    }
    out <- out[, c(base_cols, other_cols)]
  } else if (is.null(contrast)) {
    out <- out |>
      dplyr::left_join(ifrac_wide, by = "trial")
  } else if (!is.null(contrast)) {
    ifrac_k_cols <- grep("^ifrac_", names(out), value = TRUE)
    ifrac_k_cols <- ifrac_k_cols[!grepl("^ifrac2_", ifrac_k_cols)]
    base_cols <- ifrac_output_base_cols(out)
    other_cols <- c("precision", ifrac_k_cols, "reference")
    other_cols <- other_cols[other_cols %in% names(out)]
    ref_pos <- match("reference", other_cols)
    if (!is.na(ref_pos)) {
      other_cols <- c(other_cols[-ref_pos], "reference")
    }
    out <- out[, c(base_cols, other_cols)]
  }

  out <- set_subgroup_factor(out, ref_label, all_subgroups, contrast_sep)
  out <- out[, !grepl("^ifrac2_", names(out)), drop = FALSE]
  out <- interactmeta_rename_modifier_col_to_level(out)
  out <- interactmeta_rename_trial_col_to_study(out)
  coerce_ifrac_numeric_cols(out)
}

ifrac_combined_label <- function() {
  "Combined"
}

#' @noRd
is_ifrac_contrast_row <- function(subgroup, contrast_sep = "-") {
  grepl(contrast_sep, as.character(subgroup), fixed = TRUE)
}

#' @noRd
is_ifrac_combined_row <- function(subgroup, combined_label = ifrac_combined_label()) {
  as.character(subgroup) == combined_label
}

#' @noRd
is_ifrac_component_row <- function(subgroup, contrast_sep = "-", combined_label = ifrac_combined_label()) {
  !is_ifrac_contrast_row(subgroup, contrast_sep) &
    !is_ifrac_combined_row(subgroup, combined_label)
}

#' @noRd
build_ifrac_combined_rows <- function(dat_complete, combined_label = ifrac_combined_label()) {
  has_n <- "n" %in% names(dat_complete)
  ref <- if ("reference" %in% names(dat_complete)) {
    unique(as.character(dat_complete$reference))[[1L]]
  } else {
    as.character(dat_complete$subgroup)[[1L]]
  }
  level_col <- interactmeta_resolve_level_col(dat_complete)
  blocks <- comparison_trial_effect_blocks(
    data = dat_complete,
    trial_var = interactmeta_resolve_study_col(dat_complete),
    subgroup_var = level_col,
    reference = ref
  )
  out <- blocks$overall
  out$subgroup <- combined_label
  out$precision <- 1 / pmax(out$vi, 1e-12)
  if (has_n) {
    n_by_trial <- dat_complete |>
      dplyr::group_by(.data$trial) |>
      dplyr::summarise(n = sum(.data$n, na.rm = TRUE), .groups = "drop")
    out <- dplyr::left_join(out, n_by_trial, by = "trial")
  }
  out
}

ifrac_output_base_cols <- function(dat) {
  cols <- c(interactmeta_resolve_study_col(dat), interactmeta_resolve_level_col(dat), "yi", "vi")
  if ("n" %in% names(dat)) {
    cols <- c(cols, "n")
  }
  cols
}

coerce_ifrac_numeric_cols <- function(dat) {
  num_cols <- names(dat)[
    names(dat) %in% c("yi", "vi", "precision", "n") |
      grepl("^ifrac2?_", names(dat))
  ]
  if (!length(num_cols)) {
    return(dat)
  }
  for (col in num_cols) {
    dat[[col]] <- as.numeric(dat[[col]])
  }
  dat
}

set_subgroup_factor <- function(out, ref_label, all_subgroups, contrast_sep) {
  col <- interactmeta_resolve_level_col(out)
  if (!(col %in% names(out)) && "subgroup" %in% names(out)) {
    col <- "subgroup"
  }
  if (!(col %in% names(out))) {
    return(out)
  }

  ref <- ref_label
  if (is.null(ref) && "reference" %in% names(out)) {
    ref <- unique(out$reference)
    if (length(ref) != 1L || is.na(ref) || !nzchar(ref)) {
      stop("`reference` must be a single non-missing value when factoring `subgroup`.", call. = FALSE)
    }
    ref <- ref[[1]]
  }

  if (is.null(ref)) {
    return(out)
  }

  if (!(ref %in% all_subgroups)) {
    stop(
      sprintf("`reference` (%s) must match a subgroup label in the data.", ref),
      call. = FALSE
    )
  }

  sg_levels <- c(ref, all_subgroups[all_subgroups != ref])
  is_subgroup_row <- rep(TRUE, nrow(out))

  if (!is.null(contrast_sep) && nzchar(contrast_sep)) {
    is_subgroup_row <- is_ifrac_component_row(out[[col]], contrast_sep = contrast_sep)
  }

  if (all(is_subgroup_row)) {
    out[[col]] <- factor(out[[col]], levels = sg_levels)
    stopifnot(identical(levels(out[[col]])[[1]], ref))
  } else if (any(is_subgroup_row)) {
    out[[col]] <- as.character(out[[col]])
    out[[col]][is_subgroup_row] <- as.character(
      factor(out[[col]][is_subgroup_row], levels = sg_levels)
    )
  }

  out
}

check_include_combined <- function(include_combined) {
  if (!is.logical(include_combined) || length(include_combined) != 1L || is.na(include_combined)) {
    stop("`include_combined` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  invisible(TRUE)
}

check_include_interaction <- function(include_interaction) {
  if (!is.logical(include_interaction) || length(include_interaction) != 1L || is.na(include_interaction)) {
    stop("`include_interaction` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  invisible(TRUE)
}

check_contrast <- function(contrast) {
  if (is.null(contrast)) {
    return(invisible(TRUE))
  }

  if (identical(contrast, "all")) {
    return(invisible(TRUE))
  }

  if (!is.character(contrast) || length(contrast) != 2L) {
    stop("`contrast` must be NULL, \"all\", or a character vector of length 2.", call. = FALSE)
  }

  if (any(!nzchar(contrast))) {
    stop("`contrast` must not contain empty strings.", call. = FALSE)
  }

  invisible(TRUE)
}

check_contrast_sep <- function(contrast_sep) {
  if (!is.character(contrast_sep) || length(contrast_sep) != 1L) {
    stop("`contrast_sep` must be a single character string.", call. = FALSE)
  }

  invisible(TRUE)
}

check_reference <- function(reference) {
  if (is.null(reference)) {
    return(invisible(TRUE))
  }

  if (is.character(reference) && length(reference) == 1L && nzchar(reference)) {
    return(invisible(TRUE))
  }

  if (is.numeric(reference) && length(reference) == 1L && is.finite(reference)) {
    return(invisible(TRUE))
  }

  stop("`reference` must be NULL, a single subgroup label string, or a single numeric index.", call. = FALSE)
}

resolve_reference <- function(reference, subgroup_labels) {
  if (is.numeric(reference)) {
    idx <- as.integer(reference)
    if (idx < 1L || idx > length(subgroup_labels)) {
      stop("`reference` index is out of bounds for the available subgroup labels.", call. = FALSE)
    }
    return(subgroup_labels[[idx]])
  }

  ref <- as.character(reference)[[1]]
  if (!(ref %in% subgroup_labels)) {
    stop("`reference` must match an existing subgroup label.", call. = FALSE)
  }

  ref
}

is_escalc_placeholder_vi <- function(yi, vi) {
  is.finite(yi) & is.finite(vi) & yi == 0 & vi == 8
}

coerce_ifrac_inputs <- function(yi, vi, trial, subgroup, big_var = 100^2) {
  if (is.data.frame(yi) || inherits(yi, "escalc")) {
    dat <- as.data.frame(yi)

    if (!is.null(vi)) {
      stop(
        "Do not pass `vi` separately when the first argument is a data frame or `escalc` object.",
        call. = FALSE
      )
    }

    if (is.null(trial)) {
      trial <- interactmeta_resolve_study_col(dat)
    }
    if (is.null(subgroup)) {
      subgroup <- interactmeta_resolve_level_col(dat)
    }

    if (!is.character(trial) || length(trial) != 1L || !(trial %in% names(dat))) {
      trial <- interactmeta_resolve_study_col(dat, fallback = trial)
    }
    if (!(trial %in% names(dat))) {
      stop(
        sprintf("`trial` must be a single column name in the data. Got: %s.", trial),
        call. = FALSE
      )
    }

    if (!is.character(subgroup) || length(subgroup) != 1L || !(subgroup %in% names(dat))) {
      stop(
        sprintf("`subgroup` must be a single column name in the data. Got: %s.", subgroup),
        call. = FALSE
      )
    }

    if (!all(c("yi", "vi") %in% names(dat))) {
      stop("Data must contain columns `yi` and `vi`.", call. = FALSE)
    }

    ok <- is.finite(dat$yi) & is.finite(dat$vi) & dat$vi > 0
    if (!all(ok)) {
      n_drop <- sum(!ok)
      warning(
        sprintf(
          "Dropping %d row(s) with non-finite `yi`, non-finite `vi`, or non-positive `vi`.",
          n_drop
        ),
        call. = FALSE
      )
      dat <- dat[ok, , drop = FALSE]
    }

    if (nrow(dat) == 0L) {
      stop("No rows remain after removing invalid `yi` / `vi` values.", call. = FALSE)
    }

    out <- list(
      yi = as.numeric(dat$yi),
      vi = as.numeric(dat$vi),
      trial = dat[[trial]],
      subgroup = as.character(dat[[subgroup]]),
      subgroup_levels = if (is.factor(dat[[subgroup]])) {
        levels(dat[[subgroup]])
      } else {
        NULL
      }
    )
    if (inherits(yi, "escalc")) {
      placeholder <- is_escalc_placeholder_vi(out$yi, out$vi)
      if (any(placeholder)) {
        out$vi[placeholder] <- big_var
      }
    }
    if ("n" %in% names(dat)) {
      out$n <- as.numeric(dat$n)
    }
    return(out)
  }

  if (is.null(vi) || is.null(trial) || is.null(subgroup)) {
    stop(
      "`yi`, `vi`, `trial`, and `subgroup` are required when the first argument is not a data frame.",
      call. = FALSE
    )
  }

  list(
    yi = as.numeric(yi),
    vi = as.numeric(vi),
    trial = trial,
    subgroup = as.character(subgroup),
    subgroup_levels = if (is.factor(subgroup)) levels(subgroup) else NULL
  )
}

check_yi <- function(yi) {
  if (!is.numeric(yi)) {
    stop("`yi` must be numeric.", call. = FALSE)
  }

  if (any(!is.finite(yi))) {
    stop("`yi` must contain only finite values.", call. = FALSE)
  }

  invisible(TRUE)
}

check_lengths <- function(...) {
  args <- list(...)
  n <- length(args[[1]])

  bad <- vapply(args, length, integer(1L)) != n
  if (any(bad)) {
    stop(
      sprintf(
        "All inputs must have the same length. Got: %s.",
        paste(sprintf("%s=%d", names(args), vapply(args, length, integer(1L))), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

check_big_var <- function(big_var) {
  if (!is.numeric(big_var) || length(big_var) != 1L) {
    stop("`big_var` must be a single numeric value.", call. = FALSE)
  }

  if (!is.finite(big_var)) {
    stop("`big_var` must be finite.", call. = FALSE)
  }

  if (big_var <= 0) {
    stop("`big_var` must be positive.", call. = FALSE)
  }

  invisible(TRUE)
}

check_vi <- function(vi) {
  if (!is.numeric(vi)) {
    stop("`vi` must be numeric.", call. = FALSE)
  }

  if (any(!is.finite(vi))) {
    stop("`vi` must contain only finite values.", call. = FALSE)
  }

  if (any(vi <= 0)) {
    stop("`vi` must contain only positive values.", call. = FALSE)
  }

  invisible(TRUE)
}

# ---- overall_ifrac.R ----
#' @rdname ifrac
#' @export
ifrac.interactmeta <- function(x,
                               ifrac = c("interaction", "overall"),
                               ifrac_var = NULL,
                               summary = FALSE,
                               prob = 0.95,
                               ...) {
  overall_ifrac_from_fit(
    x,
    ifrac_var = ifrac_var,
    ifrac = match.arg(ifrac),
    summary = summary,
    prob = prob
  )
}


#' @rdname ifrac
#' @export
ifrac.interactmeta_overall_ifrac <- function(x, ...) {
  x
}

#' @noRd
overall_ifrac_from_fit <- function(fit,
                                   ifrac_var = NULL,
                                   ifrac = c("interaction", "overall"),
                                   summary = FALSE,
                                   prob = 0.95) {
  ifrac <- match.arg(ifrac)
  check_comparison_fit(fit)
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st)) {
    stop(
      "Fit has no `\"interactmeta\"` attribute; refit with `cams()`.",
      call. = FALSE
    )
  }

  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
  ref <- st$reference
  if (is.null(ref) || !nzchar(ref)) {
    stop("Fit has no `reference` in attribute `\"interactmeta\"`.", call. = FALSE)
  }

  d <- st$data_fit
  if (is.null(d) && !is.null(st$ifrac_data)) {
    d <- st$ifrac_data
    if (is.null(ifrac_var)) {
      ifrac_var <- st$ifrac_var
      if (is.null(ifrac_var) || !ifrac_var %in% names(d)) {
        vars <- detect_ifrac_vars(d)
        ifrac_var <- vars[[1]]
      }
    }
  } else if (is.null(ifrac_var)) {
    ifrac_var <- st$ifrac_vars_std[[1]]
    if (is.null(ifrac_var) && !is.null(st$ifrac_var)) {
      ifrac_var <- st$ifrac_var
    }
    if (is.null(ifrac_var) || !length(ifrac_var) || !ifrac_var %in% names(d)) {
      if (!is.null(d)) {
        vars <- detect_ifrac_vars(d)
        ifrac_var <- vars[[1]]
      }
    }
  }

  if (is.null(d)) {
    stop("Fit has no `data_fit` or `ifrac_data`; refit with `cams()`.", call. = FALSE)
  }

  draws <- if (identical(ifrac, "interaction")) {
    tau_gamma <- comparison_tau_gamma_samples(fit, ifrac_var = ifrac_var, trial_var = trial_var)
    interaction_ifrac_compute(
      data = d,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = ref,
      ifrac_var = ifrac_var,
      tau_gamma = tau_gamma,
      contrast_sep = "-"
    )
  } else {
    tau_hyp <- sprintf("%s__Intercept = 0", interactmeta_sd_re_group())
    tau <- comparison_sd_samples(fit, tau_hyp)
    overall_ifrac_compute(
      data = d,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = ref,
      ifrac_var = ifrac_var,
      tau = tau,
      contrast_sep = "-"
    )
  }

  trial_pi <- overall_ifrac_trial_pi(
    d,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = ref,
    ifrac_var = ifrac_var
  )
  bundle <- list(
    support = draws,
    trial_pi = trial_pi,
    reference = ref,
    ifrac_var = ifrac_var,
    subgroup_label = overall_ifrac_subgroup_label(ifrac_var, interactmeta = st),
    n_trials = nrow(trial_pi),
    summary = overall_ifrac_summary_df(draws, prob = prob, weight = interactmeta_ifrac_draw_weights(fit, draws)),
    prob = prob,
    ifrac = ifrac,
    fit = fit
  )
  obj <- overall_ifrac_result(draws, summary = summary, prob = prob, bundle = bundle)
  if (inherits(obj, "interactmeta_overall_ifrac")) {
    w <- interactmeta_ifrac_draw_weights(fit, draws)
    if (!is.null(w)) {
      attr(obj, ".weight") <- w
    }
    obj <- interactmeta_attach_origin(obj, fit, quantity = "overall_ifrac")
  }
  obj
}

#' @noRd
overall_ifrac_from_data <- function(data,
                                    fit = NULL,
                                    tau = NULL,
                                    tau_gamma = NULL,
                                    ifrac = c("interaction", "overall"),
                                    reference = NULL,
                                    trial = "study",
                                    subgroup = "level",
                                    ifrac_var = NULL,
                                    contrast_sep = "-",
                                    summary = FALSE,
                                    prob = 0.95) {
  ifrac <- match.arg(ifrac)
  if (!is.null(fit)) {
    check_comparison_fit(fit)
    if (!is_comparison_fit(fit)) {
      stop("`fit` must be an `interactmeta` from `cams()`.", call. = FALSE)
    }
    st <- attr(fit, "interactmeta", exact = TRUE)
    trial_var <- if (!is.null(st) && !is.null(st$trial_var)) st$trial_var else trial
    if (identical(ifrac, "interaction")) {
      if (is.null(ifrac_var)) {
        ifrac_var <- detect_ifrac_vars(data)[[1]]
      }
      tau_gamma <- comparison_tau_gamma_samples(fit, ifrac_var = ifrac_var, trial_var = trial_var)
    } else {
      tau_hyp <- sprintf("%s__Intercept = 0", interactmeta_sd_re_group())
      tau <- comparison_sd_samples(fit, tau_hyp)
    }
  } else if (identical(ifrac, "interaction")) {
    if (is.null(tau_gamma)) {
      stop(
        "Supply `fit` or `tau_gamma` when `object` is a data frame and `ifrac = \"interaction\"`.",
        call. = FALSE
      )
    }
  } else if (is.null(tau)) {
    stop(
      "Supply `fit` or `tau` when `object` is a data frame and `ifrac = \"overall\"`.",
      call. = FALSE
    )
  }

  if (identical(ifrac, "interaction")) {
    tau_gamma <- as.numeric(tau_gamma)
    if (!length(tau_gamma) || !all(is.finite(tau_gamma))) {
      stop("`tau_gamma` must be one or more finite numeric values.", call. = FALSE)
    }
  } else {
    tau <- as.numeric(tau)
    if (!length(tau) || !all(is.finite(tau))) {
      stop("`tau` must be one or more finite numeric values.", call. = FALSE)
    }
  }

  ref <- reference
  if (is.null(ref)) {
    if ("reference" %in% names(data)) {
      ref <- unique(data$reference)
      ref <- ref[!is.na(ref) & nzchar(ref)]
    }
    if (length(ref) != 1L) {
      stop(
        "Supply `reference` or include a single `reference` column in `object`.",
        call. = FALSE
      )
    }
    ref <- ref[[1]]
  }

  if (is.null(ifrac_var)) {
    ifrac_var <- detect_ifrac_vars(data)[[1]]
  }

  draws <- if (identical(ifrac, "interaction")) {
    interaction_ifrac_compute(
      data = data,
      trial_var = trial,
      subgroup_var = subgroup,
      reference = ref,
      ifrac_var = ifrac_var,
      tau_gamma = tau_gamma,
      contrast_sep = contrast_sep
    )
  } else {
    overall_ifrac_compute(
      data = data,
      trial_var = trial,
      subgroup_var = subgroup,
      reference = ref,
      ifrac_var = ifrac_var,
      tau = tau,
      contrast_sep = contrast_sep
    )
  }

  trial_pi <- overall_ifrac_trial_pi(
    data,
    trial_var = trial,
    subgroup_var = subgroup,
    reference = ref,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  bundle <- list(
    support = draws,
    trial_pi = trial_pi,
    reference = ref,
    ifrac_var = ifrac_var,
    subgroup_label = overall_ifrac_subgroup_label(ifrac_var),
    n_trials = nrow(trial_pi),
    summary = overall_ifrac_summary_df(draws, prob = prob),
    prob = prob,
    ifrac = ifrac
  )
  overall_ifrac_result(draws, summary = summary, prob = prob, bundle = bundle)
}

#' @noRd
overall_ifrac_prepare_trials <- function(data,
                                         trial_var,
                                         subgroup_var,
                                         reference,
                                         ifrac_var,
                                         contrast_sep = "-") {
  if (!ifrac_var %in% names(data)) {
    stop(sprintf("Column `%s` not found in data.", ifrac_var), call. = FALSE)
  }

  d <- data[
    is_ifrac_component_row(data[[subgroup_var]], contrast_sep = contrast_sep),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) {
    stop("No subgroup rows found for overall IF computation.", call. = FALSE)
  }

  sg <- as.character(d[[subgroup_var]])
  reference <- overall_ifrac_resolve_reference(sg, reference)
  d <- d[order(d[[trial_var]], sg != reference, sg), , drop = FALSE]
  d$se <- sqrt(d$vi)

  trials <- unique(as.character(d[[trial_var]]))
  n_trials <- length(trials)
  n_per <- nrow(d) / n_trials
  if (abs(n_per - round(n_per)) > 1e-8 || round(n_per) != n_per) {
    stop("Each trial must have the same number of subgroup rows.", call. = FALSE)
  }

  pi <- vapply(
    trials,
    function(tr) {
      d[[ifrac_var]][d[[trial_var]] == tr][[1]]
    },
    numeric(1)
  )

  list(data = d, trials = trials, pi = pi, reference = reference)
}

#' Match `reference` to subgroup labels in the IF table
#'
#' Fitted `data_fit` recodes subgroups to `interactmeta_sg_1` (reference)
#' / `interactmeta_sg_2` while `st$reference` stays the original label.
#' @noRd
overall_ifrac_resolve_reference <- function(sg, reference) {
  sg <- unique(as.character(sg))
  ref <- as.character(reference)[[1L]]
  if (length(ref) && ref %in% sg) {
    return(ref)
  }
  canon <- "interactmeta_sg_1"
  if (canon %in% sg) {
    return(canon)
  }
  ref
}

#' @noRd
overall_ifrac_weighted_pool <- function(pi, den, hetero_sd) {
  hetero_sd <- as.numeric(hetero_sd)
  wi <- vapply(hetero_sd, function(h) 1 / (h^2 + den), numeric(length(pi)))
  if (!is.matrix(wi)) {
    wi <- matrix(wi, nrow = length(pi), ncol = length(hetero_sd))
  }
  vapply(
    seq_len(ncol(wi)),
    function(j) {
      w <- wi[, j]
      sum(pi * w) / sum(w)
    },
    numeric(1)
  )
}

#' Trial-level IF coordinates and the two pooling denominators
#'
#' `pi` is the trial information fraction (same under both weighting schemes).
#' Combined denominator is \eqn{\mathrm{Var}(m_j)=1/\sum_g 1/\sigma_{gj}^2},
#' the overall-block sampling variance
#' \eqn{(1-\pi_j)^2\sigma_{Aj}^2+\pi_j^2\sigma_{Bj}^2}.
#' Interaction denominator is \eqn{\sigma_{Aj}^2+\sigma_{Bj}^2=\mathrm{Var}(g_j)}.
#' Weights are \eqn{w_j=1/(\tau^2+v_j)} at a supplied heterogeneity SD.
#' @noRd
overall_ifrac_trial_components <- function(data,
                                           trial_var,
                                           subgroup_var,
                                           reference,
                                           ifrac_var,
                                           contrast_sep = "-") {
  prep <- overall_ifrac_prepare_trials(
    data = data,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = reference,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  d <- prep$data
  trials <- prep$trials
  pi <- unname(prep$pi)
  reference <- prep$reference

  den_overall <- vapply(
    trials,
    function(tr) {
      ix <- d[[trial_var]] == tr
      vi <- pmax(as.numeric(d$vi[ix]), 1e-12)
      1 / sum(1 / vi)
    },
    numeric(1)
  )
  den_interaction <- vapply(
    trials,
    function(tr) {
      ix <- d[[trial_var]] == tr
      sum(d$vi[ix])
    },
    numeric(1)
  )
  max_vi <- vapply(
    trials,
    function(tr) {
      max(d$vi[d[[trial_var]] == tr], na.rm = TRUE)
    },
    numeric(1)
  )

  data.frame(
    trial = as.character(trials),
    pi = as.numeric(pi),
    den_overall = unname(as.numeric(den_overall)),
    den_interaction = unname(as.numeric(den_interaction)),
    max_vi = unname(as.numeric(max_vi)),
    stringsAsFactors = FALSE
  )
}

#' @noRd
overall_ifrac_compute <- function(data,
                                  trial_var,
                                  subgroup_var,
                                  reference,
                                  ifrac_var,
                                  tau,
                                  contrast_sep = "-") {
  comp <- overall_ifrac_trial_components(
    data = data,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = reference,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  overall_ifrac_weighted_pool(
    pi = comp$pi,
    den = comp$den_overall,
    hetero_sd = tau
  )
}

#' Interaction-IF weights: inverse of tau_gamma^2 + subgroup sampling variances
#' @noRd
interaction_ifrac_compute <- function(data,
                                      trial_var,
                                      subgroup_var,
                                      reference,
                                      ifrac_var,
                                      tau_gamma,
                                      contrast_sep = "-") {
  comp <- overall_ifrac_trial_components(
    data = data,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = reference,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  overall_ifrac_weighted_pool(
    pi = comp$pi,
    den = comp$den_interaction,
    hetero_sd = tau_gamma
  )
}

#' Posterior draws of interaction heterogeneity (trial ifrac2 SD)
#' @noRd
comparison_tau_gamma_samples <- function(fit, ifrac_var, trial_var = NULL) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(trial_var) || !nzchar(trial_var)) {
    trial_var <- if (!is.null(st) && !is.null(st$trial_var)) st$trial_var else "trial"
  }
  subgroup <- sub("^ifrac_", "", ifrac_var)
  interaction_param <- paste0("Interaction_", subgroup)
  legacy_ifrac2_param <- paste0("ifrac2_", subgroup)
  specs <- comparison_heterogeneity_specs(fit, trial_var = trial_var, component = "slope")
  if (!nrow(specs)) {
    stop(
      "No interaction heterogeneity (`Interaction_*` / `ifrac2_*` slope SD) found in the fit.",
      call. = FALSE
    )
  }
  hit <- specs$parameter %in% c(interaction_param, legacy_ifrac2_param)
  if (!any(hit)) {
    hit <- specs$type == "slope"
    if (!any(hit)) {
      stop(
        sprintf(
          "No slope heterogeneity parameter for `%s`; available: %s.",
          interaction_param,
          paste(specs$parameter, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    ix <- which(hit)[[1L]]
  } else {
    ix <- which(hit)[[1L]]
  }
  comparison_sd_samples(fit, specs$hypothesis[[ix]])
}

#' @noRd
overall_ifrac_weight_label <- function(ifrac = c("interaction", "overall")) {
  ifrac <- match.arg(ifrac)
  if (identical(ifrac, "interaction")) {
    "Interaction-weighted information fraction"
  } else {
    "Combined-weighted information fraction"
  }
}

#' Short row label for [print.interactmeta_overall_ifrac()] tables
#' @noRd
overall_ifrac_print_row_label <- function(ifrac = c("interaction", "overall")) {
  ifrac <- match.arg(ifrac)
  if (identical(ifrac, "interaction")) {
    "interaction-weighted"
  } else {
    "combined-weighted"
  }
}

#' Trial-level information fractions used in the overall IF weights
#' @noRd
overall_ifrac_trial_pi <- function(data,
                                   trial_var,
                                   subgroup_var,
                                   reference,
                                   ifrac_var,
                                   contrast_sep = "-") {
  d <- data[
    is_ifrac_component_row(data[[subgroup_var]], contrast_sep = contrast_sep),
    ,
    drop = FALSE
  ]
  trials <- unique(as.character(d[[trial_var]]))
  data.frame(
    trial = trials,
    pi = vapply(
      trials,
      function(tr) {
        d[[ifrac_var]][d[[trial_var]] == tr][[1]]
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

#' Draws plus metadata for plotting
#' @noRd
overall_ifrac_bundle <- function(object,
                                 fit = NULL,
                                 tau = NULL,
                                 tau_gamma = NULL,
                                 ifrac = c("interaction", "overall"),
                                 reference = NULL,
                                 trial = "study",
                                 subgroup = "level",
                                 ifrac_var = NULL,
                                 contrast_sep = "-",
                                 prob = 0.95,
                                 ...) {
  ifrac <- match.arg(ifrac)
  if (inherits(object, "interactmeta_overall_ifrac")) {
    return(list(
      draws = as.numeric(object),
      trial_pi = attr(object, "trial_pi", exact = TRUE),
      reference = attr(object, "reference", exact = TRUE),
      ifrac_var = attr(object, "ifrac_var", exact = TRUE),
      subgroup_label = attr(object, "subgroup_label", exact = TRUE),
      n_trials = attr(object, "n_trials", exact = TRUE),
      summary = attr(object, "summary", exact = TRUE),
      prob = attr(object, "prob", exact = TRUE) %||% prob,
      ifrac = attr(object, "ifrac", exact = TRUE) %||% "interaction",
      fit = attr(object, "interactmeta", exact = TRUE)
    ))
  }

  if (is.numeric(object) && !is.data.frame(object) && !inherits(object, "interactmeta_overall_ifrac")) {
    draws <- as.numeric(object)
    if (!length(draws) || !all(is.finite(draws))) {
      stop("`object` must be a finite numeric vector of overall IF draws.", call. = FALSE)
    }
    summ <- overall_ifrac_summary_df(draws, prob = prob)
    return(list(
      support = draws,
      trial_pi = NULL,
      reference = reference,
      ifrac_var = ifrac_var,
      subgroup_label = overall_ifrac_subgroup_label(ifrac_var),
      n_trials = NA_integer_,
      summary = summ,
      prob = prob,
      ifrac = ifrac
    ))
  }

  if (is_comparison_fit(object)) {
    check_comparison_fit(fit)
    st <- attr(object, "interactmeta", exact = TRUE)
    trial_var <- st$trial_var %||% interactmeta_study_column()
    subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
    ref <- st$reference
    d <- st$data_fit
    if (is.null(d) && !is.null(st$ifrac_data)) {
      d <- st$ifrac_data
      if (is.null(ifrac_var)) {
        ifrac_var <- st$ifrac_var
        if (is.null(ifrac_var) || !ifrac_var %in% names(d)) {
          ifrac_var <- detect_ifrac_vars(d)[[1]]
        }
      }
    } else if (is.null(ifrac_var)) {
      ifrac_var <- st$ifrac_vars_std[[1]]
      if (is.null(ifrac_var) && !is.null(st$ifrac_var)) {
        ifrac_var <- st$ifrac_var
      }
      if (is.null(ifrac_var) || !ifrac_var %in% names(d)) {
        ifrac_var <- detect_ifrac_vars(d)[[1]]
      }
    }
    draws <- overall_ifrac_from_fit(
      object,
      ifrac_var = ifrac_var,
      ifrac = ifrac,
      summary = FALSE,
      prob = prob
    )
    # Prefer attributes already attached by overall_ifrac_from_fit
    trial_pi <- attr(draws, "trial_pi", exact = TRUE)
    if (is.null(trial_pi)) {
      trial_pi <- overall_ifrac_trial_pi(
        d,
        trial_var = trial_var,
        subgroup_var = subgroup_var,
        reference = ref,
        ifrac_var = ifrac_var
      )
    }
    summ <- attr(draws, "summary", exact = TRUE)
    if (is.null(summ)) {
      summ <- overall_ifrac_summary_df(draws, prob = prob)
    }
    return(list(
      support = draws,
      trial_pi = trial_pi,
      reference = ref,
      ifrac_var = ifrac_var,
      subgroup_label = overall_ifrac_subgroup_label(ifrac_var, interactmeta = st),
      n_trials = if (is.data.frame(trial_pi)) nrow(trial_pi) else length(trial_pi),
      summary = summ,
      prob = prob,
      ifrac = ifrac,
      fit = object
    ))
  }

  ref <- reference
  if (is.null(ref) && "reference" %in% names(object)) {
    ref <- unique(object$reference)
    ref <- ref[!is.na(ref) & nzchar(ref)]
    if (length(ref) == 1L) {
      ref <- ref[[1]]
    }
  }
  if (is.null(ifrac_var)) {
    ifrac_var <- detect_ifrac_vars(object)[[1]]
  }
  draws <- overall_ifrac_from_data(
    object,
    fit = fit,
    tau = tau,
    tau_gamma = tau_gamma,
    ifrac = ifrac,
    reference = ref,
    trial = trial,
    subgroup = subgroup,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep,
    summary = FALSE,
    prob = prob
  )
  trial_pi <- overall_ifrac_trial_pi(
    object,
    trial_var = trial,
    subgroup_var = subgroup,
    reference = ref,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  list(
    support = draws,
    trial_pi = trial_pi,
    reference = ref,
    ifrac_var = ifrac_var,
    subgroup_label = overall_ifrac_subgroup_label(ifrac_var),
    n_trials = nrow(trial_pi),
    summary = overall_ifrac_summary_df(draws, prob = prob),
    prob = prob,
    ifrac = ifrac
  )
}

#' @noRd
overall_ifrac_summary_df <- function(draws, prob = 0.95, weight = NULL) {
  if (is.null(weight)) {
    weight <- attr(draws, ".weight", exact = TRUE)
  }
  mom <- comparison_draw_moments(as.numeric(draws), prob = prob, weight = weight)
  data.frame(
    estimate = mom$estimate,
    median = mom$median,
    mean = mom$mean,
    sd = mom$sd,
    lci = mom$lci,
    uci = mom$uci,
    prob = mom$prob,
    stringsAsFactors = FALSE
  )
}

#' Align quadrature node weights with IF draws from the same fit
#' @noRd
interactmeta_ifrac_draw_weights <- function(fit, draws) {
  w <- as.data.frame(fit$support)$.weight
  if (is.null(w)) {
    return(NULL)
  }
  if (length(w) == length(draws)) {
    return(as.numeric(w))
  }
  NULL
}

#' One-line summary for overall-IF plot captions (IF %, interval, study count)
#' @noRd
overall_ifrac_footnote <- function(summary,
                                   prob = NULL,
                                   n_trials = NULL,
                                   digits = 2,
                                   plugin = FALSE,
                                   ifrac = c("interaction", "overall")) {
  ifrac <- match.arg(ifrac)
  if (is.null(prob)) {
    prob <- if (!is.null(summary$prob)) summary$prob[[1L]] else 0.95
  }
  pct <- function(x) {
    format(
      round(100 * as.numeric(x), digits),
      nsmall = digits,
      trim = TRUE
    )
  }
  est_pct <- pct(summary$estimate)
  has_int <- !isTRUE(plugin) &&
    "lci" %in% names(summary) &&
    "uci" %in% names(summary) &&
    is.finite(as.numeric(summary$lci)[[1L]]) &&
    is.finite(as.numeric(summary$uci)[[1L]])
  out <- paste0(overall_ifrac_weight_label(ifrac), ": ", est_pct, "%")
  if (has_int) {
    out <- paste0(
      out,
      " [",
      pct(summary$lci),
      "%, ",
      pct(summary$uci),
      "%]"
    )
  } else {
    out <- paste0(out, " (at fitted ", "\u03c4", ")")
  }
  if (!is.null(n_trials) && !is.na(n_trials)) {
    out <- paste0(out, " \u00b7 k = ", n_trials, " studies")
  }
  out
}

#' @noRd
overall_ifrac_footnote_bundle <- function(bundle, digits = 2) {
  plugin <- FALSE
  if (!is.null(bundle$fit)) {
    plugin <- isTRUE(interactmeta_ifrac_is_plugin(bundle$fit))
  }
  overall_ifrac_footnote(
    summary = bundle$summary,
    prob = bundle$prob %||% 0.95,
    n_trials = bundle$n_trials,
    digits = digits,
    plugin = plugin,
    ifrac = bundle$ifrac %||% "interaction"
  )
}

#' @noRd
new_interactmeta_overall_ifrac <- function(draws, bundle, prob = 0.95) {
  obj <- draws
  class(obj) <- c("interactmeta_overall_ifrac", "numeric")
  attr(obj, "trial_pi") <- bundle$trial_pi
  attr(obj, "summary") <- bundle$summary
  attr(obj, "reference") <- bundle$reference
  attr(obj, "ifrac_var") <- bundle$ifrac_var
  attr(obj, "subgroup_label") <- bundle$subgroup_label
  attr(obj, "n_trials") <- bundle$n_trials
  attr(obj, "prob") <- prob
  attr(obj, "ifrac") <- bundle$ifrac %||% "interaction"
  attr(obj, "interactmeta") <- bundle$fit
  obj
}

#' @noRd
overall_ifrac_result <- function(draws, summary = FALSE, prob = 0.95, bundle = NULL) {
  if (summary) {
    return(overall_ifrac_summary_df(draws, prob = prob))
  }
  if (is.null(bundle)) {
    bundle <- list(
      support = draws,
      trial_pi = NULL,
      reference = NULL,
      ifrac_var = NULL,
      subgroup_label = NULL,
      n_trials = NA_integer_,
      summary = overall_ifrac_summary_df(draws, prob = prob),
      prob = prob
    )
  }
  new_interactmeta_overall_ifrac(draws, bundle = bundle, prob = prob)
}

#' Print overall information fraction summary
#'
#' @param x A `interactmeta_overall_ifrac` object.
#' @param digits Decimal places (default `2`).
#' @param ... Ignored.
#'
#' @export
print.interactmeta_overall_ifrac <- function(x, digits = 2, ...) {
  summ <- attr(x, "summary", exact = TRUE)
  prob <- attr(x, "prob", exact = TRUE) %||% 0.95
  if (is.null(summ)) {
    summ <- overall_ifrac_summary_df(x, prob = prob)
  }
  scheme <- attr(x, "ifrac", exact = TRUE) %||% "interaction"
  style <- interactmeta_interval_style(x)
  plugin_fit <- identical(style, "frequentist")
  summ <- as.data.frame(summ, stringsAsFactors = FALSE)
  summ$IF <- overall_ifrac_print_row_label(scheme)
  comparison_print_effects_table(
    summ,
    digits = digits,
    prob = prob,
    id_cols = "IF",
    title = NULL,
    interval_scope = "single",
    style = style
  )
  if (plugin_fit) {
    tau_lab <- if (identical(scheme, "interaction")) {
      paste0(
        "at fitted interaction heterogeneity ",
        "\u03c4",
        "_",
        "\u03b3"
      )
    } else {
      paste0("at fitted heterogeneity ", "\u03c4")
    }
    cat("(", tau_lab, ")\n", sep = "")
  }
  invisible(x)
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- prepare_ifrac_input.R ----
is_ifrac_output <- function(x) {
  if (!is.data.frame(x)) {
    return(FALSE)
  }
  cols <- grep("^ifrac_", names(x), value = TRUE)
  cols <- cols[!grepl("^ifrac2_", cols)]
  length(cols) > 0L
}

prepare_ifrac_input <- function(x,
                                reference = NULL,
                                trial = "study",
                                subgroup = "level",
                                contrast = "all",
                                contrast_sep = "-",
                                big_var = 100^2,
                                include_combined = FALSE,
                                include_interaction = FALSE) {
  if (is_ifrac_output(x)) {
    out <- tibble::as_tibble(x)
    # Already designed: reference must agree with the stored design (ifrac_<nonref>).
    ref <- interactmeta_reference(out, reference = reference, subgroup_var = subgroup)
    out <- interactmeta_ensure_reference_column(out, ref)
    # Reject ifrac_<reference> columns (pi must be non-reference only).
    ifrac_cols <- grep("^ifrac_", names(out), value = TRUE)
    ifrac_cols <- ifrac_cols[!grepl("^ifrac2_", ifrac_cols)]
    bad <- ifrac_cols[sub("^ifrac_", "", ifrac_cols) == ref]
    if (length(bad)) {
      stop(
        "Design column(s) ",
        paste(bad, collapse = ", "),
        " name the reference subgroup `",
        ref,
        "`. CAMS `ifrac_*` columns must be non-reference only ",
        "(pi = IF_nonref, reference weight = 1 - pi).",
        call. = FALSE
      )
    }
    level_col <- interactmeta_resolve_level_col(out, fallback = subgroup)
    return(add_ifrac2_design(out, subgroup_col = level_col))
  }

  if (is.null(reference)) {
    stop(
      "`reference` is required when `x` is not already output from `ifrac()`.",
      call. = FALSE
    )
  }

  designed <- ifrac(
    x,
    trial = trial,
    subgroup = subgroup,
    reference = reference,
    contrast = contrast,
    contrast_sep = contrast_sep,
    big_var = big_var,
    include_combined = include_combined,
    include_interaction = include_interaction
  )
  level_col <- interactmeta_resolve_level_col(designed, fallback = subgroup)
  add_ifrac2_design(designed, subgroup_col = level_col)
}
