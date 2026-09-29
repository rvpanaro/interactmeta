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
                             trial_col = "trial",
                             subgroup_col = "subgroup",
                             vi_col = "vi",
                             reference_col = "reference",
                             include_contrast_rows = FALSE,
                             contrast_sep = "-",
                             subgroup_levels = NULL) {
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
add_ifrac2_design <- function(dat, subgroup_col = "subgroup") {
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
