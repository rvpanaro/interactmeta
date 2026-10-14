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
#' pooling is a plug-in at the fitted \eqn{\tau} or \eqn{\tau_\gamma}.
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
#' fit <- cam(d, reference = "Lower")
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
                             trial = "trial",
                             modifier = NULL,
                             subgroup = "subgroup",
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
  blocks <- comparison_trial_effect_blocks(
    data = dat_complete,
    trial_var = "trial",
    subgroup_var = "subgroup",
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
  cols <- c("trial", "subgroup", "yi", "vi")
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
  if (!("subgroup" %in% names(out))) {
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
    is_subgroup_row <- is_ifrac_component_row(out$subgroup, contrast_sep = contrast_sep)
  }

  if (all(is_subgroup_row)) {
    out$subgroup <- factor(out$subgroup, levels = sg_levels)
    stopifnot(identical(levels(out$subgroup)[[1]], ref))
  } else if (any(is_subgroup_row)) {
    out$subgroup <- as.character(out$subgroup)
    out$subgroup[is_subgroup_row] <- as.character(
      factor(out$subgroup[is_subgroup_row], levels = sg_levels)
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
      trial <- "trial"
    }
    if (is.null(subgroup)) {
      subgroup <- "subgroup"
    }

    if (!is.character(trial) || length(trial) != 1L || !(trial %in% names(dat))) {
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
