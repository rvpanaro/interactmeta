#' Trial-level IV overall and interaction contrast blocks
#'
#' Single source for quoted (observed) trial-level CAMS quantities.
#'
#' CAMS overall rotation (IV pool):
#' `m_j = sum_g w_g y_gj` with `w_g = 1 / v_gj`, `Var(m_j) = 1 / sum_g 1/v_gj`.
#'
#' CAMS interaction contrast:
#' `g_j = y_Bj - y_Aj`.
#' Variance uses `v_Aj + v_Bj` (within-trial subgroup covariance **assumed zero**).
#'
#' Information fraction of B:
#' `pi_j = (1/v_Bj) / (1/v_Aj + 1/v_Bj)`.
#'
#' CAMS orthogonality: when `pi_j` is that IV weight,
#' `Cov(m_j, g_j) = -(1-pi_j) v_Aj + pi_j v_Bj = 0`.
#'
#' @noRd
comparison_trial_effect_blocks <- function(data,
                                     trial_var = "trial",
                                     subgroup_var = "subgroup",
                                     reference,
                                     nonref = NULL,
                                     ifrac_var = NULL,
                                     min_vi = 1e-12) {
  dat <- as.data.frame(data)
  dat[[trial_var]] <- as.character(dat[[trial_var]])
  dat[[subgroup_var]] <- as.character(dat[[subgroup_var]])
  ref <- as.character(reference)[[1L]]
  trials <- unique(dat[[trial_var]])

  if (is.null(nonref)) {
    sg_all <- unique(dat[[subgroup_var]])
    nonrefs <- setdiff(sg_all, ref)
    if (!length(nonrefs)) {
      stop("No non-reference subgroup found for trial effect blocks.", call. = FALSE)
    }
    nonref <- nonrefs[[1L]]
  }
  nonref <- as.character(nonref)[[1L]]

  overall_rows <- vector("list", length(trials))
  interaction_rows <- vector("list", length(trials))
  pi <- stats::setNames(rep(NA_real_, length(trials)), trials)

  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    d <- dat[dat[[trial_var]] == tr, , drop = FALSE]
    if (!nrow(d)) {
      next
    }
    yi <- as.numeric(d$yi)
    vi_raw <- pmax(as.numeric(d$vi), 0)
    if (all(vi_raw <= 0)) {
      ov_yi <- mean(yi)
      ov_vi <- 0
    } else {
      # CAMS overall rotation: precision-weighted mean of subgroup estimates.
      vi <- pmax(vi_raw, min_vi)
      w <- 1 / vi
      ov_yi <- sum(w * yi) / sum(w)
      ov_vi <- 1 / sum(w)
      if (all(vi_raw == 0)) {
        ov_vi <- 0
      }
    }
    ov <- data.frame(
      trial = tr,
      yi = ov_yi,
      vi = ov_vi,
      stringsAsFactors = FALSE
    )
    if (!is.null(ifrac_var) && ifrac_var %in% names(d)) {
      pi_tr <- mean(as.numeric(d[[ifrac_var]]), na.rm = TRUE)
      ov$ifrac <- pi_tr
      pi[[tr]] <- pi_tr
    }
    overall_rows[[j]] <- ov

    i_ref <- which(d[[subgroup_var]] == ref)
    i_non <- which(d[[subgroup_var]] == nonref)
    if (length(i_ref) != 1L || length(i_non) != 1L) {
      stop(
        "Each trial must have exactly one row for reference and non-reference ",
        "subgroups; trial `",
        tr,
        "` does not.",
        call. = FALSE
      )
    }
    # CAMS interaction: g_j = y_Bj - y_Aj
    # Var(g_j) = v_Aj + v_Bj when Cov(y_Aj, y_Bj) is assumed zero.
    interaction_rows[[j]] <- data.frame(
      trial = tr,
      yi = as.numeric(d$yi[[i_non]]) - as.numeric(d$yi[[i_ref]]),
      vi = pmax(as.numeric(d$vi[[i_non]]), 0) + pmax(as.numeric(d$vi[[i_ref]]), 0),
      stringsAsFactors = FALSE
    )
  }

  list(
    overall = dplyr::bind_rows(overall_rows),
    interaction = dplyr::bind_rows(interaction_rows),
    pi = pi,
    reference = ref,
    nonref = nonref,
    ifrac_var = ifrac_var
  )
}
