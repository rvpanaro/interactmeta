#' Parameter table as a data frame
#'
#' Uses stored `$draws` (DIRECT nodes for bayesmeta; one point-estimate row
#' for metafor). Adds `.draw` / `.chain` / `.iteration` when missing so
#' grouped summaries can index rows.
#' @noRd
comparison_tidy_draws <- function(fit) {
  post <- if (inherits(fit, "interactmeta")) {
    as.data.frame(fit$draws)
  } else {
    as.data.frame(fit)
  }
  n <- nrow(post)
  if (!".draw" %in% names(post)) {
    post$.draw <- seq_len(n)
  }
  if (!".chain" %in% names(post)) {
    post$.chain <- 1L
  }
  if (!".iteration" %in% names(post)) {
    post$.iteration <- seq_len(n)
  }
  post
}

#' Median + equal-tail interval for a numeric vector
#'
#' Canonical moments (median / mean / sd / equal-tail QI).
#' @noRd
comparison_draw_moments <- function(x, prob = 0.95, weight = NULL, cond_sd = NULL,
                                    exponentiate = FALSE) {
  x <- as.numeric(x)
  n_raw <- length(x)
  has_w <- !is.null(weight) && length(weight) == n_raw &&
    length(unique(as.numeric(weight)[is.finite(as.numeric(weight))])) > 1L
  has_s <- !is.null(cond_sd) && any(as.numeric(cond_sd) > 0, na.rm = TRUE)
  if (has_w || has_s) {
    return(interactmeta_mixture_moments(
      mean = x,
      weight = weight,
      cond_sd = cond_sd,
      prob = prob,
      exponentiate = exponentiate
    ))
  }
  if (isTRUE(exponentiate)) {
    x <- exp(x)
  }
  x <- x[is.finite(x)]
  empty <- list(
    median = NA_real_,
    mean = NA_real_,
    sd = NA_real_,
    lci = NA_real_,
    uci = NA_real_,
    estimate = NA_real_,
    prob = prob
  )
  if (!length(x)) {
    return(empty)
  }
  if (length(x) == 1L) {
    v <- x[[1L]]
    return(list(
      median = v,
      mean = v,
      sd = NA_real_,
      lci = NA_real_,
      uci = NA_real_,
      estimate = v,
      prob = prob
    ))
  }
  qi <- comparison_median_qi_core(x, prob = prob)
  list(
    median = qi$estimate,
    mean = mean(x),
    sd = stats::sd(x),
    lci = qi$lci,
    uci = qi$uci,
    estimate = qi$estimate,
    prob = prob
  )
}

#' Equal-tail median + QI
#' @noRd
comparison_median_qi_core <- function(x, prob = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  alpha <- (1 - prob) / 2
  qs <- stats::quantile(x, probs = c(0.5, alpha, 1 - alpha), names = FALSE, type = 7)
  list(
    estimate = qs[[1L]],
    lci = qs[[2L]],
    uci = qs[[3L]],
    prob = prob
  )
}

#' @return list(estimate, lci, uci, prob)
#' @noRd
comparison_median_qi_vec <- function(x, prob = 0.95) {
  comparison_median_qi_core(x, prob = prob)
}

#' Grouped median / mean / sd / equal-tail interval summary
#' @noRd
comparison_median_qi_grouped <- function(draws, value_col, group_cols, prob = 0.95,
                                         exponentiate = FALSE) {
  draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(d, ...) {
      w <- if (".weight" %in% names(d)) d$.weight else NULL
      s <- if (".cond_sd" %in% names(d)) d$.cond_sd else NULL
      mom <- comparison_draw_moments(
        d[[value_col]],
        prob = prob,
        weight = w,
        cond_sd = s,
        exponentiate = exponentiate
      )
      tibble::tibble(
        estimate = mom$estimate,
        mean = mom$mean,
        sd = mom$sd,
        lci = mom$lci,
        uci = mom$uci,
        prob = prob
      )
    }) |>
    dplyr::ungroup()
}
