#' Linear predictor table for a CAMS fit
#'
#' Returns the stored `$mu` matrix (quadrature nodes or a point-estimate row
#' × observations). This is the fixed-effect linear predictor, not trial-level
#' reconstruction; see [shrinkage()].
#'
#' @param object An `interactmeta` from [cam()].
#' @param summary If `TRUE`, return a data frame of medians and SDs by
#'   observation. If `FALSE`, return the numeric matrix.
#' @param ... Ignored.
#' @return A data frame when `summary = TRUE`, otherwise a numeric matrix.
#' @export
fitted.interactmeta <- function(object, summary = TRUE, ...) {
  mu <- object$mu
  if (is.null(mu)) {
    stop("`interactmeta` is missing linear-predictor values (`mu`).", call. = FALSE)
  }
  mu <- matrix(as.numeric(mu), nrow = nrow(mu), ncol = ncol(mu))
  if (isTRUE(summary)) {
    est <- apply(mu, 2L, function(col) comparison_median_qi_vec(col, prob = 0.95)$estimate)
    return(data.frame(
      Estimate = est,
      Est.Error = apply(mu, 2L, stats::sd),
      stringsAsFactors = FALSE
    ))
  }
  mu
}

#' Parameter table of a CAMS fit
#'
#' Returns the stored `$draws` data frame: DIRECT quadrature nodes for
#' bayesmeta, or a single point-estimate row for metafor.
#'
#' @param x An `interactmeta` from [cam()].
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#' @return A data frame.
#' @export
as.data.frame.interactmeta <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$draws
}

#' @noRd
comparison_fit_variables <- function(fit) {
  if (inherits(fit, "interactmeta")) {
    return(names(fit$draws))
  }
  names(as.data.frame(fit))
}
