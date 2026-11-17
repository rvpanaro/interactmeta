# Fast Bayesian cams() for testthat (quadrature, not MCMC).
# Override grid coarseness: Sys.setenv(INTERACTMETA_TEST_N_SUPPORT = "40").

fast_bayes_n_support <- function() {
  raw <- Sys.getenv("INTERACTMETA_TEST_N_SUPPORT", unset = "25")
  n <- suppressWarnings(as.integer(raw)[[1L]])
  if (!is.finite(n) || n < 8L) {
    25L
  } else {
    n
  }
}

fast_bayes_fit_args <- function() {
  list(n_support = fast_bayes_n_support())
}

fast_bayes_cams <- function(x, ..., reference, silent = 2L, fit_args = NULL) {
  base <- fast_bayes_fit_args()
  if (!is.null(fit_args)) {
    base <- utils::modifyList(base, fit_args, keep.null = TRUE)
  }
  dots <- list(...)
  if ("fit_args" %in% names(dots)) {
    base <- utils::modifyList(base, dots$fit_args, keep.null = TRUE)
    dots$fit_args <- NULL
  }
  args <- c(
    list(x = x, reference = reference, silent = silent, fit_args = base),
    dots
  )
  suppressWarnings(do.call(cams, args))
}