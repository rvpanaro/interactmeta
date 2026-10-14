#' Independent-block quadrature helpers for binary bayesmeta CAMS
#'
#' DIRECT support from `bmr` / `bayesmeta` is a mixture of normals (plus a
#' discrete distribution on tau). Independent blocks combine by a product of
#' mixture weights. That Gaussian-mixture product is an exact check on sums of
#' independent blocks. The published 1-D law of those sums is
#' `bayesmeta::convolve()` of his density and cdf, and the published interval
#' is his shortest or central rule on that convolution (or `$post.interval`
#' when the margin lives on one stored fit).
#'
#' These objects are **marginal** posteriors. Trials that share a block also
#' share hyperparameters and are not mutually independent.
#'
#' @noRd
NULL

interactmeta_log_sum_exp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(-Inf)
  }
  m <- max(x)
  m + log(sum(exp(x - m)))
}

interactmeta_gm <- function(mean, sd, weight) {
  mean <- as.numeric(mean)
  sd <- pmax(as.numeric(sd), 0)
  weight <- as.numeric(weight)
  n <- length(mean)
  if (length(sd) != n || length(weight) != n) {
    stop("Gaussian mixture `mean`, `sd`, and `weight` must have equal length.", call. = FALSE)
  }
  ok <- is.finite(mean) & is.finite(sd) & is.finite(weight) & weight > 0
  if (!any(ok)) {
    stop("Gaussian mixture has no finite positive-weight components.", call. = FALSE)
  }
  weight <- weight[ok]
  weight <- weight / sum(weight)
  structure(
    list(mean = mean[ok], sd = sd[ok], weight = weight),
    class = "interactmeta_gm"
  )
}

interactmeta_gm_sd_eval <- function(sd, min_sd = 1e-12) {
  pmax(as.numeric(sd), min_sd)
}

interactmeta_gm_density <- function(gm, x, log = FALSE) {
  x <- as.numeric(x)
  sd <- interactmeta_gm_sd_eval(gm$sd)
  logw <- log(gm$weight)
  out <- vapply(x, function(z) {
    interactmeta_log_sum_exp(logw + stats::dnorm(z, mean = gm$mean, sd = sd, log = TRUE))
  }, numeric(1))
  if (isTRUE(log)) {
    out
  } else {
    exp(out)
  }
}

interactmeta_gm_cdf <- function(gm, x) {
  x <- as.numeric(x)
  sd <- interactmeta_gm_sd_eval(gm$sd)
  n <- length(x)
  k <- length(gm$weight)
  cdf <- matrix(NA_real_, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    cdf[, j] <- stats::pnorm(x, mean = gm$mean[[j]], sd = sd[[j]])
  }
  as.numeric(cdf %*% gm$weight)
}

interactmeta_gm_quantile <- function(gm, probs) {
  probs <- as.numeric(probs)
  mu <- sum(gm$weight * gm$mean)
  v <- sum(gm$weight * (gm$sd^2 + (gm$mean - mu)^2))
  s <- sqrt(max(v, 0))
  span <- max(s, 1e-8)
  vapply(probs, function(p) {
    if (!is.finite(p)) {
      return(NA_real_)
    }
    if (p <= 0) {
      return(-Inf)
    }
    if (p >= 1) {
      return(Inf)
    }
    lo <- mu - 8 * span
    hi <- mu + 8 * span
    f_lo <- interactmeta_gm_cdf(gm, lo) - p
    f_hi <- interactmeta_gm_cdf(gm, hi) - p
    tries <- 0L
    while (f_lo > 0 && tries < 20L) {
      lo <- lo - (hi - lo + span)
      f_lo <- interactmeta_gm_cdf(gm, lo) - p
      tries <- tries + 1L
    }
    tries <- 0L
    while (f_hi < 0 && tries < 20L) {
      hi <- hi + (hi - lo + span)
      f_hi <- interactmeta_gm_cdf(gm, hi) - p
      tries <- tries + 1L
    }
    stats::uniroot(
      function(z) interactmeta_gm_cdf(gm, z) - p,
      lower = lo,
      upper = hi,
      tol = 1e-8
    )$root
  }, numeric(1))
}

interactmeta_gm_mean <- function(gm) {
  sum(gm$weight * gm$mean)
}

interactmeta_gm_var <- function(gm) {
  mu <- interactmeta_gm_mean(gm)
  sum(gm$weight * (gm$sd^2 + (gm$mean - mu)^2))
}

interactmeta_gm_summary <- function(gm, prob = 0.95,
                                    interval.type = "shortest") {
  interval.type <- interactmeta_match_interval_type(interval.type)
  mu <- interactmeta_gm_mean(gm)
  s <- sqrt(max(interactmeta_gm_var(gm), 0))
  median <- interactmeta_gm_quantile(gm, 0.5)
  bounds <- interactmeta_interval_from_distribution(
    quantile = function(p) interactmeta_gm_quantile(gm, p),
    cdf = function(x) interactmeta_gm_cdf(gm, x),
    prob = prob,
    interval.type = interval.type
  )
  qs <- c(bounds[[1L]], median[[1L]], bounds[[2L]])
  list(
    mean = mu,
    sd = s,
    median = unname(median[[1L]]),
    lower = unname(qs[[1L]]),
    upper = unname(qs[[3L]]),
    q = stats::setNames(qs, c("lower", "median", "upper")),
    interval.type = interval.type
  )
}

#' Mean and SD of exp(X) when X is a Gaussian mixture (mixture of lognormals)
#' @noRd
interactmeta_gm_lognormal_moments <- function(gm) {
  mean_y <- sum(gm$weight * exp(gm$mean + 0.5 * gm$sd^2))
  mean_y2 <- sum(gm$weight * exp(2 * gm$mean + 2 * gm$sd^2))
  list(
    mean = mean_y,
    sd = sqrt(max(mean_y2 - mean_y^2, 0))
  )
}

#' Linear combination aX + bY of **independent** Gaussian mixtures
#'
#' The product of mixture weights assumes `gm_x` and `gm_y` come from
#' independent blocks. Valid for cross-block CAMS quantities such as
#' `theta - gamma` and trial-level `m_j ± c g_j`.
#'
#' Do **not** use this to form `M_star = alpha + pi_star * theta`. `alpha`
#' and `theta` share the overall-block DIRECT support and can be
#' correlated; that linear combination must keep the paired support
#' (and `Cov(alpha, theta)`), not a product of marginal mixtures.
#' @noRd
interactmeta_gm_linear_sum <- function(gm_x, a, gm_y, b) {
  interactmeta_gm_assert_independent(gm_x, gm_y)
  a <- as.numeric(a)[[1L]]
  b <- as.numeric(b)[[1L]]
  if (!is.finite(a) || !is.finite(b)) {
    stop("Linear combination coefficients must be finite.", call. = FALSE)
  }
  ix <- seq_along(gm_x$weight)
  iy <- seq_along(gm_y$weight)
  i <- rep(ix, each = length(iy))
  j <- rep(iy, times = length(ix))
  mean <- a * gm_x$mean[i] + b * gm_y$mean[j]
  sd <- sqrt((a * gm_x$sd[i])^2 + (b * gm_y$sd[j])^2)
  weight <- gm_x$weight[i] * gm_y$weight[j]
  interactmeta_gm(mean = mean, sd = sd, weight = weight)
}

#' Tag a Gaussian mixture with its CAMS block (`"overall"` / `"interaction"`)
#' @noRd
interactmeta_gm_tag_block <- function(gm, block) {
  attr(gm, "block") <- block
  gm
}

#' Stop if both mixtures are tagged as the same CAMS block
#' @noRd
interactmeta_gm_assert_independent <- function(gm_x, gm_y) {
  bx <- attr(gm_x, "block", exact = TRUE)
  by <- attr(gm_y, "block", exact = TRUE)
  if (is.null(bx) || is.null(by) || !identical(bx, by)) {
    return(invisible(NULL))
  }
  stop(
    "`interactmeta_gm_linear_sum()` assumes independent mixtures. ",
    "Do not combine two components from the same CAMS block ",
    "(e.g. `alpha` and `theta` from the overall `bmr` support). ",
    "Form `M* = alpha + pi* theta` from the joint support, not by ",
    "convolving the two marginals.",
    call. = FALSE
  )
}

#' Discrete (point-mass) distribution on a support, e.g. tau nodes
#' @noRd
interactmeta_discrete <- function(x, weight) {
  x <- as.numeric(x)
  weight <- as.numeric(weight)
  ok <- is.finite(x) & is.finite(weight) & weight > 0
  weight <- weight[ok]
  weight <- weight / sum(weight)
  list(x = x[ok], weight = weight)
}

interactmeta_discrete_summary <- function(dist, prob = 0.95) {
  mu <- sum(dist$weight * dist$x)
  s <- sqrt(max(sum(dist$weight * (dist$x - mu)^2), 0))
  qs <- interactmeta_weighted_quantile(
    dist$x,
    dist$weight,
    probs = c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2)
  )
  list(
    mean = mu,
    sd = s,
    median = unname(qs[[2L]]),
    lower = unname(qs[[1L]]),
    upper = unname(qs[[3L]]),
    q = qs
  )
}

#' Density / CDF of Z = c * X, with Jacobian. c = 0 is degenerate at 0.
#' @noRd
interactmeta_call_density <- function(dens, x, log = FALSE) {
  if ("log" %in% names(formals(dens))) {
    dens(x, log = log)
  } else if (isTRUE(log)) {
    log(dens(x))
  } else {
    dens(x)
  }
}

interactmeta_scale_density <- function(dens, c) {
  c <- as.numeric(c)[[1L]]
  if (!is.finite(c)) {
    stop("Scale coefficient must be finite.", call. = FALSE)
  }
  if (abs(c) < 1e-15) {
    stop("Scale coefficient is zero; use a degenerate point mass at 0.", call. = FALSE)
  }
  function(x, log = FALSE) {
    z <- as.numeric(x) / c
    if (isTRUE(log)) {
      as.numeric(interactmeta_call_density(dens, z, log = TRUE) - log(abs(c)))
    } else {
      as.numeric(interactmeta_call_density(dens, z, log = FALSE) / abs(c))
    }
  }
}

interactmeta_scale_cdf <- function(cdf, c) {
  c <- as.numeric(c)[[1L]]
  if (!is.finite(c)) {
    stop("Scale coefficient must be finite.", call. = FALSE)
  }
  if (abs(c) < 1e-15) {
    stop("Scale coefficient is zero; use a degenerate point mass at 0.", call. = FALSE)
  }
  if (c > 0) {
    function(x) as.numeric(cdf(as.numeric(x) / c))
  } else {
    function(x) as.numeric(1 - cdf(as.numeric(x) / c))
  }
}

interactmeta_negate_density <- function(dens) {
  interactmeta_scale_density(dens, -1)
}

interactmeta_negate_cdf <- function(cdf) {
  interactmeta_scale_cdf(cdf, -1)
}

#' Normal-normal shrinkage conditional on hyperparameters
#'
#' tau = 0 is handled without 1/tau^2: the latent equals mu0 with variance 0
#' given the hyperparameters.
#' @noRd
interactmeta_normal_normal_shrink <- function(y, sigma2, mu0, tau) {
  y <- as.numeric(y)
  sigma2 <- as.numeric(sigma2)
  mu0 <- as.numeric(mu0)
  tau <- as.numeric(tau)
  n <- max(length(y), length(sigma2), length(mu0), length(tau))
  y <- rep_len(y, n)
  sigma2 <- rep_len(sigma2, n)
  mu0 <- rep_len(mu0, n)
  tau <- rep_len(tau, n)
  mean <- mu0
  sd <- rep(0, n)
  pos <- is.finite(tau) & tau > 0 & is.finite(sigma2) & sigma2 > 0
  if (any(pos)) {
    prec_y <- 1 / sigma2[pos]
    prec_0 <- 1 / (tau[pos]^2)
    post_prec <- prec_y + prec_0
    sd[pos] <- sqrt(1 / post_prec)
    mean[pos] <- (prec_y * y[pos] + prec_0 * mu0[pos]) / post_prec
  }
  zero_sig <- is.finite(sigma2) & sigma2 <= 0 & is.finite(tau) & tau > 0
  if (any(zero_sig)) {
    mean[zero_sig] <- y[zero_sig]
    sd[zero_sig] <- 0
  }
  cbind(mean = mean, sd = sd)
}

#' Extract mixture nodes from a bayesmeta() support matrix (mu | tau)
#'
#' Untagged on purpose: independent univariate `bayesmeta` fits (DA
#' `mu_B - mu_A`) may be combined with `interactmeta_gm_linear_sum()`.
#' The CAMS interaction mean is tagged `"interaction"` where it is
#' constructed in `interactmeta_bayesmeta_latent()`.
#' @noRd
interactmeta_gm_from_bayesmeta <- function(fit) {
  sup <- as.data.frame(fit$support)
  if (!"mu" %in% names(sup) && "mean" %in% names(sup)) {
    sup$mu <- sup$mean
  }
  interactmeta_gm(
    mean = sup$mu,
    sd = sup$sd,
    weight = sup$weight
  )
}

#' Extract mixture nodes for one bmr() coefficient
#' @noRd
interactmeta_bmr_coef_name <- function(fit, role = c("intercept", "slope")) {
  role <- match.arg(role)
  nms <- colnames(fit$support$mean)
  if (is.null(nms) || !length(nms)) {
    stop("`bmr` support$mean has no column names.", call. = FALSE)
  }
  if (identical(role, "intercept")) {
    hit <- intersect(c("intercept", "Intercept", "(Intercept)"), nms)
    if (length(hit)) {
      return(hit[[1L]])
    }
    return(nms[[1L]])
  }
  hit <- intersect(c("ifrac", "slope"), nms)
  if (length(hit)) {
    return(hit[[1L]])
  }
  if (length(nms) < 2L) {
    stop("No slope column in `bmr` support$mean.", call. = FALSE)
  }
  nms[[2L]]
}

interactmeta_gm_from_bmr <- function(fit, which) {
  mu <- as.numeric(fit$support$mean[, which])
  var <- as.numeric(fit$support$covariance[, which, which])
  interactmeta_gm_tag_block(
    interactmeta_gm(
      mean = mu,
      sd = sqrt(pmax(var, 0)),
      weight = fit$support$weight
    ),
    "overall"
  )
}

#' Mixture for latent m_j from bmr$shrink.moments
#' @noRd
interactmeta_gm_from_bmr_shrink <- function(fit, which) {
  mom <- fit$shrink.moments(which = which)
  interactmeta_gm_tag_block(
    interactmeta_gm(
      mean = mom[, "mean"],
      sd = mom[, "sd"],
      weight = fit$support$weight
    ),
    "overall"
  )
}

#' Mixture for latent g_j from bayesmeta cond.moment(individual=)
#' @noRd
interactmeta_gm_from_bayesmeta_individual <- function(fit, which) {
  tau <- as.numeric(fit$support[, "tau"])
  w <- as.numeric(fit$support[, "weight"])
  cm <- fit$cond.moment(tau = tau, individual = which)
  if (is.null(dim(cm))) {
    stop("Unexpected `cond.moment()` result for an individual study.", call. = FALSE)
  }
  if (length(dim(cm)) == 3L) {
    cm <- cm[, , 1L]
  }
  interactmeta_gm_tag_block(
    interactmeta_gm(
      mean = cm[, "mean"],
      sd = cm[, "sd"],
      weight = w
    ),
    "interaction"
  )
}

#' Unary density/CDF wrappers for convolve()
#' @noRd
interactmeta_bmr_coef_density <- function(fit, which) {
  force(fit)
  force(which)
  function(x, log = FALSE) {
    as.numeric(fit$dposterior(beta = x, which.beta = which, log = log))
  }
}

interactmeta_bmr_coef_cdf <- function(fit, which) {
  force(fit)
  force(which)
  function(x) {
    as.numeric(fit$pposterior(beta = x, which.beta = which))
  }
}

interactmeta_bayesmeta_mu_density <- function(fit) {
  force(fit)
  function(x, log = FALSE) {
    as.numeric(fit$dposterior(mu = x, log = log))
  }
}

interactmeta_bayesmeta_mu_cdf <- function(fit) {
  force(fit)
  function(x) {
    as.numeric(fit$pposterior(mu = x))
  }
}

interactmeta_bmr_shrink_density <- function(fit, which) {
  force(fit)
  force(which)
  function(x, log = FALSE) {
    as.numeric(fit$dshrink(theta = x, which = which, log = log))
  }
}

interactmeta_bmr_shrink_cdf <- function(fit, which) {
  force(fit)
  force(which)
  function(x) {
    as.numeric(fit$pshrink(theta = x, which = which))
  }
}

interactmeta_bayesmeta_individual_density <- function(fit, which) {
  force(fit)
  force(which)
  function(x, log = FALSE) {
    as.numeric(fit$dposterior(mu = x, individual = which, log = log))
  }
}

interactmeta_bayesmeta_individual_cdf <- function(fit, which) {
  force(fit)
  force(which)
  function(x) {
    as.numeric(fit$pposterior(mu = x, individual = which))
  }
}

#' Convolution of two independent continuous distributions via bayesmeta::convolve
#'
#' Inputs are density/CDF of independent random variables X and Y.
#' Output is a DIRECT approximation of the law of X+Y. Settings `delta`
#' (grid spacing) and `epsilon` (tail mass) are numerical accuracy controls,
#' not model parameters.
#' @noRd
interactmeta_convolve_sum <- function(dens1,
                                      dens2,
                                      cdf1,
                                      cdf2,
                                      delta = 0.01,
                                      epsilon = 1e-4) {
  bayesmeta::convolve(
    dens1 = dens1,
    dens2 = dens2,
    cdf1 = cdf1,
    cdf2 = cdf2,
    delta = delta,
    epsilon = epsilon
  )
}

#' `"shortest"` or `"central"`. `"evidentiary"` is a single-model bayesmeta rule.
#' @noRd
interactmeta_match_interval_type <- function(interval.type = NULL) {
  if (is.null(interval.type) || !length(interval.type) ||
        !nzchar(as.character(interval.type)[[1L]])) {
    return("shortest")
  }
  interval.type <- as.character(interval.type)[[1L]]
  if (identical(interval.type, "evidentiary")) {
    stop(
      "`interval.type = \"evidentiary\"` is defined for a single bayesmeta model, not for a derived CAMS sum.",
      call. = FALSE
    )
  }
  if (!interval.type %in% c("shortest", "central")) {
    stop(
      "`interval.type` must be \"shortest\" or \"central\".",
      call. = FALSE
    )
  }
  interval.type
}

#' Stop before ifrac() replaces a non-positive variance with a placeholder.
#'
#' A supplied `vi <= 0` is not a missing subgroup. Metafor still uses the
#' existing drop-and-placeholder path.
#' @noRd
interactmeta_reject_nonpositive_sampling_vi <- function(x) {
  if (!is.data.frame(x) && !inherits(x, "escalc")) {
    return(invisible(NULL))
  }
  if (!"vi" %in% names(x)) {
    return(invisible(NULL))
  }
  interactmeta_sigma_from_vi(x$vi)
  invisible(NULL)
}

#' Sampling standard deviations passed to bmr() / bayesmeta().
#'
#' Variances must be finite and strictly positive. No `1e-12` floor.
#' @noRd
interactmeta_sigma_from_vi <- function(vi) {
  vi <- as.numeric(vi)
  if (!length(vi) || any(!is.finite(vi)) || any(vi <= 0)) {
    stop(
      "Sampling variances passed to bayesmeta must be finite and strictly positive.",
      call. = FALSE
    )
  }
  sqrt(vi)
}

#' bmr `$post.interval` / `$qposterior` index a coefficient by position.
#'
#' A character name reaches `sum(which.beta)` inside his argument check
#' (`&` does not short-circuit), so the name is resolved here first.
#' @noRd
interactmeta_bmr_which_index <- function(object, which) {
  if (is.null(which) || is.numeric(which)) {
    return(which)
  }
  nms <- colnames(object$support$mean)
  if (is.null(nms) || !length(nms)) {
    nms <- setdiff(colnames(object$summary), "tau")
  }
  hit <- match(as.character(which)[[1L]], nms)
  if (is.na(hit)) {
    stop(
      "No bmr coefficient named ",
      deparse(as.character(which)[[1L]]),
      ".",
      call. = FALSE
    )
  }
  hit
}

#' Call `$post.interval` with that object's own formals.
#'
#' `bmr` uses `beta.level` / `which.beta` (and `shrink.interval` for a
#' study margin). `bayesmeta` uses `mu.level` / `individual`. Tau uses
#' `tau.level` on either object. Effect intervals are not pinned at 0;
#' his tau method may be, and that stays inside his function.
#' @noRd
interactmeta_post_interval <- function(object,
                                      prob,
                                      interval.type = NULL,
                                      margin = c("effect", "tau", "shrink"),
                                      which = NULL,
                                      predict = FALSE,
                                      individual = FALSE) {
  margin <- match.arg(margin)
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% object$interval.type
  )
  prob <- as.numeric(prob)[[1L]]
  if (!is.finite(prob) || prob <= 0 || prob >= 1) {
    stop("`prob` must lie in (0, 1).", call. = FALSE)
  }
  if (identical(margin, "shrink") && is.function(object$shrink.interval)) {
    fn <- object$shrink.interval
    fml <- names(formals(fn))
    args <- list()
    if ("level" %in% fml) {
      args$level <- prob
    }
    if ("which" %in% fml) {
      args$which <- which
    }
    if ("method" %in% fml) {
      args$method <- interval.type
    }
    res <- do.call(fn, args)
    return(unname(as.numeric(res)[1:2]))
  }
  fn <- object$post.interval
  if (!is.function(fn)) {
    stop("Stored fit has no `$post.interval`.", call. = FALSE)
  }
  fml <- names(formals(fn))
  args <- list()
  if ("method" %in% fml) {
    args$method <- interval.type
  }
  if (identical(margin, "tau")) {
    args$tau.level <- prob
  } else if ("beta.level" %in% fml) {
    args$beta.level <- prob
    args$which.beta <- interactmeta_bmr_which_index(object, which)
  } else {
    args$mu.level <- prob
    if (isTRUE(predict) && "predict" %in% fml) {
      args$predict <- TRUE
    }
    if (!identical(individual, FALSE) && "individual" %in% fml) {
      args$individual <- individual
    }
  }
  res <- do.call(fn, args)
  unname(as.numeric(res)[1:2])
}

#' His quantile at one or more probabilities, using that object's formals.
#' @noRd
interactmeta_native_quantile <- function(object,
                                        p,
                                        margin = c("effect", "tau", "shrink"),
                                        which = NULL,
                                        predict = FALSE,
                                        individual = FALSE) {
  margin <- match.arg(margin)
  p <- as.numeric(p)
  if (identical(margin, "shrink") && is.function(object$qshrink)) {
    return(as.numeric(vapply(p, function(pp) {
      object$qshrink(p = pp, which = which)
    }, numeric(1))))
  }
  fn <- object$qposterior
  fml <- names(formals(fn))
  vapply(p, function(pp) {
    args <- list()
    if (identical(margin, "tau")) {
      args$tau.p <- pp
    } else if ("beta.p" %in% fml) {
      args$beta.p <- pp
      args$which.beta <- interactmeta_bmr_which_index(object, which)
    } else {
      args$mu.p <- pp
      if (isTRUE(predict) && "predict" %in% fml) {
        args$predict <- TRUE
      }
      if (!identical(individual, FALSE) && "individual" %in% fml) {
        args$individual <- individual
      }
    }
    as.numeric(do.call(fn, args))[[1L]]
  }, numeric(1))
}

#' Central or shortest interval from a posterior quantile function.
#'
#' Central is `quantile(c((1-p)/2, 1-(1-p)/2))`. Shortest minimizes
#' `q(a+p)-q(a)` over the lower-tail probability `a`. Effect parameters are
#' not forced to a boundary at 0 (his tau `$post.interval` may be; call that
#' instead of this helper).
#' @noRd
interactmeta_interval_from_distribution <- function(quantile,
                                                   cdf,
                                                   prob,
                                                   interval.type = "shortest",
                                                   tol = NULL) {
  interval.type <- interactmeta_match_interval_type(interval.type)
  prob <- as.numeric(prob)[[1L]]
  if (is.null(tol) || !is.finite(tol) || tol <= 0) {
    tol <- 2^16 * .Machine$double.eps
  }
  if (identical(interval.type, "central")) {
    return(unname(as.numeric(quantile(c((1 - prob) / 2, 1 - (1 - prob) / 2)))))
  }
  force(cdf)
  tail_max <- 1 - prob
  tail_min <- tail_max / 50
  width <- function(a) {
    bounds <- as.numeric(quantile(c(a, a + prob)))
    bounds[[2L]] - bounds[[1L]]
  }
  tail <- stats::optimize(
    width,
    lower = tail_min,
    upper = tail_max,
    tol = tol
  )$minimum
  unname(as.numeric(quantile(c(tail, tail + prob))))
}

#' One-dimensional posterior (mixture check and/or published convolve / native)
#'
#' `gm` is the exact Gaussian-mixture check. It does not supply the published
#' interval when `convolve` or `published` is set. `published` is a list
#' `kind` (`"native"` or `"convolve"`), plus the stored fit and margin
#' selectors for `"native"`.
#' @noRd
interactmeta_posterior_1d <- function(gm = NULL,
                                      convolve = NULL,
                                      density = NULL,
                                      cdf = NULL,
                                      quantile = NULL,
                                      name = NULL,
                                      note = NULL,
                                      published = NULL,
                                      interval.type = NULL,
                                      tol.uniroot = NULL) {
  if (is.null(density)) {
    if (!is.null(convolve)) {
      density <- convolve$density
    } else if (!is.null(gm)) {
      density <- function(x, log = FALSE) interactmeta_gm_density(gm, x, log = log)
    }
  }
  if (is.null(cdf)) {
    if (!is.null(convolve)) {
      cdf <- convolve$cdf
    } else if (!is.null(gm)) {
      cdf <- function(x) interactmeta_gm_cdf(gm, x)
    }
  }
  if (is.null(quantile)) {
    if (!is.null(convolve)) {
      quantile <- convolve$quantile
    } else if (!is.null(published) && identical(published$kind, "native")) {
      quantile <- function(p) {
        interactmeta_native_quantile(
          published$object,
          p,
          margin = published$margin %||% "effect",
          which = published$which,
          predict = isTRUE(published$predict),
          individual = published$individual %||% FALSE
        )
      }
    } else if (!is.null(gm)) {
      quantile <- function(p) interactmeta_gm_quantile(gm, p)
    }
  }
  structure(
    list(
      name = name,
      gm = gm,
      convolve = convolve,
      density = density,
      cdf = cdf,
      quantile = quantile,
      note = note,
      published = published,
      interval.type = if (is.null(interval.type)) {
        NULL
      } else {
        interactmeta_match_interval_type(interval.type)
      },
      tol.uniroot = tol.uniroot
    ),
    class = "interactmeta_posterior_1d"
  )
}

#' Published summary of a 1-D posterior.
#'
#' Bounds come from `$post.interval` or from the shortest/central rule on
#' `bayesmeta::convolve()`. The Gaussian mixture, when stored, supplies the
#' mean and sd only.
#' @noRd
interactmeta_posterior_1d_summary <- function(x, prob = 0.95, interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(
    interval.type %||% x$interval.type %||% "shortest"
  )
  tol <- x$tol.uniroot %||% (2^16 * .Machine$double.eps)
  use_published <- !is.null(x$published) || !is.null(x$convolve)
  if (!use_published) {
    if (!is.null(x$gm)) {
      return(interactmeta_gm_summary(
        x$gm,
        prob = prob,
        interval.type = interval.type
      ))
    }
    if (is.null(x$quantile)) {
      stop("1-D posterior has neither a Gaussian mixture nor a quantile function.", call. = FALSE)
    }
  }
  med <- NA_real_
  bounds <- c(NA_real_, NA_real_)
  if (!is.null(x$published) && identical(x$published$kind, "native")) {
    pub <- x$published
    med <- interactmeta_native_quantile(
      pub$object,
      0.5,
      margin = pub$margin %||% "effect",
      which = pub$which,
      predict = isTRUE(pub$predict),
      individual = pub$individual %||% FALSE
    )[[1L]]
    bounds <- interactmeta_post_interval(
      pub$object,
      prob = prob,
      interval.type = interval.type,
      margin = pub$margin %||% "effect",
      which = pub$which,
      predict = isTRUE(pub$predict),
      individual = pub$individual %||% FALSE
    )
  } else if (!is.null(x$convolve)) {
    med <- as.numeric(x$convolve$quantile(0.5))[[1L]]
    bounds <- interactmeta_interval_from_distribution(
      quantile = x$convolve$quantile,
      cdf = x$convolve$cdf,
      prob = prob,
      interval.type = interval.type,
      tol = tol
    )
  } else if (!is.null(x$quantile)) {
    qs <- x$quantile(c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2))
    med <- unname(qs[[2L]])
    bounds <- unname(qs[c(1L, 3L)])
  }
  mean_v <- NA_real_
  sd_v <- NA_real_
  if (!is.null(x$gm)) {
    mean_v <- interactmeta_gm_mean(x$gm)
    sd_v <- sqrt(max(interactmeta_gm_var(x$gm), 0))
  }
  list(
    mean = mean_v,
    sd = sd_v,
    median = med,
    lower = bounds[[1L]],
    upper = bounds[[2L]],
    q = stats::setNames(
      c(bounds[[1L]], med, bounds[[2L]]),
      c((1 - prob) / 2, 0.5, 1 - (1 - prob) / 2)
    )
  )
}

#' Weighted Pearson correlation of two series with the same weights
#' @noRd
interactmeta_weighted_cor <- function(x, y, w) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  w <- as.numeric(w)
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  cxy <- sum(w * (x - mx) * (y - my))
  sx <- sqrt(max(sum(w * (x - mx)^2), 0))
  sy <- sqrt(max(sum(w * (y - my)^2), 0))
  if (sx <= 0 || sy <= 0) {
    return(NA_real_)
  }
  cxy / (sx * sy)
}

#' Inverse CAMS rotation (trial-level)
#'
#' `mu_Aj = m_j - pi_j * g_j`
#' `mu_Bj = m_j + (1 - pi_j) * g_j`
#' These are study-specific latent subgroup effects at that trial's
#' `pi_j`, not population `mu_A(pi*)` / `mu_B(pi*)` from `subgroups()`.
#' @noRd
interactmeta_muA_from_mg <- function(m, g, pi) {
  as.numeric(m) - as.numeric(pi) * as.numeric(g)
}

interactmeta_muB_from_mg <- function(m, g, pi) {
  as.numeric(m) + (1 - as.numeric(pi)) * as.numeric(g)
}
