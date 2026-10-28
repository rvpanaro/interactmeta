# ---- cam_prior.R ----
#' CAMS prior specification
#'
#' Shared CAMS prior used by [cams()]. Character shorthands
#' such as `"uniform"` and `"sqrt"` are accepted and normalized internally.
#'
#' `"uniform"` with no finite `upper` is the **improper** density
#' \eqn{p(\tau)\propto 1} for \eqn{\tau \ge 0}, which is what
#' [bayesmeta::bayesmeta()] means by `tau.prior = "uniform"`. It is **not**
#' \eqn{U(0, U)} unless `upper` is set. `"sqrt"` is the improper
#' \eqn{p(\tau)\propto \tau^{-1/2}} density. These are bayesmeta
#' densities; frequentist CAMS (`engine = "metafor"`) ignores Bayesian priors.
#'
#' @param family Prior family. One of `"uniform"`, `"half-normal"`,
#'   `"half-cauchy"`, `"sqrt"`, `"normal"`.
#' @param scale Scale of a half-normal or half-Cauchy.
#' @param mean,sd Location and scale of a normal prior (fixed effects).
#' @param lower,upper Support for a uniform prior. `upper = Inf` (default)
#'   is the improper uniform on \eqn{[lower, \infty)}.
#' @return An object of class `"cam_prior"`.
#' @noRd
cam_prior <- function(family,
                       scale = NULL,
                       mean = 0,
                       sd = NULL,
                       lower = 0,
                       upper = Inf) {
  family <- cam_prior_canonical_family(family)
  out <- list(
    family = family,
    scale = if (is.null(scale)) NULL else as.numeric(scale)[[1L]],
    mean = as.numeric(mean)[[1L]],
    sd = if (is.null(sd)) NULL else as.numeric(sd)[[1L]],
    lower = as.numeric(lower)[[1L]],
    upper = as.numeric(upper)[[1L]]
  )
  class(out) <- "cam_prior"
  cam_prior_validate(out)
}

#' @export
#' @noRd
print.cam_prior <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
#' @noRd
format.cam_prior <- function(x, ...) {
  switch(
    x$family,
    "half-normal" = sprintf("Half-Normal(0, %s)", format(x$scale, digits = 6)),
    "half-cauchy" = sprintf("Half-Cauchy(0, %s)", format(x$scale, digits = 6)),
    "normal" = sprintf("Normal(%s, %s)", format(x$mean, digits = 6), format(x$sd, digits = 6)),
    "sqrt" = "sqrt (p(tau) proportional to tau^(-1/2), tau >= 0)",
    "uniform" = if (cam_prior_is_improper_uniform(x)) {
      sprintf(
        "improper uniform (p proportional to 1 on [%s, Inf); bayesmeta default)",
        format(x$lower, digits = 6)
      )
    } else {
      sprintf("Uniform(%s, %s)", format(x$lower, digits = 6), format(x$upper, digits = 6))
    },
    "function" = "user density function",
    x$family
  )
}

#' Documented CAMS Bayesian defaults
#'
#' Heterogeneity: bayesmeta's improper `"uniform"` on overall and
#' interaction tau when both are omitted. Fixed effects: independent
#' Normal(0, 10) — bayesmeta's improper flat mu is not a reliable default
#' for the overall `bmr` block together with improper tau. Empirical /
#' WIP half-normal scales are documented and used in vignettes
#' (arXiv:2606.23968).
#' @noRd
cam_default_tau_prior <- function() {
  cam_prior("uniform")
}

#' @noRd
cam_default_tau_interaction_prior <- function() {
  cam_prior("uniform")
}

#' @noRd
cam_default_fixed_prior <- function() {
  cam_prior("normal", mean = 0, sd = 10)
}

#' Normalize a user prior argument to a `cam_prior`
#' @noRd
normalize_cam_prior <- function(x, role = c("tau", "fixed")) {
  role <- match.arg(role)
  if (inherits(x, "cam_prior")) {
    return(cam_prior_validate(x, role = role))
  }
  if (is.null(x)) {
    return(NULL)
  }
  if (is.function(x)) {
    out <- list(family = "function", fun = x)
    class(out) <- "cam_prior"
    return(out)
  }
  if (is.numeric(x) && length(x) == 2L && !is.null(names(x))) {
    if (!all(c("mean", "sd") %in% names(x))) {
      stop(
        "A numeric prior needs names `mean` and `sd`, as in `c(mean = 0, sd = 10)`.",
        call. = FALSE
      )
    }
    return(cam_prior_validate(
      cam_prior("normal", mean = x[["mean"]], sd = x[["sd"]]),
      role = role
    ))
  }
  if (is.character(x) && length(x) == 1L) {
    fam <- cam_prior_canonical_family(x)
    if (identical(fam, "half-normal")) {
      stop(
        "Character shorthand `\"half-normal\"` needs a scale. Pass a bayesmeta density, ",
        "for example `function(tau) bayesmeta::dhalfnormal(tau, scale = 0.5)`.",
        call. = FALSE
      )
    }
    if (identical(fam, "half-cauchy")) {
      stop(
        "Character shorthand `\"half-cauchy\"` needs a scale. Pass a bayesmeta density, ",
        "for example `function(tau) bayesmeta::dhalfcauchy(tau, scale = 1)`.",
        call. = FALSE
      )
    }
    if (identical(fam, "normal")) {
      stop(
        "Character shorthand `\"normal\"` needs a mean and an sd. ",
        "Pass `fixed.prior = c(mean = 0, sd = 10)`.",
        call. = FALSE
      )
    }
    return(cam_prior(fam))
  }
  stop(
    "Cannot interpret prior of class `",
    paste(class(x), collapse = "/"),
    "`. Use `\"uniform\"`, `\"sqrt\"`, a bayesmeta density, or `c(mean = , sd = )`.",
    call. = FALSE
  )
}

#' Resolve CAMS-level Bayesian priors (defaults once; interaction may inherit)
#' @noRd
cam_priors_resolve <- function(tau.prior = NULL,
                                tau.interaction.prior = NULL,
                                fixed.prior = NULL,
                                user_tau = FALSE) {
  tau <- normalize_cam_prior(
    if (is.null(tau.prior)) cam_default_tau_prior() else tau.prior,
    role = "tau"
  )
  if (is.null(tau.interaction.prior)) {
    tau_g <- if (isTRUE(user_tau)) {
      tau
    } else {
      normalize_cam_prior(cam_default_tau_interaction_prior(), role = "tau")
    }
  } else {
    tau_g <- normalize_cam_prior(tau.interaction.prior, role = "tau")
  }
  fixed <- normalize_cam_prior(
    if (is.null(fixed.prior)) cam_default_fixed_prior() else fixed.prior,
    role = "fixed"
  )
  structure(
    list(
      tau = tau,
      tau_interaction = tau_g,
      fixed = fixed
    ),
    class = "cam_priors"
  )
}

#' Error if the selected engine cannot represent a CAMS prior equivalently
#' @noRd
cam_priors_check_engine <- function(priors, engine) {
  invisible(priors)
}

#' @export
#' @noRd
print.cam_priors <- function(x, ...) {
  cat("CAMS priors:\n")
  cat("  tau (overall)     : ", format(x$tau), "\n", sep = "")
  cat("  tau (interaction) : ", format(x$tau_interaction), "\n", sep = "")
  cat("  population-level means: ", format(x$fixed), "\n", sep = "")
  invisible(x)
}

cam_prior_canonical_family <- function(family) {
  fam <- tolower(as.character(family)[[1L]])
  fam <- gsub("[-_ ]+", "", fam)
  aliases <- c(
    uniform = "uniform",
    flat = "uniform",
    halfnormal = "half-normal",
    hn = "half-normal",
    halfn = "half-normal",
    halfcauchy = "half-cauchy",
    hc = "half-cauchy",
    sqrt = "sqrt",
    normal = "normal",
    gaussian = "normal"
  )
  if (!fam %in% names(aliases)) {
    stop(
      "Unknown CAMS prior family `",
      family,
      "`. Supported: uniform, half-normal, half-cauchy, sqrt, normal.",
      call. = FALSE
    )
  }
  unname(aliases[[fam]])
}

cam_prior_is_improper_uniform <- function(x) {
  identical(x$family, "uniform") && (!is.finite(x$upper) || is.infinite(x$upper))
}

cam_prior_is_improper <- function(x) {
  identical(x$family, "sqrt") ||
    cam_prior_is_improper_uniform(x) ||
    identical(x$family, "function")
}

cam_prior_validate <- function(x, role = c("tau", "fixed")) {
  role <- match.arg(role)
  stopifnot(inherits(x, "cam_prior"))
  if (identical(x$family, "half-normal") || identical(x$family, "half-cauchy")) {
    if (is.null(x$scale) || !is.finite(x$scale) || x$scale <= 0) {
      stop(x$family, " priors require a finite positive `scale`.", call. = FALSE)
    }
  }
  if (identical(x$family, "normal")) {
    if (is.null(x$sd) || !is.finite(x$sd) || x$sd <= 0) {
      stop("Normal priors require a finite positive `sd`.", call. = FALSE)
    }
    if (!is.finite(x$mean)) {
      stop("Normal priors require a finite `mean`.", call. = FALSE)
    }
  }
  if (identical(x$family, "uniform")) {
    if (!is.finite(x$lower)) {
      stop("Uniform priors require a finite `lower`.", call. = FALSE)
    }
    if (is.finite(x$upper) && x$upper <= x$lower) {
      stop("Uniform `upper` must exceed `lower`.", call. = FALSE)
    }
  }
  if (identical(role, "fixed") && x$family %in% c("sqrt", "half-normal", "half-cauchy")) {
    stop("Fixed-effect priors cannot be `", x$family, "`.", call. = FALSE)
  }
  x
}

#' Evaluate a CAMS tau prior density on a grid (for plots and checks)
#' @noRd
cam_prior_tau_density <- function(prior, tau) {
  prior <- normalize_cam_prior(prior, role = "tau")
  tau <- as.numeric(tau)
  if (!length(tau)) {
    return(numeric())
  }
  if (identical(prior$family, "function")) {
    return(as.numeric(prior$fun(tau)))
  }
  bm <- prior_to_bayesmeta(prior)
  if (is.function(bm)) {
    return(as.numeric(bm(tau)))
  }
  if (identical(bm, "uniform")) {
    out <- rep(1, length(tau))
    out[tau < 0] <- 0
    return(out)
  }
  if (identical(bm, "sqrt")) {
    check_bayesmeta()
    return(as.numeric(bayesmeta::dprior(tau = tau, tau.prior = "sqrt")))
  }
  stop("No density translation for prior family `", prior$family, "`.", call. = FALSE)
}

#' Translate a CAMS tau prior into a bayesmeta / bmr `tau.prior` argument
#' @noRd
prior_to_bayesmeta <- function(prior) {
  prior <- normalize_cam_prior(prior, role = "tau")
  if (identical(prior$family, "function")) {
    return(prior$fun)
  }
  if (identical(prior$family, "sqrt")) {
    return("sqrt")
  }
  if (cam_prior_is_improper_uniform(prior)) {
    if (isTRUE(prior$lower == 0)) {
      return("uniform")
    }
    lo <- prior$lower
    return(function(tau) {
      ifelse(tau < lo, 0, 1)
    })
  }
  if (identical(prior$family, "uniform")) {
    lo <- prior$lower
    hi <- prior$upper
    return(function(tau) {
      stats::dunif(tau, min = lo, max = hi)
    })
  }
  if (identical(prior$family, "half-normal")) {
    sc <- prior$scale
    return(function(tau) {
      bayesmeta::dhalfnormal(tau, scale = sc)
    })
  }
  if (identical(prior$family, "half-cauchy")) {
    sc <- prior$scale
    return(function(tau) {
      bayesmeta::dhalfcauchy(tau, scale = sc)
    })
  }
  if (identical(prior$family, "normal")) {
    stop("A normal prior is not a valid heterogeneity prior; use half-normal.", call. = FALSE)
  }
  stop("No bayesmeta translation for prior family `", prior$family, "`.", call. = FALSE)
}

#' Translate a CAMS fixed-effect prior for bayesmeta / bmr
#'
#' The numbers are not broadcast as one normal on the CAMS coefficients
#' \eqn{(\alpha, \beta, \gamma)}. `bmr` is told the prior on \eqn{(\alpha, \theta)}
#' explicitly (`beta.prior.mean`, `beta.prior.sd`, and `beta.prior.cov`).
#' The covariance kept here is `diag(sd^2)`: independent normals on
#' \eqn{(\alpha, \theta)}. `bayesmeta` is told `mu.prior` for \eqn{\gamma}
#' separately.
#'
#' Independent normals on \eqn{(\alpha, \theta)} are not independent normals on
#' \eqn{(\alpha, \beta, \gamma)}, because \eqn{\theta = \beta + \gamma}.
#' @noRd
prior_to_bayesmeta_fixed <- function(prior) {
  prior <- normalize_cam_prior(prior, role = "fixed")
  note <- paste(
    "bmr receives the prior on (alpha, theta), including beta.prior.cov;",
    "bayesmeta receives mu.prior for gamma separately.",
    "Independent normals on (alpha, theta) are not independent normals on",
    "(alpha, beta, gamma), because theta = beta + gamma."
  )
  if (identical(prior$family, "uniform") || cam_prior_is_improper_uniform(prior)) {
    return(list(
      mu.prior = c(mean = NA_real_, sd = NA_real_),
      beta.prior.mean = NULL,
      beta.prior.sd = NULL,
      beta.prior.cov = NULL,
      note = note
    ))
  }
  if (!identical(prior$family, "normal")) {
    stop(
      "bayesmeta fixed-effect priors currently support `normal` or improper uniform.",
      call. = FALSE
    )
  }
  sd <- prior$sd
  mean <- prior$mean
  cov <- diag(c(sd, sd)^2, nrow = 2L, ncol = 2L)
  dimnames(cov) <- list(c("alpha", "theta"), c("alpha", "theta"))
  list(
    mu.prior = c(mean = mean, sd = sd),
    beta.prior.mean = c(alpha = mean, theta = mean),
    beta.prior.sd = c(alpha = sd, theta = sd),
    beta.prior.cov = cov,
    note = note
  )
}
