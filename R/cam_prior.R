#' CAMS prior specification
#'
#' Shared CAMS prior used by [cam()]. Character shorthands
#' such as `"uniform"` and `"sqrt"` are accepted and normalized internally.
#'
#' `"uniform"` with no finite `upper` is the **improper** density
#' \eqn{p(\tau)\propto 1} for \eqn{\tau \ge 0}, which is what
#' [bayesmeta::bayesmeta()] means by `tau.prior = "uniform"`. It is **not**
#' \eqn{U(0, U)} unless `upper` is set. `"sqrt"` is the improper
#' \eqn{p(\tau)\propto \tau^{-1/2}} density. Neither improper form can be
#' represented equivalently as a proper Stan density. They are bayesmeta
#' densities; frequentist CAMS (`engine = "metafor"`) ignores Bayesian priors.
#'
#' @param family Prior family. One of `"uniform"`, `"half-normal"`,
#'   `"half-cauchy"`, `"sqrt"`, `"normal"`.
#' @param scale Scale of a half-normal or half-Cauchy (Stan parameterization).
#' @param mean,sd Location and scale of a normal prior (fixed effects).
#' @param lower,upper Support for a uniform prior. `upper = Inf` (default)
#'   is the improper uniform on \eqn{[lower, \infty)}.
#' @return An object of class `"cam_prior"`.
#' @export
#' @examples
#' cam_prior("half-normal", scale = 0.5)
#' cam_prior("uniform", lower = 0, upper = 2)
#' @aliases cams_prior
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
  class(out) <- c("cam_prior", "cams_prior")
  cam_prior_validate(out)
}

#' @rdname cam_prior
#' @export
cams_prior <- cam_prior

#' @export
print.cam_prior <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.cam_prior <- function(x, ...) {
  switch(
    x$family,
    "half-normal" = sprintf("Half-Normal(0, %s)", format(x$scale, digits = 6)),
    "half-cauchy" = sprintf("Half-Cauchy(0, %s)", format(x$scale, digits = 6)),
    "normal" = sprintf("Normal(%s, %s)", format(x$mean, digits = 6), format(x$sd, digits = 6)),
    "sqrt" = "sqrt (p(tau) proportional to tau^(-1/2), tau >= 0)",
    "uniform" = if (cam_prior_is_improper_uniform(x)) {
      sprintf("improper uniform (p(tau) proportional to 1, tau >= %s)", format(x$lower, digits = 6))
    } else {
      sprintf("Uniform(%s, %s)", format(x$lower, digits = 6), format(x$upper, digits = 6))
    },
    "function" = "user density function",
    x$family
  )
}

#' @export
print.cams_prior <- print.cam_prior

#' @export
format.cams_prior <- format.cam_prior

#' Documented CAMS Bayesian defaults
#' @noRd
cam_default_tau_prior <- function() {
  cam_prior("half-normal", scale = 0.5)
}

#' @noRd
cam_default_tau_interaction_prior <- function() {
  cam_prior("half-normal", scale = 0.1)
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
  if (is.character(x) && length(x) == 1L) {
    fam <- cam_prior_canonical_family(x)
    if (identical(fam, "half-normal")) {
      stop(
        "Character shorthand `\"half-normal\"` needs a scale; use ",
        "`cam_prior(\"half-normal\", scale = ...)`.",
        call. = FALSE
      )
    }
    if (identical(fam, "half-cauchy")) {
      stop(
        "Character shorthand `\"half-cauchy\"` needs a scale; use ",
        "`cam_prior(\"half-cauchy\", scale = ...)`.",
        call. = FALSE
      )
    }
    if (identical(fam, "normal")) {
      stop(
        "Character shorthand `\"normal\"` needs `mean` and `sd`; use ",
        "`cam_prior(\"normal\", mean = ..., sd = ...)`.",
        call. = FALSE
      )
    }
    return(cam_prior(fam))
  }
  stop(
    "Cannot interpret prior of class `",
    paste(class(x), collapse = "/"),
    "`. Use a character shorthand or `cam_prior()`.",
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
    class = c("cam_priors", "cams_priors")
  )
}

#' Error if the selected engine cannot represent a CAMS prior equivalently
#' @noRd
cam_priors_check_engine <- function(priors, engine) {
  invisible(priors)
}

#' @export
print.cam_priors <- function(x, ...) {
  cat("CAMS priors:\n")
  cat("  tau (overall)     : ", format(x$tau), "\n", sep = "")
  cat("  tau (interaction) : ", format(x$tau_interaction), "\n", sep = "")
  cat("  fixed effects     : ", format(x$fixed), "\n", sep = "")
  invisible(x)
}

#' @export
print.cams_priors <- print.cam_priors

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

#' Engine-independent brms prior rows (inspectable without brms)
#' @noRd
prior_to_brms_rows <- function(priors, trial_var = "trial", ifrac2_vars = "ifrac2_1") {
  if (!inherits(priors, "cam_priors")) {
    stop("`priors` must be a `cam_priors` object.", call. = FALSE)
  }
  fe <- prior_to_brms_one(priors$fixed, role = "fixed")
  tau <- prior_to_brms_one(priors$tau, role = "tau")
  tau_g <- prior_to_brms_one(priors$tau_interaction, role = "tau")
  rows <- list(
    data.frame(
      prior = fe$prior,
      class = "Intercept",
      coef = "",
      group = "",
      lb = NA_real_,
      ub = NA_real_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      prior = fe$prior,
      class = "b",
      coef = "",
      group = "",
      lb = NA_real_,
      ub = NA_real_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      prior = tau$prior,
      class = "sd",
      coef = "Intercept",
      group = trial_var,
      lb = tau$lb,
      ub = tau$ub,
      stringsAsFactors = FALSE
    )
  )
  for (coef in as.character(ifrac2_vars)) {
    rows[[length(rows) + 1L]] <- data.frame(
      prior = tau_g$prior,
      class = "sd",
      coef = coef,
      group = trial_var,
      lb = tau_g$lb,
      ub = tau_g$ub,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

prior_to_brms_one <- function(prior, role = c("tau", "fixed")) {
  role <- match.arg(role)
  prior <- normalize_cam_prior(prior, role = role)
  if (identical(prior$family, "sqrt")) {
    stop(
      "The `sqrt` prior p(tau) proportional to tau^(-1/2) is improper and ",
      "cannot be represented equivalently in brms/Stan. ",
      "Use engine = \"bayesmeta\" or choose a proper prior.",
      call. = FALSE
    )
  }
  if (cam_prior_is_improper_uniform(prior)) {
    stop(
      "`tau.prior = \"uniform\"` is the improper density p(tau) proportional to 1 ",
      "for tau >= 0 (bayesmeta). brms/Stan cannot use that prior. ",
      "Supply an explicit bound, e.g. cam_prior(\"uniform\", lower = 0, upper = 2), ",
      "or use a half-normal.",
      call. = FALSE
    )
  }
  if (identical(prior$family, "function")) {
    stop(
      "Function priors are supported for bayesmeta only; brms needs a named family.",
      call. = FALSE
    )
  }
  if (identical(prior$family, "half-normal")) {
    return(list(prior = sprintf("normal(0, %s)", prior$scale), lb = 0, ub = NA_real_))
  }
  if (identical(prior$family, "half-cauchy")) {
    return(list(prior = sprintf("cauchy(0, %s)", prior$scale), lb = 0, ub = NA_real_))
  }
  if (identical(prior$family, "normal")) {
    return(list(
      prior = sprintf("normal(%s, %s)", prior$mean, prior$sd),
      lb = NA_real_,
      ub = NA_real_
    ))
  }
  if (identical(prior$family, "uniform")) {
    return(list(
      prior = sprintf("uniform(%s, %s)", prior$lower, prior$upper),
      lb = prior$lower,
      ub = prior$upper
    ))
  }
  stop("No brms translation for prior family `", prior$family, "`.", call. = FALSE)
}

#' Translate CAMS priors into a `brms::prior` specification
#' @noRd
prior_to_brms <- function(priors, trial_var = "trial", ifrac2_vars = "ifrac2_1") {
  check_brms()
  rows <- prior_to_brms_rows(priors, trial_var = trial_var, ifrac2_vars = ifrac2_vars)
  priors_list <- vector("list", nrow(rows))
  for (i in seq_len(nrow(rows))) {
    args <- list(
      prior = rows$prior[[i]],
      class = rows$class[[i]]
    )
    if (nzchar(rows$coef[[i]])) {
      args$coef <- rows$coef[[i]]
    }
    if (nzchar(rows$group[[i]])) {
      args$group <- rows$group[[i]]
    }
    # brms forbids `coef` together with `lb`/`ub`. SD parameters already have
    # lb = 0; a uniform(L, U) prior density is zero outside [L, U], so the
    # bound is in the density rather than a silent change of support.
    if (!nzchar(rows$coef[[i]])) {
      if (is.finite(rows$lb[[i]])) {
        args$lb <- rows$lb[[i]]
      }
      if (is.finite(rows$ub[[i]])) {
        args$ub <- rows$ub[[i]]
      }
    }
    priors_list[[i]] <- do.call(brms::set_prior, args)
  }
  do.call(c, priors_list)
}

check_brms <- function() {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package `brms` is required for `engine = \"brms\"`. ",
      "Install it with install.packages(\"brms\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
