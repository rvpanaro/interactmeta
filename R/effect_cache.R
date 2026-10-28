.interactmeta_fit_sidecar <- new.env(parent = emptyenv())

#' Sidecar store keyed by `attr(fit, "interactmeta_sidecar")` (reference semantics).
#' @noRd
interactmeta_fit_sidecar_key <- function(fit) {
  key <- attr(fit, "interactmeta_sidecar", exact = TRUE)
  if (!is.null(key) && !nzchar(key)) {
    return(NULL)
  }
  key
}

#' @noRd
interactmeta_fit_sidecar <- function(fit, create = FALSE) {
  key <- interactmeta_fit_sidecar_key(fit)
  if (is.null(key)) {
    if (!isTRUE(create)) {
      return(NULL)
    }
    return(NULL)
  }
  if (!exists(key, envir = .interactmeta_fit_sidecar, inherits = FALSE)) {
    return(NULL)
  }
  get(key, envir = .interactmeta_fit_sidecar, inherits = FALSE)
}

#' Attach sidecar key at end of [cams()] (same evaluation frame as the fit).
#' @noRd
interactmeta_init_fit_sidecar <- function(fit) {
  if (!inherits(fit, "interactmeta")) {
    return(fit)
  }
  key <- interactmeta_fit_sidecar_key(fit)
  if (is.null(key)) {
    key <- paste0(
      "s",
      paste(format(Sys.time(), "%Y%m%d%H%M%OS3"), sample.int(1e6, 1L), sep = "_")
    )
    attr(fit, "interactmeta_sidecar") <- key
    assign(key, new.env(parent = emptyenv()), envir = .interactmeta_fit_sidecar)
  }
  fit
}

#' In-session cache for marginal effect summaries (shortest-interval work)
#'
#' Repeated [modifiers()], [combined()], [interactions()], [forestplot()],
#' and [summary()] share the same DIRECT mixture intervals. Results are stored
#' on `fit$effect_cache` (an environment) keyed by request arguments.
#' @noRd
interactmeta_effect_cache_env <- function(fit) {
  if (!is.list(fit)) {
    return(NULL)
  }
  sc <- interactmeta_fit_sidecar(fit, create = TRUE)
  if (!is.null(sc)) {
    cache <- sc$effect_cache
    if (is.null(cache) || !is.environment(cache)) {
      cache <- new.env(parent = emptyenv())
      sc$effect_cache <- cache
    }
    return(cache)
  }
  cache <- fit$effect_cache
  if (is.null(cache) || !is.environment(cache)) {
    cache <- new.env(parent = emptyenv())
    fit$effect_cache <- cache
  }
  cache
}

#' @noRd
interactmeta_sync_reporting_from_sidecar <- function(fit) {
  sc <- interactmeta_fit_sidecar(fit, create = FALSE)
  if (is.null(sc) || !is.list(sc$reporting)) {
    return(invisible(fit))
  }
  fit$reporting <- sc$reporting
  if (is.environment(sc$effect_cache)) {
    fit$effect_cache <- sc$effect_cache
  }
  invisible(fit)
}

#' @noRd
interactmeta_reporting_sidecar <- function(fit) {
  sc <- interactmeta_fit_sidecar(fit, create = TRUE)
  if (is.null(sc)) {
    return(fit$reporting)
  }
  sc$reporting
}

#' @noRd
interactmeta_set_fit_reporting <- function(fit, reporting) {
  sc <- interactmeta_fit_sidecar(fit, create = TRUE)
  if (is.null(sc)) {
    fit$reporting <- reporting
    return(invisible(fit))
  }
  sc$reporting <- reporting
  fit$reporting <- reporting
  invisible(fit)
}

#' @noRd
interactmeta_effect_cache_at_token <- function(at, ifrac) {
  ifrac <- match.arg(ifrac, c("interaction", "overall"))
  if (is.null(at)) {
    return(paste0("infer:", ifrac))
  }
  if (inherits(at, "interactmeta_overall_ifrac")) {
    summ <- attr(at, "summary", exact = TRUE)
    if (is.data.frame(summ) && nrow(summ)) {
      est <- as.numeric(summ$estimate[[1L]])
      if (is.finite(est)) {
        return(sprintf("oi:%.10g", est))
      }
    }
  }
  if (is.numeric(at) && length(at) == 1L) {
    return(sprintf("plug:%.10g", as.numeric(at)[[1L]]))
  }
  if (is.numeric(at) && length(at) > 1L) {
    return(paste0("infer:", ifrac))
  }
  "custom"
}

#' @noRd
interactmeta_effect_cache_key <- function(quantity,
                                          at = NULL,
                                          ifrac = c("interaction", "overall"),
                                          type = NULL,
                                          prob = 0.95,
                                          max_support = NULL,
                                          contrast_sep = "-",
                                          exponentiate = FALSE,
                                          summary = TRUE) {
  ifrac <- match.arg(ifrac)
  ty <- if (is.null(type)) {
    ""
  } else {
    paste(sort(as.character(type)), collapse = "+")
  }
  md <- if (is.null(max_support)) {
    "all"
  } else {
    as.character(as.integer(max_support)[[1L]])
  }
  paste(
    quantity,
    interactmeta_effect_cache_at_token(at, ifrac),
    paste0("type=", ty),
    paste0("prob=", format(as.numeric(prob)[[1L]], digits = 17)),
    paste0("md=", md),
    paste0("sep=", contrast_sep),
    paste0("exp=", as.integer(isTRUE(exponentiate))),
    paste0("sum=", as.integer(isTRUE(summary))),
    sep = "|"
  )
}

#' @noRd
interactmeta_effect_cache_get <- function(fit, key) {
  cache <- interactmeta_effect_cache_env(fit)
  if (is.null(cache) || !nzchar(key)) {
    return(NULL)
  }
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  NULL
}

#' @noRd
interactmeta_effect_cache_set <- function(fit, key, value) {
  cache <- interactmeta_effect_cache_env(fit)
  if (is.null(cache) || !nzchar(key) || is.null(value)) {
    return(invisible(NULL))
  }
  assign(key, value, envir = cache)
  invisible(NULL)
}

#' Coerce a coverage probability (handles partial `fit$reporting` lists).
#' @noRd
interactmeta_default_prob <- function(prob, default = 0.95) {
  fallback <- as.numeric(default %||% 0.95)
  if (length(fallback) < 1L || !is.finite(fallback[[1L]])) {
    fallback <- 0.95
  } else {
    fallback <- fallback[[1L]]
  }
  if (is.null(prob) || !length(prob)) {
    return(fallback)
  }
  p <- suppressWarnings(as.numeric(prob)[[1L]])
  if (!length(p) || !is.finite(p) || p <= 0 || p >= 1) {
    return(fallback)
  }
  p
}

#' Whether default reporting (parameters + effect cache) is ready
#' @noRd
interactmeta_reporting_is_ready <- function(fit, prob = 0.95) {
  if (!inherits(fit, "interactmeta")) {
    return(FALSE)
  }
  rep <- interactmeta_reporting_sidecar(fit) %||% fit$reporting
  if (!is.list(rep)) {
    return(FALSE)
  }
  prob <- interactmeta_default_prob(prob, default = rep$prob %||% 0.95)
  cache <- interactmeta_effect_cache_env(fit)
  is.matrix(rep$parameters) &&
    isTRUE(all.equal(as.numeric(rep$prob), prob)) &&
    is.environment(cache)
}

#' Ensure reporting bundle (lazy after [cams()] with `precompute = FALSE`).
#' @noRd
interactmeta_ensure_reporting <- function(fit,
                                          prob = 0.95,
                                          warm_marginals = TRUE) {
  if (!inherits(fit, "interactmeta")) {
    return(invisible(fit))
  }
  if (interactmeta_reporting_is_ready(fit, prob = prob)) {
    return(invisible(fit))
  }
  fit <- interactmeta_init_fit_sidecar(fit)
  fit <- interactmeta_precompute_reporting(
    fit,
    prob = prob,
    warm_marginals = warm_marginals
  )
  invisible(fit)
}

#' Default reporting bundle at end of [cams()] (95% intervals, inferred IF).
#' @noRd
interactmeta_precompute_reporting <- function(fit,
                                              prob = 0.95,
                                              warm_marginals = TRUE) {
  if (!inherits(fit, "interactmeta")) {
    return(invisible(fit))
  }
  rep_existing <- interactmeta_reporting_sidecar(fit) %||% fit$reporting
  prob <- interactmeta_default_prob(
    prob,
    default = if (is.list(rep_existing)) rep_existing$prob %||% 0.95 else 0.95
  )
  st <- attr(fit, "interactmeta", exact = TRUE)
  fit <- interactmeta_init_fit_sidecar(fit)
  reporting <- list(
    prob = prob,
    ifrac = "interaction",
    contrast_sep = "-",
    exponentiate = FALSE,
    interval.type = interactmeta_match_interval_type(fit$interval.type),
    reference = if (is.list(st)) st$reference else NULL,
    method = fit$method,
    inference = if (interactmeta_is_bayesian_fit(fit)) {
      "Bayesian"
    } else {
      "frequentist"
    }
  )
  interactmeta_set_fit_reporting(fit, reporting)
  interactmeta_effect_cache_env(fit)
  mat <- tryCatch(
    if (interactmeta_is_bayesian_fit(fit)) {
      interactmeta_bayesmeta_param_summary(fit, prob = prob)
    } else {
      interactmeta_metafor_param_summary(fit, prob = prob)
    },
    error = function(e) NULL
  )
  if (is.matrix(mat)) {
    rep <- interactmeta_reporting_sidecar(fit) %||% fit$reporting
    rep$parameters <- mat
    interactmeta_set_fit_reporting(fit, rep)
  }
  interactmeta_sync_reporting_from_sidecar(fit)
  if (isTRUE(warm_marginals)) {
    interactmeta_warm_effect_cache(fit, prob = prob)
  }
  fit
}

#' Precompute default marginal summaries at the reporting IF.
#'
#' [combined()] and [interactions()] use engine-native bayesmeta intervals when
#' eligible; [modifiers()] still uses the DIRECT mixture grid.
#' @noRd
interactmeta_warm_effect_cache <- function(fit, prob = 0.95) {
  if (!inherits(fit, "interactmeta")) {
    return(invisible(fit))
  }
  prob <- as.numeric(prob)[[1L]]
  tryCatch(modifiers(fit, prob = prob, summary = TRUE), error = function(e) NULL)
  tryCatch(
    combined(fit, prob = prob, summary = TRUE, type = "Overall"),
    error = function(e) NULL
  )
  tryCatch(
    interactions(fit, prob = prob, summary = TRUE, type = "Overall"),
    error = function(e) NULL
  )
  invisible(fit)
}
