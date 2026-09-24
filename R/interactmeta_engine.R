#' Engine / method / engine_args helpers
#' @noRd
NULL

#' Public CAMS backends: Bayesian DIRECT (bayesmeta) and frequentist (metafor).
#' @noRd
interactmeta_public_engines <- function() {
  c("bayesmeta", "metafor")
}

#' Engines advertised to users and listed in [interactmeta_capability_matrix()].
#' @noRd
interactmeta_engines <- function() {
  interactmeta_public_engines()
}

#' Undocumented CAMS backend kept for internal equivalence tests.
#' Not part of the public engine API.
#' @noRd
interactmeta_internal_engines <- function() {
  "brms"
}

interactmeta_removed_engine_error <- function(engine) {
  engine <- as.character(engine)[[1L]]
  if (tolower(engine) %in% c("cmdstanr", "stan", "rstan", "rstanarm")) {
    stop(
      "`engine = \"cmdstanr\"` has been removed. ",
      "Bayesian CAMS uses `engine = \"bayesmeta\"`; ",
      "frequentist CAMS uses `engine = \"metafor\"`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

interactmeta_match_engine <- function(engine) {
  engine <- as.character(engine)[[1L]]
  if (identical(tolower(engine), "bayes")) {
    engine <- "bayesmeta"
  }
  interactmeta_removed_engine_error(engine)
  public <- interactmeta_public_engines()
  allowed <- c(public, interactmeta_internal_engines())
  if (!engine %in% allowed) {
    stop(
      "`engine` must be one of ",
      paste(paste0('"', public, '"'), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  engine
}

#' Whether the user left `engine` at its default vector / unspecified
#' @noRd
interactmeta_engine_unspecified <- function(engine) {
  is.null(engine) || (is.character(engine) && length(engine) > 1L)
}

#' Normalize `approach = "bayes"` / `"frequentist"` / `"REML"` / `"ML"`
#' @noRd
interactmeta_normalize_approach <- function(approach) {
  if (is.null(approach) || (is.character(approach) && length(approach) > 1L)) {
    return(NA_character_)
  }
  a <- tolower(as.character(approach)[[1L]])
  if (a %in% c("bayes", "bayesian")) {
    return("bayes")
  }
  if (a %in% c("frequentist", "freq")) {
    return("frequentist")
  }
  if (identical(a, "reml")) {
    return("REML")
  }
  if (identical(a, "ml")) {
    return("ML")
  }
  stop(
    "`approach` must be \"bayes\", \"frequentist\", \"REML\", or \"ML\".",
    call. = FALSE
  )
}

#' Public engine that implements a method
#' @noRd
interactmeta_engine_for_method <- function(method, allowed = interactmeta_public_engines()) {
  method <- as.character(method)[[1L]]
  engine <- switch(
    method,
    bayes = "bayesmeta",
    REML = "metafor",
    ML = "metafor",
    NA_character_
  )
  if (is.na(engine) || !engine %in% allowed) {
    return(NA_character_)
  }
  engine
}

#' Pick engine from inferential approach / method when `engine` is omitted
#' @noRd
interactmeta_engine_from_approach <- function(engine,
                                             method,
                                             approach = NULL,
                                             allowed = interactmeta_public_engines()) {
  approach <- interactmeta_normalize_approach(approach)
  specified_engine <- !interactmeta_engine_unspecified(engine)
  method_chr <- if (is.null(method) || (is.character(method) && length(method) > 1L)) {
    NA_character_
  } else {
    as.character(method)[[1L]]
  }
  if (!specified_engine) {
    from_method <- if (!is.na(method_chr)) {
      interactmeta_engine_for_method(method_chr, allowed)
    } else {
      NA_character_
    }
    engine <- if (!is.na(approach) && identical(approach, "bayes") &&
      "bayesmeta" %in% allowed) {
      "bayesmeta"
    } else if (!is.na(approach) &&
      approach %in% c("frequentist", "REML", "ML") &&
      "metafor" %in% allowed) {
      "metafor"
    } else if (!is.na(from_method)) {
      from_method
    } else {
      allowed[[1L]]
    }
  } else {
    engine <- as.character(engine)[[1L]]
  }
  engine <- interactmeta_match_engine(engine)
  if (!is.na(approach)) {
    if (identical(approach, "bayes") && !identical(engine, "bayesmeta") &&
      !engine %in% interactmeta_internal_engines()) {
      stop(
        "`approach = \"bayes\"` uses engine = \"bayesmeta\" when `engine` is ",
        "omitted. `engine = \"",
        engine,
        "\"` is not a Bayesian DIRECT backend.",
        call. = FALSE
      )
    }
    if (approach %in% c("frequentist", "REML", "ML") &&
      identical(engine, "bayesmeta")) {
      stop(
        "`approach = \"",
        approach,
        "\"` uses engine = \"metafor\" when `engine` is omitted. ",
        "It cannot be combined with `engine = \"bayesmeta\"`.",
        call. = FALSE
      )
    }
  }
  engine
}

#' Default method once an engine (and optional approach) is known
#' @noRd
interactmeta_default_method <- function(engine, approach = NULL) {
  approach <- interactmeta_normalize_approach(approach)
  if (identical(approach, "REML") || identical(approach, "ML")) {
    return(approach)
  }
  if (identical(engine, "metafor")) {
    return("REML")
  }
  "bayes"
}

interactmeta_resolve_method <- function(engine, method) {
  allowed <- c("bayes", "REML", "ML")
  if (is.null(method) || (is.character(method) && length(method) > 1L)) {
    method <- interactmeta_default_method(engine)
  } else {
    method <- as.character(method)[[1L]]
    if (identical(method, "mle")) {
      stop(
        "`method = \"mle\"` has been removed. ",
        "Use `engine = \"metafor\"` for frequentist CAMS (`method = \"REML\"` or `\"ML\"`).",
        call. = FALSE
      )
    }
    if (!method %in% allowed) {
      stop(
        "`method` must be one of ",
        paste(paste0('"', allowed, '"'), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  if (identical(engine, "bayesmeta") && !identical(method, "bayes")) {
    stop(
      "`engine = \"",
      engine,
      "\"` requires `method = \"bayes\"`.",
      call. = FALSE
    )
  }
  if (identical(engine, "metafor") && !method %in% c("REML", "ML")) {
    stop(
      "`engine = \"metafor\"` requires `method = \"REML\"` or `\"ML\"`.",
      call. = FALSE
    )
  }
  method
}

interactmeta_algorithm_for <- function(engine, method) {
  if (identical(engine, "bayesmeta")) {
    return("quadrature")
  }
  if (identical(engine, "metafor")) {
    return("optimization")
  }
  if (identical(method, "mle")) {
    return("optimization")
  }
  "MCMC"
}

#' Merge `...` computational knobs into `engine_args` (engine_args wins)
#' @noRd
interactmeta_merge_engine_args <- function(engine_args, dots) {
  if (is.null(engine_args)) {
    engine_args <- list()
  }
  if (!is.list(engine_args)) {
    stop("`engine_args` must be a list.", call. = FALSE)
  }
  extra <- intersect(
    names(dots),
    c(
      "chains", "iter", "warmup", "thin", "refresh", "cores", "control",
      "iter_warmup", "iter_sampling", "parallel_chains", "backend",
      "delta", "epsilon", "interval.type", "convolve_derived", "test", "verbose",
      "optimizer", "opt_method", "init", "adapt_delta"
    )
  )
  for (nm in extra) {
    if (!nm %in% names(engine_args)) {
      engine_args[[nm]] <- dots[[nm]]
    }
  }
  engine_args
}

interactmeta_computational_dot_names <- function() {
  c(
    "chains", "iter", "warmup", "thin", "refresh", "cores", "control",
    "iter_warmup", "iter_sampling", "parallel_chains", "backend",
    "delta", "epsilon", "interval.type", "convolve_derived", "test", "verbose",
    "optimizer", "opt_method", "init", "adapt_delta", "n_draws"
  )
}

.interactmeta_warn_unused_priors <- function(user_supplied) {
  if (isTRUE(user_supplied)) {
    warning(
      "`engine = \"metafor\"` is frequentist; Bayesian prior arguments are ignored.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
