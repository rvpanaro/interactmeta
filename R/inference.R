# ---- interactmeta_engine.R ----
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

interactmeta_match_engine <- function(engine) {
  engine <- as.character(engine)[[1L]]
  if (identical(tolower(engine), "bayes")) {
    engine <- "bayesmeta"
  }
  public <- interactmeta_public_engines()
  if (!engine %in% public) {
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

#' `metafor::rma()` method values accepted by [cams()] (see `?metafor::rma`)
#' @noRd
interactmeta_metafor_rma_methods <- function() {
  c(
    "EE", "DL", "HE", "CO", "HS", "HSk", "SJ", "ML", "REML", "EB", "PM",
    "GENQ", "PMM", "GENQM"
  )
}

#' Whether a `cams()` method selects bayesmeta
#' @noRd
interactmeta_cams_method_is_bayesian <- function(method) {
  tolower(as.character(method)[[1L]]) %in% c("bayes", "bayesian")
}

#' Resolve `cams(..., method = )` to backend + stored method
#' @noRd
interactmeta_resolve_cams_method <- function(method) {
  if (is.null(method) || (is.character(method) && length(method) > 1L)) {
    method <- "Bayesian"
  }
  method <- as.character(method)[[1L]]
  if (identical(method, "mle")) {
    stop(
      "`method = \"mle\"` has been removed. ",
      "Use `method = \"ML\"` or `\"REML\"` for frequentist CAMS.",
      call. = FALSE
    )
  }
  if (interactmeta_cams_method_is_bayesian(method)) {
    return(list(
      backend = "bayesmeta",
      method_key = "bayes",
      method = "Bayesian"
    ))
  }
  if (identical(tolower(method), "frequentist")) {
    return(list(
      backend = "metafor",
      method_key = "REML",
      method = "REML"
    ))
  }
  rma <- interactmeta_metafor_rma_methods()
  if (!method %in% rma) {
    stop(
      "`method` must be \"Bayesian\" or a `metafor::rma()` method: ",
      paste(paste0('"', rma, '"'), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  list(backend = "metafor", method_key = method, method = method)
}

#' Internal capability-matrix method key from user-facing `method`
#' @noRd
interactmeta_method_capability_key <- function(method) {
  method <- as.character(method)[[1L]]
  if (interactmeta_cams_method_is_bayesian(method) || identical(method, "Bayesian")) {
    return("bayes")
  }
  method
}

#' Whether a fit uses Bayesian DIRECT / mixture inference
#' @noRd
interactmeta_is_bayesian_fit <- function(x) {
  if (inherits(x, "summary.interactmeta") && is.list(x$fit)) {
    return(interactmeta_is_bayesian_fit(x$fit))
  }
  if (inherits(x, "interactmeta_overall_ifrac")) {
    fit <- attr(x, "interactmeta", exact = TRUE)
    if (!is.null(fit)) {
      return(interactmeta_is_bayesian_fit(fit))
    }
  }
  inf <- if (is.list(x)) {
    as.character(x$inference %||% NA_character_)[[1L]]
  } else {
    NA_character_
  }
  if (identical(inf, "bayesian")) {
    return(TRUE)
  }
  if (identical(inf, "frequentist")) {
    return(FALSE)
  }
  method <- if (is.list(x)) {
    as.character(x$method %||% NA_character_)[[1L]]
  } else {
    NA_character_
  }
  if (interactmeta_cams_method_is_bayesian(method) || identical(method, "Bayesian")) {
    return(TRUE)
  }
  if (nzchar(method) && !identical(method, "Bayesian")) {
    return(FALSE)
  }
  identical(interactmeta_fit_backend(x), "bayesmeta")
}

#' Backend package slot for a CAMS / construction fit
#' @noRd
interactmeta_fit_backend <- function(x) {
  if (is.list(x) && !is.null(x$bayesmeta)) {
    return("bayesmeta")
  }
  if (is.list(x) && !is.null(x$metafor)) {
    return("metafor")
  }
  org <- attr(x, "interactmeta", exact = TRUE)
  if (is.list(org) && !is.null(org$backend)) {
    return(as.character(org$backend)[[1L]])
  }
  org2 <- attr(x, "interactmeta_origin", exact = TRUE)
  if (is.list(org2) && !is.null(org2$fit_backend)) {
    return(as.character(org2$fit_backend)[[1L]])
  }
  if (is.list(x) && !is.null(x$fit) && inherits(x$fit, "interactmeta")) {
    return(interactmeta_fit_backend(x$fit))
  }
  method <- as.character(x$method %||% NA_character_)[[1L]]
  if (nzchar(method)) {
    eng <- interactmeta_engine_for_method(method)
    if (!is.na(eng)) {
      return(eng)
    }
  }
  NA_character_
}

#' User-facing `method` label (matches `cams(..., method = )` and compare tables)
#' @noRd
interactmeta_method_display <- function(x) {
  method <- as.character(x$method %||% NA_character_)[[1L]]
  if (identical(method, "Bayesian") ||
      interactmeta_cams_method_is_bayesian(method)) {
    return("Bayesian")
  }
  if (nzchar(method) && !identical(method, "bayes")) {
    return(method)
  }
  if (interactmeta_is_bayesian_fit(x)) {
    return("Bayesian")
  }
  if (identical(interactmeta_fit_backend(x), "metafor")) {
    return(as.character(x$method %||% "REML")[[1L]])
  }
  if (nzchar(method)) {
    return(method)
  }
  NA_character_
}

#' Public engine that implements a method
#' @noRd
interactmeta_engine_for_method <- function(method, allowed = interactmeta_public_engines()) {
  method <- as.character(method)[[1L]]
  if (interactmeta_cams_method_is_bayesian(method) || identical(method, "bayes")) {
    engine <- "bayesmeta"
  } else if (method %in% interactmeta_metafor_rma_methods()) {
    engine <- "metafor"
  } else {
    engine <- switch(
      method,
      bayes = "bayesmeta",
      REML = "metafor",
      ML = "metafor",
      NA_character_
    )
  }
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
    if (identical(approach, "bayes") && !identical(engine, "bayesmeta")) {
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
  if (identical(engine, "metafor") && !method %in% interactmeta_metafor_rma_methods()) {
    stop(
      "`engine = \"metafor\"` requires a `metafor::rma()` method ",
      "(e.g. `\"REML\"` or `\"ML\"`).",
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
  stop("Unknown CAMS engine `", engine, "`.", call. = FALSE)
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
    c("delta", "epsilon", "interval.type", "convolve_derived", "test", "verbose")
  )
  for (nm in extra) {
    if (!nm %in% names(engine_args)) {
      engine_args[[nm]] <- dots[[nm]]
    }
  }
  engine_args
}

#' Merge `...` into metafor `rma()` args (`engine_args` wins on conflicts)
#' @noRd
interactmeta_merge_rma_args <- function(engine_args, dots, exclude = character()) {
  if (is.null(engine_args)) {
    engine_args <- list()
  }
  if (!is.list(engine_args)) {
    stop("`engine_args` must be a list.", call. = FALSE)
  }
  exclude <- unique(c(
    exclude,
    "method", "engine", "approach", "engine_args"
  ))
  for (nm in setdiff(names(dots), exclude)) {
    if (!nm %in% names(engine_args)) {
      engine_args[[nm]] <- dots[[nm]]
    }
  }
  engine_args
}

interactmeta_computational_dot_names <- function() {
  c("delta", "epsilon", "interval.type", "convolve_derived", "test", "verbose")
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

# ---- interactmeta_capabilities.R ----
#' Analysis / inference / implementation capabilities
#'
#' Three distinct axes:
#' \itemize{
#'   \item **analysis** — DA, AD, SWADA, or CAMS
#'   \item **inference** — Bayesian or frequentist
#'   \item **engine** — computational backend that implements that inference
#' }
#' CAMS is one statistical model. Bayesian CAMS is implemented using
#' `bayesmeta` (DIRECT); frequentist CAMS is implemented using `metafor`
#' (ML/REML). DA and AD are not inherently frequentist. Bayesian inference
#' does not imply CAMS information-fraction functionality. IF quantities and
#' IF plots are CAMS-only.
#'
#' @return A data frame with one row per **public** analysis × engine ×
#'   method combination, including heterogeneity, posterior, IF functionality,
#'   and IF plots. Combinations that are theoretically possible but not part
#'   of the public API are omitted.
#' @export
#' @examples
#' interactmeta_capability_matrix()
interactmeta_capability_matrix <- function() {
  rbind(
    interactmeta_capability_row("cam", "metafor", "REML", "optimization", "frequentist", TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    interactmeta_capability_row("cam", "metafor", "ML", "optimization", "frequentist", TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    interactmeta_capability_row("cam", "bayesmeta", "bayes", "quadrature", "bayesian", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    interactmeta_capability_row("da", "metafor", "REML", "optimization", "frequentist", TRUE, FALSE, FALSE, FALSE, TRUE, TRUE),
    interactmeta_capability_row("da", "metafor", "ML", "optimization", "frequentist", TRUE, FALSE, FALSE, FALSE, TRUE, TRUE),
    interactmeta_capability_row("da", "bayesmeta", "bayes", "quadrature", "bayesian", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
    interactmeta_capability_row("ad", "metafor", "REML", "optimization", "frequentist", TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
    interactmeta_capability_row("ad", "metafor", "ML", "optimization", "frequentist", TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
    interactmeta_capability_row("ad", "bayesmeta", "bayes", "quadrature", "bayesian", TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
    interactmeta_capability_row("swada", "metafor", "REML", "optimization", "frequentist", TRUE, FALSE, FALSE, FALSE, TRUE, TRUE),
    interactmeta_capability_row("swada", "metafor", "ML", "optimization", "frequentist", TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
}

#' @noRd
interactmeta_capability_row <- function(analysis,
                                        engine,
                                        method,
                                        algorithm,
                                        inference,
                                        heterogeneity,
                                        posterior,
                                        ifrac,
                                        ifrac_plots,
                                        interaction = TRUE,
                                        subgroups = TRUE) {
  data.frame(
    analysis = analysis,
    engine = engine,
    method = method,
    algorithm = algorithm,
    inference = inference,
    heterogeneity = isTRUE(heterogeneity),
    posterior = isTRUE(posterior),
    ifrac = isTRUE(ifrac),
    ifrac_plots = isTRUE(ifrac_plots),
    interaction = isTRUE(interaction),
    subgroups = isTRUE(subgroups),
    stringsAsFactors = FALSE
  )
}

#' Engines implemented for an analysis family
#' @noRd
interactmeta_analysis_engines <- function(analysis) {
  analysis <- as.character(analysis)[[1L]]
  switch(
    analysis,
    cam = interactmeta_engines(),
    da = c("metafor", "bayesmeta"),
    ad = c("metafor", "bayesmeta"),
    swada = "metafor",
    conventional = c("metafor", "bayesmeta"),
    stop("Unknown analysis `", analysis, "`.", call. = FALSE)
  )
}

#' Resolve `method` for DA / AD / SWADA (same vocabulary as [cams()])
#' @noRd
interactmeta_resolve_construction_method <- function(analysis, method) {
  analysis <- as.character(analysis)[[1L]]
  if (is.null(method) || (is.character(method) && length(method) > 1L)) {
    method <- if (identical(analysis, "swada")) "REML" else "Bayesian"
  }
  resolved <- interactmeta_resolve_cams_method(method)
  if (identical(analysis, "swada") && !identical(resolved$backend, "metafor")) {
    stop(
      "SWADA supports frequentist inference only (`method = \"REML\"` or ",
      "`\"ML\"`). A Bayesian `method` (e.g. `method = \"Bayesian\"`) is ",
      "not implemented for [swada()].",
      call. = FALSE
    )
  }
  list(
    analysis = analysis,
    backend = resolved$backend,
    method = resolved$method,
    method_key = resolved$method_key,
    algorithm = interactmeta_algorithm_for(resolved$backend, resolved$method_key),
    inference = if (identical(resolved$method_key, "bayes")) "bayesian" else "frequentist"
  )
}

#' @noRd
interactmeta_resolve_construction_engine_method <- function(analysis,
                                                            engine,
                                                            method,
                                                            approach = NULL) {
  if (!is.null(engine) && !interactmeta_engine_unspecified(engine)) {
    stop(
      "`engine` is no longer used; set `method = \"Bayesian\"` or a ",
      "`metafor::rma()` method.",
      call. = FALSE
    )
  }
  if (!is.null(approach)) {
    stop("`approach` is no longer used; set `method` instead.", call. = FALSE)
  }
  interactmeta_resolve_construction_method(analysis, method)
}

#' Capability flags for one implemented combination
#' @noRd
interactmeta_capabilities_for <- function(analysis, engine, method) {
  mat <- interactmeta_capability_matrix()
  hit <- mat$analysis == analysis & mat$engine == engine & mat$method == method
  if (!any(hit) && identical(analysis, "cam") && identical(engine, "metafor") &&
      method %in% interactmeta_metafor_rma_methods()) {
    hit <- mat$analysis == analysis & mat$engine == engine & mat$method == "REML"
  }
  if (!any(hit)) {
    stop(
      "Unsupported combination: analysis = ",
      analysis,
      ", engine = ",
      engine,
      ", method = ",
      method,
      ".",
      call. = FALSE
    )
  }
  row <- mat[which(hit)[[1L]], , drop = FALSE]
  list(
    analysis = row$analysis[[1L]],
    engine = row$engine[[1L]],
    method = row$method[[1L]],
    algorithm = row$algorithm[[1L]],
    inference = row$inference[[1L]],
    heterogeneity = isTRUE(row$heterogeneity[[1L]]),
    posterior = isTRUE(row$posterior[[1L]]),
    ifrac = isTRUE(row$ifrac[[1L]]),
    ifrac_plots = isTRUE(row$ifrac_plots[[1L]]),
    interaction = isTRUE(row$interaction[[1L]]),
    subgroups = isTRUE(row$subgroups[[1L]])
  )
}

#' Attach analysis / method / capability metadata
#' @noRd
interactmeta_attach_capabilities <- function(x, analysis, method_display) {
  backend <- interactmeta_fit_backend(x)
  if (is.na(backend)) {
    backend <- interactmeta_engine_for_method(method_display)
  }
  method_key <- interactmeta_method_capability_key(method_display)
  cap <- interactmeta_capabilities_for(analysis, backend, method_key)
  x$analysis <- cap$analysis
  x$method <- method_display
  x$algorithm <- cap$algorithm
  x$inference <- cap$inference
  x$capabilities <- cap
  x
}

#' Analysis family of a fitted object
#' @noRd
interactmeta_fit_analysis <- function(x) {
  if (!is.null(x$analysis) && nzchar(as.character(x$analysis)[[1L]])) {
    return(as.character(x$analysis)[[1L]])
  }
  if (inherits(x, "interactmeta_da")) {
    return("da")
  }
  if (inherits(x, "interactmeta_ad")) {
    return("ad")
  }
  if (inherits(x, "interactmeta_swada")) {
    return("swada")
  }
  if (inherits(x, "interactmeta_conventional")) {
    return("conventional")
  }
  if (inherits(x, "interactmeta")) {
    return("cam")
  }
  NA_character_
}

#' @noRd
interactmeta_fit_capabilities <- function(x) {
  stored <- x$capabilities
  if (is.list(stored) && !is.null(stored$ifrac)) {
    return(stored)
  }
  analysis <- interactmeta_fit_analysis(x)
  if (is.na(analysis) || identical(analysis, "conventional")) {
    return(list(
      analysis = analysis,
      ifrac = FALSE,
      ifrac_plots = FALSE,
      posterior = FALSE,
      heterogeneity = FALSE,
      interaction = !is.na(analysis),
      subgroups = identical(analysis, "conventional")
    ))
  }
  meth <- as.character(x$method %||% NA_character_)[[1L]]
  backend <- interactmeta_fit_backend(x)
  method_key <- interactmeta_method_capability_key(meth)
  if (!is.na(backend) && nzchar(meth)) {
    return(tryCatch(
      interactmeta_capabilities_for(analysis, backend, method_key),
      error = function(e) {
        list(
          analysis = analysis,
          ifrac = identical(analysis, "cam"),
          ifrac_plots = identical(analysis, "cam"),
          interaction = TRUE,
          subgroups = !identical(analysis, "ad")
        )
      }
    ))
  }
  list(
    analysis = analysis,
    ifrac = identical(analysis, "cam"),
    ifrac_plots = identical(analysis, "cam"),
    interaction = TRUE,
    subgroups = !identical(analysis, "ad")
  )
}

#' @noRd
interactmeta_is_non_cams_analysis <- function(x) {
  analysis <- interactmeta_fit_analysis(x)
  analysis %in% c("da", "ad", "swada", "conventional") ||
    inherits(x, interactmeta_non_cams_classes())
}

#' Error when a quantity is not defined for this analysis
#' @noRd
interactmeta_require_capability <- function(x, capability, feature = NULL) {
  cap <- interactmeta_fit_capabilities(x)
  if (isTRUE(cap[[capability]])) {
    return(invisible(TRUE))
  }
  analysis <- cap$analysis %||% interactmeta_fit_analysis(x)
  feature <- feature %||% paste0(capability, "()")
  inference <- cap$inference %||% x$inference %||% NA_character_
  bayes_note <- if (identical(inference, "bayesian") &&
                      capability %in% c("ifrac", "ifrac_plots")) {
    " A Bayesian posterior does not add CAMS information-fraction functionality."
  } else {
    ""
  }
  if (capability %in% c("ifrac", "ifrac_plots")) {
    stop(
      feature,
      " requires the CAMS information-fraction construction. ",
      if (identical(analysis, "da")) {
        "DA fits do not contain CAMS information fractions."
      } else if (identical(analysis, "ad")) {
        "AD fits do not contain CAMS information fractions."
      } else if (identical(analysis, "swada")) {
        "SWADA does not use the CAMS information-fraction construction, so CAMS-specific IF quantities and IF-based visualizations are not available."
      } else if (identical(analysis, "conventional")) {
        "conventional() is a DA/AD bundle, not a CAMS fit, and does not contain CAMS information fractions."
      } else {
        paste0(
          "This object is not a CAMS fit from cams()."
        )
      },
      bayes_note,
      call. = FALSE
    )
  }
  if (identical(capability, "subgroups") && identical(analysis, "ad")) {
    stop(
      feature,
      " is not defined for AD. Average Difference (AD) estimates only the pooled within-trial interaction; use `interactions(fit)`.",
      call. = FALSE
    )
  }
  stop(
    feature,
    " is not defined for a ",
    interactmeta_analysis_acronym(analysis),
    " analysis.",
    call. = FALSE
  )
}

#' Marginal `type` for DA / AD / SWADA (default API vector means Overall only)
#' @noRd
interactmeta_construction_marginal_type <- function(type, analysis = NULL) {
  type <- as.character(type)
  if (length(type) == 2L && identical(type, c("Overall", "Prediction"))) {
    if (identical(analysis, "swada")) {
      return(type)
    }
    return("Overall")
  }
  match.arg(type, c("Overall", "Prediction"), several.ok = TRUE)
}

#' Reject CAMS-only marginal arguments on DA / AD / SWADA fits
#' @noRd
interactmeta_validate_construction_marginal_args <- function(
    fit,
    feature,
    type = NULL,
    shrinkage = FALSE,
    at = NULL,
    ifrac = NULL) {
  if (!interactmeta_is_non_cams_analysis(fit)) {
    return(invisible(TRUE))
  }
  analysis <- interactmeta_fit_analysis(fit)
  acro <- interactmeta_analysis_acronym(analysis)

  if (isTRUE(shrinkage)) {
    stop(
      feature,
      " with `shrinkage = TRUE` requires a CAMS fit from `cams()`. ",
      acro,
      " constructions do not support CAMS shrinkage.",
      call. = FALSE
    )
  }
  if (!is.null(at)) {
    stop(
      feature,
      " does not accept `at` (or `contribution`) for ",
      acro,
      " fits; reporting information fractions are CAMS-only.",
      call. = FALSE
    )
  }
  if (!is.null(ifrac) && length(ifrac) && !identical(as.character(ifrac)[[1L]], "interaction")) {
    stop(
      feature,
      " does not accept `ifrac = \"overall\"` for ",
      acro,
      " fits.",
      call. = FALSE
    )
  }
  if (!is.null(type)) {
    type <- as.character(type)
    bad <- setdiff(type, c("Overall", "Prediction"))
    if (length(bad)) {
      stop(
        feature,
        " `type` must be `\"Overall\"` and/or `\"Prediction\"`.",
        call. = FALSE
      )
    }
    bad <- setdiff(type, "Overall")
    if (length(bad) && !identical(analysis, "swada")) {
      if (identical(analysis, "ad") && identical(feature, "modifiers()")) {
        stop(
          "modifiers() is not defined for AD. Average Difference (AD) estimates only the pooled within-trial interaction; use `interactions(fit)`.",
          call. = FALSE
        )
      }
      hint <- if (identical(analysis, "ad") && identical(feature, "interactions()")) {
        " AD reports a single pooled interaction; use the default `type = \"Overall\"` only."
      } else {
        " Use `type = \"Overall\"` only, or fit `cams()` or `swada()` for SWADA prediction intervals."
      }
      stop(
        feature,
        " does not support `type = \"",
        paste(bad, collapse = "\", \""),
        "\"` for ",
        acro,
        " fits.",
        hint,
        call. = FALSE
      )
    }
    if (length(bad) && identical(analysis, "swada") &&
        identical(feature, "modifiers()")) {
      stop(
        "modifiers() for SWADA supports `type = \"Overall\"` only; ",
        "use `interactions(fit, type = \"Prediction\")` for the ",
        "SWADA interaction prediction interval.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

# ---- interactmeta_origin.R ----
#' Engine tagging and derived-effect helpers for one CAMS fit
#'
#' CAMS downstream quantities are transformations of one selected engine's
#' fit. SWADA / DA / AD are not CAMS fits and must not borrow IF quantities.
#' @noRd
NULL

#' Classes that are not CAMS fits
#' @noRd
interactmeta_non_cams_classes <- function() {
  c(
    "interactmeta_swada",
    "interactmeta_da",
    "interactmeta_ad",
    "interactmeta_conventional"
  )
}

#' Informative error when a CAMS/IF method is applied to SWADA / DA / AD
#' @noRd
interactmeta_abort_non_cams <- function(object, feature) {
  if (!interactmeta_is_non_cams_analysis(object)) {
    return(invisible(NULL))
  }
  interactmeta_require_capability(object, "ifrac", feature = feature)
}

#' Native backend object stored on a fit (no refit)
#'
#' For Bayesian CAMS (`engine = "bayesmeta"`), returns the list of stored
#' univariate fits: `$overall` ([bayesmeta::bmr()] of \eqn{m_j} on \eqn{\pi_j}),
#' `$interaction` ([bayesmeta::bayesmeta()] of \eqn{g_j}), and `$trial_blocks`.
#' Those objects are the `bayesmeta` / DIRECT backend, not CAMS estimands.
#' Plot them with the native `bayesmeta` methods
#' (`vignette("cams")`). For frequentist CAMS the
#' analogue is `$metafor`.
#'
#' @param fit A [cams()] / [da()] / [ad()] / [swada()] / [conventional()] object.
#' @return For Bayesian CAMS, the list with `$overall`, `$interaction`, and
#'   `$trial_blocks`. For frequentist CAMS, `fit$metafor`. For DA / AD /
#'   SWADA / conventional, the stored `$fits` list. `NULL` if none is stored.
#' @seealso [traceplot.interactmeta()], `vignette("cams")`
#' @keywords internal
#' @export
native_fit <- function(fit) {
  analysis <- interactmeta_fit_analysis(fit)
  if (!is.na(analysis) && analysis %in% c("da", "ad", "swada", "conventional")) {
    return(fit$fits)
  }
  check_comparison_fit(fit)
  backend <- interactmeta_fit_backend(fit)
  if (identical(backend, "bayesmeta")) {
    return(fit$bayesmeta)
  }
  if (identical(backend, "metafor")) {
    return(fit$metafor)
  }
  NULL
}

#' @noRd
interactmeta_attach_origin <- function(x,
                                      fit,
                                      quantity,
                                      fitted_or_derived = "derived",
                                      uncertainty = NULL) {
  backend <- interactmeta_fit_backend(fit)
  meth <- as.character(fit$method %||% NA_character_)[[1L]]
  algo <- as.character(fit$algorithm %||% NA_character_)[[1L]]
  method_key <- interactmeta_method_capability_key(meth)
  if (is.null(uncertainty)) {
    uncertainty <- interactmeta_default_uncertainty(backend, method_key)
  }
  attr(x, "interactmeta_origin") <- list(
    fit_backend = backend,
    fit_method = meth,
    fit_algorithm = algo,
    quantity = quantity,
    fitted_or_derived = fitted_or_derived,
    uncertainty = uncertainty,
    call = fit$call
  )
  x
}

#' Inferential interval style for print / summary tables
#'
#' `"frequentist"` is `engine = "metafor"` only. Bayesian engines keep
#' `"bayesian"` even when `$support` happens to have one row.
#' @noRd
interactmeta_interval_style <- function(x) {
  if (interactmeta_is_bayesian_fit(x)) {
    "bayesian"
  } else {
    "frequentist"
  }
}

#' @rdname interactmeta_fit_backend
#' @noRd
interactmeta_origin_engine <- function(x) {
  interactmeta_fit_backend(x)
}

#' Whether IF pooling uses fitted heterogeneity (not a posterior)
#' @noRd
interactmeta_ifrac_is_plugin <- function(x) {
  !interactmeta_is_bayesian_fit(x)
}

#' @noRd
interactmeta_default_uncertainty <- function(engine, method) {
  if (identical(engine, "bayesmeta")) {
    return(paste(
      "central interval from the DIRECT-node Gaussian mixture",
      "(.weight and within-node .cond_*)"
    ))
  }
  if (identical(engine, "metafor")) {
    return("Delta method from the selected metafor::rma coefficient covariance")
  }
  stop("Unknown CAMS engine `", engine, "`.", call. = FALSE)
}

#' Mixture / weighted moments for a scalar derived from one fit
#'
#' `cond_sd` is the within-node standard deviation.
#' @noRd
interactmeta_mixture_moments <- function(mean,
                                         weight = NULL,
                                         cond_sd = NULL,
                                         prob = 0.95,
                                         exponentiate = FALSE,
                                         interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(interval.type)
  mean <- as.numeric(mean)
  n <- length(mean)
  empty <- list(
    median = NA_real_,
    mean = NA_real_,
    sd = NA_real_,
    lci = NA_real_,
    uci = NA_real_,
    estimate = NA_real_,
    prob = prob
  )
  if (!n) {
    return(empty)
  }
  if (is.null(weight) || length(weight) != n) {
    weight <- rep(1 / n, n)
  }
  w <- as.numeric(weight)
  w[!is.finite(w) | w < 0] <- 0
  if (!sum(w)) {
    w <- rep(1 / n, n)
  }
  w <- w / sum(w)
  if (is.null(cond_sd) || length(cond_sd) != n) {
    cond_sd <- rep(0, n)
  }
  s <- pmax(as.numeric(cond_sd), 0)
  s[!is.finite(s)] <- 0

  ok <- is.finite(mean)
  if (!any(ok)) {
    return(empty)
  }
  mean <- mean[ok]
  w <- w[ok]
  s <- s[ok]
  w <- w / sum(w)

  if (length(mean) == 1L && s[[1L]] <= 0) {
    v <- mean[[1L]]
    if (isTRUE(exponentiate)) {
      v <- exp(v)
    }
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

  if (all(s <= 0) && length(unique(w)) == 1L) {
    if (isTRUE(exponentiate)) {
      mean <- exp(mean)
    }
    return(comparison_draw_moments(
      mean,
      prob = prob,
      interval.type = interval.type
    ))
  }

  gm <- interactmeta_gm(mean = mean, sd = s, weight = w)
  sm <- interactmeta_gm_summary(
    gm,
    prob = prob,
    interval.type = interval.type
  )
  if (isTRUE(exponentiate)) {
    ln <- interactmeta_gm_lognormal_moments(gm)
    return(list(
      median = exp(sm$median),
      mean = ln$mean,
      sd = ln$sd,
      lci = exp(sm$lower),
      uci = exp(sm$upper),
      estimate = exp(sm$median),
      prob = prob
    ))
  }
  list(
    median = sm$median,
    mean = sm$mean,
    sd = sm$sd,
    lci = sm$lower,
    uci = sm$upper,
    estimate = sm$median,
    prob = prob
  )
}

#' Conditional SD of a linear combination of CAMS primitives on one draws table
#'
#' Uses the stored overall-block joint moments `Var(alpha)`, `Var(theta)`,
#' and `Cov(alpha, theta)`, plus independent `Var(gamma)`. This is how
#' population `M* = alpha + pi* theta` keeps the covariance; do not replace
#' it with an independent convolution of the `alpha` and `theta` mixture
#' marginals.
#' Coefficients \eqn{(a_\alpha,a_\theta,a_\gamma)} apply to
#' \eqn{(\alpha,\theta,\gamma)}, **not** \eqn{(\alpha, b_{\mathrm{ifrac}}, \gamma)}.
#' Stored `.cond_var_theta` is \eqn{\mathrm{Var}(\theta)} from the overall
#' `bmr` block. Use \eqn{(1, \pi, 0)} for \eqn{M=\alpha+\pi\theta},
#' \eqn{(1, \pi, -\pi)} for \eqn{\mu_A}, and \eqn{(1, \pi, 1-\pi)} for \eqn{\mu_B}.
#' @noRd
interactmeta_lincomb_cond_sd <- function(post,
                                         a_alpha,
                                         a_theta,
                                         a_gamma) {
  n <- nrow(post)
  rec <- function(a) {
    a <- as.numeric(a)
    if (length(a) == 1L) {
      a <- rep(a, n)
    }
    if (length(a) != n) {
      stop("Linear-combination coefficients must recycle to the draws table.", call. = FALSE)
    }
    a
  }
  a_alpha <- rec(a_alpha)
  a_theta <- rec(a_theta)
  a_gamma <- rec(a_gamma)
  va <- if (".cond_var_alpha" %in% names(post)) {
    pmax(as.numeric(post$.cond_var_alpha), 0)
  } else {
    rep(0, n)
  }
  vt <- if (".cond_var_theta" %in% names(post)) {
    pmax(as.numeric(post$.cond_var_theta), 0)
  } else {
    rep(0, n)
  }
  cat <- if (".cond_cov_alpha_theta" %in% names(post)) {
    as.numeric(post$.cond_cov_alpha_theta)
  } else {
    rep(0, n)
  }
  vg <- if (".cond_sd_gamma" %in% names(post)) {
    pmax(as.numeric(post$.cond_sd_gamma), 0)^2
  } else {
    rep(0, n)
  }
  va[!is.finite(va)] <- 0
  vt[!is.finite(vt)] <- 0
  cat[!is.finite(cat)] <- 0
  vg[!is.finite(vg)] <- 0
  sqrt(pmax(
    a_alpha^2 * va + a_theta^2 * vt + 2 * a_alpha * a_theta * cat + a_gamma^2 * vg,
    0
  ))
}

#' Extra predictive SD from new-trial random effects (same fit's tau)
#' @noRd
interactmeta_prediction_re_sd <- function(post, omega, x = NULL) {
  n <- nrow(post)
  omega <- as.numeric(omega)
  if (length(omega) == 1L) {
    omega <- rep(omega, n)
  }
  tau_m <- if ("sd_trial__Intercept" %in% names(post)) {
    pmax(as.numeric(post$sd_trial__Intercept), 0)
  } else {
    rep(0, n)
  }
  slope <- grep("^sd_trial__ifrac2_", names(post), value = TRUE)
  tau_g <- if (length(slope)) {
    pmax(as.numeric(post[[slope[[1L]]]]), 0)
  } else {
    rep(0, n)
  }
  if (is.null(x)) {
    sqrt(pmax(tau_m^2 + (omega * tau_g)^2, 0))
  } else {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, n)
    }
    sqrt(pmax(tau_m^2 + ((x - omega) * tau_g)^2, 0))
  }
}

#' Whether this posterior table carries quadrature conditional moments
#' @noRd
interactmeta_has_cond_moments <- function(post) {
  is.data.frame(post) &&
    ".cond_var_alpha" %in% names(post) &&
    ".cond_sd_gamma" %in% names(post)
}

#' Fill metafor Wald intervals from the stored rma objects (no refit)
#' @noRd
interactmeta_metafor_lincomb <- function(fit,
                                         L_ov,
                                         L_int,
                                         pred_var = 0,
                                         exponentiate = FALSE,
                                         prob = 0.95) {
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  if (is.null(ov) || is.null(int)) {
    return(NULL)
  }
  b_ov <- as.numeric(ov$beta)
  b_int <- as.numeric(int$beta)[[1L]]
  L_ov <- as.numeric(L_ov)
  L_int <- as.numeric(L_int)[[1L]]
  if (length(b_ov) != length(L_ov)) {
    return(NULL)
  }

  est <- sum(L_ov * b_ov) + L_int * b_int
  vb_ov <- as.matrix(ov$vb)
  vb_int <- as.matrix(int$vb)
  var <- as.numeric(t(L_ov) %*% vb_ov %*% L_ov) +
    (L_int^2) * as.numeric(vb_int[[1L]]) +
    as.numeric(pred_var)[[1L]]
  se <- sqrt(max(var, 0))
  z <- stats::qnorm(1 - (1 - prob) / 2)
  lo <- est - z * se
  hi <- est + z * se
  if (isTRUE(exponentiate)) {
    est_y <- exp(est)
    return(list(
      estimate = est_y,
      mean = est_y,
      sd = est_y * se,
      lci = exp(lo),
      uci = exp(hi),
      prob = prob
    ))
  }
  list(
    estimate = est,
    mean = est,
    sd = se,
    lci = lo,
    uci = hi,
    prob = prob
  )
}

#' @noRd
interactmeta_metafor_omega <- function(fit, overall_ifrac = NULL) {
  if (inherits(overall_ifrac, "interactmeta_overall_ifrac")) {
    return(stats::median(as.numeric(overall_ifrac), na.rm = TRUE))
  }
  if (is.numeric(overall_ifrac) && length(overall_ifrac)) {
    return(stats::median(as.numeric(overall_ifrac), na.rm = TRUE))
  }
  oi <- tryCatch(ifrac(fit, summary = FALSE), error = function(e) NULL)
  if (is.null(oi)) {
    return(NA_real_)
  }
  stats::median(as.numeric(oi), na.rm = TRUE)
}

#' Replace NA intervals on a metafor effect summary using the stored rma vcov
#' @noRd
interactmeta_metafor_fill_effects <- function(fit,
                                              summary,
                                              quantity,
                                              omega,
                                              exponentiate = FALSE,
                                              prob = 0.95,
                                              delta_ifrac2 = 1) {
  if (!!interactmeta_is_bayesian_fit(fit) || is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  if (all(is.finite(summary$lci)) && all(is.finite(summary$uci))) {
    return(summary)
  }
  omega <- as.numeric(omega)[[1L]]
  if (!is.finite(omega)) {
    return(summary)
  }
  tau_m2 <- max(as.numeric(fit$metafor$overall$tau2), 0)
  tau_g2 <- max(as.numeric(fit$metafor$interaction$tau2), 0)
  for (i in seq_len(nrow(summary))) {
    ty <- if ("type" %in% names(summary)) {
      as.character(summary$type[[i]])
    } else {
      "Overall"
    }
    pred <- identical(ty, "Prediction")
    lin <- NULL
    if (identical(quantity, "interaction")) {
      L_ov <- c(0, 0)
      L_int <- 1
      pvar <- if (pred) tau_g2 * (as.numeric(delta_ifrac2)[[1L]]^2) else 0
      lin <- interactmeta_metafor_lincomb(
        fit, L_ov, L_int, pred_var = pvar,
        exponentiate = exponentiate, prob = prob
      )
    } else if (identical(quantity, "combined")) {
      L_ov <- c(1, omega)
      L_int <- 0
      pvar <- if (pred) tau_m2 + (omega^2) * tau_g2 else 0
      lin <- interactmeta_metafor_lincomb(
        fit, L_ov, L_int, pred_var = pvar,
        exponentiate = exponentiate, prob = prob
      )
    } else if (identical(quantity, "subgroups") &&
               interactmeta_modifier_summary_col(summary) %in% names(summary)) {
      st <- attr(fit, "interactmeta", exact = TRUE)
      ref <- as.character(st$reference %||% "")[[1L]]
      lvl_col <- interactmeta_modifier_summary_col(summary)
      sg <- as.character(summary[[lvl_col]][[i]])
      if (identical(sg, ref)) {
        L_ov <- c(1, omega)
        L_int <- -omega
        pvar <- if (pred) tau_m2 + (omega^2) * tau_g2 else 0
      } else {
        L_ov <- c(1, omega)
        L_int <- 1 - omega
        pvar <- if (pred) tau_m2 + ((1 - omega)^2) * tau_g2 else 0
      }
      lin <- interactmeta_metafor_lincomb(
        fit, L_ov, L_int, pred_var = pvar,
        exponentiate = exponentiate, prob = prob
      )
    }
    if (is.null(lin)) {
      next
    }
    summary$estimate[[i]] <- lin$estimate
    summary$mean[[i]] <- lin$mean
    summary$sd[[i]] <- lin$sd
    summary$lci[[i]] <- lin$lci
    summary$uci[[i]] <- lin$uci
    summary$prob[[i]] <- lin$prob
  }
  summary
}

#' Attach quadrature weights / within-tau SDs onto effect-draw rows
#' @noRd
interactmeta_decorate_effect_draws <- function(df, post, cond_sd = NULL) {
  if (is.data.frame(post) && ".weight" %in% names(post)) {
    if (nrow(df) == nrow(post)) {
      df$.weight <- post$.weight
    } else if (nrow(df) > 0L && nrow(post) > 0L && (nrow(df) %% nrow(post) == 0L)) {
      df$.weight <- rep(post$.weight, length.out = nrow(df))
    }
  }
  if (!is.null(cond_sd)) {
    df$.cond_sd <- as.numeric(cond_sd)
  }
  df
}

#' Summarize derived estimands from the same fit (no refit)
#' @noRd
interactmeta_summarize_effects <- function(draws,
                                           group_cols,
                                           prob,
                                           exponentiate,
                                           fit,
                                           quantity,
                                           omega,
                                           delta_ifrac2 = 1) {
  mix <- is.data.frame(draws) && ".cond_sd" %in% names(draws)
  if (mix) {
    # Public Bayesian intervals: mixture interval on DIRECT nodes (fit rule).
    it <- interactmeta_match_interval_type(fit$interval.type)
    summ <- comparison_effects_summary(
      draws,
      group_cols = group_cols,
      prob = prob,
      exponentiate = exponentiate,
      interval.type = it
    )
    draws <- comparison_transform_effect_draws(draws, exponentiate)
  } else {
    draws_out <- comparison_transform_effect_draws(draws, exponentiate)
    summ <- comparison_effects_summary(draws_out, group_cols = group_cols, prob = prob)
    summ <- interactmeta_metafor_fill_effects(
      fit = fit,
      summary = summ,
      quantity = quantity,
      omega = omega,
      exponentiate = exponentiate,
      prob = prob,
      delta_ifrac2 = delta_ifrac2
    )
    draws <- draws_out
  }
  list(support = draws, summary = summ)
}

# ---- fit_finalize.R ----
#' CAMS internal parameter names
#'
#' Developer glossary. Public summaries still use descriptive labels
#' (`Ecological_*`, `Interaction_*`, `tau`, …).
#'
#' \describe{
#'   \item{`alpha`}{Combined intercept; population mean at information fraction 0 for the reference subgroup.}
#'   \item{`beta`}{Ecological / composition slope (`b_ifrac`).}
#'   \item{`gamma`}{Interaction / subgroup offset (`b_subgroup`).}
#'   \item{`theta`}{Combined information-fraction slope; CAMS identity `theta = beta + gamma`.}
#'   \item{`tau_m`}{Combined-block heterogeneity SD (`sd_trial__Intercept`).}
#'   \item{`tau_g`}{Interaction-block heterogeneity SD (`sd_trial__ifrac2_*`).}
#'   \item{`pi_j`}{Trial information fraction of the non-reference subgroup.}
#'   \item{`m_j`}{Trial overall (IV-pooled) effect.}
#'   \item{`g_j`}{Trial interaction contrast `y_B - y_A`.}
#' }
#'
#' @keywords internal
#' @name cam-parameter-names
NULL

#' Attach call, model spec, priors, and trial blocks to a CAMS fit
#' @noRd
interactmeta_finalize_fit <- function(fit,
                                      call,
                                      design,
                                      data,
                                      data_fit,
                                      backend,
                                      method_display,
                                      method_key,
                                      priors,
                                      fit_args,
                                      seed,
                                      ecological_slope,
                                      big_var,
                                      contrast,
                                      trial_var,
                                      subgroup_var) {
  algorithm <- fit$algorithm %||% interactmeta_algorithm_for(backend, method_key)
  fit$call <- call
  fit$model_spec <- interactmeta_build_model_spec(
    design = design,
    engine = backend,
    method = method_key,
    algorithm = algorithm,
    ecological_slope = ecological_slope,
    contrast = contrast,
    big_var = big_var,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  fit$priors <- interactmeta_wrap_priors(priors, backend, design)
  fit$computation <- interactmeta_build_computation(
    engine = backend,
    method = method_key,
    engine_args = fit_args,
    seed = seed
  )
  fit$transformed_data <- interactmeta_build_transformed_data(
    fit = fit,
    design = design,
    data = data,
    data_fit = data_fit,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  fit$method <- method_display
  interactmeta_attach_capabilities(fit, "cam", method_display)
  invisible(fit)
}

interactmeta_build_model_spec <- function(design,
                                          engine,
                                          method,
                                          algorithm,
                                          ecological_slope,
                                          contrast,
                                          big_var,
                                          trial_var,
                                          subgroup_var) {
  parameterization <- switch(
    engine,
    bayesmeta = "orthogonal_blocks (overall bmr + interaction bayesmeta)",
    metafor = "orthogonal_blocks (overall rma + interaction rma)",
    engine
  )
  heterogeneity <- list(
    tau_m = "overall block / trial intercept SD",
    tau_g = "interaction block / ifrac2 slope SD",
    correlated = FALSE
  )
  structure(
    list(
      family = "CAMS",
      n_subgroups = design$n_subgroups,
      n_nonref = design$n_nonref,
      reference = design$reference,
      subgroup_levels = design$subgroup_levels,
      non_ref_labels = design$non_ref_labels,
      information_fraction = list(
        definition = "precision weight: I_g = 1/vi; pi_g = I_g / sum_h I_h",
        contrast = contrast,
        ifrac_vars = design$ifrac_vars,
        ifrac2_vars = design$ifrac2_vars
      ),
      ecological_slope = isTRUE(ecological_slope),
      random_effects = "independent overall-block and interaction-block trial effects",
      random_effects_covariance = "diagonal",
      heterogeneity = heterogeneity,
      parameterization = parameterization,
      method = method,
      engine = engine,
      algorithm = algorithm,
      likelihood = "yi ~ Normal(mu, se = sqrt(vi)); vi treated as known",
      missing_cell_variance = big_var,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      formula = design$formula
    ),
    class = "interactmeta_model_spec"
  )
}

#' @export
print.interactmeta_model_spec <- function(x, ...) {
  cat("CAMS model specification\n")
  cat("  subgroups          : ", x$n_subgroups, " (reference = ", x$reference, ")\n", sep = "")
  cat("  ecological slope   : ", x$ecological_slope, "\n", sep = "")
  cat("  parameterization   : ", x$parameterization, "\n", sep = "")
  cat("  method             : ", x$method, "\n", sep = "")
  cat("  likelihood         : ", x$likelihood, "\n", sep = "")
  invisible(x)
}

interactmeta_build_computation <- function(engine, method, engine_args, seed) {
  args <- if (is.null(engine_args)) list() else engine_args
  if ("convolve_derived" %in% names(args)) {
    args$convolve_derived <- NULL
  }
  if (identical(engine, "bayesmeta")) {
    if (is.null(args$delta)) args$delta <- 0.01
    if (is.null(args$epsilon)) args$epsilon <- 1e-4
  }
  list(
    engine_args = args,
    seed = if (is.na(seed)) NA else seed,
    stochastic = FALSE,
    deterministic = TRUE,
    bitwise_reproducible = FALSE
  )
}

interactmeta_build_transformed_data <- function(fit,
                                                design,
                                                data,
                                                data_fit,
                                                trial_var,
                                                subgroup_var) {
  observed <- comparison_fit_display_data(data, design)
  blocks <- NULL
  if (isTRUE(design$n_nonref == 1L)) {
    blocks <- comparison_trial_effect_blocks(
      data = observed,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = design$reference,
      nonref = design$non_ref_labels[[1L]],
      ifrac_var = design$ifrac_vars[[1L]]
    )
  }
  list(
    observed_effects = observed[, intersect(
      c(trial_var, subgroup_var, "yi", "vi", design$ifrac_vars, design$ifrac2_vars),
      names(observed)
    ), drop = FALSE],
    trial_blocks = blocks
  )
}

interactmeta_wrap_priors <- function(priors, engine, design) {
  if (!inherits(priors, "cam_priors")) {
    return(priors)
  }
  engine_rep <- interactmeta_priors_engine_representation(priors, engine, design)
  structure(
    list(
      common = priors,
      engine = engine_rep,
      tau = priors$tau,
      tau_interaction = priors$tau_interaction,
      fixed = priors$fixed
    ),
    class = "interactmeta_priors"
  )
}

#' @export
print.interactmeta_priors <- function(x, ...) {
  if (inherits(x$common, "cam_priors")) {
    print(x$common)
  } else {
    print.cam_priors(x)
  }
  invisible(x)
}

interactmeta_priors_engine_representation <- function(priors, engine, design) {
  if (identical(engine, "bayesmeta")) {
    tau_bm <- prior_to_bayesmeta(priors$tau)
    tau_g_bm <- prior_to_bayesmeta(priors$tau_interaction)
    fe <- prior_to_bayesmeta_fixed(priors$fixed)
    return(list(
      applied = TRUE,
      note = paste(
        "bmr prior is on (alpha, theta), with beta.prior.cov = diag(sd^2);",
        "mu.prior is for gamma only.",
        "Independent normals on (alpha, theta) are not independent normals on",
        "(alpha, beta, gamma), because theta = beta + gamma."
      ),
      tau = if (is.character(tau_bm)) tau_bm else format(priors$tau),
      tau_interaction = if (is.character(tau_g_bm)) tau_g_bm else format(priors$tau_interaction),
      fixed = fe,
      tau_density = if (is.function(tau_bm)) tau_bm else NULL,
      tau_interaction_density = if (is.function(tau_g_bm)) tau_g_bm else NULL
    ))
  }
  list(
    applied = FALSE,
    note = "engine = metafor is frequentist; Bayesian priors are not used"
  )
}

#' Linear predictor table for a CAMS fit
#'
#' Returns the stored `$mu` matrix (quadrature nodes or a point-estimate row
#' times observations). This is the fixed-effect linear predictor, not
#' trial-level reconstruction; see [shrinkage()].
#'
#' @param object An `interactmeta` from [cams()].
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
#' Returns the stored `$support` data frame: DIRECT quadrature nodes for
#' bayesmeta, or a single point-estimate row for metafor.
#'
#' @param x An `interactmeta` from [cams()].
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#' @return A data frame.
#' @export
as.data.frame.interactmeta <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$support
}

#' @noRd
comparison_fit_variables <- function(fit) {
  if (inherits(fit, "interactmeta")) {
    return(names(fit$support))
  }
  names(as.data.frame(fit))
}

# ---- draw_summaries.R ----
#' Parameter table as a data frame
#'
#' Uses stored `$support` (DIRECT nodes for bayesmeta; one point-estimate row
#' for metafor). Adds `.draw` / `.chain` / `.iteration` when missing so
#' grouped summaries can index rows.
#' @noRd
comparison_tidy_support <- function(fit) {
  post <- if (inherits(fit, "interactmeta")) {
    as.data.frame(fit$support)
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

#' Median + posterior interval for a numeric vector
#'
#' Canonical moments (median / mean / sd / shortest or central interval).
#' @noRd
comparison_draw_moments <- function(x, prob = 0.95, weight = NULL, cond_sd = NULL,
                                    exponentiate = FALSE,
                                    interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(interval.type)
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
      exponentiate = exponentiate,
      interval.type = interval.type
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
  qi <- comparison_median_qi_core(
    x,
    prob = prob,
    interval.type = interval.type
  )
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

#' Median with a shortest or central interval
#' @noRd
comparison_median_qi_core <- function(x, prob = 0.95,
                                      interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(interval.type)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  estimate <- stats::median(x)
  if (identical(interval.type, "central")) {
    alpha <- (1 - prob) / 2
    bounds <- stats::quantile(
      x,
      probs = c(alpha, 1 - alpha),
      names = FALSE,
      type = 7
    )
  } else {
    sx <- sort(x)
    width_n <- max(1L, ceiling(prob * length(sx)))
    starts <- seq_len(length(sx) - width_n + 1L)
    ends <- starts + width_n - 1L
    best <- which.min(sx[ends] - sx[starts])
    bounds <- c(sx[starts[[best]]], sx[ends[[best]]])
  }
  list(
    estimate = estimate,
    lci = bounds[[1L]],
    uci = bounds[[2L]],
    prob = prob,
    interval.type = interval.type
  )
}

#' @return list(estimate, lci, uci, prob)
#' @noRd
comparison_median_qi_vec <- function(x, prob = 0.95,
                                     interval.type = NULL) {
  comparison_median_qi_core(
    x,
    prob = prob,
    interval.type = interval.type
  )
}

#' Grouped median / mean / sd / posterior interval summary
#' @noRd
comparison_median_qi_grouped <- function(draws, value_col, group_cols, prob = 0.95,
                                         exponentiate = FALSE,
                                         interval.type = NULL) {
  interval.type <- interactmeta_match_interval_type(interval.type)
  out <- draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(d, ...) {
      w <- if (".weight" %in% names(d)) d$.weight else NULL
      s <- if (".cond_sd" %in% names(d)) d$.cond_sd else NULL
      mom <- comparison_draw_moments(
        d[[value_col]],
        prob = prob,
        weight = w,
        cond_sd = s,
        exponentiate = exponentiate,
        interval.type = interval.type
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
  attr(out, "interval.type") <- interval.type
  out
}
