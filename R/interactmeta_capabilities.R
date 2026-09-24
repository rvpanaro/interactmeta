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

#' Capability rows for undocumented CAMS backends (not user-facing).
#' @noRd
interactmeta_capability_matrix_internal <- function() {
  rbind(
    interactmeta_capability_matrix(),
    interactmeta_capability_row("cam", "brms", "bayes", "MCMC", "bayesian", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
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

#' Resolve engine + method for DA / AD / SWADA
#' @noRd
interactmeta_resolve_construction_engine_method <- function(analysis,
                                                            engine,
                                                            method,
                                                            approach = NULL) {
  analysis <- as.character(analysis)[[1L]]
  allowed <- interactmeta_analysis_engines(analysis)
  engine <- interactmeta_engine_from_approach(
    engine,
    method,
    approach = approach,
    allowed = allowed
  )
  if (!engine %in% allowed) {
    stop(
      interactmeta_analysis_acronym(analysis),
      " supports engine = ",
      paste(paste0('"', allowed, '"'), collapse = " or "),
      ". `engine = \"",
      engine,
      "\"` is not a ",
      interactmeta_analysis_acronym(analysis),
      " backend",
      if (engine %in% interactmeta_internal_engines()) {
        " (CAMS is fitted with bayesmeta or metafor, not as a DA/AD engine)."
      } else {
        "."
      },
      call. = FALSE
    )
  }
  if (is.null(method) || (is.character(method) && length(method) > 1L)) {
    method <- interactmeta_default_method(engine, approach)
  } else {
    method <- as.character(method)[[1L]]
  }
  if (identical(engine, "metafor") && !method %in% c("REML", "ML")) {
    stop(
      "`engine = \"metafor\"` requires `method = \"REML\"` or `\"ML\"`.",
      call. = FALSE
    )
  }
  if (identical(engine, "bayesmeta") && !identical(method, "bayes")) {
    stop("`engine = \"bayesmeta\"` requires `method = \"bayes\"`.", call. = FALSE)
  }
  list(
    analysis = analysis,
    engine = engine,
    method = method,
    algorithm = interactmeta_algorithm_for(engine, method),
    inference = if (identical(method, "bayes")) "bayesian" else "frequentist"
  )
}

#' Capability flags for one implemented combination
#' @noRd
interactmeta_capabilities_for <- function(analysis, engine, method) {
  mat <- interactmeta_capability_matrix_internal()
  hit <- mat$analysis == analysis & mat$engine == engine & mat$method == method
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

#' Attach analysis / engine / capability metadata
#' @noRd
interactmeta_attach_capabilities <- function(x, analysis, engine, method) {
  cap <- interactmeta_capabilities_for(analysis, engine, method)
  x$analysis <- cap$analysis
  x$engine <- cap$engine
  x$method <- cap$method
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
  if (inherits(x, "interactmeta") || inherits(x, "cam") || inherits(x, "cams")) {
    return("cam")
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
  eng <- as.character(x$engine %||% NA_character_)[[1L]]
  meth <- as.character(x$method %||% NA_character_)[[1L]]
  if (!is.na(eng) && !is.na(meth)) {
    return(tryCatch(
      interactmeta_capabilities_for(analysis, eng, meth),
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
          "This object is not a CAMS fit from cam()."
        )
      },
      bayes_note,
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
