# ---- results.R ----
#' Resolve marginal options from `fit$reporting` when arguments are unset
#' @noRd
interactmeta_resolve_reporting_opts <- function(fit,
                                                prob = NULL,
                                                ifrac = NULL,
                                                contrast_sep = NULL,
                                                exponentiate = NULL) {
  rep <- interactmeta_reporting_sidecar(fit) %||% fit$reporting
  if (!is.list(rep)) {
    return(list(
      prob = prob %||% 0.95,
      ifrac = ifrac %||% "interaction",
      contrast_sep = contrast_sep %||% "-",
      exponentiate = exponentiate %||% FALSE
    ))
  }
  list(
    prob = if (is.null(prob)) rep$prob %||% 0.95 else prob,
    ifrac = if (is.null(ifrac)) rep$ifrac %||% "interaction" else ifrac,
    contrast_sep = if (is.null(contrast_sep)) rep$contrast_sep %||% "-" else contrast_sep,
    exponentiate = if (is.null(exponentiate)) isTRUE(rep$exponentiate) else exponentiate
  )
}

#' Call a CAMS marginal quantity extractor
#' @noRd
interactmeta_marginal <- function(fit,
                                  quantity = c("interaction", "modifiers", "combined"),
                                  ...) {
  quantity <- match.arg(quantity)
  switch(
    quantity,
    interaction = interactions(fit, ...),
    modifiers = modifiers(fit, ...),
    combined = combined(fit, ...)
  )
}

#' Tidy all CAMS estimands from a fitted model
#'
#' Returns one long [tibble::tibble] suitable for export, QC, or downstream
#' tables. Uses the same canonical settings as the default [cams()] reporting
#' bundle (`summary = TRUE`, `type = "Overall"` for pooled rows).
#'
#' @param x An [`interactmeta`] object from [cams()].
#' @param estimands Character vector: `"interaction"`, `"modifiers"`,
#'   `"combined"`, `"heterogeneity"`.
#' @param type Passed to marginal extractors (`"Overall"`, `"Prediction"`, or both).
#' @param at Reporting information fraction; `NULL` uses the pooled [ifrac()].
#' @param ifrac Weighting for pooled IF (`"interaction"` or `"overall"`).
#' @param prob Interval probability.
#' @param exponentiate If `NULL`, uses `x$reporting$exponentiate` when set.
#' @param contrast_sep Contrast separator for interaction labels.
#' @param summary Must be `TRUE` for tidy export (pooled summary rows).
#' @param ... Ignored.
#' @return A tibble with `estimand` plus columns from [generics::tidy()].
#' @export
#' @importFrom generics tidy
tidy.interactmeta <- function(x,
                              estimands = c(
                                "interaction",
                                "modifiers",
                                "combined",
                                "heterogeneity"
                              ),
                              type = "Overall",
                              at = NULL,
                              ifrac = NULL,
                              prob = NULL,
                              exponentiate = NULL,
                              contrast_sep = NULL,
                              summary = TRUE,
                              ...) {
  if (!inherits(x, "interactmeta")) {
    stop("`x` must be an interactmeta object from `cams()`.", call. = FALSE)
  }
  if (!isTRUE(summary)) {
    stop(
      "`tidy()` on a CAMS fit requires `summary = TRUE` for export rows.",
      call. = FALSE
    )
  }
  rep <- interactmeta_reporting_sidecar(x) %||% x$reporting
  prob_use <- if (is.null(prob)) {
    if (is.list(rep)) rep$prob %||% 0.95 else 0.95
  } else {
    interactmeta_default_prob(prob)
  }
  interactmeta_ensure_reporting(x, prob = prob_use, warm_marginals = TRUE)
  opts <- interactmeta_resolve_reporting_opts(
    x,
    prob = prob,
    ifrac = ifrac,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate
  )
  ifrac_val <- match.arg(opts$ifrac, c("interaction", "overall"))
  estimands <- match.arg(estimands, several.ok = TRUE)
  st <- attr(x, "interactmeta", exact = TRUE)
  meta <- tibble::tibble(
    reference = st$reference %||% NA_character_,
    method = "CAMS",
    inference = if (interactmeta_is_bayesian_fit(x)) {
      "Bayesian"
    } else {
      "frequentist"
    },
    .method = x$method %||% NA_character_
  )
  pieces <- list()
  marg_args <- list(
    at = at,
    ifrac = ifrac_val,
    type = type,
    prob = opts$prob,
    summary = TRUE,
    contrast_sep = opts$contrast_sep,
    exponentiate = opts$exponentiate
  )
  if ("interaction" %in% estimands) {
    obj <- do.call(interactions, c(list(fit = x), marg_args))
    td <- tidy(obj)
    if (nrow(td)) {
      td$estimand <- "interaction"
      if (!"level" %in% names(td) && "interaction" %in% names(td)) {
        td$level <- td$interaction
      }
      pieces$interaction <- td
    }
  }
  if ("modifiers" %in% estimands) {
    obj <- do.call(modifiers, c(list(fit = x), marg_args))
    td <- tidy(obj)
    if (nrow(td)) {
      td$estimand <- "modifier"
      lvl <- interactmeta_modifier_summary_col(obj$summary)
      if (lvl %in% names(td) && !"level" %in% names(td)) {
        td$level <- td[[lvl]]
      }
      pieces$modifiers <- td
    }
  }
  if ("combined" %in% estimands) {
    obj <- do.call(combined, c(list(fit = x), marg_args))
    td <- tidy(obj)
    if (nrow(td)) {
      td$estimand <- "combined"
      if (!"level" %in% names(td)) {
        td$level <- NA_character_
      }
      pieces$combined <- td
    }
  }
  if ("heterogeneity" %in% estimands) {
    obj <- heterogeneity(x)
    td <- tidy(obj)
    if (nrow(td)) {
      td$estimand <- "heterogeneity"
      pieces$heterogeneity <- td
    }
  }
  if (!length(pieces)) {
    return(meta[0, , drop = FALSE])
  }
  out <- dplyr::bind_rows(pieces)
  out <- dplyr::bind_cols(out, meta[rep(1L, nrow(out)), , drop = FALSE])
  tibble::as_tibble(out)
}
