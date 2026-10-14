#' Extract DA / AD / SWADA estimands from the selected fit only
#'
#' These constructions do not use CAMS information fractions. Extractors
#' return quantities stored on the fit; they never call [cam()] or another
#' analysis to fill missing functionality.
#' @noRd
NULL

#' Wald summary on the ratio scale: exp(point and CI), delta-method SE
#' @noRd
interactmeta_exponentiate_wald <- function(out) {
  out$se <- as.numeric(out$se) * exp(as.numeric(out$estimate))
  out$estimate <- exp(as.numeric(out$estimate))
  out$ci.lb <- exp(as.numeric(out$ci.lb))
  out$ci.ub <- exp(as.numeric(out$ci.ub))
  out
}

#' @noRd
interactmeta_construction_row_to_summary <- function(row, exponentiate = FALSE) {
  row <- as.data.frame(row, stringsAsFactors = FALSE)
  est <- as.numeric(row$estimate)[[1L]]
  se <- as.numeric(row$se)[[1L]]
  lo <- as.numeric(row$ci.lb)[[1L]]
  hi <- as.numeric(row$ci.ub)[[1L]]
  if (isTRUE(exponentiate)) {
    se <- se * exp(est)
    est <- exp(est)
    lo <- exp(lo)
    hi <- exp(hi)
  }
  data.frame(
    type = "Overall",
    estimate = est,
    mean = est,
    std.error = se,
    sd = se,
    lci = lo,
    uci = hi,
    conf.low = lo,
    conf.high = hi,
    stringsAsFactors = FALSE
  )
}

#' @noRd
interactmeta_extract_interaction <- function(fit, exponentiate = FALSE, ...) {
  interactmeta_require_capability(fit, "interaction", feature = "interactions()")
  analysis <- interactmeta_fit_analysis(fit)
  if (identical(analysis, "conventional")) {
    stop(
      "conventional() is a DA/AD bundle. Call interactions() on `$da` or `$ad`.",
      call. = FALSE
    )
  }
  sm <- interactmeta_construction_row_to_summary(fit$interaction, exponentiate = exponentiate)
  sm$interaction <- as.character(fit$interaction$label)[[1L]]
  out <- list(
    summary = sm,
    analysis = analysis,
    engine = fit$engine,
    method = fit$method,
    inference = fit$inference
  )
  class(out) <- c("interactmeta_construction_interaction", "list")
  interactmeta_attach_origin(out, fit, quantity = "interaction", fitted_or_derived = "fitted")
}

#' @noRd
interactmeta_extract_subgroups <- function(fit, exponentiate = FALSE, ...) {
  interactmeta_require_capability(fit, "subgroups", feature = "subgroups()")
  analysis <- interactmeta_fit_analysis(fit)
  if (identical(analysis, "conventional")) {
    stop(
      "conventional() is a DA/AD bundle. Call subgroups() on `$da`.",
      call. = FALSE
    )
  }
  sg <- as.data.frame(fit$subgroups, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(sg)), function(i) {
    sm <- interactmeta_construction_row_to_summary(sg[i, , drop = FALSE], exponentiate = exponentiate)
    sm$subgroup <- as.character(sg$label[[i]])
    sm
  })
  out <- list(
    summary = dplyr::bind_rows(rows),
    analysis = analysis,
    engine = fit$engine,
    method = fit$method,
    inference = fit$inference
  )
  class(out) <- c("interactmeta_construction_subgroups", "list")
  interactmeta_attach_origin(out, fit, quantity = "subgroups", fitted_or_derived = "fitted")
}

#' @noRd
interactmeta_extract_heterogeneity <- function(fit, ...) {
  interactmeta_require_capability(fit, "heterogeneity", feature = "heterogeneity()")
  analysis <- interactmeta_fit_analysis(fit)
  rows <- list()
  add_tau <- function(label, tau2, tau = NA_real_, lo = NA_real_, hi = NA_real_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      parameter = label,
      type = "tau",
      estimate = if (is.finite(tau)) tau else sqrt(max(as.numeric(tau2), 0)),
      tau2 = as.numeric(tau2),
      lci = lo,
      uci = hi,
      stringsAsFactors = FALSE
    )
  }
  if (identical(analysis, "da") || identical(analysis, "conventional")) {
    sg <- as.data.frame(fit$subgroups, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(sg))) {
      add_tau(
        paste0("tau[", sg$label[[i]], "]"),
        sg$tau2[[i]],
        tau = if ("tau" %in% names(sg)) sg$tau[[i]] else NA_real_
      )
    }
  }
  if (identical(analysis, "ad") || identical(analysis, "swada") ||
        identical(analysis, "conventional")) {
    int <- as.data.frame(fit$interaction, stringsAsFactors = FALSE)
    add_tau(
      "tau[interaction]",
      int$tau2[[1L]],
      tau = if ("tau" %in% names(int)) int$tau[[1L]] else NA_real_
    )
  }
  out <- list(
    summary = dplyr::bind_rows(rows),
    analysis = analysis,
    engine = fit$engine,
    method = fit$method,
    inference = fit$inference
  )
  if (identical(fit$engine, "bayesmeta")) {
    sm <- out$summary
    if (!is.null(fit$fits$reference)) {
      tau_ref <- interactmeta_bayesmeta_tau_summary(fit$fits$reference)
      i <- which(sm$parameter == paste0("tau[", fit$reference, "]"))
      if (length(i) == 1L) {
        sm$lci[i] <- tau_ref$lower
        sm$uci[i] <- tau_ref$upper
        sm$estimate[i] <- tau_ref$median
      }
    }
    if (!is.null(fit$fits$nonref)) {
      tau_non <- interactmeta_bayesmeta_tau_summary(fit$fits$nonref)
      i <- which(sm$parameter == paste0("tau[", fit$nonref, "]"))
      if (length(i) == 1L) {
        sm$lci[i] <- tau_non$lower
        sm$uci[i] <- tau_non$upper
        sm$estimate[i] <- tau_non$median
      }
    }
    if (!is.null(fit$fits$interaction) && inherits(fit$fits$interaction, "bayesmeta")) {
      tau_g <- interactmeta_bayesmeta_tau_summary(fit$fits$interaction)
      i <- which(sm$parameter == "tau[interaction]")
      if (length(i) == 1L) {
        sm$lci[i] <- tau_g$lower
        sm$uci[i] <- tau_g$upper
        sm$estimate[i] <- tau_g$median
      }
    }
    out$summary <- sm
  }
  class(out) <- c("interactmeta_construction_heterogeneity", "list")
  interactmeta_attach_origin(out, fit, quantity = "heterogeneity", fitted_or_derived = "fitted")
}

#' @export
print.interactmeta_construction_interaction <- function(x, digits = 2, ...) {
  cat(
    interactmeta_analysis_acronym(x$analysis), " interaction (", x$engine, " / ", x$method, ")\n",
    sep = ""
  )
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
print.interactmeta_construction_subgroups <- function(x, digits = 2, ...) {
  cat(
    interactmeta_analysis_acronym(x$analysis), " ",
    interactmeta_modifier_heading(x$modifier %||% "subgroup"),
    " effects (", x$engine, " / ", x$method, ")\n",
    sep = ""
  )
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
print.interactmeta_construction_heterogeneity <- function(x, digits = 2, ...) {
  cat(
    interactmeta_analysis_acronym(x$analysis), " heterogeneity (", x$engine, " / ", x$method, ")\n",
    sep = ""
  )
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
#' @importFrom generics tidy
tidy.interactmeta_construction_interaction <- function(x, ...) {
  tibble::as_tibble(x$summary)
}

#' @export
#' @importFrom generics tidy
tidy.interactmeta_construction_subgroups <- function(x, ...) {
  tibble::as_tibble(x$summary)
}

#' @export
#' @importFrom generics tidy
tidy.interactmeta_construction_heterogeneity <- function(x, ...) {
  tibble::as_tibble(x$summary)
}
