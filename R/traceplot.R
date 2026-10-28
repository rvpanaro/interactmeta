# ---- traceplot.R ----
#' @noRd
check_bayesmeta <- function() {
  if (!requireNamespace("bayesmeta", quietly = TRUE)) {
    stop(
      "Package `bayesmeta` is required for forest-plot helpers and `traceplot()`. ",
      "Install it with install.packages(\"bayesmeta\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Rubin-style shrinkage trace plots for a CAMS fit
#'
#' Draws [bayesmeta::traceplot()] diagnostics from the **native**
#' `bayesmeta` / `bmr` objects stored on a Bayesian CAMS fit
#' (trial-specific shrinkage means vs heterogeneity τ, with the τ posterior
#' underneath). Available `type` values are **interaction** (the
#' interaction-block `bayesmeta` fit) and **combined** (the combined-population
#' effects represented by the overall-block `bmr` fit).
#'
#' This method does **not** fit auxiliary [bayesmeta::bayesmeta()] models.
#' Frequentist CAMS fits have no native bayesmeta objects; `traceplot()`
#' returns an informative error rather than substituting another
#' implementation. Per-subgroup traceplots are not stored by the orthogonal
#' CAMS factorization (`type = "modifiers"` or `"subgroups"` errors).
#'
#' @param x A Bayesian [`interactmeta`] from [cams()]
#'   (`engine = "bayesmeta"`).
#' @param type `"interaction"` or `"combined"`. `"modifiers"` / `"subgroups"` are rejected:
#'   Bayesian CAMS does not store
#'   separate per-subgroup meta-analyses.
#' @param ylab y-axis label passed to [bayesmeta::traceplot()]. Default is
#'   `"Interaction"` for `type = "interaction"` and `"Effect"` for
#'   `"combined"`.
#' @param xlab x-axis label used in place of bayesmeta's hardcoded
#'   `"heterogeneity τ"`. Default is `"heterogeneity τ_γ"` for
#'   interaction and `"heterogeneity τ"` for combined.
#' @param tau.prior Ignored. Retained for compatibility; the stored native
#'   objects already use the CAMS priors from the original [cams()] fit.
#' @param X Optional design row for [bayesmeta::traceplot.bmr]. Default is
#'   intercept plus mean trial information fraction from the **stored**
#'   CAMS trial blocks.
#' @param ... Further arguments passed to [bayesmeta::traceplot()]
#'   (`mulim`, `taulim`, `ci`, `prior`, `infinity`, colours, …).
#'
#' @return Invisibly, a named list of the native [bayesmeta::bayesmeta()] or
#'   [bayesmeta::bmr()] objects that were plotted.
#'
#' @examples
#' \dontrun{
#' library(interactmeta)
#' es <- escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
#' fit <- cams(es, reference = "Lower")
#' traceplot(fit, type = "interaction")
#' traceplot(fit, type = "combined", prior = TRUE)
#' }
#'
#' @seealso [forestplot()], [bayesmeta::traceplot()], [plot_heterogeneity()],
#'   [native_fit()], `vignette("cams")`
#' @export
#' @exportS3Method bayesmeta::traceplot interactmeta
traceplot.interactmeta <- function(x,
                               type = c("subgroups", "interaction", "combined"),
                               ylab = NULL,
                               xlab = NULL,
                               tau.prior = "sqrt",
                               X = NULL,
                               ...) {
  check_comparison_fit(x)
  check_bayesmeta()
  type <- comparison_resolve_effect_type(match.arg(type))
  force(tau.prior)
  labs <- comparison_traceplot_axis_labels(type, ylab = ylab, xlab = xlab)

  objs <- comparison_traceplot_objects(x, type = type)
  if (!length(objs)) {
    stop("No native bayesmeta objects available for `traceplot()`.", call. = FALSE)
  }

  n <- length(objs)
  if (n > 1L) {
    old_par <- graphics::par(mfrow = n2mfrow_traceplot(n), oma = c(0, 0, 0, 0))
    on.exit(graphics::par(old_par), add = TRUE)
  }

  for (nm in names(objs)) {
    obj <- objs[[nm]]
    args <- list(x = obj, ylab = labs$ylab, ...)
    if (inherits(obj, "bmr")) {
      args$X <- X %||% comparison_traceplot_default_bmr_X(x, obj)
    }
    if ("X" %in% names(args) && is.null(args$X)) {
      args$X <- NULL
    }
    comparison_traceplot_native(args, xlab = labs$xlab)
    if (n > 1L && nzchar(nm)) {
      graphics::title(main = nm, line = 0.6, cex.main = 0.95)
    }
  }

  invisible(objs)
}

#' @noRd
n2mfrow_traceplot <- function(n) {
  n <- max(1L, as.integer(n)[[1L]])
  if (n <= 1L) {
    return(c(1L, 1L))
  }
  if (n == 2L) {
    return(c(1L, 2L))
  }
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)
  c(nr, nc)
}

#' Type-specific traceplot axis labels
#' @noRd
comparison_traceplot_axis_labels <- function(type, ylab = NULL, xlab = NULL) {
  type <- comparison_resolve_effect_type(type)
  if (is.null(ylab)) {
    ylab <- if (identical(type, "interaction")) {
      "Interaction"
    } else {
      "Effect"
    }
  }
  ylab <- as.character(ylab)[[1L]]
  if (is.null(xlab)) {
    xlab <- if (identical(type, "interaction")) {
      heterogeneity_tau_xlab("Interaction")
    } else {
      heterogeneity_tau_xlab("Combined")
    }
  }
  list(ylab = ylab, xlab = xlab)
}

#' Run native [bayesmeta::traceplot()] with a CAMS τ / τ_γ xlab
#'
#' `bayesmeta::traceplot()` hardcodes `heterogeneity τ` via `mtext()` in
#' the density panel, then resets `layout()` before returning. Overlaying
#' a second xlab after that lands in the wrong panel. Substitute `mtext`
#' for the duration of the native call instead.
#' @noRd
comparison_traceplot_native <- function(args, xlab) {
  ns <- asNamespace("bayesmeta")
  imports <- parent.env(ns)
  if (!exists("mtext", envir = imports, inherits = FALSE)) {
    do.call(bayesmeta::traceplot, args)
    return(invisible(NULL))
  }
  orig <- get("mtext", envir = imports, inherits = FALSE)
  wrapped <- function(text, side = 3, line = 0, outer = FALSE, ...) {
    if (identical(as.integer(side)[[1L]], 1L) && !isTRUE(outer)) {
      text <- xlab
    }
    orig(text = text, side = side, line = line, outer = outer, ...)
  }
  was_locked <- bindingIsLocked("mtext", imports)
  if (was_locked) {
    unlockBinding("mtext", imports)
  }
  assign("mtext", wrapped, envir = imports)
  on.exit(
    {
      assign("mtext", orig, envir = imports)
      if (was_locked) {
        lockBinding("mtext", imports)
      }
    },
    add = TRUE
  )
  do.call(bayesmeta::traceplot, args)
}

#' Resolve native bayesmeta / bmr objects for CAMS traceplots (no refit)
#' @noRd
comparison_traceplot_objects <- function(fit, type, tau.prior = "sqrt") {
  force(tau.prior)
  type <- comparison_resolve_effect_type(type)
  if (!interactmeta_is_bayesian_fit(fit) || is.null(fit$bayesmeta)) {
    stop(
      "traceplot() uses the native bayesmeta/bmr objects stored on ",
      "Bayesian CAMS fits. This fit used method = \"",
      interactmeta_method_display(fit) %||% "unknown",
      "\". Auxiliary bayesmeta models are not fitted.",
      call. = FALSE
    )
  }
  if (identical(type, "subgroups")) {
    stop(
      "traceplot(type = \"subgroups\") is not available. The bayesmeta CAMS ",
      "engine stores an overall bmr fit and an interaction bayesmeta fit, ",
      "not separate per-subgroup meta-analyses. Use type = \"combined\" or ",
      "type = \"interaction\".",
      call. = FALSE
    )
  }
  if (identical(type, "interaction")) {
    obj <- fit$bayesmeta$interaction
    if (is.null(obj)) {
      stop("This CAMS fit has no stored native interaction bayesmeta object.", call. = FALSE)
    }
    return(list(Interaction = obj))
  }
  obj <- fit$bayesmeta$overall
  if (is.null(obj)) {
    stop("This CAMS fit has no stored native overall bmr object.", call. = FALSE)
  }
  list(Combined = obj)
}

#' Default X row for overall bmr: intercept + mean trial IF from the stored fit
#' @noRd
comparison_traceplot_default_bmr_X <- function(fit, overall_fit = NULL) {
  force(overall_fit)
  ifrac <- NULL
  blocks <- fit$bayesmeta$trial_blocks %||%
    fit$transformed_data$trial_blocks
  if (is.data.frame(blocks$overall) && "ifrac" %in% names(blocks$overall)) {
    ifrac <- as.numeric(blocks$overall$ifrac)
  }
  if (is.null(ifrac) || !length(ifrac) || all(!is.finite(ifrac))) {
    return(NULL)
  }
  pi <- mean(ifrac[is.finite(ifrac)])
  matrix(
    c(1, pi),
    nrow = 1L,
    dimnames = list(sprintf("Combined (IF = %.2f)", pi), c("intercept", "ifrac"))
  )
}
