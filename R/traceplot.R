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
#' (study-specific shrinkage means vs heterogeneity τ, with the τ posterior
#' underneath). Available `type` values are **interaction** (the
#' interaction-block `bayesmeta` fit) and **combined** (the combined-population
#' effects represented by the overall-block `bmr` fit).
#'
#' This method does **not** fit auxiliary [bayesmeta::bayesmeta()] models.
#' Frequentist CAMS fits have no native bayesmeta objects; `traceplot()`
#' returns an informative error rather than substituting another
#' implementation. Per-subgroup traceplots are not stored by the orthogonal
#' CAMS factorization (`type = "subgroups"` errors).
#'
#' @param x A Bayesian [`interactmeta`] from [cam()]
#'   (`engine = "bayesmeta"`).
#' @param type `"interaction"` or `"combined"`. `"subgroups"` is rejected:
#'   Bayesian CAMS does not store
#'   separate per-subgroup meta-analyses.
#' @param ylab y-axis label passed to [bayesmeta::traceplot()]. Default is
#'   `"Interaction"` for `type = "interaction"` and `"Effect"` for
#'   `"combined"`.
#' @param xlab x-axis label used in place of bayesmeta's hardcoded
#'   `"heterogeneity τ"`. Default is `"heterogeneity τ_γ"` for
#'   interaction and `"heterogeneity τ"` for combined.
#' @param tau.prior Ignored. Retained for compatibility; the stored native
#'   objects already use the CAMS priors from the original [cam()] fit.
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
#' library(interact)
#' es <- escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
#' fit <- cam(es, reference = "Younger")
#' traceplot(fit, type = "interaction")
#' traceplot(fit, type = "combined", prior = TRUE)
#' }
#'
#' @seealso [forestplot()], [bayesmeta::traceplot()], [plot_heterogeneity()],
#'   [native_fit()], `vignette("introduction-to-interact")`
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
  eng <- as.character(fit$engine %||% "")[[1L]]
  if (!identical(eng, "bayesmeta") || is.null(fit$bayesmeta)) {
    stop(
      "traceplot() uses the native bayesmeta/bmr objects stored on ",
      "engine = \"bayesmeta\" CAMS fits. This fit used engine = \"",
      if (nzchar(eng)) eng else "unknown",
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

#' Study-level component rows from a CAMS fit
#' @noRd
comparison_traceplot_study_data <- function(fit) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  dat <- if (!is.null(st$ifrac_data)) {
    as.data.frame(st$ifrac_data)
  } else {
    as.data.frame(fit$data)
  }
  sg_col <- comparison_traceplot_subgroup_var(fit)
  if (!sg_col %in% names(dat)) {
    stop("Fit data are missing the subgroup column `", sg_col, "`.", call. = FALSE)
  }
  keep <- is_ifrac_component_row(dat[[sg_col]])
  out <- dat[keep, , drop = FALSE]
  if (!nrow(out)) {
    stop("No subgroup study rows found in the CAMS fit data.", call. = FALSE)
  }
  need <- c("yi", "vi")
  missing <- setdiff(need, names(out))
  if (length(missing)) {
    stop(
      "Fit data are missing required columns: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  out
}

#' @noRd
comparison_traceplot_subgroup_var <- function(fit) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  as.character(st$subgroup_var %||% "subgroup")[[1L]]
}

#' @noRd
comparison_traceplot_trial_var <- function(fit) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  as.character(st$trial_var %||% "trial")[[1L]]
}

#' Trial-level overall / interaction blocks rebuilt from fit data
#' @noRd
comparison_traceplot_trial_blocks <- function(fit) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st)) {
    stop("Fit is missing the `\"interactmeta\"` attribute; refit with `cam()`.", call. = FALSE)
  }
  dat <- comparison_traceplot_study_data(fit)
  trial_col <- comparison_traceplot_trial_var(fit)
  sg_col <- comparison_traceplot_subgroup_var(fit)
  ref <- as.character(st$reference %||% "")[[1L]]
  if (!nzchar(ref)) {
    stop("Fit has no reference subgroup label.", call. = FALSE)
  }
  nonref <- st$coef_map$ifrac_orig %||% st$ifrac_vars
  if (!length(nonref)) {
    stop("Fit has no non-reference `ifrac_*` column.", call. = FALSE)
  }
  nonref_lab <- sub("^ifrac_", "", as.character(nonref)[[1L]])
  ifrac_col <- as.character(nonref)[[1L]]
  if (!ifrac_col %in% names(dat)) {
    # Fall back to any ifrac_* column.
    ifrac_cols <- grep("^ifrac_", names(dat), value = TRUE)
    ifrac_cols <- ifrac_cols[!grepl("^ifrac2_", ifrac_cols)]
    if (!length(ifrac_cols)) {
      stop("Fit data have no `ifrac_*` columns for overall IF.", call. = FALSE)
    }
    ifrac_col <- ifrac_cols[[1L]]
  }

  blocks <- comparison_trial_effect_blocks(
    data = dat,
    trial_var = trial_col,
    subgroup_var = sg_col,
    reference = ref,
    nonref = nonref_lab,
    ifrac_var = ifrac_col
  )
  if (!nrow(blocks$interaction)) {
    stop(
      "Could not build interaction contrasts for `traceplot()` ",
      "(need one row each for reference and non-reference per trial).",
      call. = FALSE
    )
  }
  list(
    overall = blocks$overall,
    interaction = blocks$interaction,
    reference = ref,
    nonref = nonref_lab,
    ifrac_var = ifrac_col
  )
}
