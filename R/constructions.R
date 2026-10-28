# ---- conventional.R ----
#' Pair binary-subgroup rows into trial-level estimates
#'
#' Returns one row per trial with reference / non-reference `yi`/`vi` and the
#' CAMS interaction contrast `g = y_nonref - y_ref`. Interaction variances use
#' `v_ref + v_nonref` (within-trial subgroup covariance assumed zero).
#' @noRd
interactmeta_binary_trial_table <- function(data,
                                            trial_var = "study",
                                            subgroup_var = "level",
                                            reference,
                                            yi_var = "yi",
                                            vi_var = "vi") {
  dat <- as.data.frame(data)
  ref <- as.character(reference)[[1L]]
  sg <- as.character(dat[[subgroup_var]])
  levels_sg <- unique(sg)
  nonrefs <- setdiff(levels_sg, ref)
  if (length(nonrefs) != 1L) {
    stop(
      "DA, AD, and SWADA currently require exactly two subgroup levels ",
      "(one reference). Found: ",
      paste(levels_sg, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  nonref <- nonrefs[[1L]]
  trials <- unique(as.character(dat[[trial_var]]))
  rows <- lapply(trials, function(tr) {
    d <- dat[as.character(dat[[trial_var]]) == tr, , drop = FALSE]
    i_ref <- which(as.character(d[[subgroup_var]]) == ref)
    i_non <- which(as.character(d[[subgroup_var]]) == nonref)
    if (length(i_ref) != 1L || length(i_non) != 1L) {
      stop(
        "Trial `",
        tr,
        "` must have exactly one row for `",
        ref,
        "` and `",
        nonref,
        "`.",
        call. = FALSE
      )
    }
    yA <- as.numeric(d[[yi_var]][[i_ref]])
    yB <- as.numeric(d[[yi_var]][[i_non]])
    vA <- as.numeric(d[[vi_var]][[i_ref]])
    vB <- as.numeric(d[[vi_var]][[i_non]])
    nA <- if ("n" %in% names(d)) as.numeric(d$n[[i_ref]]) else NA_real_
    nB <- if ("n" %in% names(d)) as.numeric(d$n[[i_non]]) else NA_real_
    data.frame(
      trial = tr,
      y_ref = yA,
      v_ref = vA,
      n_ref = nA,
      y_nonref = yB,
      v_nonref = vB,
      n_nonref = nB,
      g = yB - yA,
      vg = vA + vB,
      pi = (1 / vB) / (1 / vA + 1 / vB),
      m = (vB * yA + vA * yB) / (vA + vB),
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  attr(out, "reference") <- ref
  attr(out, "nonref") <- nonref
  out
}

interactmeta_rma_summary <- function(fit, label) {
  z <- stats::qnorm(0.975)
  est <- as.numeric(fit$beta)[[1L]]
  se <- as.numeric(fit$se)[[1L]]
  data.frame(
    label = label,
    estimate = est,
    se = se,
    ci.lb = est - z * se,
    ci.ub = est + z * se,
    tau2 = as.numeric(fit$tau2)[[1L]],
    method = fit$method,
    stringsAsFactors = FALSE
  )
}

interactmeta_da_from_pairs <- function(pairs, method, rma_args, call) {
  ref <- attr(pairs, "reference")
  nonref <- attr(pairs, "nonref")
  fit_ref <- do.call(
    metafor::rma,
    c(
      list(
        yi = pairs$y_ref,
        vi = pairs$v_ref,
        slab = pairs$trial,
        method = method
      ),
      rma_args
    )
  )
  fit_non <- do.call(
    metafor::rma,
    c(
      list(
        yi = pairs$y_nonref,
        vi = pairs$v_nonref,
        slab = pairs$trial,
        method = method
      ),
      rma_args
    )
  )
  sm_ref <- interactmeta_rma_summary(fit_ref, ref)
  sm_non <- interactmeta_rma_summary(fit_non, nonref)
  z <- stats::qnorm(0.975)
  da <- sm_non$estimate - sm_ref$estimate
  da_se <- sqrt(sm_non$se^2 + sm_ref$se^2)
  list(
    reference = ref,
    nonref = nonref,
    pairs = pairs,
    subgroups = rbind(sm_ref, sm_non),
    interaction = data.frame(
      label = paste0(nonref, " - ", ref, " (DA)"),
      estimate = da,
      se = da_se,
      ci.lb = da - z * da_se,
      ci.ub = da + z * da_se,
      tau2 = NA_real_,
      method = method,
      stringsAsFactors = FALSE
    ),
    weights = data.frame(
      trial = pairs$trial,
      weight_ref = as.numeric(stats::weights(fit_ref)),
      weight_nonref = as.numeric(stats::weights(fit_non)),
      stringsAsFactors = FALSE
    ),
    fits = list(reference = fit_ref, nonref = fit_non),
    method = method,
    call = call
  )
}

interactmeta_ad_from_pairs <- function(pairs, method, rma_args, call) {
  ref <- attr(pairs, "reference")
  nonref <- attr(pairs, "nonref")
  fit_g <- do.call(
    metafor::rma,
    c(
      list(
        yi = pairs$g,
        vi = pairs$vg,
        slab = pairs$trial,
        method = method
      ),
      rma_args
    )
  )
  list(
    reference = ref,
    nonref = nonref,
    pairs = pairs,
    interaction = interactmeta_rma_summary(fit_g, paste0(nonref, " - ", ref)),
    weights = data.frame(
      trial = pairs$trial,
      weight = as.numeric(stats::weights(fit_g)),
      stringsAsFactors = FALSE
    ),
    fits = list(interaction = fit_g),
    method = method,
    call = call
  )
}

#' Estimate a treatment-by-subgroup interaction using the Difference-of-Averages construction
#'
#' **Difference of Averages (DA)** synthesizes each subgroup series separately
#' and takes the difference of those pooled treatment effects. Study
#' inverse-variance weights generally differ across the two subgroup
#' analyses, so DA need not equal [ad()]. Panaro, Röver and Friede (2026)
#' treat DA and AD as **different constructions**.
#'
#' For binary subgroups with reference \eqn{A} and non-reference \eqn{B},
#' \eqn{\mu_A} and \eqn{\mu_B} are independent random-effects meta-analyses
#' of the subgroup series, and
#' \eqn{\gamma_{\mathrm{DA}} = \mu_B - \mu_A}.
#'
#' DA is an **analysis construction**, not an inference framework. It may be
#' fitted with frequentist inference (`engine = "metafor"`, REML/ML) or
#' Bayesian inference (`engine = "bayesmeta"`). Bayesian DA fits two independent
#' univariate [bayesmeta::bayesmeta()] models and forms \eqn{\mu_B - \mu_A} by
#' Cartesian-product Gaussian-mixture (and optional convolution). It does
#' **not** use quantile matching, and it does **not** acquire CAMS
#' information-fraction functionality merely because a posterior exists.
#'
#' Fisher et al. (2017) discussed a related across-trial contrast of
#' subgroup averages as the “deluded” approach. That paper does not define
#' the DA/AD pair used here. DA is **not** the joint bivariate mixed model
#' of van Houwelingen, Arends and Stijnen (2002). That joint model remains a
#' related reference for modelling two subgroup series together with a
#' heterogeneity matrix \eqn{\Sigma}; [da()] uses separate univariate
#' pooling instead.
#'
#' [ad()] meta-analyses the within-trial differences. [swada()] imposes the
#' same study weights on both constructions. [cams()] is a different object:
#' the CAMS statistical model, which is the analysis family that uses the
#' information-fraction construction.
#'
#' @param data A data frame with `trial`, `subgroup`, `yi`, and `vi`.
#' @param reference Reference subgroup label.
#' @param trial,subgroup,yi,vi Column names.
#' @param engine Computational backend. `"metafor"` (default) or `"bayesmeta"`.
#'   CAMS backends are not DA backends. `"bayes"` is accepted as a synonym of
#'   `"bayesmeta"`.
#' @param approach Inferential approach. `"bayes"` uses `engine = "bayesmeta"`
#'   when `engine` is omitted; `"frequentist"`, `"REML"`, or `"ML"` use
#'   `engine = "metafor"`. `NULL` (default) follows `engine` / `method`.
#' @param method Inferential method. `NULL` (default) is `"REML"` for metafor
#'   and `"bayes"` for bayesmeta. metafor also accepts `"ML"`. If `engine` is
#'   omitted, `"bayes"` uses bayesmeta and `"REML"` / `"ML"` use metafor.
#' @param tau.prior,mu.prior Bayesian priors when `engine = "bayesmeta"`
#'   (defaults: improper uniform on \eqn{\tau}, matching bayesmeta;
#'   \eqn{N(0,10)} on \eqn{\mu}. AD uses the same tau default as CAMS
#'   interaction). Ignored by metafor.
#' @param ... Computational options passed to [metafor::rma()] or
#'   [bayesmeta::bayesmeta()] (including `engine_args`).
#' @return An object of class `"interactmeta_da"`.
#' @references
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \url{https://doi.org/10.48550/arXiv.2508.15531}
#'
#' Fisher DJ, Carpenter JR, Morris TP, Freeman SC, Tierney JF (2017).
#' Meta-analytical methods to identify who benefits most from treatments:
#' daft, deluded, or deft approach? \emph{BMJ}, 356, j573.
#' \url{https://doi.org/10.1136/bmj.j573}
#'
#' van Houwelingen HC, Arends LR, Stijnen T (2002).
#' Advanced methods in meta-analysis: multivariate approach and
#' meta-regression. \emph{Statistics in Medicine}, 21(4), 589--624.
#' \url{https://doi.org/10.1002/sim.1040}
#'
#' van Houwelingen HC, Zwinderman KH, Stijnen T (1993).
#' A bivariate approach to meta-analysis. \emph{Statistics in Medicine},
#' 12(24), 2273--2284. \url{https://doi.org/10.1002/sim.4780122405}
#' @param modifier Column of the two levels. Overrides `subgroup` when set.
#'   The default column remains `"subgroup"`.
#' @seealso [ad()], [swada()], [cams()], [conventional()], [interactmeta_capability_matrix()]
#' @export
da <- function(data,
               reference,
               modifier = NULL,
               trial = "study",
               subgroup = "level",
               yi = "yi",
               vi = "vi",
               method = "Bayesian",
               tau.prior = NULL,
               mu.prior = NULL,
               ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup, dat = data)
  interactmeta_fit_da_ad(
    analysis = "da",
    data = data,
    reference = reference,
    trial = trial,
    subgroup = subgroup,
    yi = yi,
    vi = vi,
    method = method,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    call = match.call(),
    ...
  )
}

#' Estimate a treatment-by-subgroup interaction using the Average-Difference construction
#'
#' **Average Difference (AD)** first forms the within-trial
#' treatment-by-subgroup difference and then meta-analyzes those
#' differences. For binary subgroups with reference \eqn{A} and non-reference
#' \eqn{B},
#' \eqn{g_j = y_{Bj} - y_{Aj}} with
#' \eqn{\sigma_{gj}^2 = \sigma_{Aj}^2 + \sigma_{Bj}^2} (within-trial subgroup
#' covariance assumed zero), and \eqn{\gamma_{\mathrm{AD}}} is a random-effects
#' meta-analysis of the \eqn{g_j}.
#'
#' AD is a **different construction** from [da()] (Panaro, Röver & Friede
#' 2026), not an inference framework. It may be fitted with
#' `engine = "metafor"` (REML/ML) or `engine = "bayesmeta"` (Bayesian RE-MA
#' of \eqn{g_j}). The CAMS manuscript’s BIM (Bayesian interaction
#' meta-analysis) is the Bayesian analogue of this construction. A Bayesian
#' AD fit is internally complete: it does not borrow uncertainty or
#' heterogeneity from Bayesian DA or CAMS. It still does **not** use CAMS
#' information fractions.
#'
#' Study weights are those of the interaction series,
#' \eqn{w_j \propto 1/(\sigma_{gj}^2 + \tau_g^2)}. They generally differ from
#' the subgroup-specific weights used by [da()]. [swada()] replaces those
#' analysis-specific weights with one common vector. [cams()] is a
#' different object: the CAMS statistical model.
#'
#' @inheritParams da
#' @return An object of class `"interactmeta_ad"`.
#' @references
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \url{https://doi.org/10.48550/arXiv.2508.15531}
#'
#' Fisher DJ, Carpenter JR, Morris TP, Freeman SC, Tierney JF (2017).
#' Meta-analytical methods to identify who benefits most from treatments:
#' daft, deluded, or deft approach? \emph{BMJ}, 356, j573.
#' \url{https://doi.org/10.1136/bmj.j573}
#'
#' Panaro RV, Röver C, Friede T (2025).
#' Consistent Bayesian meta-analysis on subgroup specific effects and
#' interactions. arXiv:2512.18785. \url{https://doi.org/10.48550/arXiv.2512.18785}
#' @param modifier Column of the two levels. Overrides `subgroup` when set.
#'   The default column remains `"subgroup"`.
#' @seealso [da()], [swada()], [cams()], [conventional()], [interactmeta_capability_matrix()]
#' @export
ad <- function(data,
               reference,
               modifier = NULL,
               trial = "study",
               subgroup = "level",
               yi = "yi",
               vi = "vi",
               method = "Bayesian",
               tau.prior = NULL,
               mu.prior = NULL,
               ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup, dat = data)
  interactmeta_fit_da_ad(
    analysis = "ad",
    data = data,
    reference = reference,
    trial = trial,
    subgroup = subgroup,
    yi = yi,
    vi = vi,
    method = method,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    call = match.call(),
    ...
  )
}

#' @noRd
interactmeta_fit_da_ad <- function(analysis,
                                   data,
                                   reference,
                                   trial,
                                   subgroup,
                                   yi,
                                   vi,
                                   method,
                                   tau.prior,
                                   mu.prior,
                                   call,
                                   ...) {
  resolved <- interactmeta_resolve_construction_method(analysis, method)
  if (identical(resolved$backend, "metafor") &&
      (!is.null(tau.prior) || !is.null(mu.prior))) {
    .interactmeta_warn_unused_priors(TRUE)
  }
  dat <- as.data.frame(data)
  pairs <- interactmeta_binary_trial_table(
    dat,
    trial_var = trial,
    subgroup_var = subgroup,
    reference = reference,
    yi_var = yi,
    vi_var = vi
  )
  dots <- list(...)
  engine_args <- dots$engine_args
  if (is.null(engine_args)) {
    engine_args <- list()
  }
  dots$engine_args <- NULL
  if (identical(resolved$backend, "bayesmeta")) {
    tau.prior <- interactmeta_construction_tau_prior(analysis, tau.prior)
    mu.prior <- interactmeta_construction_mu_prior(mu.prior)
    if ("convolve_derived" %in% names(engine_args) || "convolve_derived" %in% names(dots)) {
      stop(
        "`convolve_derived` cannot change reported intervals. Independent sums use bayesmeta::convolve().",
        call. = FALSE
      )
    }
    engine_args <- interactmeta_merge_engine_args(engine_args, dots)
    if (identical(analysis, "da")) {
      built <- interactmeta_da_from_pairs_bayesmeta(
        pairs,
        tau.prior = tau.prior,
        mu.prior = mu.prior,
        engine_args = engine_args,
        call = call
      )
    } else {
      built <- interactmeta_ad_from_pairs_bayesmeta(
        pairs,
        tau.prior = tau.prior,
        mu.prior = mu.prior,
        engine_args = engine_args,
        call = call
      )
    }
  } else {
    dots[intersect(names(dots), interactmeta_computational_dot_names())] <- NULL
    if (identical(analysis, "da")) {
      built <- interactmeta_da_from_pairs(
        pairs,
        method = resolved$method_key,
        rma_args = dots,
        call = call
      )
    } else {
      built <- interactmeta_ad_from_pairs(
        pairs,
        method = resolved$method_key,
        rma_args = dots,
        call = call
      )
    }
  }
  notes <- if (identical(analysis, "da")) {
    "Difference of Averages (DA): independent RE-MA per subgroup series; interaction is the difference of the two pooled means. Not a CAMS information-fraction analysis."
  } else {
    "Average Difference (AD): RE-MA of within-trial contrasts g_j = y_B - y_A. Not a CAMS information-fraction analysis."
  }
  out <- list(
    approach = analysis,
    call = built$call,
    reference = built$reference,
    nonref = built$nonref,
    pairs = built$pairs,
    subgroups = built$subgroups,
    interaction = built$interaction,
    weights = built$weights,
    fits = built$fits,
    posterior = built$posterior,
    priors = built$priors,
    data = dat,
    modifier = subgroup,
    assumptions = list(
      subgroup_covariance = "assumed_zero",
      notes = notes
    )
  )
  if (identical(analysis, "ad")) {
    out$subgroups <- NULL
  }
  class(out) <- c(paste0("interactmeta_", analysis), "list")
  interactmeta_attach_capabilities(out, analysis, resolved$method)
}

#' Bundle of DA and AD with analysis-specific weights
#'
#' Runs [da()] and [ad()] on the same binary subgroup data. These are
#' **different constructions** (Panaro, Röver & Friede 2026). Inverse-variance
#' weights generally **differ** across the two subgroup series and the
#' within-trial interaction, so the DA interaction need not equal AD.
#'
#' This is a convenience bundle, not a fourth statistical model. Prefer
#' [da()] and [ad()] when you want the constructions separately. [swada()]
#' applies one common weight vector. [cams()] fits the CAMS model.
#'
#' The function is **not** the joint bivariate mixed model of van Houwelingen,
#' Arends and Stijnen (2002). That paper (and van Houwelingen, Zwinderman and
#' Stijnen 1993) remains the reference for modelling \eqn{(y_{Aj}, y_{Bj})}
#' together with a between-trial heterogeneity matrix \eqn{\Sigma}.
#'
#' @inheritParams da
#' @return An object of class `"interactmeta_conventional"` containing both
#'   DA and AD results from the selected inferential backend.
#' @references
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \url{https://doi.org/10.48550/arXiv.2508.15531}
#' @param modifier Column of the two levels. Overrides `subgroup` when set.
#'   The default column remains `"subgroup"`.
#' @seealso [da()], [ad()], [swada()], [cams()]
#' @export
conventional <- function(data,
                         reference,
                         modifier = NULL,
                         trial = "study",
                         subgroup = "level",
                         yi = "yi",
                         vi = "vi",
                         engine = c("metafor", "bayesmeta"),
                         method = NULL,
                         approach = NULL,
                         tau.prior = NULL,
                         mu.prior = NULL,
                         ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup, dat = data)
  da_fit <- da(
    data,
    reference = reference,
    trial = trial,
    subgroup = subgroup,
    yi = yi,
    vi = vi,
    engine = engine,
    method = method,
    approach = approach,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    ...
  )
  ad_fit <- ad(
    data,
    reference = reference,
    trial = trial,
    subgroup = subgroup,
    yi = yi,
    vi = vi,
    engine = engine,
    method = method,
    approach = approach,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    ...
  )
  out <- list(
    approach = "conventional",
    analysis = "conventional",
    engine = da_fit$engine,
    method = da_fit$method,
    algorithm = da_fit$algorithm,
    inference = da_fit$inference,
    call = match.call(),
    reference = da_fit$reference,
    nonref = da_fit$nonref,
    pairs = da_fit$pairs,
    subgroups = da_fit$subgroups,
    interaction = ad_fit$interaction,
    interaction_da = da_fit$interaction,
    weights = data.frame(
      trial = da_fit$weights$trial,
      weight_ref = da_fit$weights$weight_ref,
      weight_nonref = da_fit$weights$weight_nonref,
      weight_interaction = ad_fit$weights$weight,
      stringsAsFactors = FALSE
    ),
    fits = list(
      reference = da_fit$fits$reference,
      nonref = da_fit$fits$nonref,
      interaction = ad_fit$fits$interaction
    ),
    data = da_fit$data,
    da = da_fit,
    ad = ad_fit,
    assumptions = list(
      subgroup_covariance = "assumed_zero",
      notes = "Bundle of Difference of Averages (DA) and Average Difference (AD) with analysis-specific weights. Not the 2002 joint bivariate mixed model. Not a CAMS information-fraction analysis."
    ),
    capabilities = list(
      analysis = "conventional",
      engine = da_fit$engine,
      method = da_fit$method,
      inference = da_fit$inference,
      heterogeneity = TRUE,
      posterior = isTRUE(da_fit$capabilities$posterior),
      ifrac = FALSE,
      ifrac_plots = FALSE,
      interaction = TRUE,
      subgroups = TRUE
    )
  )
  class(out) <- c("interactmeta_conventional", "list")
  out
}

interactmeta_print_effect_rows <- function(df, digits) {
  keep <- intersect(c("label", "subgroup", "estimate", "se", "ci.lb", "ci.ub"), names(df))
  out <- as.data.frame(df[, keep, drop = FALSE], stringsAsFactors = FALSE)
  num <- vapply(out, is.numeric, logical(1L))
  out[num] <- lapply(out[num], function(z) round(z, digits))
  print(out, row.names = FALSE)
  invisible(out)
}

#' Heterogeneity block for DA / AD / SWADA console output
#' @noRd
interactmeta_print_fit_heterogeneity <- function(x, digits = 2, prob = NULL) {
  if (!isTRUE(interactmeta_fit_capabilities(x)$heterogeneity)) {
    return(invisible(NULL))
  }
  opts <- comparison_print_digits_prob(x, digits = digits, prob = prob)
  het <- tryCatch(
    interactmeta_extract_heterogeneity(x, prob = opts$prob),
    error = function(e) NULL
  )
  if (is.null(het) || is.null(het$summary) || !nrow(het$summary)) {
    return(invisible(NULL))
  }
  cat("\nBetween-trial heterogeneity:\n", sep = "")
  comparison_print_effects_table(
    het$summary,
    digits = opts$digits,
    prob = opts$prob,
    id_cols = c("parameter"),
    title = NULL,
    interval_scope = "single",
    style = interactmeta_interval_style(x)
  )
  invisible(het)
}

#' @export
print.interactmeta_da <- function(x, digits = 2, ...) {
  cat("Difference of Averages (DA)\n")
  cat("  analysis  : DA\n")
  cat("  method    : ", x$method, "\n", sep = "")
  cat("  inference : ", x$inference %||% NA_character_, "\n", sep = "")
  cat("  reference : ", x$reference, "\n", sep = "")
  cat("\nWeighted averages (log scale):\n")
  interactmeta_print_effect_rows(x$subgroups, digits)
  cat("\nDifference of (weighted) averages:\n")
  interactmeta_print_effect_rows(x$interaction, digits)
  interactmeta_print_fit_heterogeneity(x, digits = digits, ...)
  invisible(x)
}

#' @export
print.interactmeta_ad <- function(x, digits = 2, ...) {
  cat("Average Difference (AD)\n")
  cat("  analysis  : AD\n")
  cat("  method    : ", x$method, "\n", sep = "")
  cat("  inference : ", x$inference %||% NA_character_, "\n", sep = "")
  cat("  reference : ", x$reference, "\n", sep = "")
  cat("\nWeighted average difference:\n")
  interactmeta_print_effect_rows(x$interaction, digits)
  interactmeta_print_fit_heterogeneity(x, digits = digits, ...)
  invisible(x)
}

#' @export
print.interactmeta_conventional <- function(x, digits = 2, ...) {
  cat("DA and AD with analysis-specific weights\n")
  cat("  method    : ", x$method, "\n", sep = "")
  cat("  inference : ", x$inference %||% NA_character_, "\n", sep = "")
  cat("  reference : ", x$reference, "\n", sep = "")
  cat("\nSubgroup pooled effects (DA; log scale):\n")
  interactmeta_print_effect_rows(x$subgroups, digits)
  cat("\nAverage Difference (AD; within-trial contrasts):\n")
  interactmeta_print_effect_rows(x$interaction, digits)
  cat("\nDifference of Averages (DA; not collapsible with AD):\n")
  interactmeta_print_effect_rows(x$interaction_da, digits)
  interactmeta_print_fit_heterogeneity(x, digits = digits, ...)
  invisible(x)
}

#' @export
#' @importFrom generics tidy
tidy.interactmeta_da <- function(x, exponentiate = FALSE, ...) {
  sg <- x$subgroups
  names(sg)[names(sg) == "label"] <- "subgroup"
  sg$component <- "subgroup"
  int <- x$interaction
  names(int)[names(int) == "label"] <- "subgroup"
  int$component <- "interaction"
  out <- dplyr::bind_rows(sg, int)
  if (isTRUE(exponentiate)) {
    out <- interactmeta_exponentiate_wald(out)
  }
  tibble::as_tibble(out)
}

#' @export
#' @importFrom generics tidy
tidy.interactmeta_ad <- function(x, exponentiate = FALSE, ...) {
  int <- x$interaction
  names(int)[names(int) == "label"] <- "subgroup"
  int$component <- "interaction"
  if (isTRUE(exponentiate)) {
    int <- interactmeta_exponentiate_wald(int)
  }
  tibble::as_tibble(int)
}

#' @export
#' @importFrom generics tidy
tidy.interactmeta_conventional <- function(x, exponentiate = FALSE, ...) {
  sg <- x$subgroups
  names(sg)[names(sg) == "label"] <- "subgroup"
  sg$component <- "subgroup"
  ad_row <- x$interaction
  names(ad_row)[names(ad_row) == "label"] <- "subgroup"
  ad_row$component <- "interaction"
  ad_row$construction <- "AD"
  da_row <- x$interaction_da
  names(da_row)[names(da_row) == "label"] <- "subgroup"
  da_row$component <- "interaction"
  da_row$construction <- "DA"
  out <- dplyr::bind_rows(sg, ad_row, da_row)
  if (isTRUE(exponentiate)) {
    out <- interactmeta_exponentiate_wald(out)
  }
  tibble::as_tibble(out)
}

#' @name ifrac
#' @rdname ifrac
#' @export
ifrac.interactmeta_da <- function(x, ...) {
  interactmeta_abort_non_cams(x, "ifrac()")
}

#' @rdname ifrac
#' @export
ifrac.interactmeta_ad <- ifrac.interactmeta_da

#' @rdname ifrac
#' @export
ifrac.interactmeta_conventional <- ifrac.interactmeta_da

#' @export
plot.interactmeta_da <- function(x, y, ...) {
  interactmeta_plot_construction(x, ...)
}

#' @export
plot.interactmeta_ad <- plot.interactmeta_da

#' @export
plot.interactmeta_conventional <- plot.interactmeta_da

#' Non-IF forest of stored DA / AD / SWADA estimates
#' @noRd
interactmeta_plot_construction <- function(x, ...) {
  dots <- list(...)
  if (any(c("ifrac", "ifrac_var", "overall_ifrac", "type") %in% names(dots))) {
    interactmeta_require_capability(x, "ifrac_plots", feature = "This visualization")
  }
  td <- tidy(x)
  ycol <- if ("subgroup" %in% names(td)) "subgroup" else "label"
  td$.row <- as.character(td[[ycol]])
  if ("component" %in% names(td)) {
    td$.row <- paste(td$component, td$.row, sep = ": ")
  }
  xlab <- "Effect (log scale)"
  ggplot2::ggplot(td, ggplot2::aes(x = .data$estimate, y = .data$.row)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9A9A9A", linewidth = 0.3) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$ci.lb, xmax = .data$ci.ub), linewidth = 0.45) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::labs(x = xlab, y = NULL, title = interactmeta_analysis_display(x)) +
    interactmeta_theme_report()
}


# ---- construction_extractors.R ----
#' Extract DA / AD / SWADA estimands from the selected fit only
#'
#' These constructions do not use CAMS information fractions. Extractors
#' return quantities stored on the fit; they never call [cams()] or another
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
interactmeta_construction_row_to_summary <- function(row,
                                                     exponentiate = FALSE,
                                                     type = "Overall") {
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
    type = type,
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
interactmeta_swada_prediction_summary_row <- function(pr, type, exponentiate) {
  row <- data.frame(
    label = type,
    estimate = pr$estimate,
    se = pr$se,
    ci.lb = pr$ci.lb,
    ci.ub = pr$ci.ub,
    stringsAsFactors = FALSE
  )
  interactmeta_construction_row_to_summary(
    row,
    exponentiate = exponentiate,
    type = type
  )
}

#' @noRd
interactmeta_extract_interaction <- function(fit,
                                             type = "Overall",
                                             prob = 0.95,
                                             exponentiate = FALSE,
                                             ...) {
  interactmeta_require_capability(fit, "interaction", feature = "interactions()")
  analysis <- interactmeta_fit_analysis(fit)
  if (identical(analysis, "conventional")) {
    stop(
      "conventional() is a DA/AD bundle. Call interactions() on `$da` or `$ad`.",
      call. = FALSE
    )
  }
  type <- interactmeta_construction_marginal_type(type, analysis = analysis)
  int_lab <- as.character(fit$interaction$label)[[1L]]
  rows <- lapply(type, function(ty) {
    if (identical(ty, "Overall")) {
      sm <- interactmeta_construction_row_to_summary(
        fit$interaction,
        exponentiate = exponentiate,
        type = "Overall"
      )
    } else if (identical(analysis, "swada")) {
      pr <- interactmeta_swada_prediction(
        fit,
        target = "interaction",
        prob = prob
      )
      sm <- interactmeta_swada_prediction_summary_row(
        pr,
        type = "Prediction",
        exponentiate = exponentiate
      )
    } else {
      stop("Internal: unsupported construction Prediction.", call. = FALSE)
    }
    sm$interaction <- int_lab
    sm
  })
  out <- list(
    summary = dplyr::bind_rows(rows),
    analysis = analysis,
    engine = fit$engine,
    method = fit$method,
    inference = fit$inference,
    prob = prob
  )
  class(out) <- c("interactmeta_construction_interaction", "list")
  interactmeta_attach_origin(out, fit, quantity = "interaction", fitted_or_derived = "fitted")
}

#' @noRd
interactmeta_extract_subgroups <- function(fit, exponentiate = FALSE, ...) {
  interactmeta_require_capability(fit, "subgroups", feature = "modifiers()")
  analysis <- interactmeta_fit_analysis(fit)
  if (identical(analysis, "conventional")) {
    stop(
      "conventional() is a DA/AD bundle. Call modifiers() on `$da`.",
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

#' Profile / likelihood CIs on tau for DA / AD / SWADA metafor fits
#' @noRd
interactmeta_construction_metafor_fill_heterogeneity <- function(fit, summary, prob = 0.95) {
  if (interactmeta_is_bayesian_fit(fit) || is.null(summary) || !nrow(summary)) {
    return(summary)
  }
  fits <- fit$fits
  if (is.null(fits) || !length(fits)) {
    return(summary)
  }
  ref <- as.character(fit$reference)[[1L]]
  nonref <- as.character(fit$nonref)[[1L]]
  map <- list()
  if (!is.null(fits$reference) && inherits(fits$reference, "rma")) {
    map[[paste0("tau[", ref, "]")]] <- fits$reference
  }
  if (!is.null(fits$nonref) && inherits(fits$nonref, "rma")) {
    map[[paste0("tau[", nonref, "]")]] <- fits$nonref
  }
  if (!is.null(fits$interaction) && inherits(fits$interaction, "rma")) {
    map[["tau[interaction]"]] <- fits$interaction
  } else if (identical(interactmeta_fit_analysis(fit), "da")) {
    int_rma <- construction_forest_interaction_rma(fit)
    if (!is.null(int_rma) && inherits(int_rma, "rma")) {
      map[["tau[interaction]"]] <- int_rma
    }
  }
  prob <- as.numeric(prob)[[1L]]
  for (i in seq_len(nrow(summary))) {
    par <- as.character(summary$parameter[[i]])
    model <- map[[par]]
    if (is.null(model)) {
      next
    }
    tau2 <- as.numeric(model$tau2)[[1L]]
    tau <- sqrt(max(tau2, 0))
    summary$estimate[[i]] <- tau
    summary$tau2[[i]] <- tau2
    interval <- interactmeta_metafor_rma_tau_ci(model, prob = prob)
    summary$lci[[i]] <- interval[["lci"]]
    summary$uci[[i]] <- interval[["uci"]]
  }
  summary
}

#' @noRd
interactmeta_extract_heterogeneity <- function(fit, prob = 0.95, ...) {
  interactmeta_require_capability(fit, "heterogeneity", feature = "heterogeneity()")
  prob <- as.numeric(prob)[[1L]]
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
  if (identical(analysis, "da") || identical(analysis, "ad") ||
        identical(analysis, "swada") ||
        identical(analysis, "conventional")) {
    int <- as.data.frame(fit$interaction, stringsAsFactors = FALSE)
    tau2 <- int$tau2[[1L]]
    tau <- if ("tau" %in% names(int)) int$tau[[1L]] else NA_real_
    if (identical(analysis, "da") &&
        (!is.finite(tau2) || is.na(tau2))) {
      int_rma <- construction_forest_interaction_rma(fit)
      if (!is.null(int_rma) && inherits(int_rma, "rma")) {
        tau2 <- as.numeric(int_rma$tau2)[[1L]]
        tau <- sqrt(max(tau2, 0))
      }
    }
    add_tau(
      "tau[interaction]",
      tau2,
      tau = tau
    )
  }
  out <- list(
    summary = dplyr::bind_rows(rows),
    analysis = analysis,
    engine = fit$engine,
    method = fit$method,
    inference = fit$inference
  )
  if (interactmeta_is_bayesian_fit(fit)) {
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
  } else {
    out$summary <- interactmeta_construction_metafor_fill_heterogeneity(
      fit,
      out$summary,
      prob = prob
    )
  }
  attr(out, "prob") <- prob
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
tidy.interactmeta_construction_subgroups <- tidy.interactmeta_construction_interaction

#' @export
#' @importFrom generics tidy
tidy.interactmeta_construction_heterogeneity <- tidy.interactmeta_construction_interaction

# ---- da_ad_bayes.R ----
#' Bayesian DA and AD via independent bayesmeta random-effects models
#'
#' DA: two independent univariate Bayesian MAs, then μ_B − μ_A by
#' bayesmeta::convolve(). The Gaussian mixture is a check. AD: one Bayesian
#' MA of within-trial g_j, with intervals from $post.interval. Neither uses CAMS IF.
#' @noRd
NULL

#' @noRd
interactmeta_construction_tau_prior <- function(analysis, tau.prior = NULL) {
  if (!is.null(tau.prior)) {
    return(normalize_cam_prior(tau.prior, role = "tau"))
  }
  if (identical(analysis, "ad")) {
    cam_default_tau_interaction_prior()
  } else {
    cam_default_tau_prior()
  }
}

#' @noRd
interactmeta_construction_mu_prior <- function(mu.prior = NULL) {
  if (is.null(mu.prior)) {
    return(cam_default_fixed_prior())
  }
  normalize_cam_prior(mu.prior, role = "fixed")
}

#' @noRd
interactmeta_bayesmeta_univariate <- function(y,
                                              vi,
                                              labels,
                                              tau.prior,
                                              mu.prior,
                                              engine_args = list()) {
  check_bayesmeta()
  tau_bm <- prior_to_bayesmeta(tau.prior)
  fe <- prior_to_bayesmeta_fixed(mu.prior)
  args <- if (is.null(engine_args)) list() else engine_args
  if ("convolve_derived" %in% names(args)) {
    stop(
      "`convolve_derived` cannot change reported intervals. Independent sums use bayesmeta::convolve().",
      call. = FALSE
    )
  }
  interval.type <- interactmeta_match_interval_type(args$interval.type)
  args$interval.type <- NULL
  args[c("method", "engine", "chains", "iter", "warmup")] <- NULL
  if (!anyNA(fe$mu.prior)) {
    args$mu.prior <- unname(fe$mu.prior)
  }
  do.call(
    bayesmeta::bayesmeta,
    c(
      list(
        y = as.numeric(y),
        sigma = interactmeta_sigma_from_vi(vi),
        labels = as.character(labels),
        tau.prior = tau_bm,
        interval.type = interval.type
      ),
      args
    )
  )
}

#' @noRd
interactmeta_bayesmeta_tau_summary <- function(fit, prob = 0.95) {
  tau <- as.numeric(fit$support[, "tau"])
  w <- as.numeric(fit$support[, "weight"])
  interactmeta_discrete_summary(interactmeta_discrete(tau, w), prob = prob)
}

#' @noRd
interactmeta_bayesmeta_tau2_mean <- function(fit) {
  tau <- as.numeric(fit$support[, "tau"])
  w <- as.numeric(fit$support[, "weight"])
  w <- w / sum(w)
  sum(w * tau^2)
}

#' @noRd
interactmeta_bayesmeta_re_weights <- function(fit) {
  tau2 <- interactmeta_bayesmeta_tau2_mean(fit)
  sig2 <- as.numeric(fit$sigma)^2
  w <- 1 / (sig2 + tau2)
  w / sum(w)
}

#' @noRd
interactmeta_bayesmeta_mu_row <- function(fit, label, method = "bayes", prob = 0.95) {
  it <- interactmeta_match_interval_type(fit$interval.type)
  loc <- interactmeta_native_location(fit, prob = prob, interval.type = it, margin = "effect")
  tau_loc <- interactmeta_native_location(fit, prob = prob, interval.type = it, margin = "tau")
  data.frame(
    label = label,
    estimate = loc$median,
    se = as.numeric(fit$summary["sd", "mu"])[[1L]],
    ci.lb = loc$lower,
    ci.ub = loc$upper,
    tau2 = interactmeta_bayesmeta_tau2_mean(fit),
    tau = tau_loc$median,
    method = method,
    stringsAsFactors = FALSE
  )
}

#' Bayesian DA: independent RE-MA per subgroup, then μ_B − μ_A
#' @noRd
interactmeta_da_from_pairs_bayesmeta <- function(pairs,
                                                 tau.prior,
                                                 mu.prior,
                                                 engine_args,
                                                 call) {
  ref <- attr(pairs, "reference")
  nonref <- attr(pairs, "nonref")
  fit_ref <- interactmeta_bayesmeta_univariate(
    y = pairs$y_ref,
    vi = pairs$v_ref,
    labels = pairs$trial,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    engine_args = engine_args
  )
  fit_non <- interactmeta_bayesmeta_univariate(
    y = pairs$y_nonref,
    vi = pairs$v_nonref,
    labels = pairs$trial,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    engine_args = engine_args
  )
  gm_ref <- interactmeta_gm_from_bayesmeta(fit_ref)
  gm_non <- interactmeta_gm_from_bayesmeta(fit_non)
  # Exact mixture is a check. The published difference is convolve().
  gm_da <- interactmeta_gm_linear_sum(gm_non, 1, gm_ref, -1)
  grid_delta <- fit_non$delta %||% fit_ref$delta
  grid_eps <- fit_non$epsilon %||% fit_ref$epsilon
  conv <- interactmeta_convolve_sum(
    dens1 = interactmeta_bayesmeta_mu_density(fit_non),
    dens2 = interactmeta_negate_density(interactmeta_bayesmeta_mu_density(fit_ref)),
    cdf1 = interactmeta_bayesmeta_mu_cdf(fit_non),
    cdf2 = interactmeta_negate_cdf(interactmeta_bayesmeta_mu_cdf(fit_ref)),
    delta = grid_delta,
    epsilon = grid_eps
  )
  it <- interactmeta_match_interval_type(fit_non$interval.type)
  tol <- fit_non$tol.uniroot %||% fit_ref$tol.uniroot
  bounds <- interactmeta_interval_from_distribution(
    conv$quantile, conv$cdf, prob = 0.95, interval.type = it, tol = tol
  )
  sm_da <- list(
    median = as.numeric(conv$quantile(0.5))[[1L]],
    mean = interactmeta_gm_mean(gm_da),
    sd = sqrt(max(interactmeta_gm_var(gm_da), 0)),
    lower = bounds[[1L]],
    upper = bounds[[2L]]
  )
  sm_ref <- interactmeta_bayesmeta_mu_row(fit_ref, ref)
  sm_non <- interactmeta_bayesmeta_mu_row(fit_non, nonref)
  list(
    reference = ref,
    nonref = nonref,
    pairs = pairs,
    subgroups = rbind(sm_ref, sm_non),
    interaction = data.frame(
      label = paste0(nonref, " - ", ref, " (DA)"),
      estimate = sm_da$median,
      se = sm_da$sd,
      ci.lb = sm_da$lower,
      ci.ub = sm_da$upper,
      tau2 = NA_real_,
      tau = NA_real_,
      method = "bayes",
      stringsAsFactors = FALSE
    ),
    weights = data.frame(
      trial = pairs$trial,
      weight_ref = interactmeta_bayesmeta_re_weights(fit_ref),
      weight_nonref = interactmeta_bayesmeta_re_weights(fit_non),
      stringsAsFactors = FALSE
    ),
    fits = list(reference = fit_ref, nonref = fit_non),
    posterior = list(
      mu_reference = interactmeta_posterior_1d(
        gm = gm_ref,
        name = "mu_reference",
        published = interactmeta_native_published(fit_ref, "effect"),
        interval.type = it,
        tol.uniroot = tol
      ),
      mu_nonref = interactmeta_posterior_1d(
        gm = gm_non,
        name = "mu_nonref",
        published = interactmeta_native_published(fit_non, "effect"),
        interval.type = it,
        tol.uniroot = tol
      ),
      interaction = interactmeta_posterior_1d(
        gm = gm_da,
        convolve = conv,
        name = "gamma_DA",
        note = "Independent Bayesian DA: published posterior is convolve(mu_B, -mu_A); gm is the mixture check.",
        interval.type = it,
        tol.uniroot = tol
      )
    ),
    priors = list(tau = tau.prior, fixed = mu.prior),
    method = "bayes",
    call = call
  )
}

#' Bayesian AD: RE-MA of g_j = y_B − y_A
#' @noRd
interactmeta_ad_from_pairs_bayesmeta <- function(pairs,
                                                 tau.prior,
                                                 mu.prior,
                                                 engine_args,
                                                 call) {
  ref <- attr(pairs, "reference")
  nonref <- attr(pairs, "nonref")
  fit_g <- interactmeta_bayesmeta_univariate(
    y = pairs$g,
    vi = pairs$vg,
    labels = pairs$trial,
    tau.prior = tau.prior,
    mu.prior = mu.prior,
    engine_args = engine_args
  )
  sm <- interactmeta_bayesmeta_mu_row(fit_g, paste0(nonref, " - ", ref))
  gm <- interactmeta_gm_from_bayesmeta(fit_g)
  list(
    reference = ref,
    nonref = nonref,
    pairs = pairs,
    interaction = sm,
    weights = data.frame(
      trial = pairs$trial,
      weight = interactmeta_bayesmeta_re_weights(fit_g),
      stringsAsFactors = FALSE
    ),
    fits = list(interaction = fit_g),
    posterior = list(
      interaction = interactmeta_posterior_1d(
        gm = gm,
        name = "gamma_AD",
        published = interactmeta_native_published(fit_g, "effect"),
        interval.type = fit_g$interval.type
      )
    ),
    priors = list(tau = tau.prior, fixed = mu.prior),
    method = "bayes",
    call = call
  )
}

# ---- swada.R ----
#' Apply the Same Weighting Across Different Analyses principle
#'
#' **SWADA** (Same Weighting Across Different Analyses) is a *consistency
#' principle*, not a statistical model parallel to [da()], [ad()], or [cams()].
#' [da()] (Difference of Averages) and [ad()] (Average Difference) are
#' **different constructions** (Panaro, Röver & Friede 2026) and can assign
#' studies different effective contributions.
#' SWADA requires the same study weights for the subgroup-mean analyses and
#' for the within-trial interaction, so that
#' \eqn{\mu_B - \mu_A = \gamma}.
#'
#' The default (`weights = "interaction"`) is Interaction RE-weights SWADA
#' (Panaro, Röver & Friede 2026): take the [ad()] weights
#' \eqn{w_j \propto 1/(v_{g j}+\hat\tau_g^2)} and apply them to both
#' subgroup series. Under that common vector the DA interaction is
#' algebraically identical to AD: \eqn{\sum_j w_j(y_{Bj}-y_{Aj})=\sum_j w_j g_j}.
#' The interaction estimate itself is the AD/`rma` of \eqn{g_j}; SWADA does not
#' introduce a third interaction model.
#'
#' @param data A data frame with `trial`, `subgroup`, `yi`, and `vi`.
#' @param reference Reference subgroup label.
#' @param trial,subgroup,yi,vi Column names. `subgroup` remains accepted.
#' @param modifier Column of the two levels. Overrides `subgroup` when set.
#'   The default column remains `"subgroup"`.
#' @param weights Weighting scheme. Currently `"interaction"` (recommended).
#' @param engine Computational backend. Currently `"metafor"` only. Bayesian
#'   SWADA is not implemented; Bayesian inference does not add CAMS IF
#'   functionality.
#' @param approach Inferential approach. SWADA is frequentist only; `"bayes"`
#'   is an error.
#' @param method Estimator for the interaction RE-MA (`"REML"` default, or `"ML"`).
#' @param ... Passed to [ad()] / [metafor::rma()] for the interaction fit.
#' @return An object of class `"interactmeta_swada"`.
#' @references
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \url{https://doi.org/10.48550/arXiv.2508.15531}
#' @seealso [da()], [ad()], [cams()], [conventional()]
#' @export
swada <- function(data,
                  reference,
                  modifier = NULL,
                  trial = "study",
                  subgroup = "level",
                  yi = "yi",
                  vi = "vi",
                  weights = c("interaction"),
                  method = "REML",
                  ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup, dat = data)
  weights <- match.arg(weights)
  dots <- list(...)
  resolved <- interactmeta_resolve_construction_method("swada", method)
  ad_fit <- do.call(
    ad,
    c(
      list(
        data = data,
        reference = reference,
        trial = trial,
        subgroup = subgroup,
        yi = yi,
        vi = vi,
        method = resolved$method
      ),
      dots
    )
  )
  dat <- ad_fit$data
  pairs <- ad_fit$pairs
  ref <- ad_fit$reference
  nonref <- ad_fit$nonref
  fit_g <- ad_fit$fits$interaction
  w_raw <- as.numeric(ad_fit$weights$weight)
  w <- w_raw / sum(w_raw)
  sm_g <- ad_fit$interaction
  sm_ref <- interactmeta_swada_weighted_mean(
    y = pairs$y_ref,
    vi = pairs$v_ref,
    w = w,
    label = ref
  )
  sm_non <- interactmeta_swada_weighted_mean(
    y = pairs$y_nonref,
    vi = pairs$v_nonref,
    w = w,
    label = nonref
  )
  g_from_sub <- sm_non$estimate - sm_ref$estimate
  out <- list(
      approach = "swada",
      call = match.call(),
      reference = ref,
      nonref = nonref,
      weights_scheme = weights,
      pairs = pairs,
      subgroups = rbind(sm_ref, sm_non),
      interaction = sm_g,
      collapsibility = list(
        mu_nonref_minus_mu_ref = g_from_sub,
        interaction = sm_g$estimate,
        equal = isTRUE(all.equal(g_from_sub, sm_g$estimate, tolerance = 1e-8))
      ),
      weights = data.frame(
        trial = pairs$trial,
        weight = w,
        weight_ref = w,
        weight_nonref = w,
        weight_percent = 100 * w,
        stringsAsFactors = FALSE
      ),
      fits = list(interaction = fit_g),
      ad = ad_fit,
      data = dat,
      modifier = subgroup,
      assumptions = list(
        subgroup_covariance = "assumed_zero",
        notes = "SWADA principle: one weight vector for DA and AD. Default Interaction RE-weights make the DA interaction identical to AD. Not a CAMS information-fraction analysis."
      )
    )
  class(out) <- c("interactmeta_swada", "list")
  interactmeta_attach_capabilities(out, "swada", resolved$method)
}

interactmeta_swada_weighted_mean <- function(y, vi, w, label) {
  y <- as.numeric(y)
  vi <- as.numeric(vi)
  w <- as.numeric(w)
  w <- w / sum(w)
  est <- sum(w * y)
  se <- sqrt(sum(w^2 * vi))
  z <- stats::qnorm(0.975)
  data.frame(
    label = label,
    estimate = est,
    se = se,
    ci.lb = est - z * se,
    ci.ub = est + z * se,
    tau2 = NA_real_,
    method = "SWADA-WLS",
    stringsAsFactors = FALSE
  )
}

#' SWADA predictive interval: common weights + interaction heterogeneity
#'
#' Point estimate is the SWADA weighted mean; the interval adds
#' \eqn{\hat\tau_g^2} from the interaction RE-MA (same rule as forest plots).
#' @noRd
interactmeta_swada_wls_prediction <- function(y, vi, w, tau2, prob = 0.95) {
  y <- as.numeric(y)
  vi <- as.numeric(vi)
  w <- as.numeric(w)
  w <- w / sum(w)
  est <- sum(w * y)
  se_fe <- sqrt(sum(w^2 * vi))
  tau2 <- as.numeric(tau2)[[1L]]
  se <- sqrt(se_fe^2 + max(tau2, 0, na.rm = TRUE))
  z <- stats::qnorm(1 - (1 - prob) / 2)
  list(
    estimate = est,
    se = se,
    ci.lb = est - z * se,
    ci.ub = est + z * se
  )
}

#' @noRd
interactmeta_swada_interaction_tau2 <- function(fit) {
  tau2 <- fit$interaction$tau2
  if (length(tau2) && is.finite(tau2[[1L]])) {
    return(as.numeric(tau2)[[1L]])
  }
  fit_g <- fit$fits$interaction
  if (!is.null(fit_g) && inherits(fit_g, "rma")) {
    return(as.numeric(fit_g$tau2)[[1L]])
  }
  0
}

#' SWADA interaction or subgroup prediction (forest plots and `interactions()`)
#' @noRd
interactmeta_swada_prediction <- function(fit,
                                           target = c("interaction", "subgroup"),
                                           label = NULL,
                                           prob = 0.95) {
  target <- match.arg(target)
  pairs <- fit$pairs
  w <- fit$weights$weight
  tau2 <- interactmeta_swada_interaction_tau2(fit)
  if (identical(target, "interaction")) {
    interactmeta_swada_wls_prediction(
      y = pairs$g,
      vi = pairs$vg,
      w = w,
      tau2 = tau2,
      prob = prob
    )
  } else {
    if (is.null(label)) {
      stop("`label` is required for subgroup SWADA prediction.", call. = FALSE)
    }
    y <- if (identical(label, fit$reference)) pairs$y_ref else pairs$y_nonref
    vi <- if (identical(label, fit$reference)) pairs$v_ref else pairs$v_nonref
    interactmeta_swada_wls_prediction(y, vi, w, tau2, prob = prob)
  }
}

#' @export
print.interactmeta_swada <- function(x, digits = 2, ...) {
  cat("SWADA principle (same weighting across DA and AD)\n")
  cat("  scheme    : ", x$weights_scheme, "\n", sep = "")
  cat("  method    : ", x$method, "\n", sep = "")
  cat("  inference : ", x$inference %||% NA_character_, "\n", sep = "")
  cat("  reference : ", x$reference, "\n", sep = "")
  cat("\nSubgroup-specific treatment effects (common weights; log scale):\n")
  interactmeta_print_effect_rows(x$subgroups, digits)
  cat("\nSWADA-consistent interaction (equals AD under interaction weights):\n")
  interactmeta_print_effect_rows(x$interaction, digits)
  cat(
    "\nCollapsible: DA interaction ",
    if (isTRUE(x$collapsibility$equal)) "==" else "!=",
    " AD interaction.\n",
    sep = ""
  )
  interactmeta_print_fit_heterogeneity(x, digits = digits, ...)
  invisible(x)
}

#' @export
#' @importFrom generics tidy
tidy.interactmeta_swada <- function(x, exponentiate = FALSE, ...) {
  sg <- x$subgroups
  names(sg)[names(sg) == "label"] <- "subgroup"
  sg$component <- "subgroup"
  int <- x$interaction
  names(int)[names(int) == "label"] <- "subgroup"
  int$component <- "interaction"
  out <- dplyr::bind_rows(sg, int)
  if (isTRUE(exponentiate)) {
    out <- interactmeta_exponentiate_wald(out)
  }
  tibble::as_tibble(out)
}

#' Coerce a conventional bundle to its DA component
#' @noRd
interactmeta_as_da <- function(x) {
  if (inherits(x, "interactmeta_da")) {
    return(x)
  }
  if (inherits(x, "interactmeta_conventional")) {
    if (!is.null(x$da) && inherits(x$da, "interactmeta_da")) {
      return(x$da)
    }
    x$approach <- "da"
    class(x) <- c("interactmeta_da", "list")
    x$interaction <- x$interaction_da
    return(x)
  }
  stop("`da` must be a `da()` or `conventional()` fit.", call. = FALSE)
}

#' Coerce a conventional bundle or dedicated AD object
#' @noRd
interactmeta_as_ad <- function(x, from_conventional = NULL) {
  if (inherits(x, "interactmeta_ad")) {
    return(x)
  }
  if (inherits(x, "interactmeta_conventional")) {
    if (!is.null(x$ad) && inherits(x$ad, "interactmeta_ad")) {
      return(x$ad)
    }
    x$approach <- "ad"
    class(x) <- c("interactmeta_ad", "list")
    return(x)
  }
  if (!is.null(from_conventional) && inherits(from_conventional, "interactmeta_conventional")) {
    return(interactmeta_as_ad(from_conventional))
  }
  stop("`ad` must be an `ad()` or `conventional()` fit.", call. = FALSE)
}

#' Side-by-side modifier-level treatment effects from DA, SWADA, and CAMS
#'
#' Compares pooled **modifier-level** effects from `DA`, `SWADA`, and `CAMS`
#' ([da()] or the DA part of [conventional()], [swada()], [cams()]).
#' These are not three equivalent statistical models: DA pools each modifier
#' series separately, SWADA applies one common weight vector, and CAMS is the
#' contribution-adjusted model. Average Difference ([ad()]) does not pool
#' modifier-specific treatment effects, so it is omitted here.
#'
#' Warns if mixed constructions do not share the same `engine` and `method`.
#' A like-for-like frequentist comparison uses `metafor` / `REML` throughout
#' (`cams(..., method = "REML")`). Two CAMS fits of the same model with
#' different methods are labelled (`CAMS — Bayesian` vs `CAMS — REML`) and
#' are not treated as a mismatch.
#'
#' Two or more [cams()] fits may be passed instead of DA / SWADA / CAMS, in
#' which case [plot()] draws the compared subgroup estimates.
#'
#' @param da A [da()] or [conventional()] fit, or a [cams()] fit when comparing
#'   two CAMS approaches.
#' @param swada A [swada()] fit, or a second [cams()] fit.
#' @param cam A [cams()] CAMS fit. Optional when `da` and `swada`
#'   are already CAMS fits.
#' @param exponentiate If `TRUE` (default), report on the exponentiated
#'   (ratio) scale. Use `FALSE` for the fitted log scale.
#' @param type CAMS layer (`"Overall"` default).
#' @return A tibble with one row per approach × subgroup and class
#'   `interactmeta_compare_subgroups` (plottable).
#' @name compare_modifiers
#' @aliases compare_subgroups
#' @keywords internal
compare_modifiers <- function(da,
                                        swada,
                                        cam = NULL,
                                        exponentiate = TRUE,
                                        type = "Overall") {
  cams_fits <- interactmeta_compact_fits(da, swada, cam)
  if (interactmeta_all_cams(cams_fits)) {
    return(interactmeta_compare_subgroups_fits(
      cams_fits,
      exponentiate = exponentiate,
      type = type
    ))
  }
  da <- interactmeta_as_da(da)
  stopifnot(inherits(swada, "interactmeta_swada"))
  if (is.null(cam) || !inherits(cam, "interactmeta")) {
    stop(
      "`compare_modifiers()` needs DA, SWADA, and CAMS fits, ",
      "or two or more CAMS fits.",
      call. = FALSE
    )
  }
  interactmeta_warn_incomparable_fits(
    list(da, swada, cam),
    c("DA", "SWADA", "CAMS"),
    "compare_modifiers"
  )
  da_tab <- tidy(da, exponentiate = exponentiate)
  da_tab <- da_tab[da_tab$component == "subgroup", , drop = FALSE]
  da_tab$approach <- "DA"
  sw <- tidy(swada, exponentiate = exponentiate)
  sw <- sw[sw$component == "subgroup", , drop = FALSE]
  sw$approach <- "SWADA"
  sg <- modifiers(cam, type = type, exponentiate = exponentiate, summary = TRUE)
  td <- tidy(sg)
  td <- td[td$type == type, , drop = FALSE]
  lvl_col <- interactmeta_modifier_summary_col(td)
  cams_tab <- data.frame(
    subgroup = as.character(td[[lvl_col]]),
    estimate = td$estimate,
    se = td$std.error,
    ci.lb = td$conf.low,
    ci.ub = td$conf.high,
    approach = "CAMS",
    stringsAsFactors = FALSE
  )
  out <- dplyr::bind_rows(
    da_tab[, c("approach", "subgroup", "estimate", "se", "ci.lb", "ci.ub")],
    sw[, c("approach", "subgroup", "estimate", "se", "ci.lb", "ci.ub")],
    cams_tab
  )
  sg_levels <- c(da$reference, da$nonref)
  out$subgroup <- factor(as.character(out$subgroup), levels = sg_levels)
  out$approach <- factor(
    out$approach,
    levels = c("DA", "SWADA", "CAMS")
  )
  out <- out[order(out$subgroup, out$approach), , drop = FALSE]
  rownames(out) <- NULL
  xlab <- if (isTRUE(exponentiate)) "Treatment effect" else "Log treatment effect"
  interactmeta_compare_result(
    out,
    "interactmeta_compare_subgroups",
    exponentiate = exponentiate,
    xlab = xlab
  )
}

#' Side-by-side treatment-by-subgroup interactions from DA, AD, SWADA, and CAMS
#'
#' Under the default Interaction RE-weights, the SWADA-consistent interaction
#' equals [ad()]; it is listed once as `SWADA` rather than as two duplicate
#' SWADA methods. [da()] and [ad()] are different constructions (Panaro,
#' Röver & Friede 2026) and generally differ because they use different
#' study weights.
#'
#' The first argument may be a [conventional()] bundle, in which case DA and
#' AD are read from that object and `ad` may be omitted when the remaining
#' arguments are `swada` and `cam` (historical three-argument call).
#'
#' Warns if mixed constructions do not share the same `engine` and `method`.
#' A like-for-like frequentist comparison uses `metafor` / `REML` throughout
#' (`da()` / `ad()` / `swada()` defaults, plus `cams(..., method = "REML")`).
#' Two CAMS fits of the same model with different methods are labelled
#' (`CAMS — Bayesian` vs `CAMS — REML`) and are not treated as a mismatch.
#'
#' Two or more [cams()] fits may be passed instead of DA / AD / SWADA / CAMS.
#' [plot()] draws the compared interaction estimates.
#'
#' @param da A [da()] or [conventional()] fit, or a [cams()] fit when comparing
#'   two CAMS approaches.
#' @param ad An [ad()] fit. Optional when `da` is a [conventional()] bundle, or
#'   a second [cams()] fit.
#' @param swada A [swada()] fit.
#' @param cam A [cams()] CAMS fit.
#' @inheritParams compare_modifiers
#' @param exponentiate If `TRUE`, report on the exponentiated (ratio) scale.
#'   Default `FALSE` (fitted log scale).
#' @return A tibble with class `interactmeta_compare_interactions`.
#' @keywords internal
compare_interactions <- function(da,
                                           ad = NULL,
                                           swada = NULL,
                                           cam = NULL,
                                           exponentiate = FALSE,
                                           type = "Overall") {
  cams_fits <- interactmeta_compact_fits(da, ad, swada, cam)
  if (interactmeta_all_cams(cams_fits)) {
    return(interactmeta_compare_interactions_fits(
      cams_fits,
      exponentiate = exponentiate,
      type = type
    ))
  }
  unpacked <- interactmeta_unpack_interaction_args(da, ad, swada, cam)
  da_fit <- unpacked$da
  ad_fit <- unpacked$ad
  swada <- unpacked$swada
  cam <- unpacked$cam
  stopifnot(inherits(da_fit, "interactmeta_da"))
  stopifnot(inherits(ad_fit, "interactmeta_ad"))
  stopifnot(inherits(swada, "interactmeta_swada"))
  stopifnot(inherits(cam, "interactmeta"))
  interactmeta_warn_incomparable_fits(
    list(da_fit, ad_fit, swada, cam),
    c("DA", "AD", "SWADA", "CAMS"),
    "compare_interactions"
  )
  pick <- function(approach, row, exponentiate, contrast) {
    est <- row$estimate
    se <- row$se
    lo <- row$ci.lb
    hi <- row$ci.ub
    if (isTRUE(exponentiate)) {
      se <- se * exp(est)
      est <- exp(est)
      lo <- exp(lo)
      hi <- exp(hi)
    }
    data.frame(
      approach = approach,
      contrast = contrast,
      estimate = est,
      se = se,
      ci.lb = lo,
      ci.ub = hi,
      stringsAsFactors = FALSE
    )
  }
  ticket <- interactmeta_compare_contrast_ticket(da_fit$nonref, da_fit$reference)
  da_row <- pick("DA", da_fit$interaction, exponentiate, paste0(ticket, " (DA)"))
  ad_row <- pick("AD", ad_fit$interaction, exponentiate, ticket)
  sw_row <- pick("SWADA", swada$interaction, exponentiate, ticket)
  ix <- interactions(cam, type = type, exponentiate = exponentiate, summary = TRUE)
  td <- tidy(ix)
  td <- td[td$type == type, , drop = FALSE]
  cams_tab <- data.frame(
    approach = "CAMS",
    contrast = ticket,
    estimate = td$estimate[[1L]],
    se = td$std.error[[1L]],
    ci.lb = td$conf.low[[1L]],
    ci.ub = td$conf.high[[1L]],
    stringsAsFactors = FALSE
  )
  out <- dplyr::bind_rows(da_row, ad_row, sw_row, cams_tab)
  out$approach <- factor(
    out$approach,
    levels = c("DA", "AD", "SWADA", "CAMS")
  )
  xlab <- if (isTRUE(exponentiate)) "Treatment effect" else "Log treatment effect"
  interactmeta_compare_result(
    out,
    "interactmeta_compare_interactions",
    exponentiate = exponentiate,
    xlab = xlab
  )
}

#' @noRd
interactmeta_unpack_interaction_args <- function(da, ad, swada, cam) {
  # Historical call: compare_interactions(conventional, swada, cam)
  if (inherits(da, "interactmeta_conventional") &&
      inherits(ad, "interactmeta_swada") &&
      inherits(swada, "interactmeta") &&
      is.null(cam)) {
    return(list(
      da = interactmeta_as_da(da),
      ad = interactmeta_as_ad(da),
      swada = ad,
      cam = swada
    ))
  }
  if (inherits(da, "interactmeta_conventional") && is.null(ad)) {
    ad <- interactmeta_as_ad(da)
  }
  list(
    da = interactmeta_as_da(da),
    ad = interactmeta_as_ad(ad, from_conventional = da),
    swada = swada,
    cam = cam
  )
}

#' Shared contrast ticket for AD / SWADA / CAMS; DA appends ` (DA)`
#' @noRd
interactmeta_compare_contrast_ticket <- function(nonref, ref) {
  paste(as.character(nonref)[[1L]], as.character(ref)[[1L]], sep = " - ")
}

#' Engine / method label stored on a construction or CAMS fit
#' @noRd
interactmeta_fit_engine_method <- function(x) {
  engine <- as.character(x$engine %||% NA_character_)[[1L]]
  method <- as.character(x$method %||% NA_character_)[[1L]]
  paste0(engine, " / ", method)
}

#' Warn when mixed constructions do not share engine and method
#'
#' Same-analysis comparisons (two CAMS approaches, two DA approaches, …) are the
#' point of an inference comparison and do not warn.
#' @noRd
interactmeta_warn_incomparable_fits <- function(fits, labels, caller) {
  keys <- vapply(fits, interactmeta_fit_engine_method, character(1L))
  if (length(unique(keys)) <= 1L) {
    return(invisible(NULL))
  }
  analyses <- vapply(fits, interactmeta_fit_analysis, character(1L))
  if (length(unique(analyses)) <= 1L) {
    return(invisible(NULL))
  }
  lines <- paste0("  ", labels, ": ", keys, collapse = "\n")
  warning(
    caller,
    "() mixes inferential settings:\n",
    lines,
    "\nEstimates and intervals are not like-for-like. ",
    "Fit every object with the same approach and method ",
    "(metafor / REML for a frequentist comparison).",
    call. = FALSE
  )
}

#' Trial contribution weights from DA, AD, and SWADA
#'
#' [da()] generally uses **different** inverse-variance weights for each
#' subgroup series. [ad()] uses the interaction-series weights. [swada()]
#' applies **one common vector** to both subgroup series (`swada_reference`
#' and `swada_nonref` are identical); under the default this equals the AD
#' weights.
#'
#' Warns if the supplied fits do not share the same `engine` and `method`.
#'
#' Historical two-argument call `compare_weights(conventional, swada)`
#' still works. CAMS information fractions are not SWADA weights; compare
#' pooled IF with [compare_ifrac()].
#'
#' @param da A [da()] or [conventional()] fit.
#' @param ad An [ad()] fit, or a [swada()] fit in the historical two-argument call.
#' @param swada A [swada()] fit. Optional when `ad` is already a SWADA object.
#' @return A tibble with one row per trial and class
#'   `interactmeta_compare_weights` (plottable).
#' @keywords internal
compare_weights <- function(da, ad = NULL, swada = NULL) {
  cams_fits <- interactmeta_compact_fits(da, ad, swada)
  if (interactmeta_all_cams(cams_fits)) {
    stop(
      "`compare_weights()` compares DA / AD / SWADA contribution weights. ",
      "For CAMS pooled information fractions use compare_ifrac().",
      call. = FALSE
    )
  }
  if (inherits(ad, "interactmeta_swada") && is.null(swada)) {
    swada <- ad
    ad <- NULL
  }
  stopifnot(inherits(swada, "interactmeta_swada"))
  da_fit <- interactmeta_as_da(da)
  if (is.null(ad) && inherits(da, "interactmeta_conventional")) {
    ad <- interactmeta_as_ad(da)
  }
  ad_fit <- interactmeta_as_ad(ad, from_conventional = da)
  interactmeta_warn_incomparable_fits(
    list(da_fit, ad_fit, swada),
    c("DA", "AD", "SWADA"),
    "compare_weights"
  )
  w_da <- da_fit$weights
  w_ad <- ad_fit$weights
  w_s <- swada$weights
  w_ad_vec <- if ("weight" %in% names(w_ad)) {
    as.numeric(w_ad$weight)
  } else {
    as.numeric(w_ad$weight_interaction)
  }
  stopifnot(identical(as.character(w_da$trial), as.character(w_s$trial)))
  w_sw <- as.numeric(w_s$weight)
  w_sw_norm <- w_sw / sum(w_sw)
  interactmeta_compare_result(
    tibble::tibble(
      trial = as.character(w_da$trial),
      da_reference = as.numeric(w_da$weight_ref) / sum(w_da$weight_ref),
      da_nonref = as.numeric(w_da$weight_nonref) / sum(w_da$weight_nonref),
      ad = w_ad_vec / sum(w_ad_vec),
      swada_reference = w_sw_norm,
      swada_nonref = w_sw_norm
    ),
    "interactmeta_compare_weights",
    exponentiate = FALSE,
    xlab = "Weight"
  )
}

#' Trial-level CAMS quantities \eqn{\pi_j}, \eqn{m_j}, \eqn{g_j}
#'
#' Reads the transformed trial blocks stored on a [cams()] fit
#' (`fit$transformed_data$trial_blocks`). For binary CAMS,
#' \eqn{g_j = y_{B j} - y_{A j}} with `A` the reference subgroup.
#'
#' @param fit An `interactmeta` fit.
#' @return A data frame with `trial`, `pi`, `m`, `v_m`, `g`, and `v_g`.
#' @name cams_trial_quantities
#' @export
cams_trial_quantities <- function(fit) {
  stopifnot(inherits(fit, "interactmeta"))
  blocks <- fit$transformed_data$trial_blocks
  if (is.null(blocks)) {
    stop(
      "No binary trial blocks on `fit$transformed_data`. ",
      "CAMS intermediate quantities are currently stored for binary subgroups.",
      call. = FALSE
    )
  }
  ov <- blocks$overall
  int <- blocks$interaction
  pi <- blocks$pi
  out <- tibble::tibble(
    trial = as.character(ov$trial),
    pi = as.numeric(unname(pi[as.character(ov$trial)])),
    m = as.numeric(ov$yi),
    v_m = as.numeric(ov$vi),
    g = as.numeric(int$yi),
    v_g = as.numeric(int$vi)
  )
  attr(out, "reference") <- blocks$reference
  attr(out, "nonref") <- blocks$nonref
  out
}

#' @rdname cams_trial_quantities
#' @export
cam_trial_quantities <- cams_trial_quantities

#' Comparable forest of modifier-level treatment effects
#'
#' Plots DA, SWADA, and CAMS pooled modifier-level treatment effects on a
#' common axis. Observed trial estimates (shared input) can be included using
#' the original published intervals when those columns are present.
#'
#' @inheritParams compare_modifiers
#' @param include_observed If `TRUE`, add published trial rows from `da$data`.
#' @param xlim Shared x-axis limits on the displayed scale.
#' @return A ggplot object.
#' @name plot_modifiers
#' @export
plot_modifiers <- function(da,
                                     swada,
                                     cam,
                                     exponentiate = TRUE,
                                     type = "Overall",
                                     include_observed = TRUE,
                                     xlim = NULL) {
  da_fit <- interactmeta_as_da(da)
  tab <- compare_modifiers(
    da_fit,
    swada,
    cam,
    exponentiate = exponentiate,
    type = type
  )
  tab$row <- as.character(tab$approach)
  tab$trial <- NA_character_
  src <- if (inherits(da, "interactmeta_conventional")) da$data else da_fit$data
  if (isTRUE(include_observed) && all(c("hr", "ci_lower", "ci_upper") %in% names(src))) {
    obs <- as.data.frame(src, stringsAsFactors = FALSE)
    obs_level_col <- interactmeta_resolve_level_col(obs)
    obs_tab <- data.frame(
      approach = "Published",
      subgroup = as.character(obs[[obs_level_col]]),
      estimate = if (isTRUE(exponentiate)) obs$hr else obs$yi,
      se = NA_real_,
      ci.lb = if (isTRUE(exponentiate)) obs$ci_lower else log(obs$ci_lower),
      ci.ub = if (isTRUE(exponentiate)) obs$ci_upper else log(obs$ci_upper),
      row = as.character(obs$trial),
      trial = as.character(obs$trial),
      stringsAsFactors = FALSE
    )
    tab <- dplyr::bind_rows(obs_tab, tab)
  }
  sg_levels <- c(da_fit$reference, da_fit$nonref)
  tab$subgroup <- factor(as.character(tab$subgroup), levels = sg_levels)
  row_levels <- unique(as.character(tab$row))
  pooled_levels <- c("DA", "SWADA", "CAMS")
  row_levels <- c(
    setdiff(row_levels, pooled_levels),
    intersect(pooled_levels, row_levels)
  )
  tab$row <- factor(tab$row, levels = rev(row_levels))
  pal <- interactmeta_approach_colours(approaches = as.character(tab$approach))
  const <- interactmeta_plot_constants()
  xlab <- if (isTRUE(exponentiate)) "Treatment effect" else "Log treatment effect"
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$estimate, y = .data$row, colour = .data$approach)) +
    ggplot2::geom_vline(
      xintercept = if (isTRUE(exponentiate)) 1 else 0,
      colour = const$null_line,
      linewidth = 0.3
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$ci.lb, xmax = .data$ci.ub),
      linewidth = const$line_hair
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_grid(subgroup ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_manual(values = pal, breaks = names(pal)) +
    ggplot2::labs(x = xlab, y = NULL, colour = NULL) +
    interactmeta_theme_report()
  if (isTRUE(exponentiate)) {
    p <- p + ggplot2::scale_x_log10(limits = xlim)
  } else if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  p
}

#' @rdname ifrac
#' @export
ifrac.interactmeta_swada <- ifrac.interactmeta_da

#' @export
plot.interactmeta_swada <- plot.interactmeta_da

#' @noRd
fit_construction <- function(type, ...) {
  type <- match.arg(type, c("da", "ad", "swada", "cams"))
  switch(
    type,
    da = do.call(da, list(...)),
    ad = do.call(ad, list(...)),
    swada = do.call(swada, list(...)),
    cams = do.call(cams, list(...))
  )
}
