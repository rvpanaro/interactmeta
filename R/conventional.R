#' Pair binary-subgroup rows into trial-level estimates
#'
#' Returns one row per trial with reference / non-reference `yi`/`vi` and the
#' CAMS interaction contrast `g = y_nonref - y_ref`. Interaction variances use
#' `v_ref + v_nonref` (within-trial subgroup covariance assumed zero).
#' @noRd
interactmeta_binary_trial_table <- function(data,
                                            trial_var = "trial",
                                            subgroup_var = "subgroup",
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
#' same study weights on both constructions. [cam()] is a different object:
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
#' @seealso [ad()], [swada()], [cam()], [conventional()], [interactmeta_capability_matrix()]
#' @export
da <- function(data,
               reference,
               modifier = NULL,
               trial = "trial",
               subgroup = "subgroup",
               yi = "yi",
               vi = "vi",
               engine = c("metafor", "bayesmeta"),
               method = NULL,
               approach = NULL,
               tau.prior = NULL,
               mu.prior = NULL,
               ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup)
  interactmeta_fit_da_ad(
    analysis = "da",
    data = data,
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
    call = match.call(),
    ...
  )
}

#' Estimate a treatment-by-subgroup interaction using the Average-Difference construction
#'
#' **Average Difference (AD)** first forms the within-study
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
#' analysis-specific weights with one common vector. [cam()] is a
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
#' @seealso [da()], [swada()], [cam()], [conventional()], [interactmeta_capability_matrix()]
#' @export
ad <- function(data,
               reference,
               modifier = NULL,
               trial = "trial",
               subgroup = "subgroup",
               yi = "yi",
               vi = "vi",
               engine = c("metafor", "bayesmeta"),
               method = NULL,
               approach = NULL,
               tau.prior = NULL,
               mu.prior = NULL,
               ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup)
  interactmeta_fit_da_ad(
    analysis = "ad",
    data = data,
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
                                   engine,
                                   method,
                                   tau.prior,
                                   mu.prior,
                                   call,
                                   approach = NULL,
                                   ...) {
  resolved <- interactmeta_resolve_construction_engine_method(
    analysis,
    engine,
    method,
    approach = approach
  )
  if (identical(resolved$engine, "metafor") &&
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
  if (identical(resolved$engine, "bayesmeta")) {
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
        method = resolved$method,
        rma_args = dots,
        call = call
      )
    } else {
      built <- interactmeta_ad_from_pairs(
        pairs,
        method = resolved$method,
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
  interactmeta_attach_capabilities(out, analysis, resolved$engine, resolved$method)
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
#' applies one common weight vector. [cam()] fits the CAMS model.
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
#' @seealso [da()], [ad()], [swada()], [cam()]
#' @export
conventional <- function(data,
                         reference,
                         modifier = NULL,
                         trial = "trial",
                         subgroup = "subgroup",
                         yi = "yi",
                         vi = "vi",
                         engine = c("metafor", "bayesmeta"),
                         method = NULL,
                         approach = NULL,
                         tau.prior = NULL,
                         mu.prior = NULL,
                         ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup)
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

#' @export
print.interactmeta_da <- function(x, digits = 2, ...) {
  cat("Difference of Averages (DA)\n")
  cat("  analysis  : DA\n")
  cat("  engine    : ", x$engine, "\n", sep = "")
  cat("  method    : ", x$method, "\n", sep = "")
  cat("  inference : ", x$inference %||% NA_character_, "\n", sep = "")
  cat("  reference : ", x$reference, "\n", sep = "")
  cat("\nSubgroup pooled effects (log scale):\n")
  interactmeta_print_effect_rows(x$subgroups, digits)
  cat("\nTreatment-by-subgroup interaction (difference of pooled means):\n")
  interactmeta_print_effect_rows(x$interaction, digits)
  invisible(x)
}

#' @export
print.interactmeta_ad <- function(x, digits = 2, ...) {
  cat("Average Difference (AD)\n")
  cat("  analysis  : AD\n")
  cat("  engine    : ", x$engine, "\n", sep = "")
  cat("  method    : ", x$method, "\n", sep = "")
  cat("  inference : ", x$inference %||% NA_character_, "\n", sep = "")
  cat("  reference : ", x$reference, "\n", sep = "")
  cat("\nTreatment-by-subgroup interaction (RE-MA of within-trial contrasts):\n")
  interactmeta_print_effect_rows(x$interaction, digits)
  invisible(x)
}

#' @export
print.interactmeta_conventional <- function(x, digits = 2, ...) {
  cat("DA and AD with analysis-specific weights\n")
  cat("  engine    : ", x$engine, "\n", sep = "")
  cat("  method    : ", x$method, "\n", sep = "")
  cat("  inference : ", x$inference %||% NA_character_, "\n", sep = "")
  cat("  reference : ", x$reference, "\n", sep = "")
  cat("\nSubgroup pooled effects (DA; log scale):\n")
  interactmeta_print_effect_rows(x$subgroups, digits)
  cat("\nAverage Difference (AD; within-trial contrasts):\n")
  interactmeta_print_effect_rows(x$interaction, digits)
  cat("\nDifference of Averages (DA; not collapsible with AD):\n")
  interactmeta_print_effect_rows(x$interaction_da, digits)
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
ifrac.interactmeta_ad <- function(x, ...) {
  interactmeta_abort_non_cams(x, "ifrac()")
}

#' @rdname ifrac
#' @export
ifrac.interactmeta_conventional <- function(x, ...) {
  interactmeta_abort_non_cams(x, "ifrac()")
}

#' @export
plot.interactmeta_da <- function(x, y, ...) {
  interactmeta_plot_construction(x, ...)
}

#' @export
plot.interactmeta_ad <- function(x, y, ...) {
  interactmeta_plot_construction(x, ...)
}

#' @export
plot.interactmeta_conventional <- function(x, y, ...) {
  interactmeta_plot_construction(x, ...)
}

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

