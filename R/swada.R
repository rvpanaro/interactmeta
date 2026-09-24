#' Apply the Same Weighting Across Different Analyses principle
#'
#' **SWADA** (Same Weighting Across Different Analyses) is a *consistency
#' principle*, not a statistical model parallel to [da()], [ad()], or [cam()].
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
#' The interaction estimate itself is the AD/`rma` of \(g_j\); SWADA does not
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
#' arXiv:2508.15531. \doi{10.48550/arXiv.2508.15531}
#' @seealso [da()], [ad()], [cam()], [conventional()]
#' @export
swada <- function(data,
                  reference,
                  modifier = NULL,
                  trial = "trial",
                  subgroup = "subgroup",
                  yi = "yi",
                  vi = "vi",
                  weights = c("interaction"),
                  engine = "metafor",
                  method = NULL,
                  approach = NULL,
                  ...) {
  subgroup <- interactmeta_resolve_modifier_column(modifier, subgroup)
  weights <- match.arg(weights)
  dots <- list(...)
  dots$engine <- NULL
  resolved <- interactmeta_resolve_construction_engine_method(
    "swada",
    engine,
    method,
    approach = approach
  )
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
        engine = resolved$engine,
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
  interactmeta_attach_capabilities(out, "swada", resolved$engine, resolved$method)
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

#' @export
print.interactmeta_swada <- function(x, digits = 2, ...) {
  cat("SWADA principle (same weighting across DA and AD)\n")
  cat("  scheme    : ", x$weights_scheme, "\n", sep = "")
  cat("  engine    : ", x$engine, "\n", sep = "")
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

#' Side-by-side level-specific treatment effects from DA, SWADA, and CAMS
#'
#' `compare_modifiers()` is the name to use. `compare_subgroups()` remains
#' an alias. Display labels are `DA` ([da()] or the DA part of [conventional()]),
#' `SWADA` ([swada()]), and `CAMS` ([cam()]). These are not three equivalent
#' statistical models: DA is separate subgroup pooling, SWADA applies one
#' common weight vector, and CAMS is the contribution-adjusted model.
#' Average Difference ([ad()]) does not itself pool subgroup-specific
#' treatment effects, so it is omitted here.
#'
#' Warns if mixed constructions do not share the same `engine` and `method`.
#' A like-for-like frequentist comparison uses `metafor` / `REML` throughout
#' (`cam(..., engine = "metafor")`). Two CAMS fits of the same model with
#' different approaches are labelled by approach (`CAMS — bayesmeta` vs
#' `CAMS — metafor (REML)`) and are not treated as a mismatch.
#'
#' Two or more [cam()] fits may be passed instead of DA / SWADA / CAMS, in
#' which case [plot()] draws the compared subgroup estimates.
#'
#' @param da A [da()] or [conventional()] fit, or a [cam()] fit when comparing
#'   two CAMS approaches.
#' @param swada A [swada()] fit, or a second [cam()] fit.
#' @param cam A [cam()] CAMS fit. Optional when `da` and `swada`
#'   are already CAMS fits.
#' @param exponentiate If `TRUE` (default), report on the exponentiated
#'   (ratio) scale. Use `FALSE` for the fitted log scale.
#' @param type CAMS layer (`"Pooled"` default).
#' @return A tibble with one row per approach × subgroup and class
#'   `interactmeta_compare_subgroups` (plottable).
#' @export
compare_subgroups <- function(da,
                                        swada,
                                        cam = NULL,
                                        exponentiate = TRUE,
                                        type = "Pooled") {
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
      "`compare_subgroups()` needs DA, SWADA, and CAMS fits, ",
      "or two or more CAMS fits.",
      call. = FALSE
    )
  }
  interactmeta_warn_incomparable_fits(
    list(da, swada, cam),
    c("DA", "SWADA", "CAMS"),
    "compare_subgroups"
  )
  da_tab <- tidy(da, exponentiate = exponentiate)
  da_tab <- da_tab[da_tab$component == "subgroup", , drop = FALSE]
  da_tab$approach <- "DA"
  sw <- tidy(swada, exponentiate = exponentiate)
  sw <- sw[sw$component == "subgroup", , drop = FALSE]
  sw$approach <- "SWADA"
  sg <- subgroups(cam, type = type, exponentiate = exponentiate, summary = TRUE)
  td <- tidy(sg)
  td <- td[td$type == type, , drop = FALSE]
  cams_tab <- data.frame(
    subgroup = td$subgroup,
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

#' @rdname compare_subgroups
#' @export
compare_modifiers <- compare_subgroups

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
#' (`da()` / `ad()` / `swada()` defaults, plus `cam(..., engine = "metafor")`).
#' Two CAMS fits of the same model with different approaches are labelled by
#' approach (`CAMS — bayesmeta` vs `CAMS — metafor (REML)`) and are not treated
#' as a mismatch.
#'
#' Two or more [cam()] fits may be passed instead of DA / AD / SWADA / CAMS.
#' [plot()] draws the compared interaction estimates.
#'
#' @param da A [da()] or [conventional()] fit, or a [cam()] fit when comparing
#'   two CAMS approaches.
#' @param ad An [ad()] fit. Optional when `da` is a [conventional()] bundle, or
#'   a second [cam()] fit.
#' @param swada A [swada()] fit.
#' @param cam A [cam()] CAMS fit.
#' @inheritParams compare_subgroups
#' @param exponentiate If `TRUE`, report on the exponentiated (ratio) scale.
#'   Default `FALSE` (fitted log scale).
#' @return A tibble with class `interactmeta_compare_interactions`.
#' @export
compare_interactions <- function(da,
                                           ad = NULL,
                                           swada = NULL,
                                           cam = NULL,
                                           exponentiate = FALSE,
                                           type = "Pooled") {
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
#' uses one common vector; under the default this equals the AD weights.
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
#' @export
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
  interactmeta_compare_result(
    tibble::tibble(
      trial = as.character(w_da$trial),
      da_reference = as.numeric(w_da$weight_ref) / sum(w_da$weight_ref),
      da_nonref = as.numeric(w_da$weight_nonref) / sum(w_da$weight_nonref),
      ad = w_ad_vec / sum(w_ad_vec),
      swada = as.numeric(w_s$weight)
    ),
    "interactmeta_compare_weights",
    exponentiate = FALSE,
    xlab = "Weight"
  )
}

#' Trial-level CAMS quantities \eqn{\pi_j}, \eqn{m_j}, \eqn{g_j}
#'
#' Reads the transformed trial blocks stored on a [cam()] fit
#' (`fit$transformed_data$trial_blocks`). For binary CAMS,
#' \eqn{g_j = y_{B j} - y_{A j}} with `A` the reference subgroup.
#'
#' @param fit An `interactmeta` fit.
#' @return A data frame with `trial`, `pi`, `m`, `v_m`, `g`, and `v_g`.
#' @aliases cams_trial_quantities
#' @export
cam_trial_quantities <- function(fit) {
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

#' @rdname cam_trial_quantities
#' @export
cams_trial_quantities <- cam_trial_quantities

#' Comparable forest of level-specific treatment effects
#'
#' `plot_modifiers()` is the name to use. `plot_subgroups()` remains an alias.
#' Plots DA, SWADA, and CAMS pooled level-specific treatment effects on a
#' common axis. Observed trial estimates (shared input) can be included using
#' the original published intervals when those columns are present.
#'
#' @inheritParams compare_subgroups
#' @param include_observed If `TRUE`, add published trial rows from `da$data`.
#' @param xlim Shared x-axis limits on the displayed scale.
#' @return A ggplot object.
#' @export
plot_subgroups <- function(da,
                                     swada,
                                     cam,
                                     exponentiate = TRUE,
                                     type = "Pooled",
                                     include_observed = TRUE,
                                     xlim = NULL) {
  da_fit <- interactmeta_as_da(da)
  tab <- compare_subgroups(
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
    obs_tab <- data.frame(
      approach = "Published",
      subgroup = as.character(obs$subgroup),
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

#' @rdname plot_subgroups
#' @export
plot_modifiers <- plot_subgroups

#' @rdname ifrac
#' @export
ifrac.interactmeta_swada <- function(x, ...) {
  interactmeta_abort_non_cams(x, "ifrac()")
}

#' @export
plot.interactmeta_swada <- function(x, y, ...) {
  interactmeta_plot_construction(x, ...)
}

