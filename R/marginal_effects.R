#' Effects at each level of the modifier
#'
#' `subgroups()` is the name to use. `modifiers()` is another name for the same extractor.
#' The examples are subgroup analyses: each study contributes two estimates.
#' The same extractor applies to any within-study contrast stored that way
#' (dose, biomarker, time point). Level labels stay as stored. Display
#' headings capitalize the first letter of the `modifier` column name.
#'
#' Evaluates CAMS level-specific treatment effects at the pooled
#' (precision-weighted) information fraction from [ifrac()]. This is **not**
#' a marginal treatment effect in the survey-sampling sense: it is the
#' subgroup-specific linear predictor at the combined-population contribution
#' structure \eqn{\widetilde{\pi}}.
#'
#' For each posterior draw, the linear predictor is evaluated at that draw's
#' combined-population information fraction with the appropriate subgroup
#' fixed effect:
#' \deqn{\mu_A(\pi^\star)=\alpha+\pi^\star(\theta-\gamma),}
#' \deqn{\mu_B(\pi^\star)=\alpha+\pi^\star(\theta-\gamma)+\gamma.}
#' On the orthogonal CAMS factorization, \eqn{\theta} is the overall-block IF slope and
#' stored `b_ifrac` is \eqn{\beta=\theta-\gamma}. This is **not** the
#' trial-specific convolution on `fit$latent$mu_A[[j]]` /
#' `fit$latent$mu_B[[j]]` (\eqn{\mu_{Aj}=m_j-\pi_j g_j} at that study's
#' \eqn{\pi_j}). A scalar `contribution` is a plug-in \eqn{\pi^\star};
#' `NULL` mixes the [ifrac()] posterior.
#'
#' \eqn{\alpha} and \eqn{\theta} share the overall-block support. Population
#' \eqn{M^\star=\alpha+\pi^\star\theta} is formed from that joint, not by
#' convolving the two mixture marginals.
#'
#' @param fit A [`interactmeta`] from [cam()] / [cam()].
#' @param contribution Combined-population information fraction \eqn{\pi}
#'   to evaluate at. If `NULL` (default), uses the posterior draws from
#'   [ifrac()] for `fit`. A numeric vector is aligned to posterior draws;
#'   a scalar (e.g. `contribution = 0.5`) is a plug-in at every draw.
#' @param ifrac When `contribution` is `NULL`, passed to [ifrac()]: `"interaction"`
#'   (default; interaction-weighted pooling of trial-level \eqn{\pi_j}) or
#'   `"overall"` for overall-block precision weights.
#' @param type One or both of `"Pooled"` and `"Prediction"`. `"Pooled"` uses fixed
#'   effects only; `"Prediction"` adds a new-trial random effect drawn from the
#'   posterior trial SDs at the combined-population IF (effect in a new **study**).
#' @param prob Coverage probability for `lci` / `uci` when `summary = TRUE`.
#' @param max_draws Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports (weighted quadrature nodes); omit it.
#' @param summary If `TRUE` (default), include a summary data frame.
#' @param contrast_sep Rows whose `subgroup` label contains this separator are
#'   excluded when resolving subgroup levels (default `"-"`).
#' @param exponentiate If `TRUE`, transform posterior draws with `exp()` before
#'   summarizing (e.g. log rate ratios to ratios).
#' @param shrinkage If `TRUE`, delegate to [shrinkage()] with
#'   `type = "subgroups"` (trial-specific reconstruction at each study's
#'   \eqn{\pi_j}). If `FALSE` (default), return only Pooled / Prediction at the
#'   combined-population information fraction.
#'
#' @return A `interactmeta_subgroup_effects` object (list with `draws`, `summary`,
#'   `overall_ifrac`, `overall_ifrac_source` (`"fit"` or `"user"`), `reference`,
#'   `prob`, `exponentiate`). Draws are indexed by `subgroup` and `type`.
#'   `summary` includes `estimate` (median), `mean`, `sd`, `lci`, and `uci`.
#'   For Bayesian CAMS, `lci`/`uci` are equal-tail mixture quantiles from the
#'   DIRECT-node table; native `$post.interval` / `convolve()` stay on `$latent`.
#'
#' @seealso [combined()], [interactions()], [shrinkage()], [ifrac()],
#'   [forestplot.interactmeta()]
#' @export
subgroups <- function(fit,
                             contribution = NULL,
                             ifrac = c("interaction", "overall"),
                             type = c("Pooled", "Prediction"),
                             prob = 0.95,
                             max_draws = NULL,
                             summary = TRUE,
                             contrast_sep = "-",
                             exponentiate = FALSE,
                             shrinkage = FALSE) {
  st_mod <- attr(fit, "interactmeta", exact = TRUE)
  modifier_name <- st_mod$modifier %||% fit$modifier %||% st_mod$subgroup_var %||% "subgroup"
  if (interactmeta_is_non_cams_analysis(fit)) {
    out <- interactmeta_extract_subgroups(fit, exponentiate = exponentiate)
    if (is.list(out)) {
      out$modifier <- modifier_name
    }
    return(out)
  }
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cam()`.", call. = FALSE)
  }
  if (isTRUE(shrinkage)) {
    out <- shrinkage(
      fit,
      type = "subgroups",
      prob = prob,
      max_draws = max_draws,
      summary = summary,
      contrast_sep = contrast_sep,
      ifrac = match.arg(ifrac),
      contribution = contribution,
      exponentiate = exponentiate
    )
    if (is.list(out)) {
      out$modifier <- modifier_name
    }
    return(out)
  }
  type <- match.arg(type, c("Pooled", "Prediction"), several.ok = TRUE)
  ifrac <- match.arg(ifrac)

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = contribution,
    ifrac = ifrac,
    max_draws = max_draws,
    contrast_sep = contrast_sep
  )

  rows <- lapply(type, function(ty) {
    lapply(prep$design$subgroup_levels, function(lbl) {
      eff <- comparison_marginal_effect_draws(
        post = prep$post,
        omega = prep$omega,
        design = prep$design,
        subgroup_label = lbl,
        subgroup_coef = prep$sg_coef[[lbl]],
        type = ty,
        fit = fit
      )
      cond_sd <- NULL
      if (interactmeta_has_cond_moments(prep$post)) {
        x <- as.numeric(!identical(as.character(lbl), as.character(prep$design$reference)))
        a_gamma <- ifelse(x == 0, -prep$omega, 1 - prep$omega)
        cond_sd <- interactmeta_lincomb_cond_sd(prep$post, 1, prep$omega, a_gamma)
        if (identical(ty, "Prediction")) {
          re_sd <- interactmeta_prediction_re_sd(prep$post, prep$omega, x = x)
          cond_sd <- sqrt(cond_sd^2 + re_sd^2)
        }
      }
      interactmeta_decorate_effect_draws(
        data.frame(
          .draw = prep$draw_ids$.draw,
          .chain = prep$draw_ids$.chain,
          .iteration = prep$draw_ids$.iteration,
          subgroup = lbl,
          type = ty,
          effect = eff,
          stringsAsFactors = FALSE
        ),
        prep$post,
        cond_sd = cond_sd
      )
    })
  })
  draws <- dplyr::bind_rows(unlist(rows, recursive = FALSE))

  out <- list(
    draws = draws,
    overall_ifrac = prep$overall_ifrac,
    overall_ifrac_source = prep$overall_ifrac_source,
    ifrac = prep$ifrac,
    reference = prep$design$reference,
    prob = prob,
    exponentiate = isTRUE(exponentiate),
    shrinkage = FALSE,
    modifier = modifier_name
  )
  if (isTRUE(summary)) {
    fin <- interactmeta_summarize_effects(
      draws = draws,
      group_cols = c("subgroup", "type"),
      prob = prob,
      exponentiate = exponentiate,
      fit = fit,
      quantity = "subgroups",
      omega = interactmeta_metafor_omega(fit, prep$overall_ifrac)
    )
    out$draws <- fin$draws
    out$summary <- fin$summary
  } else {
    out$draws <- comparison_transform_effect_draws(draws, exponentiate)
  }
  class(out) <- c("interactmeta_subgroup_effects", class(out))
  interactmeta_attach_origin(out, fit, quantity = "subgroups")
}

#' @rdname subgroups
#' @export
modifiers <- subgroups

#' Combined-population treatment effect
#'
#' Evaluates the CAMS **combined-population treatment effect** at the pooled
#' information fraction from [ifrac()]. This is the effect corresponding to
#' the combined population at its contribution-weighted subgroup composition
#' \eqn{\widetilde{\pi}}, using the linear predictor
#' \eqn{\alpha + \widetilde{\pi}\,(\beta + \gamma)} (binary case; the
#' implementation sums all non-reference IF and subgroup coefficients).
#'
#' When subgroup variances scale with sample size (or events) and allocation is
#' similar, \eqn{\pi} is close to the trial's subgroup share, so this is
#' typically a trial-mixture effect. It is still **not** automatically a
#' target-population average: the pooled \eqn{\widetilde{\pi}} is
#' precision-weighted across trials, and on log HR / log OR a mixture of
#' subgroup effects is not collapsible. Trial-specific combined-population
#' effects are the transformed quantities
#' \eqn{m_j = (1-\pi_j)y_{Aj}+\pi_j y_{Bj}}.
#'
#' Prefer [combined()] in new code. [overall()] is a backward-compatible alias.
#'
#' @param fit A [`interactmeta`] from [cam()] / [cam()].
#' @param contribution Combined-population information fraction \eqn{\pi}
#'   to evaluate at; see [subgroups()].
#' @param ifrac When `contribution` is `NULL`, passed to [ifrac()]: `"interaction"`
#'   (default) or `"overall"` for overall-block pooling of trial
#'   information fractions.
#' @param type One or both of `"Pooled"` and `"Prediction"` (new-trial random
#'   effect; see [subgroups()]).
#' @param prob,max_draws,summary,contrast_sep,exponentiate,shrinkage As in
#'   [subgroups()].
#'
#' @return A `interactmeta_overall_effect` object (list with `draws`, `summary`,
#'   `overall_ifrac`, `overall_ifrac_source`, `ifrac`, `reference`, `prob`).
#'   Also inherits class `interactmeta_combined_effect`. Draws are indexed by
#'   `type`. `summary` includes `estimate` (median), `mean`, `sd`, `lci`, and
#'   `uci`.
#'
#' @seealso [subgroups()], [interactions()], [ifrac()]
#' @export
combined <- function(fit,
                     contribution = NULL,
                     ifrac = c("interaction", "overall"),
                     type = c("Pooled", "Prediction"),
                     prob = 0.95,
                     max_draws = NULL,
                     summary = TRUE,
                     contrast_sep = "-",
                     exponentiate = FALSE,
                     shrinkage = FALSE) {
  interactmeta_abort_non_cams(fit, "combined()")
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cam()`.", call. = FALSE)
  }
  if (isTRUE(shrinkage)) {
    return(shrinkage(
      fit,
      type = "overall",
      prob = prob,
      max_draws = max_draws,
      summary = summary,
      contrast_sep = contrast_sep,
      ifrac = match.arg(ifrac),
      contribution = contribution,
      exponentiate = exponentiate
    ))
  }
  type <- match.arg(type, c("Pooled", "Prediction"), several.ok = TRUE)
  ifrac <- match.arg(ifrac)

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = contribution,
    ifrac = ifrac,
    max_draws = max_draws,
    contrast_sep = contrast_sep
  )

  rows <- lapply(type, function(ty) {
    eff <- comparison_marginal_overall_effect_draws(
      post = prep$post,
      omega = prep$omega,
      design = prep$design,
      fit = fit,
      type = ty
    )
    cond_sd <- NULL
    if (interactmeta_has_cond_moments(prep$post)) {
      cond_sd <- interactmeta_lincomb_cond_sd(prep$post, 1, prep$omega, 0)
      if (identical(ty, "Prediction")) {
        re_sd <- interactmeta_prediction_re_sd(prep$post, prep$omega, x = NULL)
        cond_sd <- sqrt(cond_sd^2 + re_sd^2)
      }
    }
    interactmeta_decorate_effect_draws(
      data.frame(
        .draw = prep$draw_ids$.draw,
        .chain = prep$draw_ids$.chain,
        .iteration = prep$draw_ids$.iteration,
        type = ty,
        effect = eff,
        stringsAsFactors = FALSE
      ),
      prep$post,
      cond_sd = cond_sd
    )
  })
  draws <- dplyr::bind_rows(rows)

  out <- list(
    draws = draws,
    overall_ifrac = prep$overall_ifrac,
    overall_ifrac_source = prep$overall_ifrac_source,
    ifrac = prep$ifrac,
    reference = prep$design$reference,
    prob = prob,
    exponentiate = isTRUE(exponentiate),
    shrinkage = FALSE
  )
  if (isTRUE(summary)) {
    fin <- interactmeta_summarize_effects(
      draws = draws,
      group_cols = "type",
      prob = prob,
      exponentiate = exponentiate,
      fit = fit,
      quantity = "combined",
      omega = interactmeta_metafor_omega(fit, prep$overall_ifrac)
    )
    out$draws <- fin$draws
    out$summary <- fin$summary
  } else {
    out$draws <- comparison_transform_effect_draws(draws, exponentiate)
  }
  class(out) <- c("interactmeta_combined_effect", "interactmeta_overall_effect", class(out))
  interactmeta_attach_origin(out, fit, quantity = "combined")
}

#' @rdname combined
#' @export
overall <- function(fit,
                    contribution = NULL,
                    ifrac = c("interaction", "overall"),
                    type = c("Pooled", "Prediction"),
                    prob = 0.95,
                    max_draws = NULL,
                    summary = TRUE,
                    contrast_sep = "-",
                    exponentiate = FALSE,
                    shrinkage = FALSE) {
  combined(
    fit = fit,
    contribution = contribution,
    ifrac = ifrac,
    type = type,
    prob = prob,
    max_draws = max_draws,
    summary = summary,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate,
    shrinkage = shrinkage
  )
}

#' @noRd
interactmeta_subgroup_effects_shrinkage <- function(fit,
                                            prob,
                                            max_draws,
                                            summary,
                                            contrast_sep,
                                            exponentiate,
                                            overall_ifrac,
                                            ifrac = c("interaction", "overall"),
                                            bundle = NULL) {
  ifrac <- match.arg(ifrac)
  bundle <- bundle %||% comparison_shrinkage_rotation_bundle(
    fit,
    max_draws = max_draws,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    overall_ifrac = overall_ifrac
  )
  comparison_shrinkage_effect_from_bundle(
    bundle = bundle,
    type = "subgroups",
    prob = prob,
    summary = summary,
    exponentiate = exponentiate
  )
}

#' @noRd
interactmeta_overall_effects_shrinkage <- function(fit,
                                           prob,
                                           max_draws,
                                           summary,
                                           contrast_sep,
                                           exponentiate,
                                           overall_ifrac,
                                           ifrac = c("interaction", "overall"),
                                           bundle = NULL) {
  ifrac <- match.arg(ifrac)
  bundle <- bundle %||% comparison_shrinkage_rotation_bundle(
    fit,
    max_draws = max_draws,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    overall_ifrac = overall_ifrac
  )
  comparison_shrinkage_effect_from_bundle(
    bundle = bundle,
    type = "overall",
    prob = prob,
    summary = summary,
    exponentiate = exponentiate
  )
}

#' Build a typed shrinkage effect object from one rotation bundle
#' @noRd
comparison_shrinkage_effect_from_bundle <- function(bundle,
                                              type = c("subgroups", "interaction", "overall", "combined"),
                                              prob = 0.95,
                                              summary = TRUE,
                                              exponentiate = FALSE) {
  type <- comparison_resolve_effect_type(match.arg(type))
  if (identical(type, "subgroups")) {
    draws <- comparison_shrinkage_draws_from_bundle(bundle)
    fin <- comparison_shrinkage_finalize_draws(
      draws = draws,
      group_cols = c("trial", "subgroup", "type"),
      prob = prob,
      summary = summary,
      exponentiate = exponentiate
    )
    out <- list(
      draws = fin$draws,
      overall_ifrac = bundle$overall_ifrac,
      overall_ifrac_source = bundle$overall_ifrac_source,
      ifrac = bundle$ifrac,
      reference = bundle$reference,
      prob = prob,
      exponentiate = isTRUE(exponentiate),
      shrinkage = TRUE,
      components = comparison_shrinkage_components_list(bundle)
    )
    out$summary <- interactmeta_bayesmeta_maybe_fill_shrinkage(
      bundle, fin$summary, type = "subgroups", prob = prob, exponentiate = exponentiate
    )
    class(out) <- c("interactmeta_subgroup_effects", class(out))
    return(out)
  }
  if (identical(type, "interaction")) {
    draws <- comparison_rotation_matrix_to_long(
      mat = bundle$g,
      post = bundle$post,
      trials = bundle$trials,
      extra = list(interaction = bundle$interaction_label),
      type = "Shrinkage",
      sd_mat = bundle$sd_g
    )
    fin <- comparison_shrinkage_finalize_draws(
      draws = draws,
      group_cols = c("trial", "interaction", "type"),
      prob = prob,
      summary = summary,
      exponentiate = exponentiate
    )
    out <- list(
      draws = fin$draws,
      overall_ifrac = bundle$overall_ifrac,
      overall_ifrac_source = bundle$overall_ifrac_source,
      ifrac = bundle$ifrac,
      reference = bundle$reference,
      prob = prob,
      exponentiate = isTRUE(exponentiate),
      shrinkage = TRUE,
      components = comparison_shrinkage_components_list(bundle)
    )
    out$summary <- interactmeta_bayesmeta_maybe_fill_shrinkage(
      bundle, fin$summary, type = "interaction", prob = prob, exponentiate = exponentiate
    )
    class(out) <- c("interactmeta_interaction_effects", class(out))
    return(out)
  }
  # overall
  draws <- comparison_rotation_matrix_to_long(
    mat = bundle$m,
    post = bundle$post,
    trials = bundle$trials,
    type = "Shrinkage",
    sd_mat = bundle$sd_m
  )
  fin <- comparison_shrinkage_finalize_draws(
    draws = draws,
    group_cols = c("trial", "type"),
    prob = prob,
    summary = summary,
    exponentiate = exponentiate
  )
  out <- list(
    draws = fin$draws,
    overall_ifrac = bundle$overall_ifrac,
    overall_ifrac_source = bundle$overall_ifrac_source,
    ifrac = bundle$ifrac,
    reference = bundle$reference,
    prob = prob,
    exponentiate = isTRUE(exponentiate),
    shrinkage = TRUE,
    components = comparison_shrinkage_components_list(bundle)
  )
  out$summary <- interactmeta_bayesmeta_maybe_fill_shrinkage(
    bundle, fin$summary, type = "overall", prob = prob, exponentiate = exponentiate
  )
  class(out) <- c("interactmeta_overall_effect", class(out))
  out
}

#' Public shrinkage summaries keep the DIRECT-mixture equal-tail intervals.
#' Native `$latent` convolve / shrink intervals stay available on the fit.
#' @noRd
interactmeta_bayesmeta_maybe_fill_shrinkage <- function(bundle,
                                                       summary,
                                                       type,
                                                       prob,
                                                       exponentiate) {
  summary
}

#' Summarize shrinkage draws; exponentiate only after reconstruction
#'
#' Mixture / weighted nodes (bayesmeta DIRECT, metafor BLUP SEs) must be
#' summarized on the model scale, then optionally transformed. MCMC draws
#' keep the previous transform-then-summarize path (median/QI commute with
#' `exp()`).
#' @noRd
comparison_shrinkage_finalize_draws <- function(draws,
                                          group_cols,
                                          prob,
                                          summary,
                                          exponentiate) {
  mix <- comparison_draws_use_mixture_moments(draws)
  summ <- NULL
  if (isTRUE(summary)) {
    if (mix) {
      summ <- comparison_effects_summary(
        draws,
        group_cols = group_cols,
        prob = prob,
        exponentiate = exponentiate
      )
    } else {
      transformed <- comparison_transform_effect_draws(draws, exponentiate)
      summ <- comparison_effects_summary(
        transformed,
        group_cols = group_cols,
        prob = prob
      )
    }
  }
  list(
    draws = comparison_transform_effect_draws(draws, exponentiate),
    summary = summ
  )
}

#' @noRd
comparison_draws_use_mixture_moments <- function(draws) {
  if (!is.data.frame(draws) || !nrow(draws)) {
    return(FALSE)
  }
  has_w <- ".weight" %in% names(draws) &&
    length(unique(as.numeric(draws$.weight)[is.finite(as.numeric(draws$.weight))])) > 1L
  has_s <- ".cond_sd" %in% names(draws) &&
    any(as.numeric(draws$.cond_sd) > 0, na.rm = TRUE)
  has_w || has_s
}

#' Trial-level overall from subgroup draws (legacy / tests)
#'
#' Uses trial IF `pi` when supplied; otherwise a shared `omega` vector.
#' @noRd
comparison_overall_shrinkage_draws <- function(sg_draws,
                                         omega,
                                         reference,
                                         data,
                                         trial_var = "trial",
                                         subgroup_var = "subgroup",
                                         pi = NULL) {
  ref <- as.character(reference)[[1L]]
  trials <- unique(as.character(sg_draws$trial))
  subgroups <- unique(as.character(sg_draws$subgroup))
  nonrefs <- setdiff(subgroups, ref)
  if (length(nonrefs) != 1L) {
    stop(
      "`comparison_overall_shrinkage_draws()` supports exactly one non-reference subgroup.",
      call. = FALSE
    )
  }
  ref_template <- sg_draws[
    sg_draws$trial == trials[[1L]] & sg_draws$subgroup == ref,
    ,
    drop = FALSE
  ]
  n_draw <- nrow(ref_template)
  if (is.null(pi)) {
    omega <- as.numeric(omega)
    if (length(omega) == 1L) {
      omega <- rep(omega, n_draw)
    }
    pi_list <- stats::setNames(rep(list(omega), length(trials)), trials)
    per_draw <- TRUE
  } else {
    pi <- as.numeric(pi)
    if (is.null(names(pi))) {
      names(pi) <- trials
    }
    per_draw <- FALSE
  }

  rows <- lapply(seq_along(trials), function(i) {
    tr <- trials[[i]]
    rref <- sg_draws[
      sg_draws$trial == tr & sg_draws$subgroup == ref,
      ,
      drop = FALSE
    ]
    non <- sg_draws[
      sg_draws$trial == tr & sg_draws$subgroup == nonrefs[[1L]],
      ,
      drop = FALSE
    ]
    pij <- if (isTRUE(per_draw)) pi_list[[tr]] else rep(pi[[tr]], n_draw)
    data.frame(
      .draw = rref$.draw,
      .chain = rref$.chain,
      .iteration = rref$.iteration,
      row = i,
      trial = tr,
      type = "Shrinkage",
      effect = (1 - pij) * rref$effect + pij * non$effect,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' @export
print.interactmeta_subgroup_effects <- function(x, digits = 2, ...) {
  opts <- comparison_print_digits_prob(x, digits = digits, ...)
  heading <- interactmeta_modifier_heading(x$modifier %||% "subgroup")
  title <- if (isTRUE(x$shrinkage)) {
    paste0(heading, " shrinkage effects (study-specific)")
  } else if (identical(x$ifrac, "interaction")) {
    paste0(heading, " effects at interaction-weighted information fraction")
  } else {
    paste0(heading, " effects at overall information fraction")
  }
  if (!is.null(x$summary)) {
    comparison_print_effects_table(
      x$summary,
      digits = opts$digits,
      prob = opts$prob,
      id_cols = c("trial", "subgroup", "type"),
      title = title,
      style = interactmeta_interval_style(x)
    )
  } else {
    cat(title, " (", nrow(x$draws), " rows)\n", sep = "")
  }
  invisible(x)
}

#' @export
print.interactmeta_overall_effect <- function(x, digits = 2, ...) {
  opts <- comparison_print_digits_prob(x, digits = digits, ...)
  title <- if (isTRUE(x$shrinkage)) {
    "Combined-population shrinkage effects (study-specific)"
  } else {
    "Combined-population treatment effect at the combined-population information fraction"
  }
  if (!is.null(x$summary)) {
    comparison_print_effects_table(
      x$summary,
      digits = opts$digits,
      prob = opts$prob,
      id_cols = c("trial", "type"),
      title = title,
      style = interactmeta_interval_style(x)
    )
  } else {
    cat(title, " (", nrow(x$draws), " rows)\n", sep = "")
  }
  invisible(x)
}

#' @noRd
comparison_transform_effect_draws <- function(draws, exponentiate) {
  if (!isTRUE(exponentiate)) {
    return(draws)
  }
  draws$effect <- exp(draws$effect)
  draws
}

#' @noRd
comparison_prepare_marginal_draws <- function(fit,
                                        overall_ifrac = NULL,
                                        ifrac = c("interaction", "overall"),
                                        max_draws = NULL,
                                        contrast_sep = "-") {
  resolved <- comparison_resolve_marginal_omega(
    fit,
    overall_ifrac = overall_ifrac,
    ifrac = ifrac,
    max_draws = max_draws
  )
  design <- comparison_marginal_design(fit, contrast_sep = contrast_sep)
  list(
    omega = resolved$omega,
    post = resolved$post,
    draw_ids = comparison_posterior_draw_ids(resolved$post),
    design = design,
    sg_coef = comparison_subgroup_coef_terms(fit, design),
    overall_ifrac = resolved$overall_ifrac,
    overall_ifrac_source = resolved$source,
    ifrac = resolved$ifrac
  )
}

#' @noRd
comparison_resolve_marginal_omega <- function(fit,
                                        overall_ifrac = NULL,
                                        ifrac = c("interaction", "overall"),
                                        max_draws = NULL) {
  ifrac <- match.arg(ifrac)
  source <- "user"
  oi <- overall_ifrac
  if (is.null(oi)) {
    oi <- interactmeta_eval_ifrac(fit, ifrac = ifrac)
    source <- "fit"
  }
  omega <- as.numeric(oi)
  if (!length(omega) || !all(is.finite(omega))) {
    stop(
      "`contribution` must be `NULL` or a finite numeric vector (or scalar).",
      call. = FALSE
    )
  }

  post_full <- comparison_posterior_df(fit)
  n_post <- nrow(post_full)
  w_post <- if (".weight" %in% names(post_full)) post_full$.weight else NULL

  if (length(omega) == 1L) {
    idx <- comparison_marginal_draw_index(
      n_post,
      max_draws = max_draws,
      weight = w_post
    )
    post <- post_full[idx, , drop = FALSE]
    omega <- rep(omega, nrow(post))
  } else if (length(omega) == n_post) {
    idx <- comparison_marginal_draw_index(
      n_post,
      max_draws = max_draws,
      weight = w_post
    )
    omega <- omega[idx]
    post <- post_full[idx, , drop = FALSE]
  } else {
    idx <- comparison_marginal_draw_index(
      length(omega),
      max_draws = max_draws,
      weight = NULL
    )
    if (length(idx) > n_post) {
      stop(
        "`contribution` has more values than posterior draws in `fit`.",
        call. = FALSE
      )
    }
    if (!is.null(w_post) && length(idx) < n_post &&
          interactmeta_has_direct_weights(w_post, n_post)) {
      stop(
        "`max_draws` cannot thin a weighted DIRECT support. ",
        "Bayesian CAMS summaries use the full quadrature grid; omit `max_draws`.",
        call. = FALSE
      )
    }
    omega <- omega[idx]
    post <- post_full[idx, , drop = FALSE]
  }

  oi_ifrac <- if (inherits(oi, "interactmeta_overall_ifrac")) {
    attr(oi, "ifrac", exact = TRUE)
  } else {
    NA_character_
  }
  if (is.null(oi_ifrac) || !nzchar(oi_ifrac)) {
    oi_ifrac <- if (identical(source, "fit")) ifrac else NA_character_
  }

  list(
    omega = omega,
    post = post,
    overall_ifrac = oi,
    source = source,
    ifrac = oi_ifrac
  )
}

#' Design metadata for marginal CAMS effects
#' @noRd
comparison_marginal_design <- function(fit, contrast_sep = "-") {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st) || is.null(st$ifrac_data)) {
    stop(
      "Fit has no `ifrac_data` in attribute `\"interactmeta\"`; refit with `cam()`.",
      call. = FALSE
    )
  }
  subgroup_var <- if (is.null(st$subgroup_var)) "subgroup" else st$subgroup_var
  dat <- st$ifrac_data
  dat <- dat[
    is_ifrac_component_row(dat[[subgroup_var]], contrast_sep = contrast_sep),
    ,
    drop = FALSE
  ]
  subgroup_levels <- sort(unique(as.character(dat[[subgroup_var]])))
  if (!length(subgroup_levels)) {
    stop("No subgroup rows found for marginal effects.", call. = FALSE)
  }

  ifrac_vars_std <- comparison_ifrac_coef_names(fit)
  if (!length(ifrac_vars_std)) {
    stop("No fixed-effect `ifrac_*` coefficients found in the model.", call. = FALSE)
  }

  list(
    reference = st$reference,
    subgroup_levels = subgroup_levels,
    ifrac_vars_std = ifrac_vars_std,
    subgroup_var = subgroup_var
  )
}

#' @noRd
comparison_ifrac_coef_names <- function(fit) {
  vars <- comparison_fit_variables(fit)
  hits <- vars[grepl("^b_ifrac", vars) & !grepl("ifrac2", vars)]
  sub("^b_", "", hits)
}

#' @noRd
comparison_subgroup_coef_terms <- function(fit, design) {
  vars <- comparison_fit_variables(fit)
  sub_coefs <- vars[grepl("^b_subgroup", vars)]
  out <- stats::setNames(
    rep(NA_character_, length(design$subgroup_levels)),
    design$subgroup_levels
  )
  if (!length(sub_coefs)) {
    return(out)
  }
  ref <- design$reference
  for (lbl in design$subgroup_levels) {
    if (identical(lbl, ref)) {
      next
    }
    pat <- paste0("^b_subgroup", lbl, "$")
    hit <- sub_coefs[grepl(pat, sub_coefs)]
    if (length(hit) == 1L) {
      out[[lbl]] <- sub("^b_", "", hit[[1]])
    }
  }
  out
}

#' @noRd
comparison_posterior_df <- function(fit) {
  post <- comparison_tidy_draws(fit)
  raw <- as.data.frame(fit$draws)
  keep <- grep("^\\.(weight|cond_)", names(raw), value = TRUE)
  for (nm in keep) {
    if (!nm %in% names(post) && nrow(raw) == nrow(post)) {
      post[[nm]] <- raw[[nm]]
    }
  }
  post
}

#' @noRd
comparison_posterior_draw_ids <- function(post) {
  data.frame(
    .draw = post$.draw,
    .chain = post$.chain,
    .iteration = post$.iteration,
    stringsAsFactors = FALSE
  )
}

#' Index rows of a draws table, optionally capped by `max_draws`
#'
#' Weighted DIRECT supports (Bayesian CAMS) must not be thinned uniformly:
#' that over-represents low-weight quadrature nodes. Omit `max_draws` for
#' those fits.
#'
#' @param n Number of rows.
#' @param max_draws Optional maximum number of rows. Ignored when `NULL` or
#'   at least `n`.
#' @param weight Optional DIRECT node weights aligned with the `n` rows.
#' @noRd
comparison_marginal_draw_index <- function(n, max_draws = NULL, weight = NULL) {
  n <- as.integer(n)[[1L]]
  if (!is.finite(n) || n < 0L) {
    stop("`n` must be a non-negative integer.", call. = FALSE)
  }
  if (is.null(max_draws) || !is.finite(max_draws) || max_draws >= n) {
    return(seq_len(n))
  }
  if (interactmeta_has_direct_weights(weight, n)) {
    stop(
      "`max_draws` cannot thin a weighted DIRECT support. ",
      "Bayesian CAMS summaries use the full quadrature grid; omit `max_draws`.",
      call. = FALSE
    )
  }
  set.seed(1L)
  sample.int(n, as.integer(max_draws))
}

#' TRUE when `weight` is a usable DIRECT-node weight vector of length `n`
#' @noRd
interactmeta_has_direct_weights <- function(weight, n) {
  if (is.null(weight) || length(weight) != as.integer(n)[[1L]]) {
    return(FALSE)
  }
  w <- as.numeric(weight)
  any(is.finite(w) & w > 0)
}

#' Overall-block IF slope \eqn{\theta = \beta_{\mathrm{pkg}} + \gamma}
#'
#' Stored `b_ifrac_*` is the common subgroup slope
#' \eqn{\beta_{\mathrm{pkg}}=\theta-\gamma}, not the slope of \eqn{M(\pi)}.
#' This composite is \eqn{\theta}, so \eqn{M(\pi)=\alpha+\theta\pi}.
#' @noRd
comparison_marginal_overall_composite <- function(post, design, fit) {
  composite <- rep(0, nrow(post))
  for (v in design$ifrac_vars_std) {
    col <- paste0("b_", v)
    if (!col %in% names(post)) {
      stop(sprintf("Missing posterior column `%s`.", col), call. = FALSE)
    }
    composite <- composite + post[[col]]
  }
  sg_coef <- comparison_subgroup_coef_terms(fit, design)
  for (lbl in design$subgroup_levels) {
    if (identical(lbl, design$reference)) {
      next
    }
    coef <- sg_coef[[lbl]]
    if (!is.na(coef) && nzchar(coef)) {
      col <- paste0("b_", coef)
      if (!col %in% names(post)) {
        stop(sprintf("Missing posterior column `%s`.", col), call. = FALSE)
      }
      composite <- composite + post[[col]]
    }
  }
  composite
}

#' Posterior draws of \eqn{M(\pi)=\alpha+\theta\pi} with \eqn{\theta=b_{\mathrm{ifrac}}+\gamma}
#' @noRd
comparison_marginal_overall_effect_draws <- function(post,
                                               omega,
                                               design,
                                               fit,
                                               type = "Pooled") {
  composite <- comparison_marginal_overall_composite(post, design, fit)
  lp <- post[["b_Intercept"]] + omega * composite
  if (identical(type, "Prediction")) {
    lp <- lp + comparison_new_trial_re_draws(fit, omega, post = post)
  }
  lp
}

#' Canonical CAMS linear predictor at a chosen information fraction
#'
#' \deqn{\mu=\alpha+\mathrm{RE}
#'   +\beta\cdot\mathrm{evaluation\_if}
#'   +(x-\mathrm{evaluation\_if})\,v
#'   +\gamma\,1(\mathrm{nonref}).}
#'
#' @param evaluation_if Scalar or length-`nrow(post)` IF at which to evaluate
#'   (trial \eqn{\pi_j} for shrinkage; overall \eqn{\omega} for pooled /
#'   prediction / standardized trial effects).
#' @param re `"zero"` (Pooled), `"new"` (Prediction), or `"stored"` (Shrinkage /
#'   Standardized).
#' @param trial Required when `re = "stored"`.
#' @noRd
comparison_effect_linpred_draws <- function(post,
                                      evaluation_if,
                                      design,
                                      subgroup_label,
                                      subgroup_coef,
                                      re = c("zero", "new", "stored"),
                                      fit = NULL,
                                      trial = NULL) {
  re <- match.arg(re)
  evaluation_if <- as.numeric(evaluation_if)
  n <- nrow(post)
  if (length(evaluation_if) == 1L) {
    evaluation_if <- rep(evaluation_if, n)
  }
  if (length(evaluation_if) != n) {
    stop(
      "`evaluation_if` must be scalar or length nrow(post) (",
      n,
      ").",
      call. = FALSE
    )
  }

  lp <- comparison_marginal_fe_draws(
    post = post,
    evaluation_if = evaluation_if,
    design = design,
    subgroup_coef = subgroup_coef
  )
  ref <- as.character(design$reference)[[1L]]
  x <- as.numeric(!identical(as.character(subgroup_label)[[1L]], ref))

  if (identical(re, "zero")) {
    return(lp)
  }
  if (identical(re, "new")) {
    if (is.null(fit)) {
      stop("`re = \"new\"` requires `fit`.", call. = FALSE)
    }
    return(lp + comparison_new_trial_re_draws(fit, evaluation_if, post = post, x = x))
  }
  # stored
  if (is.null(fit) || is.null(trial) || !nzchar(as.character(trial)[[1L]])) {
    stop("`re = \"stored\"` requires `fit` and `trial`.", call. = FALSE)
  }
  lp + comparison_stored_trial_re(
    post = post,
    fit = fit,
    trial = as.character(trial)[[1L]],
    evaluation_if = evaluation_if,
    x = x
  )
}

#' Posterior draws of the marginalized effect (linear predictor at overall IF)
#'
#' Thin wrapper around `comparison_effect_linpred_draws()` for Pooled / Prediction /
#' Standardized (ω + stored RE). Study-specific shrinkage uses
#' `comparison_shrinkage_rotation_bundle()` with `evaluation = "trial"` instead.
#' @noRd
comparison_marginal_effect_draws <- function(post,
                                       omega,
                                       design,
                                       subgroup_label,
                                       subgroup_coef,
                                       type,
                                       fit,
                                       trial = NULL) {
  re <- switch(
    type,
    Pooled = "zero",
    Prediction = "new",
    Standardized = "stored",
    Shrinkage = stop(
      "`type = \"Shrinkage\"` is not evaluated at overall ω. Use ",
      "`comparison_shrinkage_rotation_bundle()` / `shrinkage()` (pi_j) or ",
      "`type = \"Standardized\"` (ω + stored RE).",
      call. = FALSE
    ),
    stop("Unknown marginal effect type `", type, "`.", call. = FALSE)
  )
  if (identical(type, "Standardized")) {
    if (is.null(trial) || !nzchar(as.character(trial)[[1L]])) {
      stop("`type = \"Standardized\"` requires a trial id.", call. = FALSE)
    }
  }
  comparison_effect_linpred_draws(
    post = post,
    evaluation_if = omega,
    design = design,
    subgroup_label = subgroup_label,
    subgroup_coef = subgroup_coef,
    re = re,
    fit = fit,
    trial = trial
  )
}

#' Fixed-effect CAMS linpred at a chosen IF (no trial RE)
#'
#' Reference: \eqn{\alpha + \pi^\star(\theta-\gamma)}.
#' Non-reference: \eqn{\alpha + \pi^\star(\theta-\gamma)+\gamma}.
#' Uses the joint overall-block draw of \eqn{(\alpha,\theta)} already stored
#' as `b_Intercept` and `b_ifrac` (`beta = theta - gamma` on that row).
#' Not a convolution of the `alpha` and `theta` mixture marginals, and not
#' `fit$latent$mu_A[[j]]`.
#' @noRd
comparison_marginal_fe_draws <- function(post, evaluation_if, design, subgroup_coef) {
  lp <- post[["b_Intercept"]]
  for (v in design$ifrac_vars_std) {
    col <- paste0("b_", v)
    if (!col %in% names(post)) {
      stop(sprintf("Missing posterior column `%s`.", col), call. = FALSE)
    }
    lp <- lp + evaluation_if * post[[col]]
  }
  if (!is.na(subgroup_coef) && nzchar(subgroup_coef)) {
    col <- paste0("b_", subgroup_coef)
    if (!col %in% names(post)) {
      stop(sprintf("Missing posterior column `%s`.", col), call. = FALSE)
    }
    lp <- lp + post[[col]]
  }
  lp
}

#' Stored trial RE contribution: \eqn{u_j+(x-\mathrm{if})v_j}
#' @noRd
comparison_stored_trial_re <- function(post, fit, trial, evaluation_if, x) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% "trial"
  ifrac2_var <- st$ifrac2_var
  if (is.null(ifrac2_var) || !nzchar(as.character(ifrac2_var)[[1L]])) {
    ifrac2_var <- st$ifrac2_vars_std[[1L]] %||% st$ifrac2_vars[[1L]]
  }
  if (is.null(ifrac2_var) && length(st$ifrac_vars)) {
    ifrac2_var <- paste0("ifrac2_", sub("^ifrac_", "", st$ifrac_vars[[1L]]))
  }
  n <- nrow(post)
  evaluation_if <- as.numeric(evaluation_if)
  if (length(evaluation_if) == 1L) {
    evaluation_if <- rep(evaluation_if, n)
  }
  u_col <- paste0("r_", trial_var, "[", trial, ",Intercept]")
  u <- if (u_col %in% names(post)) {
    as.numeric(post[[u_col]])
  } else {
    rep(0, n)
  }
  v <- rep(0, n)
  if (!is.null(ifrac2_var) && nzchar(as.character(ifrac2_var)[[1L]])) {
    v_col <- paste0("r_", trial_var, "[", trial, ",", ifrac2_var, "]")
    if (v_col %in% names(post)) {
      v <- as.numeric(post[[v_col]])
    }
  }
  u + (x - evaluation_if) * v
}

# Soft alias kept for any external/internal call sites
#' @noRd
comparison_stored_trial_re_at_omega <- function(post, fit, trial, omega, x) {
  comparison_stored_trial_re(
    post = post,
    fit = fit,
    trial = trial,
    evaluation_if = omega,
    x = x
  )
}

#' Random effect for a new trial at a chosen IF (one draw per MCMC iteration)
#'
#' Samples `u ~ N(0, sd_trial__Intercept)` and `v ~ N(0, sd_trial__ifrac2)`
#' independently at each posterior draw. With `x` (0 = reference, 1 = nonref),
#' returns the CAMS design contribution \eqn{u+(x-\mathrm{if})v}. If `x` is
#' `NULL` (overall-effect path), returns \eqn{u+\mathrm{if}\,v} for backward
#' compatibility.
#' @noRd
comparison_new_trial_re_draws <- function(fit, ifrac, post = NULL, x = NULL) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  trial_var <- if (is.null(st$trial_var)) "trial" else st$trial_var
  ifrac2_var <- st$ifrac2_var
  if (is.null(ifrac2_var) || !nzchar(ifrac2_var)) {
    ifrac2_var <- st$ifrac2_vars_std[[1]]
  }
  if (is.null(post)) {
    post <- as.data.frame(fit)
  }
  n <- nrow(post)
  ifrac <- as.numeric(ifrac)
  if (length(ifrac) == 1L) {
    ifrac <- rep(ifrac, n)
  }
  if (interactmeta_has_cond_moments(post)) {
    return(rep(0, n))
  }
  sd_int <- paste0("sd_", trial_var, "__Intercept")
  if (!sd_int %in% names(post)) {
    stop(
      "No posterior column `",
      sd_int,
      "` found for new-trial prediction.",
      call. = FALSE
    )
  }
  comparison_with_seed(comparison_prediction_seed(), {
    re <- stats::rnorm(n, mean = 0, sd = post[[sd_int]])
    if (!is.null(ifrac2_var) && nzchar(ifrac2_var)) {
      sd_slope <- paste0("sd_", trial_var, "__", ifrac2_var)
      if (sd_slope %in% names(post)) {
        v <- stats::rnorm(n, mean = 0, sd = post[[sd_slope]])
        if (is.null(x)) {
          re <- re + ifrac * v
        } else {
          re <- re + (x - ifrac) * v
        }
      }
    }
    re
  })
}

#' Posterior draws of fitted values per subgroup data row
#'
#' Trial-level fitted values are re-centred from each study's information
#' fraction to the precision-weighted overall information fraction from [ifrac()]
#' (same weights as pooled marginal effects).
#' @noRd
comparison_posterior_trial_ifrac2_re <- function(post, trial_var, ifrac2_var, trials) {
  n_draw <- nrow(post)
  mat <- matrix(0, nrow = n_draw, ncol = length(trials))
  colnames(mat) <- trials
  for (j in seq_along(trials)) {
    col <- paste0("r_", trial_var, "[", trials[[j]], ",", ifrac2_var, "]")
    if (col %in% names(post)) {
      mat[, j] <- post[[col]]
    }
  }
  mat
}

#' @noRd
comparison_shrinkage_draws_from_bundle <- function(bundle) {
  # Emit rows in fit-data order: scan trials as stored, ref then nonref only if
  # that matches typical cam() ordering; prefer matching observed row order
  # via trial × subgroup pairs from the fitted data when available.
  pairs <- list()
  for (tr in bundle$trials) {
    pairs[[length(pairs) + 1L]] <- list(trial = tr, subgroup = bundle$reference, mat = "y_A")
    pairs[[length(pairs) + 1L]] <- list(trial = tr, subgroup = bundle$nonref, mat = "y_B")
  }
  rows <- lapply(seq_along(pairs), function(i) {
    p <- pairs[[i]]
    mat <- if (identical(p$mat, "y_A")) bundle$y_A else bundle$y_B
    sd_mat <- if (identical(p$mat, "y_A")) bundle$sd_A else bundle$sd_B
    j <- match(p$trial, bundle$trials)
    df <- data.frame(
      .draw = bundle$post$.draw,
      .chain = bundle$post$.chain,
      .iteration = bundle$post$.iteration,
      row = i,
      trial = p$trial,
      subgroup = p$subgroup,
      type = "Shrinkage",
      effect = mat[, j],
      stringsAsFactors = FALSE
    )
    if (!is.null(sd_mat)) {
      df$.cond_sd <- as.numeric(sd_mat[, j])
    }
    df
  })
  interactmeta_decorate_effect_draws(dplyr::bind_rows(rows), bundle$post)
}

#' Subgroup shrinkage draws via CAMS m/g rotation (no global IF recentering)
#' @noRd
comparison_shrinkage_draws <- function(fit,
                                 contrast_sep = "-",
                                 max_draws = NULL,
                                 trial_var = "trial",
                                 subgroup_var = "subgroup",
                                 ifrac = c("interaction", "overall"),
                                 resolved = NULL) {
  force(trial_var)
  force(subgroup_var)
  force(resolved)
  ifrac <- match.arg(ifrac)
  bundle <- comparison_shrinkage_rotation_bundle(
    fit,
    max_draws = max_draws,
    contrast_sep = contrast_sep,
    ifrac = ifrac
  )
  comparison_shrinkage_draws_from_bundle(bundle)
}


#' @noRd
comparison_effects_summary <- function(draws, group_cols, prob = 0.95,
                                       exponentiate = FALSE) {
  comparison_median_qi_grouped(
    draws = draws,
    value_col = "effect",
    group_cols = group_cols,
    prob = prob,
    exponentiate = exponentiate
  )
}

#' Subgroup-difference (interaction) effects
#'
#' Reports the CAMS treatment-by-subgroup interaction \eqn{\gamma}.
#' For binary CAMS, \eqn{\mu_B(\pi)-\mu_A(\pi)=\gamma} for every
#' information fraction \eqn{\pi}, so this quantity is **not** evaluated
#' at the overall IF. `"Prediction"` adds new-trial interaction
#' heterogeneity. Named `interactions()` to match [subgroups()] and so it
#' does not mask `stats::interaction()`.
#'
#' @param fit A [`interactmeta`] from [cam()].
#' @param contribution,ifrac,prob,max_draws,summary,contrast_sep,exponentiate,shrinkage
#'   As in [subgroups()]. `contribution` / `ifrac` are unused for the
#'   pooled interaction \eqn{\gamma}; they are retained for API compatibility
#'   and are used when `shrinkage = TRUE`.
#' @param type One or both of `"Pooled"` and `"Prediction"` when `shrinkage = FALSE`.
#'   `"Pooled"` is the fixed interaction \eqn{\gamma}. `"Prediction"` is
#'   the interaction in a **new trial**: \eqn{\gamma} plus
#'   trial-level slope heterogeneity `v * (ifrac2_nonref - ifrac2_ref)`
#'   (intercept trial RE cancels within a study).
#'
#' @return A `interactmeta_interaction_effects` object analogous to [subgroups()],
#'   with an `interaction` column instead of `subgroup`.
#'
#' @seealso [subgroups()], [shrinkage()], [combined()],
#'   [forestplot()]
#' @export
interactions <- function(fit,
                                contribution = NULL,
                                ifrac = c("interaction", "overall"),
                                type = c("Pooled", "Prediction"),
                                prob = 0.95,
                                max_draws = NULL,
                                summary = TRUE,
                                contrast_sep = "-",
                                exponentiate = FALSE,
                                shrinkage = FALSE) {
  if (interactmeta_is_non_cams_analysis(fit)) {
    return(interactmeta_extract_interaction(fit, exponentiate = exponentiate))
  }
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cam()`.", call. = FALSE)
  }
  if (isTRUE(shrinkage)) {
    return(shrinkage(
      fit,
      type = "interaction",
      prob = prob,
      max_draws = max_draws,
      summary = summary,
      contrast_sep = contrast_sep,
      ifrac = match.arg(ifrac),
      contribution = contribution,
      exponentiate = exponentiate
    ))
  }
  type <- match.arg(type, c("Pooled", "Prediction"), several.ok = TRUE)
  ifrac <- match.arg(ifrac)

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = contribution,
    ifrac = ifrac,
    max_draws = max_draws,
    contrast_sep = contrast_sep
  )
  int_design <- comparison_interaction_design(fit, contrast_sep = contrast_sep)
  if (!length(int_design$interaction_levels)) {
    stop("No interaction contrasts found for interaction effects.", call. = FALSE)
  }

  rows <- lapply(type, function(ty) {
    lapply(int_design$interaction_levels, function(lbl) {
      coef <- comparison_interaction_coef_term(fit, lbl, int_design, contrast_sep)
      non_ref <- comparison_contrast_non_reference(lbl, int_design$reference, contrast_sep)
      delta_ifrac2 <- comparison_ifrac2_contrast_delta(
        fit,
        non_ref = non_ref,
        ref = int_design$reference,
        contrast_sep = contrast_sep
      )
      eff <- comparison_marginal_interaction_draws(
        post = prep$post,
        subgroup_coef = coef,
        type = ty,
        fit = fit,
        delta_ifrac2 = delta_ifrac2
      )
      cond_sd <- NULL
      if (interactmeta_has_cond_moments(prep$post)) {
        cond_sd <- interactmeta_lincomb_cond_sd(prep$post, 0, 0, 1)
        if (identical(ty, "Prediction")) {
          slope <- grep("^sd_trial__ifrac2_", names(prep$post), value = TRUE)
          tau_g <- if (length(slope)) {
            pmax(as.numeric(prep$post[[slope[[1L]]]]), 0)
          } else {
            rep(0, nrow(prep$post))
          }
          cond_sd <- sqrt(cond_sd^2 + (tau_g * delta_ifrac2)^2)
        }
      }
      interactmeta_decorate_effect_draws(
        data.frame(
          .draw = prep$draw_ids$.draw,
          .chain = prep$draw_ids$.chain,
          .iteration = prep$draw_ids$.iteration,
          interaction = lbl,
          type = ty,
          effect = eff,
          stringsAsFactors = FALSE
        ),
        prep$post,
        cond_sd = cond_sd
      )
    })
  })
  draws <- dplyr::bind_rows(unlist(rows, recursive = FALSE))

  out <- list(
    draws = draws,
    overall_ifrac = prep$overall_ifrac,
    overall_ifrac_source = prep$overall_ifrac_source,
    ifrac = prep$ifrac,
    reference = int_design$reference,
    prob = prob,
    exponentiate = isTRUE(exponentiate),
    shrinkage = FALSE
  )
  if (isTRUE(summary)) {
    delta <- 1
    fin <- interactmeta_summarize_effects(
      draws = draws,
      group_cols = c("interaction", "type"),
      prob = prob,
      exponentiate = exponentiate,
      fit = fit,
      quantity = "interaction",
      omega = interactmeta_metafor_omega(fit, prep$overall_ifrac),
      delta_ifrac2 = delta
    )
    out$draws <- fin$draws
    out$summary <- fin$summary
  } else {
    out$draws <- comparison_transform_effect_draws(draws, exponentiate)
  }
  class(out) <- c("interactmeta_interaction_effects", class(out))
  interactmeta_attach_origin(out, fit, quantity = "interaction")
}

#' @noRd
comparison_interaction_shrinkage_summary_once <- function(summary) {
  if (!nrow(summary) || !"interaction" %in% names(summary)) {
    return(summary)
  }
  summary[!duplicated(summary$interaction), , drop = FALSE]
}

#' @export
print.interactmeta_interaction_effects <- function(x, digits = 2, ...) {
  opts <- comparison_print_digits_prob(x, digits = digits, ...)
  title <- if (isTRUE(x$shrinkage)) {
    "Interaction shrinkage effects (study-specific)"
  } else {
    "Treatment-by-subgroup interaction"
  }
  if (!is.null(x$summary)) {
    summ <- x$summary
    if (isTRUE(x$shrinkage)) {
      summ <- comparison_interaction_shrinkage_summary_once(summ)
    }
    comparison_print_effects_table(
      summ,
      digits = opts$digits,
      prob = opts$prob,
      id_cols = c("trial", "interaction", "type"),
      title = title,
      style = interactmeta_interval_style(x)
    )
  } else {
    cat(title, " (", nrow(x$draws), " rows)\n", sep = "")
  }
  invisible(x)
}

#' @noRd
comparison_interaction_design <- function(fit, contrast_sep = "-") {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st) || is.null(st$ifrac_data)) {
    stop(
      "Fit has no `ifrac_data` in attribute `\"interactmeta\"`; refit with `cam()`.",
      call. = FALSE
    )
  }
  trial_var <- if (is.null(st$trial_var)) "trial" else st$trial_var
  subgroup_var <- if (is.null(st$subgroup_var)) "subgroup" else st$subgroup_var
  ref <- st$reference
  dat <- st$ifrac_data
  sg_chr <- as.character(dat[[subgroup_var]])
  int_dat <- dat[grepl(contrast_sep, sg_chr, fixed = TRUE), , drop = FALSE]
  if (!nrow(int_dat)) {
    sg_dat <- dat[
      is_ifrac_component_row(dat[[subgroup_var]], contrast_sep = contrast_sep),
      ,
      drop = FALSE
    ]
    int_dat <- build_contrast_rows_from_subgroups(
      sg_dat,
      contrast_sep = contrast_sep,
      reference = ref
    )
  }
  if (!nrow(int_dat)) {
    return(list(
      reference = ref,
      interaction_levels = character(0),
      interaction_data = int_dat,
      trial_var = trial_var,
      subgroup_var = subgroup_var
    ))
  }
  int_dat[[subgroup_var]] <- as.character(int_dat[[subgroup_var]])
  list(
    reference = ref,
    interaction_levels = unique(int_dat[[subgroup_var]]),
    interaction_data = int_dat,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
}

#' @noRd
comparison_interaction_coef_term <- function(fit, label, int_design, contrast_sep) {
  non_ref <- comparison_contrast_non_reference(label, int_design$reference, contrast_sep)
  if (is.na(non_ref)) {
    return(NA_character_)
  }
  sg_design <- list(
    reference = int_design$reference,
    subgroup_levels = unique(c(int_design$reference, non_ref))
  )
  sg_coef <- comparison_subgroup_coef_terms(fit, sg_design)
  sg_coef[[non_ref]]
}

#' @noRd
comparison_contrast_non_reference <- function(label, reference, contrast_sep) {
  if (!grepl(contrast_sep, label, fixed = TRUE)) {
    return(NA_character_)
  }
  parts <- strsplit(label, contrast_sep, fixed = TRUE)[[1]]
  if (length(parts) != 2L) {
    return(NA_character_)
  }
  if (!identical(parts[[2]], reference)) {
    stop(
      "Contrast `",
      label,
      "` must place the reference (`",
      reference,
      "`) after `",
      contrast_sep,
      "`.",
      call. = FALSE
    )
  }
  parts[[1]]
}

#' Mean within-trial `ifrac2` contrast (non-reference minus reference)
#' @noRd
comparison_ifrac2_contrast_delta <- function(fit,
                                       non_ref,
                                       ref,
                                       contrast_sep = "-") {
  if (is.na(non_ref) || !nzchar(non_ref)) {
    return(0)
  }
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st) || is.null(st$ifrac_data)) {
    return(1)
  }
  subgroup_var <- if (is.null(st$subgroup_var)) "subgroup" else st$subgroup_var
  trial_var <- if (is.null(st$trial_var)) "trial" else st$trial_var
  dat <- st$ifrac_data
  dat <- dat[
    is_ifrac_component_row(dat[[subgroup_var]], contrast_sep = contrast_sep),
    ,
    drop = FALSE
  ]
  ifrac2_col <- if (!is.null(st$ifrac2_var) && nzchar(st$ifrac2_var)) {
    st$ifrac2_var
  } else {
    paste0("ifrac2_", non_ref)
  }
  if (!(ifrac2_col %in% names(dat))) {
    return(1)
  }
  trials <- unique(dat[[trial_var]])
  deltas <- vapply(
    trials,
    function(tr) {
      d <- dat[dat[[trial_var]] == tr, , drop = FALSE]
      i2 <- as.numeric(d[[ifrac2_col]])
      sg_tr <- as.character(d[[subgroup_var]])
      i2_non <- i2[sg_tr == non_ref]
      i2_ref <- i2[sg_tr == ref]
      if (!length(i2_non) || !length(i2_ref)) {
        return(NA_real_)
      }
      i2_non[[1]] - i2_ref[[1]]
    },
    numeric(1)
  )
  deltas <- deltas[is.finite(deltas)]
  if (!length(deltas)) {
    1
  } else {
    mean(deltas)
  }
}

#' Posterior draws of trial slope RE (ifrac2) for prediction
#'
#' Uses a fixed seed (see `comparison_prediction_seed()`) for reproducibility.
#' @noRd
comparison_new_trial_slope_re_draws <- function(fit, post = NULL) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  trial_var <- if (is.null(st$trial_var)) "trial" else st$trial_var
  ifrac2_var <- st$ifrac2_var
  if (is.null(ifrac2_var) || !nzchar(ifrac2_var)) {
    ifrac2_var <- st$ifrac2_vars_std[[1]]
  }
  if (is.null(post)) {
    post <- as.data.frame(fit)
  }
  n <- nrow(post)
  if (is.null(ifrac2_var) || !nzchar(ifrac2_var)) {
    return(rep(0, n))
  }
  if (interactmeta_has_cond_moments(post)) {
    return(rep(0, n))
  }
  sd_slope <- paste0("sd_", trial_var, "__", ifrac2_var)
  if (!sd_slope %in% names(post)) {
    return(rep(0, n))
  }
  comparison_with_seed(
    comparison_prediction_seed(),
    stats::rnorm(n, mean = 0, sd = post[[sd_slope]])
  )
}

#' Seed used for Monte Carlo new-trial Prediction draws
#' @noRd
comparison_prediction_seed <- function() {
  getOption("interactmeta.prediction_seed", 1L)
}

#' Evaluate an expression under a temporary RNG seed
#' @noRd
comparison_with_seed <- function(seed, expr) {
  seed <- suppressWarnings(as.integer(seed)[[1L]])
  if (is.na(seed)) {
    return(force(expr))
  }
  envir <- globalenv()
  had_seed <- exists(".Random.seed", envir = envir, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = envir, inherits = FALSE)
  } else {
    NULL
  }
  set.seed(seed)
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = envir)
      } else if (exists(".Random.seed", envir = envir, inherits = FALSE)) {
        rm(".Random.seed", envir = envir)
      }
    },
    add = TRUE
  )
  force(expr)
}

#' Posterior draws of a subgroup contrast at overall IF
#'
#' `"Pooled"`: fixed subgroup contrast only. `"Prediction"`: adds slope
#' heterogeneity for a new meta-analysis, `b_contrast + v * delta_ifrac2`.
#' @noRd
comparison_marginal_interaction_draws <- function(post,
                                            subgroup_coef,
                                            type = "Pooled",
                                            fit = NULL,
                                            delta_ifrac2 = 1) {
  if (is.na(subgroup_coef) || !nzchar(subgroup_coef)) {
    lp <- rep(0, nrow(post))
  } else {
    col <- paste0("b_", subgroup_coef)
    if (!col %in% names(post)) {
      stop(sprintf("Missing posterior column `%s`.", col), call. = FALSE)
    }
    lp <- post[[col]]
  }
  if (identical(type, "Prediction")) {
    if (is.null(fit)) {
      stop("`fit` is required for interaction `Prediction` draws.", call. = FALSE)
    }
    lp <- lp + comparison_new_trial_slope_re_draws(fit, post = post) * delta_ifrac2
  }
  lp
}

#' @noRd
interactmeta_interaction_effects_shrinkage <- function(fit,
                                               prob,
                                               max_draws,
                                               summary,
                                               contrast_sep,
                                               exponentiate,
                                               overall_ifrac,
                                               ifrac = c("interaction", "overall"),
                                               bundle = NULL) {
  ifrac <- match.arg(ifrac)
  bundle <- bundle %||% comparison_shrinkage_rotation_bundle(
    fit,
    max_draws = max_draws,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    overall_ifrac = overall_ifrac
  )
  comparison_shrinkage_effect_from_bundle(
    bundle = bundle,
    type = "interaction",
    prob = prob,
    summary = summary,
    exponentiate = exponentiate
  )
}

#' Convert effect summary rows to `yi` / `vi` for forest plots
#'
#' Pooled / Prediction diamonds use the **mean** of posterior draws as `yi`
#' and `sd^2` as `vi` (bayesmeta-compatible). Shrinkage greys use the
#' **median and equal-tail QI** of the reconstructed draws via
#' `interactmeta_effects_merge_yivi()` (`use_posterior_qi = TRUE`).
#' @noRd
interactmeta_effects_summary_to_yivi <- function(summary, type = NULL) {
  s <- summary
  if (!is.null(type) && "type" %in% names(s)) {
    s <- s[s$type == type, , drop = FALSE]
  }
  if (!nrow(s)) {
    return(data.frame(yi = numeric(0), vi = numeric(0)))
  }
  data.frame(
    yi = s$mean,
    vi = pmax(s$sd^2, 1e-12, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Merge effect summaries onto a template data frame by trial and label
#'
#' @param use_posterior_qi If `TRUE` (shrinkage forest rows), set `yi` to the
#'   median of reconstructed draws and attach `ci.lb` / `ci.ub` from the
#'   equal-tail posterior QI. If `FALSE` (default), use mean / `sd^2`.
#' @noRd
interactmeta_effects_merge_yivi <- function(template,
                                    summary,
                                    trial_var,
                                    label_var,
                                    summary_label_var = label_var,
                                    type = NULL,
                                    use_posterior_qi = FALSE) {
  if (!nrow(template)) {
    return(template)
  }
  s <- summary
  if (!is.null(type) && "type" %in% names(s)) {
    s <- s[s$type == type, , drop = FALSE]
  }
  key_t <- paste(template[[trial_var]], template[[label_var]], sep = "\r")
  key_s <- paste(s$trial, s[[summary_label_var]], sep = "\r")
  idx <- match(key_t, key_s)
  if (any(is.na(idx))) {
    stop(
      "Effect summaries do not align with forest template rows.",
      call. = FALSE
    )
  }
  out <- template
  if (isTRUE(use_posterior_qi) &&
      all(c("estimate", "lci", "uci") %in% names(s))) {
    out$yi <- as.numeric(s$estimate[idx])
    out$vi <- pmax(as.numeric(s$sd[idx])^2, 1e-12, na.rm = TRUE)
    out$ci.lb <- as.numeric(s$lci[idx])
    out$ci.ub <- as.numeric(s$uci[idx])
  } else {
    out$yi <- s$mean[idx]
    out$vi <- pmax(s$sd[idx]^2, 1e-12, na.rm = TRUE)
  }
  out
}
