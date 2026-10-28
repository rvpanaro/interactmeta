# ---- marginal_effects.R ----
#' Apply `fit$reporting` defaults for marginal extractors
#' @noRd
interactmeta_finalize_marginal_args <- function(fit,
                                                prob = NULL,
                                                ifrac = NULL,
                                                contrast_sep = NULL,
                                                exponentiate = NULL) {
  ifrac_resolved <- if (is.null(ifrac)) {
    NULL
  } else {
    match.arg(ifrac, c("interaction", "overall"))
  }
  opts <- interactmeta_resolve_reporting_opts(
    fit,
    prob = prob,
    ifrac = ifrac_resolved,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate
  )
  list(
    prob = opts$prob,
    ifrac = match.arg(opts$ifrac, c("interaction", "overall")),
    contrast_sep = opts$contrast_sep,
    exponentiate = opts$exponentiate
  )
}

interactmeta_resolve_marginal_at <- function(at = NULL, contribution) {
  if (!missing(contribution)) {
    if (!is.null(at)) {
      stop("Supply only one of `at` or `contribution`.", call. = FALSE)
    }
    if (!is.null(contribution)) {
      warning(
        "`contribution` is deprecated; use `at` instead.",
        call. = FALSE
      )
    }
    return(contribution)
  }
  at
}

#' Modifier-level treatment effects
#'
#' Evaluates CAMS **modifier-specific** treatment effects at the pooled
#' (precision-weighted) information fraction from [ifrac()]. Vignettes often
#' use published **subgroup** tables (IV/NIV, age strata, …); the same
#' extractor applies whenever each trial contributes one estimate per
#' **modifier level** (long format). The modifier column defaults to
#' `"level"`; pass **`modifier =`** on [cams()] to use another column name.
#' Display headings capitalize the first letter of that column name.
#'
#' This is **not** a marginal treatment effect in the survey-sampling sense:
#' it is the modifier-level linear predictor at the combined-population
#' contribution structure \eqn{\widetilde{\pi}}.
#'
#' For each posterior draw, the linear predictor is evaluated at that draw's
#' combined-population information fraction with the appropriate subgroup
#' fixed effect:
#' \deqn{\mu_A(\pi^\star)=\alpha+\pi^\star(\theta-\gamma),}
#' \deqn{\mu_B(\pi^\star)=\alpha+\pi^\star(\theta-\gamma)+\gamma.}
#' On the orthogonal CAMS factorization, \eqn{\theta} is the overall-block IF slope and
#' stored `b_ifrac` is \eqn{\beta=\theta-\gamma}. This is **not** the
#' trial-specific convolution on `fit$latent$mu_A[[j]]` /
#' `fit$latent$mu_B[[j]]` (\eqn{\mu_{Aj}=m_j-\pi_j g_j} at that trial's
#' \eqn{\pi_j}). A scalar `at` is a plug-in \eqn{\pi^\star};
#' `NULL` mixes the [ifrac()] posterior.
#'
#' \eqn{\alpha} and \eqn{\theta} share the overall-block support. Population
#' \eqn{M^\star=\alpha+\pi^\star\theta} is formed from that joint, not by
#' convolving the two mixture marginals.
#'
#' @param fit A [`interactmeta`] from [cams()].
#' @param at Combined-population **information fraction** \eqn{\pi} at which to
#'   evaluate subgroup effects. If `NULL` (default), uses the pooled reporting
#'   IF from [ifrac()] for `fit`. A numeric vector is aligned to posterior draws;
#'   a scalar (e.g. `at = 0.5`) is a plug-in at every draw.
#' @param ifrac When `at` is `NULL`, passed to [ifrac()]: `"interaction"`
#'   (default; interaction-weighted pooling of trial-level \eqn{\pi_j}) or
#'   `"overall"` for overall-block precision weights.
#' @param type One or both of `"Overall"` and `"Prediction"`. `"Overall"` uses fixed
#'   effects only; `"Prediction"` adds a new-trial random effect drawn from the
#'   posterior trial SDs at the combined-population IF (effect in a new **trial**).
#' @param prob Coverage probability for `lci` / `uci` when `summary = TRUE`.
#' @param max_support Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports (weighted quadrature nodes); omit it.
#' @param summary If `TRUE` (default), include a summary data frame.
#' @param contrast_sep Rows whose `subgroup` label contains this separator are
#'   excluded when resolving subgroup levels (default `"-"`).
#' @param exponentiate If `TRUE`, transform posterior draws with `exp()` before
#'   summarizing (e.g. log rate ratios to ratios).
#' @param shrinkage If `TRUE`, delegate to [shrinkage()] with
#'   `type = "modifiers"` (trial-specific reconstruction at each trial's
#'   \eqn{\pi_j}). If `FALSE` (default), return only Overall / Prediction at the
#'   combined-population information fraction.
#' @param contribution Deprecated; use `\code{at}`.
#'
#' @return A `interactmeta_subgroup_effects` object (list with `draws`, `summary`,
#'   `overall_ifrac`, `overall_ifrac_source` (`"fit"` or `"user"`), `reference`,
#'   `prob`, `exponentiate`). Draws are indexed by `subgroup` and `type`.
#'   `summary` includes `estimate` (median), `mean`, `sd`, `lci`, and `uci`.
#'   For Bayesian CAMS, `lci`/`uci` follow `fit$interval.type` (default
#'   `central`) on the DIRECT-node mixture.
#'
#' @seealso [combined()], [interactions()], [shrinkage()], [ifrac()],
#'   [forestplot.interactmeta()]
#' @name modifiers
#' @aliases subgroups
#' @export
modifiers <- function(fit,
                             at = NULL,
                             ifrac = NULL,
                             type = c("Overall", "Prediction"),
                             prob = NULL,
                             max_support = NULL,
                             summary = TRUE,
                             contrast_sep = NULL,
                             exponentiate = NULL,
                             shrinkage = FALSE,
                             contribution) {
  at <- interactmeta_resolve_marginal_at(at = at, contribution = contribution)
  marg <- interactmeta_finalize_marginal_args(
    fit,
    prob = prob,
    ifrac = ifrac,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate
  )
  prob <- marg$prob
  ifrac <- marg$ifrac
  contrast_sep <- marg$contrast_sep
  exponentiate <- marg$exponentiate
  st_mod <- attr(fit, "interactmeta", exact = TRUE)
  modifier_name <- st_mod$modifier %||% fit$modifier %||% st_mod$subgroup_var %||% interactmeta_level_column()
  if (interactmeta_is_non_cams_analysis(fit)) {
    interactmeta_validate_construction_marginal_args(
      fit,
      "modifiers()",
      type = interactmeta_construction_marginal_type(type),
      shrinkage = shrinkage,
      at = at,
      ifrac = ifrac
    )
    out <- interactmeta_extract_subgroups(fit, exponentiate = exponentiate)
    if (is.list(out)) {
      out$modifier <- modifier_name
    }
    return(out)
  }
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cams()`.", call. = FALSE)
  }
  if (isTRUE(shrinkage)) {
    out <- shrinkage(
      fit,
      type = "modifiers",
      prob = prob,
      max_support = max_support,
      summary = summary,
      contrast_sep = contrast_sep,
      ifrac = ifrac,
      at = at,
      exponentiate = exponentiate
    )
    if (is.list(out)) {
      out$modifier <- modifier_name
    }
    return(out)
  }
  type <- match.arg(type, c("Overall", "Prediction"), several.ok = TRUE)

  cache_key <- interactmeta_effect_cache_key(
    "subgroups",
    at = at,
    ifrac = ifrac,
    type = type,
    prob = prob,
    max_support = max_support,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate,
    summary = summary
  )
  cached <- interactmeta_effect_cache_get(fit, cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = at,
    ifrac = ifrac,
    max_support = max_support,
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
          level = lbl,
          type = ty,
          effect = eff,
          stringsAsFactors = FALSE
        ),
        prep$post,
        cond_sd = cond_sd
      )
    })
  })
  support <- dplyr::bind_rows(unlist(rows, recursive = FALSE))

  out <- list(
    support = support,
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
      draws = support,
      group_cols = c("level", "type"),
      prob = prob,
      exponentiate = exponentiate,
      fit = fit,
      quantity = "subgroups",
      omega = interactmeta_metafor_omega(fit, prep$overall_ifrac)
    )
    out$support <- fin$support
    out$summary <- fin$summary
  } else {
    out$support <- comparison_transform_effect_draws(support, exponentiate)
  }
  class(out) <- c("interactmeta_subgroup_effects", class(out))
  out <- interactmeta_attach_origin(out, fit, quantity = "subgroups")
  interactmeta_effect_cache_set(fit, cache_key, out)
  out
}

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
#' @param fit A [`interactmeta`] from [cams()].
#' @param at Combined-population information fraction \eqn{\pi}; see [modifiers()].
#' @param ifrac When `at` is `NULL`, passed to [ifrac()]: `"interaction"`
#'   (default) or `"overall"` for overall-block pooling of trial
#'   information fractions.
#' @param type One or both of `"Overall"` and `"Prediction"` (new-trial random
#'   effect; see [modifiers()]).
#' @param prob,max_support,summary,contrast_sep,exponentiate,shrinkage,contribution
#'   As in [modifiers()].
#'
#' @return An `interactmeta_combined_effect` object (list with `draws`,
#'   `summary`, `overall_ifrac`, `overall_ifrac_source`, `ifrac`, `reference`,
#'   `prob`). Draws are indexed by `type`. `summary` includes `estimate`
#'   (median), `mean`, `sd`, `lci`, and `uci`.
#'
#' @seealso [modifiers()] for level-specific \eqn{\mu_A(\pi)} and \eqn{\mu_B(\pi)} at the
#'   same reporting IF; [interactions()], [ifrac()]
#' @export
combined <- function(fit,
                     at = NULL,
                     ifrac = NULL,
                     type = c("Overall", "Prediction"),
                     prob = NULL,
                     max_support = NULL,
                     summary = TRUE,
                     contrast_sep = NULL,
                     exponentiate = NULL,
                     shrinkage = FALSE,
                     contribution) {
  at <- interactmeta_resolve_marginal_at(at = at, contribution = contribution)
  marg <- interactmeta_finalize_marginal_args(
    fit,
    prob = prob,
    ifrac = ifrac,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate
  )
  prob <- marg$prob
  ifrac <- marg$ifrac
  contrast_sep <- marg$contrast_sep
  exponentiate <- marg$exponentiate
  interactmeta_abort_non_cams(fit, "combined()")
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cams()`.", call. = FALSE)
  }
  if (isTRUE(shrinkage)) {
    return(shrinkage(
      fit,
      type = "combined",
      prob = prob,
      max_support = max_support,
      summary = summary,
      contrast_sep = contrast_sep,
      ifrac = ifrac,
      at = at,
      exponentiate = exponentiate
    ))
  }
  type <- match.arg(type, c("Overall", "Prediction"), several.ok = TRUE)

  cache_key <- interactmeta_effect_cache_key(
    "combined",
    at = at,
    ifrac = ifrac,
    type = type,
    prob = prob,
    max_support = max_support,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate,
    summary = summary
  )
  cached <- interactmeta_effect_cache_get(fit, cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  if (length(type) == 1L && identical(type[[1L]], "Overall")) {
    native <- interactmeta_try_native_combined(
      fit,
      at = at,
      ifrac = ifrac,
      prob = prob,
      contrast_sep = contrast_sep
    )
    if (!is.null(native)) {
      interactmeta_effect_cache_set(fit, cache_key, native)
      return(native)
    }
  }

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = at,
    ifrac = ifrac,
    max_support = max_support,
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
  support <- dplyr::bind_rows(rows)

  out <- list(
    support = support,
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
      draws = support,
      group_cols = "type",
      prob = prob,
      exponentiate = exponentiate,
      fit = fit,
      quantity = "combined",
      omega = interactmeta_metafor_omega(fit, prep$overall_ifrac)
    )
    out$support <- fin$support
    out$summary <- comparison_tag_combined_population_summary(fin$summary)
  } else {
    out$support <- comparison_transform_effect_draws(support, exponentiate)
  }
  class(out) <- c("interactmeta_combined_effect", class(out))
  out <- interactmeta_attach_origin(out, fit, quantity = "combined")
  interactmeta_effect_cache_set(fit, cache_key, out)
  out
}

#' @noRd
interactmeta_subgroup_effects_shrinkage <- function(fit,
                                            prob,
                                            max_support,
                                            summary,
                                            contrast_sep,
                                            exponentiate,
                                            overall_ifrac,
                                            ifrac = c("interaction", "overall"),
                                            bundle = NULL) {
  ifrac <- match.arg(ifrac)
  bundle <- bundle %||% comparison_shrinkage_rotation_bundle(
    fit,
    max_support = max_support,
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
interactmeta_combined_effects_shrinkage <- function(fit,
                                           prob,
                                           max_support,
                                           summary,
                                           contrast_sep,
                                           exponentiate,
                                           overall_ifrac,
                                           ifrac = c("interaction", "overall"),
                                           bundle = NULL) {
  ifrac <- match.arg(ifrac)
  bundle <- bundle %||% comparison_shrinkage_rotation_bundle(
    fit,
    max_support = max_support,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    overall_ifrac = overall_ifrac
  )
  comparison_shrinkage_effect_from_bundle(
    bundle = bundle,
    type = "combined",
    prob = prob,
    summary = summary,
    exponentiate = exponentiate
  )
}

#' Build a typed shrinkage effect object from one rotation bundle
#' @noRd
comparison_shrinkage_effect_from_bundle <- function(bundle,
                                              type = c("modifiers", "subgroups", "interaction", "combined"),
                                              prob = 0.95,
                                              summary = TRUE,
                                              exponentiate = FALSE) {
  type <- comparison_resolve_effect_type(match.arg(type))
  if (identical(type, "subgroups")) {
    draws <- comparison_shrinkage_draws_from_bundle(bundle)
    fin <- comparison_shrinkage_finalize_draws(
      draws = draws,
      group_cols = c("trial", "level", "type"),
      prob = prob,
      summary = summary,
      exponentiate = exponentiate
    )
    out <- list(
      support = fin$support,
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
      support = fin$support,
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
    support = fin$support,
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
    bundle, fin$summary, type = "combined", prob = prob, exponentiate = exponentiate
  )
  class(out) <- c("interactmeta_combined_effect", class(out))
  out
}

#' Public shrinkage summaries use shortest DIRECT-mixture intervals.
#' Native `$latent` convolve / shrink intervals remain available on the fit.
#' @noRd
interactmeta_bayesmeta_maybe_fill_shrinkage <- function(bundle,
                                                       summary,
                                                       type,
                                                       prob,
                                                       exponentiate) {
  fit <- bundle$fit
  if (is.null(summary) || !nrow(summary) ||
      !interactmeta_is_bayesian_fit(fit) ||
      !identical(bundle$evaluation %||% "trial", "trial")) {
    return(summary)
  }
  fit <- interactmeta_ensure_bayesmeta_latent(fit, full_trials = TRUE)
  bundle$fit <- fit
  if (is.null(fit$latent)) {
    return(summary)
  }
  interval.type <- interactmeta_match_interval_type(fit$interval.type)
  for (i in seq_len(nrow(summary))) {
    tr <- as.character(summary$trial[[i]])
    posterior <- if (identical(type, "subgroups")) {
      lvl_col <- interactmeta_modifier_summary_col(summary)
      subgroup <- as.character(summary[[lvl_col]][[i]])
      if (identical(subgroup, bundle$reference)) {
        fit$latent$mu_A[[tr]]
      } else if (identical(subgroup, bundle$nonref)) {
        fit$latent$mu_B[[tr]]
      } else {
        NULL
      }
    } else if (identical(type, "interaction")) {
      fit$latent$g[[tr]]
    } else if (identical(type, "combined")) {
      fit$latent$m[[tr]]
    } else {
      NULL
    }
    if (is.null(posterior)) {
      next
    }
    sm <- interactmeta_posterior_1d_summary(
      posterior,
      prob = prob,
      interval.type = interval.type
    )
    if (isTRUE(exponentiate)) {
      summary$estimate[[i]] <- exp(sm$median)
      summary$lci[[i]] <- exp(sm$lower)
      summary$uci[[i]] <- exp(sm$upper)
      if (!is.null(posterior$gm)) {
        moments <- interactmeta_gm_lognormal_moments(posterior$gm)
        summary$mean[[i]] <- moments$mean
        summary$sd[[i]] <- moments$sd
      }
    } else {
      summary$estimate[[i]] <- sm$median
      summary$mean[[i]] <- sm$mean
      summary$sd[[i]] <- sm$sd
      summary$lci[[i]] <- sm$lower
      summary$uci[[i]] <- sm$upper
    }
  }
  attr(summary, "interval.type") <- interval.type
  summary
}

#' Summarize shrinkage draws; exponentiate only after reconstruction
#'
#' Mixture / weighted nodes (bayesmeta DIRECT, metafor BLUP SEs) are
#' summarized on the model scale and then optionally transformed.
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
    support = comparison_transform_effect_draws(draws, exponentiate),
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
                                         trial_var = "study",
                                         subgroup_var = "level",
                                         pi = NULL) {
  ref <- as.character(reference)[[1L]]
  trials <- unique(as.character(sg_draws$trial))
  lvl_col <- interactmeta_modifier_summary_col(sg_draws)
  subgroups <- unique(as.character(sg_draws[[lvl_col]]))
  nonrefs <- setdiff(subgroups, ref)
  if (length(nonrefs) != 1L) {
    stop(
      "`comparison_overall_shrinkage_draws()` supports exactly one non-reference subgroup.",
      call. = FALSE
    )
  }
  ref_template <- sg_draws[
    sg_draws$trial == trials[[1L]] & sg_draws[[lvl_col]] == ref,
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
      sg_draws$trial == tr & sg_draws[[lvl_col]] == ref,
      ,
      drop = FALSE
    ]
    non <- sg_draws[
      sg_draws$trial == tr & sg_draws[[lvl_col]] == nonrefs[[1L]],
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

#' Row label for population [combined()] summaries (forest `All` block)
#' @noRd
comparison_combined_population_label <- function() {
  "All"
}

#' @noRd
comparison_tag_combined_population_summary <- function(summary, shrinkage = FALSE) {
  if (is.null(summary) || !nrow(summary) || isTRUE(shrinkage)) {
    return(summary)
  }
  lab <- comparison_combined_population_label()
  if (!"combined" %in% names(summary)) {
    summary$combined <- lab
  } else {
    empty <- is.na(summary$combined) | !nzchar(as.character(summary$combined))
    summary$combined[empty] <- lab
  }
  summary
}

#' @noRd
comparison_combined_print_id_cols <- function(x) {
  if (isTRUE(x$shrinkage)) {
    c("trial", "type")
  } else {
    "combined"
  }
}

#' @export
print.interactmeta_subgroup_effects <- function(x, digits = 2, ...) {
  opts <- comparison_print_digits_prob(x, digits = digits, ...)
  heading <- interactmeta_modifier_heading(x$modifier %||% interactmeta_level_column())
  title <- if (isTRUE(x$shrinkage)) {
    paste0(heading, " shrinkage effects (trial-specific)")
  } else {
    interactmeta_modifier_specific_effects_title(x$modifier %||% interactmeta_level_column())
  }
  if (!is.null(x$summary)) {
    lvl_col <- interactmeta_modifier_summary_col(x$summary)
    id_cols <- if (isTRUE(x$shrinkage)) {
      c("trial", lvl_col, "type")
    } else {
      c(lvl_col, "type")
    }
    comparison_print_effects_table(
      x$summary,
      digits = opts$digits,
      prob = opts$prob,
      id_cols = id_cols,
      title = NULL,
      style = interactmeta_interval_style(x)
    )
  } else {
    cat(title, " (", nrow(x$support), " rows)\n", sep = "")
  }
  invisible(x)
}

#' @export
print.interactmeta_combined_effect <- function(x, digits = 2, ...) {
  opts <- comparison_print_digits_prob(x, digits = digits, ...)
  title <- if (isTRUE(x$shrinkage)) {
    "Combined-population shrinkage effects (trial-specific)"
  } else {
    "Combined-population treatment effect"
  }
  if (!is.null(x$summary)) {
    comparison_print_effects_table(
      x$summary,
      digits = opts$digits,
      prob = opts$prob,
      id_cols = comparison_combined_print_id_cols(x),
      title = NULL,
      style = interactmeta_interval_style(x)
    )
  } else {
    cat(title, " (", nrow(x$support), " rows)\n", sep = "")
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
                                        max_support = NULL,
                                        contrast_sep = "-") {
  resolved <- comparison_resolve_marginal_omega(
    fit,
    overall_ifrac = overall_ifrac,
    ifrac = ifrac,
    max_support = max_support
  )
  design <- comparison_marginal_design(fit, contrast_sep = contrast_sep)
  list(
    omega = resolved$omega,
    post = resolved$post,
    draw_ids = comparison_support_row_ids(resolved$post),
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
                                        max_support = NULL) {
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
      "`at` must be `NULL` or a finite numeric vector (or scalar).",
      call. = FALSE
    )
  }

  post_full <- comparison_support_df(fit)
  n_post <- nrow(post_full)
  w_post <- if (".weight" %in% names(post_full)) post_full$.weight else NULL

  if (length(omega) == 1L) {
    idx <- comparison_marginal_draw_index(
      n_post,
      max_support = max_support,
      weight = w_post
    )
    post <- post_full[idx, , drop = FALSE]
    omega <- rep(omega, nrow(post))
  } else if (length(omega) == n_post) {
    idx <- comparison_marginal_draw_index(
      n_post,
      max_support = max_support,
      weight = w_post
    )
    omega <- omega[idx]
    post <- post_full[idx, , drop = FALSE]
  } else {
    idx <- comparison_marginal_draw_index(
      length(omega),
      max_support = max_support,
      weight = NULL
    )
    if (length(idx) > n_post) {
      stop(
        "`at` has more values than posterior draws in `fit`.",
        call. = FALSE
      )
    }
    if (!is.null(w_post) && length(idx) < n_post &&
          interactmeta_has_direct_weights(w_post, n_post)) {
      stop(
        "`max_support` cannot thin a weighted DIRECT support. ",
        "Bayesian CAMS summaries use the full quadrature grid; omit `max_support`.",
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
      "Fit has no `ifrac_data` in attribute `\"interactmeta\"`; refit with `cams()`.",
      call. = FALSE
    )
  }
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
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
comparison_support_df <- function(fit) {
  post <- comparison_tidy_support(fit)
  raw <- as.data.frame(fit$support)
  keep <- grep("^\\.(weight|cond_)", names(raw), value = TRUE)
  for (nm in keep) {
    if (!nm %in% names(post) && nrow(raw) == nrow(post)) {
      post[[nm]] <- raw[[nm]]
    }
  }
  post
}

#' @noRd
comparison_support_row_ids <- function(post) {
  data.frame(
    .draw = post$.draw,
    .chain = post$.chain,
    .iteration = post$.iteration,
    stringsAsFactors = FALSE
  )
}

#' Index rows of a draws table, optionally capped by `max_support`
#'
#' Weighted DIRECT supports (Bayesian CAMS) must not be thinned uniformly:
#' that over-represents low-weight quadrature nodes. Omit `max_support` for
#' those fits.
#'
#' @param n Number of rows.
#' @param max_support Optional maximum number of rows. Ignored when `NULL` or
#'   at least `n`.
#' @param weight Optional DIRECT node weights aligned with the `n` rows.
#' @noRd
comparison_marginal_draw_index <- function(n, max_support = NULL, weight = NULL) {
  n <- as.integer(n)[[1L]]
  if (!is.finite(n) || n < 0L) {
    stop("`n` must be a non-negative integer.", call. = FALSE)
  }
  if (is.null(max_support) || !is.finite(max_support) || max_support >= n) {
    return(seq_len(n))
  }
  if (interactmeta_has_direct_weights(weight, n)) {
    stop(
      "`max_support` cannot thin a weighted DIRECT support. ",
      "Bayesian CAMS summaries use the full quadrature grid; omit `max_support`.",
      call. = FALSE
    )
  }
  set.seed(1L)
  sample.int(n, as.integer(max_support))
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

#' Combined-block IF slope \eqn{\theta = \beta_{\mathrm{pkg}} + \gamma}
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
                                               type = "Overall") {
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
#' @param re `"zero"` (Overall / Standardized) or `"new"` (Prediction).
#' @noRd
comparison_effect_linpred_draws <- function(post,
                                      evaluation_if,
                                      design,
                                      subgroup_label,
                                      subgroup_coef,
                                      re = c("zero", "new"),
                                      fit = NULL) {
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
  lp
}

#' Posterior draws of the marginalized effect (linear predictor at overall IF)
#'
#' Thin wrapper around `comparison_effect_linpred_draws()` for Overall / Prediction /
#' Standardized. Study-specific shrinkage uses
#' `comparison_shrinkage_rotation_bundle()` with `evaluation = "trial"` instead.
#' @noRd
comparison_marginal_effect_draws <- function(post,
                                       omega,
                                       design,
                                       subgroup_label,
                                       subgroup_coef,
                                       type,
                                       fit) {
  re <- switch(
    type,
    Overall = "zero",
    Prediction = "new",
    Standardized = "zero",
    Shrinkage = stop(
      "`type = \"Shrinkage\"` is not evaluated at overall \u03c9. Use ",
      "`comparison_shrinkage_rotation_bundle()` / `shrinkage()` (pi_j) or ",
      "`type = \"Standardized\"`.",
      call. = FALSE
    ),
    stop("Unknown marginal effect type `", type, "`.", call. = FALSE)
  )
  comparison_effect_linpred_draws(
    post = post,
    evaluation_if = omega,
    design = design,
    subgroup_label = subgroup_label,
    subgroup_coef = subgroup_coef,
    re = re,
    fit = fit
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

#' Random effect for a new trial at a chosen IF
#'
#' Samples `u ~ N(0, sd_trial__Intercept)` and `v ~ N(0, sd_trial__ifrac2)`
#' independently at each posterior draw. With `x` (0 = reference, 1 = nonref),
#' returns the CAMS design contribution \eqn{u+(x-\mathrm{if})v}. If `x` is
#' `NULL` (overall-effect path), returns \eqn{u+\mathrm{if}\,v} for backward
#' compatibility.
#' @noRd
comparison_new_trial_re_draws <- function(fit, ifrac, post = NULL, x = NULL) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% interactmeta_study_column()
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
  sd_int <- paste0("sd_", interactmeta_sd_re_group(), "__Intercept")
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
      sd_slope <- paste0("sd_", interactmeta_sd_re_group(), "__", ifrac2_var)
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

#' @noRd
comparison_shrinkage_draws_from_bundle <- function(bundle) {
  # Emit rows in fit-data order: scan trials as stored, ref then nonref only if
  # that matches typical cams() ordering; prefer matching observed row order
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
      level = p$subgroup,
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


#' @noRd
comparison_effects_summary <- function(draws, group_cols, prob = 0.95,
                                       exponentiate = FALSE,
                                       interval.type = NULL) {
  comparison_median_qi_grouped(
    draws = draws,
    value_col = "effect",
    group_cols = group_cols,
    prob = prob,
    exponentiate = exponentiate,
    interval.type = interval.type
  )
}

#' Subgroup-difference (interaction) effects
#'
#' Reports the CAMS treatment-by-subgroup interaction \eqn{\gamma}.
#' For binary CAMS, \eqn{\mu_B(\pi)-\mu_A(\pi)=\gamma} for every
#' information fraction \eqn{\pi}, so this quantity is **not** evaluated
#' at the overall IF. `"Prediction"` adds new-trial interaction
#' heterogeneity. Named `interactions()` to match [modifiers()] and so it
#' does not mask `stats::interaction()`.
#'
#' @param fit A [`interactmeta`] from [cams()].
#' @param at,ifrac,prob,max_support,summary,contrast_sep,exponentiate,shrinkage,
#'   contribution As in [modifiers()]. `at` / `ifrac` are unused for the pooled
#'   interaction \eqn{\gamma}; they are retained for API compatibility and are
#'   used when `shrinkage = TRUE`.
#' @param type One or both of `"Overall"` and `"Prediction"` when `shrinkage = FALSE`.
#'   For CAMS, `"Overall"` is the fixed interaction \eqn{\gamma}; `"Prediction"`
#'   adds new-trial random effects. For [swada()], `"Prediction"` is the
#'   SWADA predictive interval (common study weights on \eqn{g_j} plus
#'   interaction \eqn{\hat\tau_g^2}); DA / AD accept `"Overall"` only.
#'
#' @return A `interactmeta_interaction_effects` object analogous to [modifiers()],
#'   with an `interaction` column instead of `subgroup`.
#'
#' @seealso [modifiers()], [shrinkage()], [combined()],
#'   [forestplot()]
#' @export
interactions <- function(fit,
                                at = NULL,
                                ifrac = NULL,
                                type = c("Overall", "Prediction"),
                                prob = NULL,
                                max_support = NULL,
                                summary = TRUE,
                                contrast_sep = NULL,
                                exponentiate = NULL,
                                shrinkage = FALSE,
                                contribution) {
  at <- interactmeta_resolve_marginal_at(at = at, contribution = contribution)
  marg <- interactmeta_finalize_marginal_args(
    fit,
    prob = prob,
    ifrac = ifrac,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate
  )
  prob <- marg$prob
  ifrac <- marg$ifrac
  contrast_sep <- marg$contrast_sep
  exponentiate <- marg$exponentiate
  if (interactmeta_is_non_cams_analysis(fit)) {
    analysis <- interactmeta_fit_analysis(fit)
    marg_type <- interactmeta_construction_marginal_type(type, analysis = analysis)
    interactmeta_validate_construction_marginal_args(
      fit,
      "interactions()",
      type = marg_type,
      shrinkage = shrinkage,
      at = at,
      ifrac = ifrac
    )
    return(interactmeta_extract_interaction(
      fit,
      type = marg_type,
      prob = prob,
      exponentiate = exponentiate
    ))
  }
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cams()`.", call. = FALSE)
  }
  if (isTRUE(shrinkage)) {
    return(shrinkage(
      fit,
      type = "interaction",
      prob = prob,
      max_support = max_support,
      summary = summary,
      contrast_sep = contrast_sep,
      ifrac = ifrac,
      at = at,
      exponentiate = exponentiate
    ))
  }
  type <- match.arg(type, c("Overall", "Prediction"), several.ok = TRUE)

  cache_key <- interactmeta_effect_cache_key(
    "interactions",
    at = at,
    ifrac = ifrac,
    type = type,
    prob = prob,
    max_support = max_support,
    contrast_sep = contrast_sep,
    exponentiate = exponentiate,
    summary = summary
  )
  cached <- interactmeta_effect_cache_get(fit, cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  if (length(type) == 1L && identical(type[[1L]], "Overall")) {
    native <- interactmeta_try_native_interactions(
      fit,
      at = at,
      ifrac = ifrac,
      prob = prob,
      contrast_sep = contrast_sep
    )
    if (!is.null(native)) {
      interactmeta_effect_cache_set(fit, cache_key, native)
      return(native)
    }
  }

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = at,
    ifrac = ifrac,
    max_support = max_support,
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
  support <- dplyr::bind_rows(unlist(rows, recursive = FALSE))

  out <- list(
    support = support,
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
      draws = support,
      group_cols = c("interaction", "type"),
      prob = prob,
      exponentiate = exponentiate,
      fit = fit,
      quantity = "interaction",
      omega = interactmeta_metafor_omega(fit, prep$overall_ifrac),
      delta_ifrac2 = delta
    )
    out$support <- fin$support
    out$summary <- fin$summary
  } else {
    out$support <- comparison_transform_effect_draws(support, exponentiate)
  }
  class(out) <- c("interactmeta_interaction_effects", class(out))
  out <- interactmeta_attach_origin(out, fit, quantity = "interaction")
  interactmeta_effect_cache_set(fit, cache_key, out)
  out
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
    "Interaction shrinkage effects (trial-specific)"
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
      title = NULL,
      style = interactmeta_interval_style(x)
    )
  } else {
    cat(title, " (", nrow(x$support), " rows)\n", sep = "")
  }
  invisible(x)
}

#' @noRd
comparison_interaction_design <- function(fit, contrast_sep = "-") {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st) || is.null(st$ifrac_data)) {
    stop(
      "Fit has no `ifrac_data` in attribute `\"interactmeta\"`; refit with `cams()`.",
      call. = FALSE
    )
  }
  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
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
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
  trial_var <- st$trial_var %||% interactmeta_study_column()
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
  trial_var <- st$trial_var %||% interactmeta_study_column()
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
  sd_slope <- paste0("sd_", interactmeta_sd_re_group(), "__", ifrac2_var)
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
#' `"Overall"`: fixed subgroup contrast only. `"Prediction"`: adds slope
#' heterogeneity for a new meta-analysis, `b_contrast + v * delta_ifrac2`.
#' @noRd
comparison_marginal_interaction_draws <- function(post,
                                            subgroup_coef,
                                            type = "Overall",
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
                                               max_support,
                                               summary,
                                               contrast_sep,
                                               exponentiate,
                                               overall_ifrac,
                                               ifrac = c("interaction", "overall"),
                                               bundle = NULL) {
  ifrac <- match.arg(ifrac)
  bundle <- bundle %||% comparison_shrinkage_rotation_bundle(
    fit,
    max_support = max_support,
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

#' Merge effect summaries onto a template data frame by trial and label
#'
#' @param use_posterior_qi If `TRUE` (shrinkage forest rows), set `yi` to the
#'   median of reconstructed draws and attach `ci.lb` / `ci.ub` from the
#'   shortest posterior interval. If `FALSE` (default), use mean / `sd^2`.
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

# ---- observed_effects.R ----
#' Observed trial-level treatment-effect estimates and variances
#'
#' Stable extractor for the observed estimates used in CAMS post-fit analysis.
#' Prefer this over reading undocumented fields on `object$data`.
#'
#' @param object A [`interactmeta`] from [cams()].
#' @param type Estimand:
#'   - `"modifiers"` (default; `"subgroup"` / deprecated `"subgroups"` accepted): one entry per
#'     trial × modifier-level row
#'   - `"combined"`: trial-specific combined-population
#'     effect \eqn{m_j} (inverse-variance combination within each trial)
#'   - `"interaction"`: trial-specific treatment-by-subgroup interaction
#'     \eqn{g_j} (non-reference minus reference within each trial)
#' @param contrast_sep Separator for interaction labels (default `"-"`).
#'
#' @return A list with:
#'   - `estimate`: named numeric vector
#'   - `variance`: named numeric vector (sampling variances)
#'   - `type`: resolved type
#'   - `labels`: character vector of names (same order)
#'
#' For `type = "modifiers"`, names are `"<trial> / <subgroup>"`. For
#' `"combined"` and `"interaction"`, names are trial IDs.
#'
#' @seealso [shrinkage()], [fitted.interactmeta()]
#' @examples
#' \dontrun{
#' observed <- observed_effects(fit, type = "combined")
#' observed$estimate["RESPIRE1"]
#' observed$variance["RESPIRE1"]
#' }
#' @export
observed_effects <- function(object,
                             type = c("modifiers", "combined", "interaction", "subgroups"),
                             contrast_sep = "-") {
  check_comparison_fit(object)
  type_chr <- as.character(type)[[1L]]
  if (identical(type_chr, "subgroups")) {
    .Deprecated(msg = 'type = "subgroups" is deprecated; use type = "modifiers".')
  }
  type <- comparison_resolve_effect_type(match.arg(type))
  dat <- comparison_predict_data(object, newdata = NULL)
  st <- attr(object, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
  ref <- as.character(st$reference)[[1L]]

  if (identical(type, "subgroups")) {
    labs <- paste(
      as.character(dat[[trial_var]]),
      as.character(dat[[subgroup_var]]),
      sep = " / "
    )
    if (anyDuplicated(labs)) {
      labs <- make.unique(labs, sep = " / ")
    }
    est <- stats::setNames(as.numeric(dat$yi), labs)
    var <- stats::setNames(pmax(as.numeric(dat$vi), 0), labs)
    return(list(
      estimate = est,
      variance = var,
      type = type,
      labels = labs
    ))
  }

  trials <- unique(as.character(dat[[trial_var]]))
  comparison_assert_unique_names(trials, what = "trial")

  nonrefs <- setdiff(unique(as.character(dat[[subgroup_var]])), ref)
  blocks <- comparison_trial_effect_blocks(
    data = dat,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = ref,
    nonref = if (length(nonrefs)) nonrefs[[1L]] else NULL
  )

  if (identical(type, "combined")) {
    est <- stats::setNames(blocks$overall$yi, blocks$overall$trial)
    var <- stats::setNames(blocks$overall$vi, blocks$overall$trial)
    return(list(
      estimate = est[trials],
      variance = var[trials],
      type = type,
      labels = trials
    ))
  }

  # interaction
  est <- stats::setNames(blocks$interaction$yi, blocks$interaction$trial)
  var <- stats::setNames(blocks$interaction$vi, blocks$interaction$trial)
  list(
    estimate = est[trials],
    variance = var[trials],
    type = type,
    labels = trials
  )
}

#' Resolve the public estimand `type`
#' @noRd
comparison_resolve_effect_type <- function(type) {
  type <- as.character(type)[[1L]]
  if (identical(type, "subgroup")) {
    type <- "subgroups"
  }
  if (identical(type, "modifiers")) {
    type <- "subgroups"
  }
  match.arg(type, c("subgroups", "combined", "interaction"))
}

#' @noRd
comparison_predict_data <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    dat <- object$data
  } else {
    dat <- newdata
  }
  if (is.null(dat) || !nrow(dat)) {
    stop("No data available for prediction.", call. = FALSE)
  }
  need <- c("yi", "vi")
  miss <- setdiff(need, names(dat))
  if (length(miss)) {
    stop(
      "Prediction data must include columns: ",
      paste(need, collapse = ", "),
      ". Missing: ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  st <- attr(object, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
  if (!trial_var %in% names(dat) || !subgroup_var %in% names(dat)) {
    stop(
      "Prediction data must include `",
      trial_var,
      "` and `",
      subgroup_var,
      "`.",
      call. = FALSE
    )
  }
  dat
}

#' @noRd
comparison_assert_unique_names <- function(x, what = "name") {
  x <- as.character(x)
  if (!length(x)) {
    stop("No ", what, " names found.", call. = FALSE)
  }
  if (anyNA(x) || any(!nzchar(x))) {
    stop(what, " names must be non-missing and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop(what, " names must be unique.", call. = FALSE)
  }
  invisible(x)
}
