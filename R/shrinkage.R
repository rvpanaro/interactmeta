#' Shrinkage (study-level) effects from a CAMS fit
#'
#' Trial-specific subgroup effects reconstructed at each trial's **own**
#' information fraction \eqn{\pi_j} (the same estimand as the observed black
#' squares on [forestplot()]). Two latent components are shrunk and then
#' inverted:
#' \deqn{\widetilde\mu_{Aj}=\widetilde m_j-\pi_j\widetilde g_j,\qquad
#'   \widetilde\mu_{Bj}=\widetilde m_j+(1-\pi_j)\widetilde g_j.}
#' with
#' \deqn{g_j=\widetilde\mu_{Bj}-\widetilde\mu_{Aj},\qquad
#'   m_j=(1-\pi_j)\widetilde\mu_{Aj}+\pi_j\widetilde\mu_{Bj}.}
#' These identities hold **draw by draw** on the model scale; exponentiation
#' (if requested) happens only afterwards. A grey forest marker therefore need
#' not lie between the quoted estimate and the subgroup pooled diamond.
#'
#' Engine-native trial effects:
#' - **bayesmeta**: DIRECT individual posteriors
#'   (`bmr$shrink.moments` / `bayesmeta$cond.moment(individual=)`), the same
#'   quantities stored in `fit$latent$mu_A` / `$mu_B`.
#' - **metafor**: empirical-Bayes BLUPs of the orthogonal overall and
#'   interaction blocks, then the same reconstruction.
#'
#' This is **not** the fixed-effect linear predictor at \eqn{\pi_j}.
#' [fitted()] returns that FE mean (`fit$mu`); [shrinkage()] is the
#' trial-level reconstruction the forest greys use. Overall / Prediction use
#' overall IF \eqn{\omega}; see [subgroups()] and
#' [standardized_trial_effects()].
#'
#' @param fit A [`interactmeta`] from [cam()].
#' @param type Which estimand (same labels as [forestplot()]):
#'   - `"subgroups"` (default; `"subgroup"` accepted): \eqn{\mu_{A,j}, \mu_{B,j}}
#'   - `"interaction"`: \eqn{g_j}
#'   - `"combined"`: \eqn{m_j}
#' @param format Output representation:
#'   - `"list"` (default): S3 object with long `$draws`, optional `$summary`,
#'     and `$components` (observed / shrunk matrices)
#'   - `"summary"`: summary data frame only
#' @param prob Coverage probability for intervals when a summary is built.
#' @param max_draws Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports; omit it.
#' @param summary If `TRUE` (default), attach / build a summary when
#'   `format = "list"` or `"summary"`.
#' @param contrast_sep Passed to contrast resolution (default `"-"`).
#' @param ifrac How overall IF \eqn{\omega} is pooled via [ifrac()] when
#'   `contribution` is `NULL` (metadata / annotations only for shrinkage;
#'   draws themselves use \eqn{\pi_j}): `"interaction"` (default) or `"overall"`.
#' @param contribution Optional combined-population information fraction
#'   (same as [subgroups()]).
#' @param exponentiate If `TRUE`, transform posterior draws with `exp()` after
#'   the linear-scale reconstruction (for reporting only).
#' @param component Soft-deprecated alias for `type`.
#'
#' @return Depends on `format`. For `format = "list"`, see `$components` and
#'   [shrinkage_components()].
#'
#' @seealso [shrinkage_components()], [standardized_trial_effects()],
#'   [subgroups()], [interactions()], [combined()], [forestplot()]
#' @export
shrinkage <- function(fit,
                      type = c("subgroups", "interaction", "combined"),
                      format = c("list", "summary"),
                      prob = 0.95,
                      max_draws = NULL,
                      summary = TRUE,
                      contrast_sep = "-",
                      ifrac = c("interaction", "overall"),
                      contribution = NULL,
                      exponentiate = FALSE,
                      component = NULL) {
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cam()`.", call. = FALSE)
  }

  type <- comparison_shrinkage_resolve_type(type = type, component = component)
  format <- match.arg(format)
  ifrac <- match.arg(ifrac)
  need_summary <- isTRUE(summary) || identical(format, "summary")
  obj <- interactmeta_shrinkage_component(
    type = type,
    fit = fit,
    prob = prob,
    max_draws = max_draws,
    summary = need_summary,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    exponentiate = exponentiate,
    overall_ifrac = contribution
  )
  comparison_format_effect_object(
    obj = obj,
    format = format,
    type = type,
    fit = fit
  )
}

#' Audit components for CAMS shrinkage (observed, IF, m/g, subgroup draws)
#'
#' Returns the quantities used to build [shrinkage()] so rotation identities
#' can be checked draw by draw on the model scale. Draws are evaluated at each
#' trial's own \eqn{\pi_j} (same as [shrinkage()] and forest greys), not at
#' overall \eqn{\omega}.
#'
#' @param fit A [`interactmeta`], or a shrinkage return value that already
#'   contains `$components`.
#' @param ... Passed through when `fit` is an `interactmeta`: `max_draws`,
#'   `ifrac`, `contrast_sep`, `contribution`.
#' @return A list with:
#'   - `pi`: named trial information fractions
#'   - `observed`: named vectors `y_A`, `y_B`, `m`, `g` and sampling variances
#'   - `draws`: draws × trial matrices `y_A`, `y_B`, `m`, `g` at \eqn{\pi_j}
#'   - `reference`, `nonref`, `trials`, `interaction_label`
#' @keywords internal
#' @export
shrinkage_components <- function(fit, ..., contribution = NULL) {
  if (is.list(fit) && !is.null(fit$components)) {
    return(fit$components)
  }
  check_comparison_fit(fit)
  dots <- list(...)
  if (!is.null(contribution) && is.null(dots$overall_ifrac)) {
    dots$overall_ifrac <- contribution
  }
  # Always the study-specific (π_j) estimand for shrinkage_components()
  dots$evaluation <- NULL
  args <- c(list(fit = fit, evaluation = "trial"), dots)
  comparison_shrinkage_components_list(do.call(comparison_shrinkage_rotation_bundle, args))
}

#' Trial-specific effects standardized to the overall IF
#'
#' Same stored random effects as [shrinkage()], but evaluated at the overall
#' information fraction \eqn{\omega} rather than each trial's \eqn{\pi_j}:
#' \deqn{\mu_{ij}=\alpha+u_j+\beta\omega+(x_i-\omega)v_j+\gamma\,1(i=B).}
#' This is **trial-specific standardization to the target population**, not
#' shrinkage. Prefer [shrinkage()] for forest greys and study-specific
#' posterior means.
#'
#' @param fit A [`interactmeta`] from [cam()].
#' @param type Which estimand: `"subgroups"`, `"interaction"`, or `"combined"`.
#' @param format Output representation (same choices as [shrinkage()]).
#' @param prob Coverage probability for intervals when a summary is built.
#' @param max_draws Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports; omit it.
#' @param summary If `TRUE` (default), attach / build a summary when
#'   `format = "list"` or `"summary"`.
#' @param contrast_sep Passed to contrast resolution (default `"-"`).
#' @param ifrac How overall IF \eqn{\omega} is pooled via [ifrac()] when
#'   `contribution` is `NULL`: `"interaction"` (default) or `"overall"`.
#' @param contribution Optional combined-population information fraction
#'   (same as [subgroups()]).
#' @param exponentiate If `TRUE`, transform posterior draws with `exp()`.
#' @param component Soft-deprecated alias for `type`.
#' @return Same structure as [shrinkage()] for the chosen `type` / `format`.
#' @seealso [shrinkage()], [subgroups()]
#' @export
standardized_trial_effects <- function(fit,
                                       type = c("subgroups", "interaction", "combined"),
                                       format = c("list", "summary"),
                                       prob = 0.95,
                                       max_draws = NULL,
                                       summary = TRUE,
                                       contrast_sep = "-",
                                       ifrac = c("interaction", "overall"),
                                       contribution = NULL,
                                       exponentiate = FALSE,
                                       component = NULL) {
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cam()`.", call. = FALSE)
  }

  type <- comparison_shrinkage_resolve_type(type = type, component = component)
  format <- match.arg(format)
  ifrac <- match.arg(ifrac)
  need_summary <- isTRUE(summary) || identical(format, "summary")
  bundle <- comparison_shrinkage_rotation_bundle(
    fit,
    max_draws = max_draws,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    overall_ifrac = contribution,
    evaluation = "overall"
  )
  obj <- comparison_shrinkage_effect_from_bundle(
    bundle = bundle,
    type = type,
    prob = prob,
    summary = need_summary,
    exponentiate = exponentiate
  )
  obj$standardized <- TRUE
  comparison_format_effect_object(
    obj = obj,
    format = format,
    type = type,
    fit = fit
  )
}

#' Resolve shrinkage `type`, accepting forestplot labels and legacy `component`
#' @noRd
comparison_shrinkage_resolve_type <- function(type = c("subgroups", "interaction", "combined"),
                                        component = NULL) {
  if (!is.null(component)) {
    warning(
      "`component` is soft-deprecated; use `type` with the same labels as ",
      "`forestplot()` / `traceplot()` (`\"subgroups\"`, `\"interaction\"`, ",
      "or `\"combined\"`).",
      call. = FALSE
    )
    component <- as.character(component)[[1L]]
    type <- switch(
      component,
      subgroup = "subgroups",
      subgroups = "subgroups",
      interaction = "interaction",
      combined = "combined",
      all = stop(
        "`component = \"all\"` is no longer supported; call `shrinkage()` ",
        "once per `type` (`\"subgroups\"`, `\"interaction\"`, or ",
        "`\"combined\"`).",
        call. = FALSE
      ),
      stop(
        "`component` must be one of \"subgroup\", \"subgroups\", ",
        "\"interaction\", or \"combined\".",
        call. = FALSE
      )
    )
    return(comparison_resolve_effect_type(type))
  }
  type <- as.character(type)[[1L]]
  if (identical(type, "all")) {
    stop(
      "`type = \"all\"` is not supported; use `\"subgroups\"`, ",
      "`\"interaction\"`, or `\"combined\"` (same as `forestplot()`).",
      call. = FALSE
    )
  }
  type <- match.arg(type, c("subgroups", "interaction", "combined"))
  comparison_resolve_effect_type(type)
}

#' Diagnostic table for study-specific CAMS shrinkage
#'
#' Per-trial table of observed vs shrunk quantities, both at trial \eqn{\pi_j}.
#' Draw-level identity checked:
#' \deqn{\widetilde y_{A,j}=\widetilde m_j-\pi_j\widetilde g_j.}
#'
#' @param fit A [`interactmeta`].
#' @param ... Passed to [shrinkage_components()] (`max_draws`, `ifrac`, …).
#' @return A data frame with one row per trial.
#' @seealso [shrinkage()], [shrinkage_components()]
#' @export
shrinkage_diagnostics <- function(fit, ...) {
  cmp <- shrinkage_components(fit, ...)
  ref <- cmp$reference
  nonref <- cmp$nonref
  om_draws <- as.numeric(cmp$omega)
  if (!length(om_draws)) {
    stop("Shrinkage components are missing overall IF `omega`.", call. = FALSE)
  }
  if (length(om_draws) == 1L) {
    om_draws <- rep(om_draws, nrow(cmp$draws$m))
  }
  rows <- lapply(cmp$trials, function(tr) {
    pi_j <- unname(cmp$pi[[tr]])
    m_obs <- unname(cmp$observed$m[[tr]])
    g_obs <- unname(cmp$observed$g[[tr]])
    yA_obs <- unname(cmp$observed$y_A[[tr]])
    yB_obs <- unname(cmp$observed$y_B[[tr]])
    w <- cmp$weight
    sd_list <- cmp$sd
    m_s <- comparison_shrinkage_node_median(
      cmp$draws$m[, tr], w, if (!is.null(sd_list$m)) sd_list$m[, tr] else NULL
    )
    g_s <- comparison_shrinkage_node_median(
      cmp$draws$g[, tr], w, if (!is.null(sd_list$g)) sd_list$g[, tr] else NULL
    )
    yA_s <- comparison_shrinkage_node_median(
      cmp$draws$y_A[, tr], w, if (!is.null(sd_list$y_A)) sd_list$y_A[, tr] else NULL
    )
    yB_s <- comparison_shrinkage_node_median(
      cmp$draws$y_B[, tr], w, if (!is.null(sd_list$y_B)) sd_list$y_B[, tr] else NULL
    )
    id_err <- max(abs(
      cmp$draws$y_A[, tr] - (cmp$draws$m[, tr] - pi_j * cmp$draws$g[, tr])
    ))
    data.frame(
      trial = tr,
      reference = ref,
      nonref = nonref,
      g_orientation = paste(nonref, "-", ref),
      pi = pi_j,
      pi_nonref = pi_j,
      pi_reference = 1 - pi_j,
      omega = stats::median(om_draws),
      observed_m = m_obs,
      shrunk_m = m_s,
      observed_g = g_obs,
      shrunk_g = g_s,
      observed_reference = yA_obs,
      reconstructed_shrunk_reference = yA_s,
      observed_nonref = yB_obs,
      reconstructed_shrunk_nonref = yB_s,
      reference_identity_max_abs = id_err,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Mixture / weighted median of shrinkage nodes for one trial
#' @noRd
comparison_shrinkage_node_median <- function(x, weight = NULL, cond_sd = NULL) {
  interactmeta_mixture_moments(
    mean = x,
    weight = weight,
    cond_sd = cond_sd
  )$median
}

#' @noRd
interactmeta_shrinkage_component <- function(type, fit, ...) {
  out <- switch(
    type,
    subgroups = interactmeta_subgroup_effects_shrinkage(fit, ...),
    interaction = interactmeta_interaction_effects_shrinkage(fit, ...),
    combined = interactmeta_combined_effects_shrinkage(fit, ...),
    stop("Unknown shrinkage type `", type, "`.", call. = FALSE)
  )
  interactmeta_attach_origin(out, fit, quantity = paste0("shrinkage_", type))
}
