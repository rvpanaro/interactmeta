# ---- shrinkage.R ----
#' Shrinkage (trial-level) effects from a CAMS fit
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
#' @param fit A [`interactmeta`] from [cams()].
#' @param type Which estimand (same labels as [forestplot()]):
#'   - `"modifiers"` (default; `"subgroup"` / deprecated `"subgroups"` accepted): \eqn{\mu_{A,j}, \mu_{B,j}}
#'   - `"interaction"`: \eqn{g_j}
#'   - `"combined"`: \eqn{m_j}
#' @param format Output representation:
#'   - `"list"` (default): S3 object with long `$support`, optional `$summary`,
#'     and `$components` (observed / shrunk matrices)
#'   - `"summary"`: summary data frame only
#' @param prob Coverage probability for intervals when a summary is built.
#' @param max_support Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports; omit it.
#' @param summary If `TRUE` (default), attach / build a summary when
#'   `format = "list"` or `"summary"`.
#' @param contrast_sep Passed to contrast resolution (default `"-"`).
#' @param ifrac How overall IF \eqn{\omega} is pooled via [ifrac()] when
#'   `at` is `NULL` (metadata / annotations only for shrinkage; draws themselves
#'   use \eqn{\pi_j}): `"interaction"` (default) or `"overall"`.
#' @param at Optional reporting information fraction (same as [modifiers()]).
#' @param contribution Deprecated; use `\code{at}`.
#' @param exponentiate If `TRUE`, transform posterior draws with `exp()` after
#'   the linear-scale reconstruction (for reporting only).
#' @param component Soft-deprecated alias for `type`.
#'
#' @return Depends on `format`. For `format = "list"`, see `$components` and
#'   [shrinkage_components()].
#'
#' @seealso [shrinkage_components()], [standardized_trial_effects()],
#'   [modifiers()], [interactions()], [combined()], [forestplot()]
#' @export
shrinkage <- function(fit,
                      type = c("modifiers", "interaction", "combined", "subgroups"),
                      format = c("list", "summary"),
                      prob = 0.95,
                      max_support = NULL,
                      summary = TRUE,
                      contrast_sep = "-",
                      ifrac = c("interaction", "overall"),
                      at = NULL,
                      exponentiate = FALSE,
                      component = NULL,
                      contribution) {
  at <- interactmeta_resolve_marginal_at(at = at, contribution = contribution)
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cams()`.", call. = FALSE)
  }

  type <- comparison_shrinkage_resolve_type(type = type, component = component)
  format <- match.arg(format)
  ifrac <- match.arg(ifrac)
  need_summary <- isTRUE(summary) || identical(format, "summary")
  obj <- interactmeta_shrinkage_component(
    type = type,
    fit = fit,
    prob = prob,
    max_support = max_support,
    summary = need_summary,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    exponentiate = exponentiate,
    overall_ifrac = at
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
#' @param ... Passed through when `fit` is an `interactmeta`: `max_support`,
#'   `ifrac`, `contrast_sep`, `at`.
#' @return A list with:
#'   - `pi`: named trial information fractions
#'   - `observed`: named vectors `y_A`, `y_B`, `m`, `g` and sampling variances
#'   - `support`: support × trial matrices `y_A`, `y_B`, `m`, `g` at \eqn{\pi_j}
#'   - `reference`, `nonref`, `trials`, `interaction_label`
#' @keywords internal
#' @export
shrinkage_components <- function(fit, ..., at = NULL, contribution) {
  at <- interactmeta_resolve_marginal_at(at = at, contribution = contribution)
  if (is.list(fit) && !is.null(fit$components)) {
    return(fit$components)
  }
  check_comparison_fit(fit)
  dots <- list(...)
  if (!is.null(at) && is.null(dots$overall_ifrac)) {
    dots$overall_ifrac <- at
  }
  # Always the study-specific (π_j) estimand for shrinkage_components()
  dots$evaluation <- NULL
  args <- c(list(fit = fit, evaluation = "trial"), dots)
  comparison_shrinkage_components_list(do.call(comparison_shrinkage_rotation_bundle, args))
}

#' Trial-specific effects standardized to the overall IF
#'
#' Evaluates the CAMS fixed-effect reconstruction at the overall information
#' fraction \eqn{\omega}, rather than using each trial's \eqn{\pi_j} and
#' engine-native trial latent effect:
#' \deqn{\mu_i=\alpha+\beta\omega+\gamma\,1(i=B).}
#' This is standardization to the target population, not shrinkage. Prefer
#' [shrinkage()] for forest greys and trial-specific posterior means.
#'
#' @param fit A [`interactmeta`] from [cams()].
#' @param type Which estimand: `"modifiers"`, `"interaction"`, or `"combined"`.
#' @param format Output representation (same choices as [shrinkage()]).
#' @param prob Coverage probability for intervals when a summary is built.
#' @param max_support Optional cap on unweighted posterior rows. Rejected for
#'   Bayesian CAMS DIRECT supports; omit it.
#' @param summary If `TRUE` (default), attach / build a summary when
#'   `format = "list"` or `"summary"`.
#' @param contrast_sep Passed to contrast resolution (default `"-"`).
#' @param ifrac How overall IF \eqn{\omega} is pooled via [ifrac()] when
#'   `at` is `NULL`: `"interaction"` (default) or `"overall"`.
#' @param at Optional reporting information fraction (same as [modifiers()]).
#' @param contribution Deprecated; use `\code{at}`.
#' @param exponentiate If `TRUE`, transform posterior draws with `exp()`.
#' @param component Soft-deprecated alias for `type`.
#' @return Same structure as [shrinkage()] for the chosen `type` / `format`.
#' @seealso [shrinkage()], [modifiers()]
#' @export
standardized_trial_effects <- function(fit,
                                       type = c("modifiers", "interaction", "combined", "subgroups"),
                                       format = c("list", "summary"),
                                       prob = 0.95,
                                       max_support = NULL,
                                       summary = TRUE,
                                       contrast_sep = "-",
                                       ifrac = c("interaction", "overall"),
                                       at = NULL,
                                       exponentiate = FALSE,
                                       component = NULL,
                                       contribution) {
  at <- interactmeta_resolve_marginal_at(at = at, contribution = contribution)
  check_comparison_fit(fit)
  if (!inherits(fit, "interactmeta")) {
    stop("`fit` must be an interactmeta object from `cams()`.", call. = FALSE)
  }

  type <- comparison_shrinkage_resolve_type(type = type, component = component)
  format <- match.arg(format)
  ifrac <- match.arg(ifrac)
  need_summary <- isTRUE(summary) || identical(format, "summary")
  bundle <- comparison_shrinkage_rotation_bundle(
    fit,
    max_support = max_support,
    contrast_sep = contrast_sep,
    ifrac = ifrac,
    overall_ifrac = at,
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
comparison_shrinkage_resolve_type <- function(type = c("modifiers", "interaction", "combined", "subgroups"),
                                        component = NULL) {
  if (!is.null(component)) {
    warning(
      "`component` is soft-deprecated; use `type` with the same labels as ",
      "`forestplot()` / `traceplot()` (`\"modifiers\"`, `\"interaction\"`, ",
      "or `\"combined\"`).",
      call. = FALSE
    )
    component <- as.character(component)[[1L]]
    type <- switch(
      component,
      subgroup = "subgroups",
      subgroups = "subgroups",
      modifiers = "subgroups",
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
      "`type = \"all\"` is not supported; use `\"modifiers\"`, ",
      "`\"interaction\"`, or `\"combined\"`.",
      call. = FALSE
    )
  }
  if (identical(type, "subgroups")) {
    .Deprecated(msg = 'type = "subgroups" is deprecated; use type = "modifiers".')
    type <- "subgroups"
  }
  type <- match.arg(type, c("modifiers", "subgroups", "interaction", "combined"))
  comparison_resolve_effect_type(type)
}

#' Diagnostic table for trial-specific CAMS shrinkage
#'
#' Per-trial table of observed vs shrunk quantities, both at trial \eqn{\pi_j}.
#' Draw-level identity checked:
#' \deqn{\widetilde y_{A,j}=\widetilde m_j-\pi_j\widetilde g_j.}
#'
#' @param fit A [`interactmeta`].
#' @param ... Passed to [shrinkage_components()] (`max_support`, `ifrac`, …).
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
    om_draws <- rep(om_draws, nrow(cmp$support$m))
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
      cmp$support$m[, tr], w, if (!is.null(sd_list$m)) sd_list$m[, tr] else NULL
    )
    g_s <- comparison_shrinkage_node_median(
      cmp$support$g[, tr], w, if (!is.null(sd_list$g)) sd_list$g[, tr] else NULL
    )
    yA_s <- comparison_shrinkage_node_median(
      cmp$support$y_A[, tr], w, if (!is.null(sd_list$y_A)) sd_list$y_A[, tr] else NULL
    )
    yB_s <- comparison_shrinkage_node_median(
      cmp$support$y_B[, tr], w, if (!is.null(sd_list$y_B)) sd_list$y_B[, tr] else NULL
    )
    id_err <- max(abs(
      cmp$support$y_A[, tr] - (cmp$support$m[, tr] - pi_j * cmp$support$g[, tr])
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

# ---- shrinkage_format.R ----
#' Format an interactmeta effect object according to `format`
#' @noRd
comparison_format_effect_object <- function(obj,
                                      format = c("list", "summary"),
                                      type,
                                      fit = NULL) {
  format <- match.arg(format)
  if (identical(format, "list")) {
    return(obj)
  }
  if (identical(format, "summary")) {
    if (is.null(obj$summary)) {
      stop(
        "No summary available; call with `summary = TRUE` or use another format.",
        call. = FALSE
      )
    }
    return(tibble::as_tibble(obj$summary))
  }
  stop("Unknown format.", call. = FALSE)
}

# ---- shrinkage_rotation.R ----
#' CAMS study-level effects via the canonical linpred
#'
#' Default **shrinkage** (`evaluation = "trial"`) evaluates at each trial's own
#' information fraction \eqn{\pi_j} with engine-native trial latent effects
#' (same estimand as the observed black squares):
#' \deqn{\widetilde\mu_{ij}
#'   =\alpha+u_j+\beta\pi_j+(x_i-\pi_j)v_j+\gamma\,1(i=B).}
#' Optional `evaluation = "overall"` evaluates the fixed-effect reconstruction
#' at pooled IF \eqn{\omega} (standardization — not shrinkage).
#' Orthogonal summaries use the same evaluation IF.
#'
#' @param evaluation `"trial"` (default, shrinkage at \eqn{\pi_j}) or
#'   `"overall"` (standardized to \eqn{\omega}).
#' @noRd
comparison_shrinkage_rotation_bundle <- function(fit,
                                           max_support = NULL,
                                           contrast_sep = "-",
                                           ifrac = c("interaction", "overall"),
                                           overall_ifrac = NULL,
                                           evaluation = c("trial", "overall")) {
  ifrac <- match.arg(ifrac)
  evaluation <- match.arg(evaluation)
  check_comparison_fit(fit)
  if (identical(evaluation, "trial") && interactmeta_is_bayesian_fit(fit)) {
    fit <- interactmeta_ensure_bayesmeta_latent(fit, full_trials = TRUE)
  }
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st) || is.null(st$ifrac_data)) {
    stop(
      "Fit has no `ifrac_data` in attribute `\"interactmeta\"`; refit with `cams()`.",
      call. = FALSE
    )
  }
  trial_var <- st$trial_var %||% interactmeta_study_column()
  subgroup_var <- st$subgroup_var %||% interactmeta_level_column()
  ref <- interactmeta_reference(fit)

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = overall_ifrac,
    ifrac = ifrac,
    max_support = max_support,
    contrast_sep = contrast_sep
  )
  post <- prep$post
  omega <- prep$omega
  n_draw <- nrow(post)
  design <- prep$design

  data_eval <- fit$data
  if (is.null(data_eval) || !nrow(data_eval)) {
    stop("Fit has no model data; refit with `cams()`.", call. = FALSE)
  }
  ifrac_vars <- comparison_ifrac_coef_names(fit)
  if (!length(ifrac_vars)) {
    stop("No information-fraction columns found on the fit.", call. = FALSE)
  }
  if (length(ifrac_vars) > 1L) {
    stop(
      "CAMS shrinkage currently supports a single non-reference ",
      "subgroup (one IF column). Found: ",
      paste(ifrac_vars, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  ifrac_var <- ifrac_vars[[1L]]
  if (!ifrac_var %in% names(data_eval)) {
    stop("Missing IF column `", ifrac_var, "` in fit data.", call. = FALSE)
  }

  sg_levels <- as.character(design$subgroup_levels)
  nonrefs <- setdiff(sg_levels, ref)
  if (length(nonrefs) != 1L) {
    stop(
      "CAMS shrinkage requires exactly one non-reference subgroup; found: ",
      paste(nonrefs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  nonref <- nonrefs[[1L]]

  trials <- unique(as.character(data_eval[[trial_var]]))
  comparison_assert_unique_names(trials, what = "trial")

  pi <- comparison_trial_information_fractions(
    data = data_eval,
    trial_var = trial_var,
    ifrac_var = ifrac_var,
    trials = trials
  )

  # Observed (quoted) at each trial's own π_j
  y_A_obs <- stats::setNames(rep(NA_real_, length(trials)), trials)
  y_B_obs <- y_A_obs
  v_A <- y_A_obs
  v_B <- y_A_obs
  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    ii <- which(as.character(data_eval[[trial_var]]) == tr)
    sg <- as.character(data_eval[[subgroup_var]][ii])
    i_ref <- ii[sg == ref]
    i_non <- ii[sg == nonref]
    if (length(i_ref) != 1L || length(i_non) != 1L) {
      stop(
        "Trial `",
        tr,
        "` must have exactly one reference and one non-reference row.",
        call. = FALSE
      )
    }
    y_A_obs[[tr]] <- as.numeric(data_eval$yi[[i_ref]])
    y_B_obs[[tr]] <- as.numeric(data_eval$yi[[i_non]])
    v_A[[tr]] <- pmax(as.numeric(data_eval$vi[[i_ref]]), 0)
    v_B[[tr]] <- pmax(as.numeric(data_eval$vi[[i_non]]), 0)
  }

  ifrac2_var <- paste0("ifrac2_", sub("^ifrac_", "", ifrac_var))
  if (!ifrac2_var %in% names(data_eval)) {
    ifrac2_var <- paste0("ifrac2_", nonref)
  }
  if (!ifrac2_var %in% names(data_eval)) {
    stop("Missing IF design column `", ifrac2_var, "` in fit data.", call. = FALSE)
  }
  pi_A <- stats::setNames(rep(NA_real_, length(trials)), trials)
  pi_B <- pi_A
  i2_A <- pi_A
  i2_B <- pi_A
  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    ii <- which(as.character(data_eval[[trial_var]]) == tr)
    sg <- as.character(data_eval[[subgroup_var]][ii])
    i_ref <- ii[sg == ref][[1L]]
    i_non <- ii[sg == nonref][[1L]]
    pi_A[[tr]] <- as.numeric(data_eval[[ifrac_var]][[i_ref]])
    pi_B[[tr]] <- as.numeric(data_eval[[ifrac_var]][[i_non]])
    i2_A[[tr]] <- as.numeric(data_eval[[ifrac2_var]][[i_ref]])
    i2_B[[tr]] <- as.numeric(data_eval[[ifrac2_var]][[i_non]])
  }
  comparison_assert_pi_is_nonref(
    pi = pi,
    v_A = v_A,
    v_B = v_B,
    pi_A = pi_A,
    pi_B = pi_B,
    i2_A = i2_A,
    i2_B = i2_B,
    reference = ref,
    nonref = nonref,
    ifrac_var = ifrac_var
  )

  m_obs <- (1 - pi) * y_A_obs + pi * y_B_obs
  g_obs <- y_B_obs - y_A_obs
  v_m <- (1 - pi)^2 * v_A + pi^2 * v_B
  v_g <- v_A + v_B

  native <- NULL
  if (identical(evaluation, "trial")) {
    native <- comparison_shrinkage_native_trial_latents(
      fit = fit,
      post = post,
      trials = trials,
      pi = pi,
      reference = ref,
      nonref = nonref
    )
  }

  sd_A <- sd_B <- sd_m <- sd_g <- NULL
  eval_if_by_trial <- vector("list", length(trials))
  names(eval_if_by_trial) <- trials
  if (!is.null(native)) {
    y_A <- native$y_A
    y_B <- native$y_B
    m <- native$m
    g <- native$g
    sd_A <- native$sd_A
    sd_B <- native$sd_B
    sd_m <- native$sd_m
    sd_g <- native$sd_g
    for (tr in trials) {
      eval_if_by_trial[[tr]] <- unname(pi[[tr]])
    }
  } else {
    # Canonical linear predictor for standardization at ω. Native
    # bayesmeta / metafor trial latents are handled above for shrinkage.
    y_A <- comparison_named_trial_matrix(n_draw, trials)
    y_B <- comparison_named_trial_matrix(n_draw, trials)
    for (j in seq_along(trials)) {
      tr <- trials[[j]]
      evaluation_if <- if (identical(evaluation, "trial")) {
        unname(pi[[tr]])
      } else {
        omega
      }
      eval_if_by_trial[[tr]] <- evaluation_if
      y_A[, j] <- comparison_effect_linpred_draws(
        post = post,
        evaluation_if = evaluation_if,
        design = design,
        subgroup_label = ref,
        subgroup_coef = prep$sg_coef[[ref]],
        re = "zero",
        fit = fit
      )
      y_B[, j] <- comparison_effect_linpred_draws(
        post = post,
        evaluation_if = evaluation_if,
        design = design,
        subgroup_label = nonref,
        subgroup_coef = prep$sg_coef[[nonref]],
        re = "zero",
        fit = fit
      )
    }
    m <- comparison_named_trial_matrix(n_draw, trials)
    g <- comparison_named_trial_matrix(n_draw, trials)
    for (j in seq_along(trials)) {
      tr <- trials[[j]]
      eif <- eval_if_by_trial[[tr]]
      if (length(eif) == 1L) {
        eif <- rep(eif, n_draw)
      }
      m[, j] <- (1 - eif) * y_A[, j] + eif * y_B[, j]
      g[, j] <- y_B[, j] - y_A[, j]
    }
  }

  list(
    post = post,
    reference = ref,
    nonref = nonref,
    trials = trials,
    pi = pi,
    omega = omega,
    evaluation = evaluation,
    evaluation_if = eval_if_by_trial,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep,
    interaction_label = paste(nonref, ref, sep = contrast_sep),
    y_A = y_A,
    y_B = y_B,
    m = m,
    g = g,
    sd_A = sd_A,
    sd_B = sd_B,
    sd_m = sd_m,
    sd_g = sd_g,
    y_A_obs = y_A_obs,
    y_B_obs = y_B_obs,
    m_obs = m_obs,
    g_obs = g_obs,
    v_A = v_A,
    v_B = v_B,
    v_m = v_m,
    v_g = v_g,
    overall_ifrac = prep$overall_ifrac,
    overall_ifrac_source = prep$overall_ifrac_source,
    ifrac = prep$ifrac,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    fit = fit
  )
}

#' Named draws × trial matrix
#' @noRd
comparison_named_trial_matrix <- function(n_draw, trials, value = NA_real_) {
  out <- matrix(value, nrow = n_draw, ncol = length(trials))
  colnames(out) <- trials
  out
}

#' Engine-native trial-level CAMS latents at each study's own π_j
#'
#' bayesmeta: DIRECT individual posteriors
#' \eqn{\widetilde m_j} (`bmr$shrink.moments`) and \eqn{\widetilde g_j}
#' (`bayesmeta$cond.moment(individual=)`), reconstructed as
#' \eqn{\widetilde\mu_{Aj}=\widetilde m_j-\pi_j\widetilde g_j},
#' \eqn{\widetilde\mu_{Bj}=\widetilde m_j+(1-\pi_j)\widetilde g_j}.
#' metafor: empirical-Bayes BLUPs of the same orthogonal blocks.
#' @noRd
comparison_shrinkage_native_trial_latents <- function(fit,
                                                post,
                                                trials,
                                                pi,
                                                reference,
                                                nonref) {
  if (interactmeta_is_bayesian_fit(fit) &&
      !is.null(fit$latent) &&
      !is.null(fit$latent$m) &&
      !is.null(fit$latent$g)) {
    return(comparison_shrinkage_from_bayesmeta_latent(
      fit = fit,
      post = post,
      trials = trials,
      pi = pi
    ))
  }
  if (!interactmeta_is_bayesian_fit(fit) &&
      (!is.null(fit$metafor) || !is.null(fit$latent))) {
    return(comparison_shrinkage_from_metafor_blup(
      fit = fit,
      post = post,
      trials = trials,
      pi = pi
    ))
  }
  NULL
}

#' Align posterior rows to a full DIRECT product grid
#' @noRd
comparison_shrinkage_draw_index <- function(post, n_full) {
  n_draw <- nrow(post)
  draw_id <- as.integer(post$.draw)
  if (length(draw_id) == n_draw &&
      !anyNA(draw_id) &&
      min(draw_id) >= 1L &&
      max(draw_id) <= n_full) {
    return(draw_id)
  }
  if (n_draw == n_full) {
    return(seq_len(n_full))
  }
  stop(
    "Cannot align trial-level shrinkage nodes with posterior draws.",
    call. = FALSE
  )
}

#' Product-grid trial latents from bayesmeta DIRECT individual posteriors
#' @noRd
comparison_shrinkage_from_bayesmeta_latent <- function(fit, post, trials, pi) {
  lat <- fit$latent
  n_draw <- nrow(post)
  y_A <- comparison_named_trial_matrix(n_draw, trials)
  y_B <- comparison_named_trial_matrix(n_draw, trials)
  m <- comparison_named_trial_matrix(n_draw, trials)
  g <- comparison_named_trial_matrix(n_draw, trials)
  sd_A <- comparison_named_trial_matrix(n_draw, trials, 0)
  sd_B <- comparison_named_trial_matrix(n_draw, trials, 0)
  sd_m <- comparison_named_trial_matrix(n_draw, trials, 0)
  sd_g <- comparison_named_trial_matrix(n_draw, trials, 0)

  gm0_m <- lat$m[[trials[[1L]]]]$gm
  gm0_g <- lat$g[[trials[[1L]]]]$gm
  if (is.null(gm0_m) || is.null(gm0_g)) {
    stop("bayesmeta latent shrinkage is missing Gaussian-mixture nodes.", call. = FALSE)
  }
  n_m <- length(gm0_m$mean)
  n_g <- length(gm0_g$mean)
  n_prod <- n_m * n_g
  i_ov <- rep(seq_len(n_m), each = n_g)
  i_int <- rep(seq_len(n_g), times = n_m)
  idx <- comparison_shrinkage_draw_index(post, n_prod)

  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    gm_m <- lat$m[[tr]]$gm
    gm_g <- lat$g[[tr]]$gm
    if (is.null(gm_m) || is.null(gm_g)) {
      stop("Missing latent m/g mixture for trial `", tr, "`.", call. = FALSE)
    }
    if (length(gm_m$mean) != n_m || length(gm_g$mean) != n_g) {
      stop(
        "Latent mixture node counts differ across trials; cannot form a product grid.",
        call. = FALSE
      )
    }
    pi_j <- unname(pi[[tr]])
    m_full <- gm_m$mean[i_ov]
    g_full <- gm_g$mean[i_int]
    sd_m_full <- gm_m$sd[i_ov]
    sd_g_full <- gm_g$sd[i_int]
    yA_full <- interactmeta_muA_from_mg(m_full, g_full, pi_j)
    yB_full <- interactmeta_muB_from_mg(m_full, g_full, pi_j)
    m[, j] <- m_full[idx]
    g[, j] <- g_full[idx]
    y_A[, j] <- yA_full[idx]
    y_B[, j] <- yB_full[idx]
    sd_m[, j] <- sd_m_full[idx]
    sd_g[, j] <- sd_g_full[idx]
    sd_A[, j] <- sqrt(pmax(sd_m_full[idx]^2 + (pi_j * sd_g_full[idx])^2, 0))
    sd_B[, j] <- sqrt(pmax(sd_m_full[idx]^2 + ((1 - pi_j) * sd_g_full[idx])^2, 0))
  }
  list(
    y_A = y_A, y_B = y_B, m = m, g = g,
    sd_A = sd_A, sd_B = sd_B, sd_m = sd_m, sd_g = sd_g
  )
}

#' Match a BLUP row to a trial id
#' @noRd
comparison_blup_row <- function(bl, trial, slab = NULL) {
  bl <- as.data.frame(bl)
  rn <- rownames(bl)
  ix <- NA_integer_
  if (!is.null(rn) && trial %in% rn) {
    ix <- match(trial, rn)[[1L]]
  } else if (!is.null(slab)) {
    sl <- as.character(slab)
    if (trial %in% sl) {
      ix <- match(trial, sl)[[1L]]
    }
  }
  if (!is.finite(ix)) {
    stop("No metafor BLUP row for trial `", trial, "`.", call. = FALSE)
  }
  pred <- if ("pred" %in% names(bl)) {
    as.numeric(bl$pred[[ix]])
  } else {
    as.numeric(bl[[1L]][[ix]])
  }
  se <- if ("se" %in% names(bl)) {
    as.numeric(bl$se[[ix]])
  } else {
    NA_real_
  }
  list(pred = pred, se = pmax(se, 0, na.rm = TRUE))
}

#' Frequentist trial-level CAMS latents from metafor BLUPs
#' @noRd
comparison_shrinkage_from_metafor_blup <- function(fit, post, trials, pi) {
  mf <- fit$metafor
  if (is.null(mf$overall) || is.null(mf$interaction)) {
    return(NULL)
  }
  bl_m <- metafor::blup(mf$overall)
  bl_g <- metafor::blup(mf$interaction)
  n_draw <- nrow(post)
  y_A <- comparison_named_trial_matrix(n_draw, trials)
  y_B <- comparison_named_trial_matrix(n_draw, trials)
  m <- comparison_named_trial_matrix(n_draw, trials)
  g <- comparison_named_trial_matrix(n_draw, trials)
  sd_A <- comparison_named_trial_matrix(n_draw, trials, 0)
  sd_B <- comparison_named_trial_matrix(n_draw, trials, 0)
  sd_m <- comparison_named_trial_matrix(n_draw, trials, 0)
  sd_g <- comparison_named_trial_matrix(n_draw, trials, 0)
  slab_m <- mf$overall$slab
  slab_g <- mf$interaction$slab
  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    pi_j <- unname(pi[[tr]])
    bm <- comparison_blup_row(bl_m, tr, slab_m)
    bg <- comparison_blup_row(bl_g, tr, slab_g)
    m[, j] <- bm$pred
    g[, j] <- bg$pred
    y_A[, j] <- interactmeta_muA_from_mg(bm$pred, bg$pred, pi_j)
    y_B[, j] <- interactmeta_muB_from_mg(bm$pred, bg$pred, pi_j)
    sd_m[, j] <- bm$se
    sd_g[, j] <- bg$se
    sd_A[, j] <- sqrt(pmax(bm$se^2 + (pi_j * bg$se)^2, 0))
    sd_B[, j] <- sqrt(pmax(bm$se^2 + ((1 - pi_j) * bg$se)^2, 0))
  }
  list(
    y_A = y_A, y_B = y_B, m = m, g = g,
    sd_A = sd_A, sd_B = sd_B, sd_m = sd_m, sd_g = sd_g
  )
}

#' Assert trial π is the non-reference information fraction (not 1-π)
#'
#' CAMS defines
#' \eqn{\pi_j = I_B/(I_A+I_B)} for non-reference subgroup B, with reference A
#' carrying weight \eqn{1-\pi_j} in
#' \eqn{m_j=(1-\pi_j)y_{A,j}+\pi_j y_{B,j}}. Design columns must satisfy
#' `ifrac_<B>` = π on both rows and `ifrac2_<B>` = x - π with x∈{0,1}.
#'
#' If the design π equals the **reference** precision share (and not the
#' non-reference share), that is treated as a sign/orientation error.
#' Hand-crafted design columns that match neither share are allowed (tests /
#' sensitivity analyses) as long as the `ifrac` / `ifrac2` geometry is valid.
#' @noRd
comparison_assert_pi_is_nonref <- function(pi,
                                     v_A,
                                     v_B,
                                     pi_A,
                                     pi_B,
                                     i2_A,
                                     i2_B,
                                     reference,
                                     nonref,
                                     ifrac_var,
                                     tol = 1e-8) {
  trials <- names(pi)
  for (tr in trials) {
    pij <- unname(pi[[tr]])
    if (!is.finite(pij) || pij < -tol || pij > 1 + tol) {
      stop(
        "Internal error: trial `", tr, "` has non-finite or out-of-range pi.",
        call. = FALSE
      )
    }
    if (abs(pi_A[[tr]] - pij) > tol || abs(pi_B[[tr]] - pij) > tol) {
      stop(
        "Internal error: design `", ifrac_var, "` for trial `", tr,
        "` is not constant and equal to non-reference pi on both subgroup rows.",
        call. = FALSE
      )
    }
    if (abs(i2_A[[tr]] - (0 - pij)) > tol ||
        abs(i2_B[[tr]] - (1 - pij)) > tol) {
      stop(
        "Internal error: `ifrac2` for trial `", tr,
        "` is not x - pi with x = 0 on reference `", reference,
        "` and x = 1 on non-reference `", nonref,
        "` (pi must be IF_", nonref, ", not 1 - IF_", nonref, ").",
        call. = FALSE
      )
    }

    ia <- 1 / v_A[[tr]]
    ib <- 1 / v_B[[tr]]
    if (!is.finite(ia) || !is.finite(ib) || ia <= 0 || ib <= 0) {
      next
    }
    pi_nonref <- ib / (ia + ib)
    pi_ref <- ia / (ia + ib)
    matches_nonref <- abs(pij - pi_nonref) <= tol
    matches_ref <- abs(pij - pi_ref) <= tol
    # True flip: design stores the reference share as ifrac_<nonref>
    if (matches_ref && !matches_nonref) {
      stop(
        "Internal error: trial `", tr, "` design `", ifrac_var,
        "` equals the reference (`", reference, "`) information fraction ",
        format(pi_ref, digits = 6),
        " instead of the non-reference (`", nonref, "`) fraction ",
        format(pi_nonref, digits = 6),
        ". Do not use 1 - IF_", nonref, " as `", ifrac_var, "`.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' One information fraction per trial
#' @noRd
comparison_trial_information_fractions <- function(data,
                                             trial_var,
                                             ifrac_var,
                                             trials) {
  # Prefer shared overall_ifrac_trial_pi when subgroup column exists
  subgroup_var <- NULL
  for (cand in c("level", "subgroup", "group", "sg")) {
    if (cand %in% names(data)) {
      subgroup_var <- cand
      break
    }
  }
  if (!is.null(subgroup_var) && ifrac_var %in% names(data)) {
    tp <- overall_ifrac_trial_pi(
      data = data,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = NA_character_,
      ifrac_var = ifrac_var
    )
    pi <- stats::setNames(as.numeric(tp$pi), as.character(tp$trial))
    pi <- pi[trials]
    if (length(pi) == length(trials) && all(is.finite(pi))) {
      names(pi) <- trials
      return(pi)
    }
  }
  pi <- stats::setNames(rep(NA_real_, length(trials)), trials)
  for (tr in trials) {
    ii <- which(as.character(data[[trial_var]]) == tr)
    vals <- as.numeric(data[[ifrac_var]][ii])
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      stop("No finite IF values for trial `", tr, "`.", call. = FALSE)
    }
    uniq <- unique(vals)
    if (length(uniq) == 1L) {
      pi[[tr]] <- uniq[[1L]]
    } else {
      pi[[tr]] <- mean(vals)
    }
  }
  if (any(!is.finite(pi))) {
    stop("Could not resolve trial information fractions.", call. = FALSE)
  }
  pi
}

#' Long draws data frame from a draws × trial matrix
#' @noRd
comparison_rotation_matrix_to_long <- function(mat,
                                         post,
                                         trials,
                                         extra = NULL,
                                         type = "Shrinkage",
                                         sd_mat = NULL) {
  rows <- lapply(seq_along(trials), function(j) {
    tr <- trials[[j]]
    df <- data.frame(
      .draw = post$.draw,
      .chain = post$.chain,
      .iteration = post$.iteration,
      row = j,
      trial = tr,
      type = type,
      effect = mat[, j],
      stringsAsFactors = FALSE
    )
    if (!is.null(sd_mat)) {
      df$.cond_sd <- as.numeric(sd_mat[, j])
    }
    if (!is.null(extra)) {
      for (nm in names(extra)) {
        val <- extra[[nm]]
        df[[nm]] <- if (length(val) == 1L) {
          val
        } else if (!is.null(names(val)) && tr %in% names(val)) {
          val[[tr]]
        } else if (length(val) == length(trials)) {
          val[[j]]
        } else {
          val
        }
      }
    }
    df
  })
  interactmeta_decorate_effect_draws(dplyr::bind_rows(rows), post)
}

#' Component bundle attached to shrinkage objects
#' @noRd
comparison_shrinkage_components_list <- function(bundle) {
  list(
    pi = bundle$pi,
    omega = bundle$omega,
    evaluation = bundle$evaluation %||% "trial",
    evaluation_if = bundle$evaluation_if,
    reference = bundle$reference,
    nonref = bundle$nonref,
    trials = bundle$trials,
    observed = list(
      y_A = bundle$y_A_obs,
      y_B = bundle$y_B_obs,
      m = bundle$m_obs,
      g = bundle$g_obs,
      variance_A = bundle$v_A,
      variance_B = bundle$v_B,
      variance_m = bundle$v_m,
      variance_g = bundle$v_g
    ),
    # support × trial matrices on the model scale (π_j + native trial latents)
    support = list(
      y_A = bundle$y_A,
      y_B = bundle$y_B,
      m = bundle$m,
      g = bundle$g
    ),
    sd = list(
      y_A = bundle$sd_A,
      y_B = bundle$sd_B,
      m = bundle$sd_m,
      g = bundle$sd_g
    ),
    weight = if (is.data.frame(bundle$post) && ".weight" %in% names(bundle$post)) {
      bundle$post$.weight
    } else {
      NULL
    },
    ifrac_var = bundle$ifrac_var,
    interaction_label = bundle$interaction_label
  )
}
