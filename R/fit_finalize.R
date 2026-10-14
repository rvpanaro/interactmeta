#' CAMS internal parameter names
#'
#' Developer glossary. Public summaries still use descriptive labels
#' (`Ecological_*`, `Interaction_*`, `tau`, …).
#'
#' \describe{
#'   \item{`alpha`}{Combined intercept; population mean at information fraction 0 for the reference subgroup.}
#'   \item{`beta`}{Ecological / composition slope (`b_ifrac`).}
#'   \item{`gamma`}{Interaction / subgroup offset (`b_subgroup`).}
#'   \item{`theta`}{Combined information-fraction slope; CAMS identity `theta = beta + gamma`.}
#'   \item{`tau_m`}{Combined-block heterogeneity SD (`sd_trial__Intercept`).}
#'   \item{`tau_g`}{Interaction-block heterogeneity SD (`sd_trial__ifrac2_*`).}
#'   \item{`pi_j`}{Trial information fraction of the non-reference subgroup.}
#'   \item{`m_j`}{Trial overall (IV-pooled) effect.}
#'   \item{`g_j`}{Trial interaction contrast `y_B - y_A`.}
#' }
#'
#' @keywords internal
#' @name cam-parameter-names
NULL

#' Attach call, model spec, priors, and trial blocks to a CAMS fit
#' @noRd
interactmeta_finalize_fit <- function(fit,
                                      call,
                                      design,
                                      data,
                                      data_fit,
                                      engine,
                                      method,
                                      priors,
                                      engine_args,
                                      seed,
                                      ecological_slope,
                                      big_var,
                                      contrast,
                                      trial_var,
                                      subgroup_var) {
  algorithm <- fit$algorithm %||% interactmeta_algorithm_for(engine, method)
  fit$call <- call
  fit$model_spec <- interactmeta_build_model_spec(
    design = design,
    engine = engine,
    method = method,
    algorithm = algorithm,
    ecological_slope = ecological_slope,
    contrast = contrast,
    big_var = big_var,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  fit$priors <- interactmeta_wrap_priors(priors, engine, design)
  fit$computation <- interactmeta_build_computation(
    engine = engine,
    method = method,
    engine_args = engine_args,
    seed = seed
  )
  fit$transformed_data <- interactmeta_build_transformed_data(
    fit = fit,
    design = design,
    data = data,
    data_fit = data_fit,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  interactmeta_attach_capabilities(fit, "cam", engine, method)
}

interactmeta_build_model_spec <- function(design,
                                          engine,
                                          method,
                                          algorithm,
                                          ecological_slope,
                                          contrast,
                                          big_var,
                                          trial_var,
                                          subgroup_var) {
  parameterization <- switch(
    engine,
    bayesmeta = "orthogonal_blocks (overall bmr + interaction bayesmeta)",
    metafor = "orthogonal_blocks (overall rma + interaction rma)",
    brms = "joint_mixed_model (brms / Stan)",
    engine
  )
  heterogeneity <- list(
    tau_m = "overall block / trial intercept SD",
    tau_g = "interaction block / ifrac2 slope SD",
    correlated = FALSE
  )
  structure(
    list(
      family = "CAMS",
      n_subgroups = design$n_subgroups,
      n_nonref = design$n_nonref,
      reference = design$reference,
      subgroup_levels = design$subgroup_levels,
      non_ref_labels = design$non_ref_labels,
      information_fraction = list(
        definition = "precision weight: I_g = 1/vi; pi_g = I_g / sum_h I_h",
        contrast = contrast,
        ifrac_vars = design$ifrac_vars,
        ifrac2_vars = design$ifrac2_vars
      ),
      ecological_slope = isTRUE(ecological_slope),
      random_effects = "independent overall-block and interaction-block trial effects",
      random_effects_covariance = "diagonal",
      heterogeneity = heterogeneity,
      parameterization = parameterization,
      method = method,
      engine = engine,
      algorithm = algorithm,
      likelihood = "yi ~ Normal(mu, se = sqrt(vi)); vi treated as known",
      missing_cell_variance = big_var,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      formula = design$formula
    ),
    class = "interactmeta_model_spec"
  )
}

#' @export
print.interactmeta_model_spec <- function(x, ...) {
  cat("CAMS model specification\n")
  cat("  subgroups          : ", x$n_subgroups, " (reference = ", x$reference, ")\n", sep = "")
  cat("  ecological slope   : ", x$ecological_slope, "\n", sep = "")
  cat("  parameterization   : ", x$parameterization, "\n", sep = "")
  cat("  method             : ", x$method, "\n", sep = "")
  cat("  likelihood         : ", x$likelihood, "\n", sep = "")
  invisible(x)
}

interactmeta_build_computation <- function(engine, method, engine_args, seed) {
  args <- if (is.null(engine_args)) list() else engine_args
  if ("convolve_derived" %in% names(args)) {
    args$convolve_derived <- NULL
  }
  if (identical(engine, "bayesmeta")) {
    if (is.null(args$delta)) args$delta <- 0.01
    if (is.null(args$epsilon)) args$epsilon <- 1e-4
  }
  stochastic <- identical(engine, "brms")
  list(
    engine_args = args,
    seed = if (is.na(seed)) NA else seed,
    stochastic = stochastic,
    deterministic = !stochastic,
    bitwise_reproducible = FALSE
  )
}

interactmeta_build_transformed_data <- function(fit,
                                                design,
                                                data,
                                                data_fit,
                                                trial_var,
                                                subgroup_var) {
  observed <- comparison_fit_display_data(data, design)
  blocks <- NULL
  if (isTRUE(design$n_nonref == 1L)) {
    blocks <- comparison_trial_effect_blocks(
      data = observed,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = design$reference,
      nonref = design$non_ref_labels[[1L]],
      ifrac_var = design$ifrac_vars[[1L]]
    )
  }
  list(
    observed_effects = observed[, intersect(
      c(trial_var, subgroup_var, "yi", "vi", design$ifrac_vars, design$ifrac2_vars),
      names(observed)
    ), drop = FALSE],
    trial_blocks = blocks
  )
}

interactmeta_wrap_priors <- function(priors, engine, design) {
  if (!inherits(priors, "cam_priors")) {
    return(priors)
  }
  engine_rep <- interactmeta_priors_engine_representation(priors, engine, design)
  structure(
    list(
      common = priors,
      engine = engine_rep,
      tau = priors$tau,
      tau_interaction = priors$tau_interaction,
      fixed = priors$fixed
    ),
    class = c("interactmeta_priors", "cam_priors")
  )
}

#' @export
print.interactmeta_priors <- function(x, ...) {
  if (inherits(x$common, "cam_priors")) {
    print(x$common)
  } else {
    print.cam_priors(x)
  }
  invisible(x)
}

interactmeta_priors_engine_representation <- function(priors, engine, design) {
  if (identical(engine, "brms")) {
    return(list(
      applied = TRUE,
      note = "brms::set_prior rows translated from the common CAMS priors",
      rows = prior_to_brms_rows(
        priors,
        trial_var = design$trial_var %||% "trial",
        ifrac2_vars = design$coef_map$ifrac2_std
      )
    ))
  }
  if (identical(engine, "bayesmeta")) {
    tau_bm <- prior_to_bayesmeta(priors$tau)
    tau_g_bm <- prior_to_bayesmeta(priors$tau_interaction)
    fe <- prior_to_bayesmeta_fixed(priors$fixed)
    return(list(
      applied = TRUE,
      note = paste(
        "bmr prior is on (alpha, theta), with beta.prior.cov = diag(sd^2);",
        "mu.prior is for gamma only.",
        "Independent normals on (alpha, theta) are not independent normals on",
        "(alpha, beta, gamma), because theta = beta + gamma."
      ),
      tau = if (is.character(tau_bm)) tau_bm else format(priors$tau),
      tau_interaction = if (is.character(tau_g_bm)) tau_g_bm else format(priors$tau_interaction),
      fixed = fe,
      tau_density = if (is.function(tau_bm)) tau_bm else NULL,
      tau_interaction_density = if (is.function(tau_g_bm)) tau_g_bm else NULL
    ))
  }
  list(
    applied = FALSE,
    note = "engine = metafor is frequentist; Bayesian priors are not used"
  )
}
