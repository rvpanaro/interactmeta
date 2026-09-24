#' CAMS internal parameter names
#'
#' Developer glossary. Public summaries still use descriptive labels
#' (`Ecological_*`, `Interaction_*`, `tau`, …).
#'
#' \describe{
#'   \item{`alpha`}{Overall intercept; population mean at information fraction 0 for the reference subgroup.}
#'   \item{`beta`}{Ecological / composition slope (`b_ifrac`).}
#'   \item{`gamma`}{Interaction / subgroup offset (`b_subgroup`).}
#'   \item{`theta`}{Overall information-fraction slope; CAMS identity `theta = beta + gamma`.}
#'   \item{`tau_m`}{Overall-block heterogeneity SD (`sd_trial__Intercept`).}
#'   \item{`tau_g`}{Interaction-block heterogeneity SD (`sd_trial__ifrac2_*`).}
#'   \item{`pi_j`}{Trial information fraction of the non-reference subgroup.}
#'   \item{`m_j`}{Trial overall (IV-pooled) effect.}
#'   \item{`g_j`}{Trial interaction contrast `y_B - y_A`.}
#' }
#'
#' @keywords internal
#' @name cam-parameter-names
NULL

#' Statistical audit of a CAMS fit
#'
#' Summarizes provenance needed to inspect or reproduce a [cam()]
#' fit: the call, model specification, priors, computational settings,
#' assumptions, transformed trial-level data, and software versions.
#' Does **not** refit the model.
#'
#' @param x An `interactmeta` from [cam()].
#' @param ... Ignored.
#' @return An object of class `"interactmeta_audit"` with components
#'   `call`, `model`, `engine`, `priors`, `computation`, `assumptions`,
#'   `transformed_data`, `draws_semantics`, `derived_inference`,
#'   `warnings`, and `software`.
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   trial = rep(c("T1", "T2"), each = 2),
#'   subgroup = rep(c("A", "B"), 2),
#'   yi = c(0.10, 0.30, 0.12, 0.28),
#'   vi = c(0.04, 0.04, 0.05, 0.03)
#' )
#' fit <- cam(dat, reference = "A", engine = "bayesmeta", silent = 2)
#' a <- audit(fit)
#' a$model
#' a$priors
#' a$transformed_data$trial_blocks
#' }
#' @seealso [cam()], [interactmeta-auditability]
#' @export
audit <- function(x, ...) {
  UseMethod("audit")
}

#' @rdname audit
#' @export
audit.interactmeta <- function(x, ...) {
  check_comparison_fit(x)
  st <- attr(x, "interactmeta", exact = TRUE)
  model <- x$model_spec %||% interactmeta_model_spec_from_fit(x)
  engine <- x$engine %||% st$engine
  method <- x$method %||% st$method
  algorithm <- x$algorithm %||% st$algorithm
  priors <- x$priors
  common <- if (is.list(priors) && !is.null(priors$common)) priors$common else priors
  warnings <- interactmeta_audit_warnings(x, engine, method, model)
  structure(
    list(
      call = x$call,
      model = model,
      engine = list(
        engine = engine,
        method = method,
        algorithm = algorithm
      ),
      priors = priors,
      computation = x$computation,
      assumptions = x$assumptions,
      transformed_data = x$transformed_data,
      draws_semantics = x$draws_semantics,
      derived_inference = interactmeta_derived_inference(engine, model),
      warnings = warnings,
      software = x$provenance
    ),
    class = "interactmeta_audit"
  )
}

#' @export
print.interactmeta_audit <- function(x, ...) {
  eng <- x$engine$engine %||% NA_character_
  meth <- x$engine$method %||% NA_character_
  algo <- x$engine$algorithm %||% NA_character_
  model <- x$model
  cat("CAMS statistical audit\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n  ")
    print(x$call)
  }
  cat("\nModel:\n")
  cat("  family           : CAMS\n")
  cat("  subgroups        : ", model$n_subgroups %||% NA, " (reference = ",
      model$reference %||% NA, ")\n", sep = "")
  cat("  ecological slope : ", model$ecological_slope %||% NA, "\n", sep = "")
  cat("  random effects   : ", model$random_effects %||% NA, "\n", sep = "")
  cat("  parameterization : ", model$parameterization %||% NA, "\n", sep = "")
  cat("\nEngine:\n")
  cat("  engine           : ", eng, "\n", sep = "")
  cat("  method           : ", meth, "\n", sep = "")
  cat("  algorithm        : ", algo, "\n", sep = "")
  if (identical(x$computation$stochastic, TRUE)) {
    cat("  stochastic       : TRUE (seed = ",
        format(x$computation$seed), ")\n", sep = "")
  } else {
    cat("  stochastic       : FALSE (deterministic)\n")
  }
  cat("\nPriors / estimator:\n")
  if (identical(eng, "metafor") || identical(meth, "REML") || identical(meth, "ML")) {
    cat("  frequentist estimator: ", meth, "\n", sep = "")
    cat("  Bayesian priors stored but not used\n")
  } else if (inherits(x$priors, "cam_priors") || inherits(x$priors$common, "cam_priors")) {
    pr <- if (inherits(x$priors$common, "cam_priors")) x$priors$common else x$priors
    cat("  tau_m (overall)     : ", format(pr$tau), "\n", sep = "")
    cat("  tau_g (interaction) : ", format(pr$tau_interaction), "\n", sep = "")
    cat("  fixed effects       : ", format(pr$fixed), "\n", sep = "")
  } else {
    cat("  (see $priors)\n")
  }
  if (!is.null(x$assumptions$subgroup_covariance)) {
    cat("\nAssumptions:\n")
    cat("  subgroup covariance : ", x$assumptions$subgroup_covariance, "\n", sep = "")
    cat("  RE covariance       : ", x$assumptions$random_effects_covariance, "\n", sep = "")
    cat("  mu includes REs     : ", x$assumptions$mu_includes_random_effects, "\n", sep = "")
  }
  if (length(x$warnings)) {
    cat("\nWarnings / limitations:\n")
    for (w in x$warnings) {
      cat("  - ", w, "\n", sep = "")
    }
  }
  if (!is.null(x$software)) {
    cat("\nSoftware:\n")
    cat("  R                   : ", x$software$r_version %||% NA, "\n", sep = "")
    pkgs <- x$software$packages
    if (length(pkgs)) {
      for (nm in names(pkgs)) {
        cat("  ", format(nm, width = 20L), ": ", pkgs[[nm]], "\n", sep = "")
      }
    }
  }
  cat("\nInspect $model, $priors, $computation, $assumptions, ",
      "$transformed_data, $draws_semantics for structured fields.\n",
      sep = "")
  invisible(x)
}

#' Attach provenance metadata to a CAMS fit (does not change estimates)
#' @noRd
interactmeta_attach_audit_metadata <- function(fit,
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
  model_spec <- interactmeta_build_model_spec(
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
  computation <- interactmeta_build_computation(
    engine = engine,
    method = method,
    engine_args = engine_args,
    seed = seed
  )
  transformed <- interactmeta_build_transformed_data(
    fit = fit,
    design = design,
    data = data,
    data_fit = data_fit,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  assumptions <- interactmeta_build_assumptions(
    engine = engine,
    method = method,
    design = design,
    big_var = big_var,
    ecological_slope = ecological_slope
  )
  draws_semantics <- interactmeta_draws_semantics(engine, method)
  priors_wrapped <- interactmeta_wrap_priors(priors, engine, design)

  fit$call <- call
  fit$model_spec <- model_spec
  fit$priors <- priors_wrapped
  fit$computation <- computation
  fit$assumptions <- assumptions
  fit$transformed_data <- transformed
  fit$draws_semantics <- draws_semantics
  fit$provenance <- interactmeta_build_provenance(engine, seed)
  fit <- interactmeta_attach_capabilities(fit, "cam", engine, method)
  fit
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
    trial_blocks = blocks,
    identities = list(
      overall = "m_j = sum_g w_g y_gj with w_g = 1/v_gj (IV pool)",
      overall_variance = "Var(m_j) = 1 / sum_g 1/v_gj",
      interaction = "g_j = y_Bj - y_Aj",
      interaction_variance = "Var(g_j) = v_Aj + v_Bj (subgroup covariance assumed zero)",
      information_fraction = "pi_j = (1/v_Bj) / (1/v_Aj + 1/v_Bj)",
      orthogonality = "Cov(m_j, g_j) = 0 when pi_j is the IV weight of B",
      slope = "theta = beta + gamma",
      inverse = "mu_Aj = m_j - pi_j * g_j; mu_Bj = m_j + (1 - pi_j) * g_j"
    )
  )
}

interactmeta_build_assumptions <- function(engine,
                                           method,
                                           design,
                                           big_var,
                                           ecological_slope) {
  list(
    known_sampling_variances = TRUE,
    subgroup_covariance = "assumed_zero",
    random_effects_covariance = "uncorrelated",
    overall_interaction_independence = engine %in% c("bayesmeta", "metafor"),
    missing_cells_imputed = TRUE,
    missing_cell_variance = big_var,
    ecological_slope = isTRUE(ecological_slope),
    mu_includes_random_effects = identical(engine, "brms"),
    priors_applied = !identical(engine, "metafor")
  )
}

interactmeta_draws_semantics <- function(engine, method) {
  if (identical(engine, "bayesmeta")) {
    return(list(
      row_type = "quadrature_node",
      iid_posterior_draws = FALSE,
      weight_column = ".weight",
      weight_meaning = "product of DIRECT tau-node weights from overall and interaction blocks",
      cell_meaning = "conditional mean of the parameter given the tau node; within-tau variance is in $latent and .cond_* columns",
      use_weight_for_moments = TRUE,
      public_summary_uses_weight = TRUE
    ))
  }
  if (identical(engine, "metafor")) {
    return(list(
      row_type = "point_estimate",
      iid_posterior_draws = FALSE,
      weight_column = NA_character_,
      cell_meaning = "point estimate (one row)",
      use_weight_for_moments = FALSE,
      public_summary_uses_weight = FALSE
    ))
  }
  list(
    row_type = "mcmc_draw",
    iid_posterior_draws = TRUE,
    weight_column = NA_character_,
    cell_meaning = "posterior draw (equal weight)",
    use_weight_for_moments = FALSE,
    public_summary_uses_weight = FALSE
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

interactmeta_build_provenance <- function(engine, seed) {
  pkgs <- c("interactmeta", "bayesmeta", "metafor", "posterior")
  if (identical(engine, "brms")) {
    pkgs <- c(pkgs, "brms")
  }
  versions <- vapply(pkgs, function(p) {
    if (!requireNamespace(p, quietly = TRUE)) {
      return(NA_character_)
    }
    as.character(utils::packageVersion(p))
  }, character(1))
  list(
    r_version = as.character(getRversion()),
    packages = as.list(versions),
    engine = engine,
    seed = if (is.na(seed)) NA else seed,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}

interactmeta_derived_inference <- function(engine, model) {
  list(
    beta = "theta - gamma (orthogonal engines) or b_ifrac (joint mixed model)",
    theta = "beta + gamma",
    subgroup_reference = "alpha + beta * pi (population); plus trial RE if mu includes REs",
    subgroup_nonreference = "alpha + beta * pi + gamma",
    inverse_cams = "mu_A = m - pi * g; mu_B = m + (1 - pi) * g",
    intervals = if (identical(engine, "metafor")) {
      "delta method from the stored rma coefficient covariance (same fit)"
    } else if (identical(engine, "bayesmeta")) {
      paste(
        "shortest or central intervals from the stored bmr/bayesmeta $post.interval,",
        "and from bayesmeta::convolve() for sums of independent blocks"
      )
    } else {
      "equal-tailed posterior intervals from this fit's draws (credible intervals)"
    },
    mu = if (identical(engine, "brms") || identical(model$engine, "brms")) {
      "conditional linear predictor including trial random effects"
    } else {
      "fixed-effect linear predictor only; shrunken trial effects in $latent (bayesmeta)"
    }
  )
}

interactmeta_audit_warnings <- function(x, engine, method, model) {
  out <- character()
  if (identical(engine, "bayesmeta")) {
    out <- c(
      out,
      paste(
        "bayesmeta $draws rows are quadrature nodes (conditional means given tau)",
        "with .weight and within-tau moments (.cond_*). Public intervals use",
        "$post.interval on one stored fit, or bayesmeta::convolve() for an",
        "independent-block sum. They are not quantiles of the conditional means.",
        "Alpha and theta stay on the joint overall-block support."
      )
    )
    out <- c(
      out,
      paste(
        "Fixed-effect priors: bmr uses independent N on (alpha, theta);",
        "the joint mixed model uses independent N on (alpha, beta, gamma).",
        "theta = beta + gamma, so the joints are not identical."
      )
    )
  }
  if (identical(engine, "metafor")) {
    out <- c(out, "metafor is frequentist; $draws is a single point-estimate row. Intervals for derived estimands use the stored rma covariance (delta method), not a second fit.")
  }
  if (!isTRUE(model$ecological_slope)) {
    out <- c(out, "ecological_slope = FALSE: beta is omitted (identically 0).")
  }
  if (isTRUE(x$assumptions$subgroup_covariance == "assumed_zero")) {
    out <- c(
      out,
      "Interaction variances use v_A + v_B (within-trial subgroup covariance assumed zero; not observed)."
    )
  }
  unique(out)
}

interactmeta_model_spec_from_fit <- function(x) {
  st <- attr(x, "interactmeta", exact = TRUE)
  list(
    family = "CAMS",
    reference = st$reference,
    ecological_slope = st$ecological_slope %||% TRUE,
    engine = x$engine %||% st$engine,
    method = x$method %||% st$method,
    n_subgroups = length(unique(as.character(x$data[[st$subgroup_var %||% "subgroup"]])))
  )
}
