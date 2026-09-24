#' Contribution-adjusted meta-analysis of treatment-by-subgroup interactions
#'
#' Fits the CAMS model (Panaro, Röver & Friede 2025) for **consistent
#' inference on treatment-by-subgroup interactions** from subgroup-level
#' aggregate data. Subgroup-specific treatment effects are inputs and
#' derived outputs; the methodological target is not a separate meta-analysis
#' within each subgroup.
#'
#' [cam()] is the preferred name for this fitter. [cams()] and
#' [subgroupmeta()] are aliases of the same function.
#'
#' By default, [ifrac()] is run inside the fitter on the input (e.g. a
#' [`metafor::escalc`] object). You can also pass an existing [ifrac()] tibble.
#'
#' The model for binary subgroups is the orthogonal CAMS factorization
#' (Panaro, Röver & Friede 2025, §3.4): an overall information-fraction
#' meta-regression of \(m_j\) plus an interaction random-effects meta-analysis
#' of \(g_j\). Bayesian CAMS is implemented using DIRECT quadrature via
#' **bayesmeta**. Frequentist CAMS is implemented using **metafor**. Markov
#' chain Monte Carlo (MCMC) is not required for CAMS inference.
#'
#' @section Model variants:
#' One call; switch with `ecological_slope`:
#' * `ecological_slope = TRUE` (default) — **comparison** model with free
#'   ecological / composition slope \eqn{\beta}
#' * `ecological_slope = FALSE` — **subgroup-specific** model with \eqn{\beta}
#'   omitted (identically 0); preferred for source-only fits under
#'   [reproducibility()] same-IF. Currently implemented for internal Stan
#'   backends only, not for public Bayesian or frequentist CAMS.
#'
#' @section Inference:
#' CAMS can be fitted using Bayesian or frequentist inference. These are not
#' two different scientific models. The approaches are **bayesmeta** (DIRECT,
#' the default) and **metafor** (REML or ML). See `vignette("engine-comparison")`.
#'
#' Computational knobs (quadrature controls, optimizer settings) go in
#' `engine_args`, not as model arguments.
#'
#' [swada()], [da()], and [ad()] are **not** alternative CAMS models. SWADA
#' does not use the CAMS information-fraction (IF) construction; IF
#' quantities and IF-based plots are not available for SWADA objects.
#' Comparison helpers may combine multiple user-supplied fits only when
#' those fits are passed explicitly.
#'
#' @section Related analyses:
#' CAMS is the **statistical model**. [da()] (Difference of Averages) and
#' [ad()] (Average Difference) are **different constructions** (Panaro,
#' Röver & Friede 2026) for a treatment-by-subgroup interaction from
#' subgroup-level aggregate data. They can use different study weights and
#' need not agree. [swada()] applies the SWADA *weighting principle* so
#' those analyses share one contribution vector; under the default
#' Interaction RE-weights, the DA interaction equals AD. [conventional()]
#' is a compatibility bundle of DA and AD. None of these is the 2002 van
#' Houwelingen joint bivariate mixed model.
#'
#' @section Method:
#' `method` is the inferential method (`"bayes"`, `"REML"` / `"ML"`), not the
#' approach name (`bayesmeta` vs `metafor`). Prefer
#' [bayescam()] for Bayesian CAMS and `cam(..., engine = "metafor")` for
#' frequentist CAMS. The fit stores `fit$inference`, `fit$engine`, `fit$method`,
#' and `fit$algorithm` (`"quadrature"` or `"optimization"`).
#'
#' @section Priors:
#' Bayesian CAMS uses one CAMS-level prior specification (see [cam_prior()]),
#' stored on `fit$priors`. Defaults are:
#' \describe{
#'   \item{fixed effects}{\eqn{\mathrm{Normal}(0, 10)} on intercept, IF slope
#'     \eqn{\beta}, and subgroup offset \eqn{\gamma}}
#'   \item{`tau.prior`}{\eqn{\mathrm{Normal}^+(0, 0.5)} on overall
#'     heterogeneity \eqn{\tau_m}}
#'   \item{`tau.interaction.prior`}{\eqn{\mathrm{Normal}^+(0, 0.1)} on
#'     interaction heterogeneity \eqn{\tau_g}; if `NULL` and the user supplied
#'     `tau.prior`, the overall prior is inherited}
#' }
#' Character shorthands `"uniform"` and `"sqrt"` are the **improper** bayesmeta
#' densities \eqn{p(\tau)\propto 1} and \eqn{p(\tau)\propto \tau^{-1/2}}
#' (\eqn{\tau \ge 0}). They are **not** \eqn{U(0, U)}.
#' `engine = "metafor"` ignores Bayesian priors (warning if the user supplied
#' them).
#'
#' **Fixed-effect caveat:** bayesmeta/bmr place independent \eqn{N(0, 10)}
#' on \eqn{(\alpha,\theta)} with \eqn{\theta = \beta + \gamma} and on
#' \eqn{\gamma}. Stored `b_ifrac` is \eqn{\beta=\theta-\gamma}.
#'
#' @param x A [`metafor::escalc`] object, a data frame with `yi` / `vi`, or output
#'   from [ifrac()] (must include `ifrac_<k>` design columns). Also accepts a
#'   data frame built from study-level columns (`trial`, `subgroup`, `yi`, `vi`).
#' @param reference Subgroup label used as reference in [ifrac()] and in the
#'   model. Required when `x` is not already an [ifrac()] result.
#' @param modifier Column that holds the two levels of the within-study
#'   contrast (default `"subgroup"`, so the shipped examples keep working).
#'   Overrides `subgroup` and `subgroup_var` when set. Level labels stay as
#'   stored. Display headings capitalize the first letter of this name
#'   (`dose` is shown as `Dose`).
#' @param trial,subgroup When `x` is `escalc` / a data frame, column names passed
#'   to [ifrac()] (defaults `"trial"` and `"subgroup"`). `subgroup` remains
#'   accepted; `modifier` is the name to use.
#' @param trial_var,subgroup_var Column names used in the model (defaults
#'   match `trial` and `subgroup`). `subgroup_var` remains accepted.
#' @param contrast,contrast_sep,big_var,include_overall,include_interaction Passed
#'   to [ifrac()] when information fractions are computed inside `cam()`.
#' @param method Inferential method: `"bayes"` (default for Bayesian CAMS),
#'   `"REML"` or `"ML"` (frequentist CAMS). Not a backend name. If `engine` is
#'   omitted, `"bayes"` uses `engine = "bayesmeta"` and `"REML"` / `"ML"` use
#'   `engine = "metafor"`. An error is thrown only when `engine` and `method`
#'   are both specified and incompatible.
#' @param ecological_slope Logical; if `TRUE` (default), include the ecological /
#'   composition slope \eqn{\beta} (`b_ifrac`) — the **comparison** model. If
#'   `FALSE`, omit \eqn{\beta} — the **subgroup-specific** model (\eqn{\beta}
#'   identically 0). Useful as a source-only fit for [reproducibility()] under
#'   the same-information-fraction assumption.
#' @param approach Inferential approach. `"bayes"` selects Bayesian CAMS
#'   (`engine = "bayesmeta"` when `engine` is omitted). `"frequentist"`,
#'   `"REML"`, or `"ML"` select frequentist CAMS (`engine = "metafor"`).
#'   `NULL` (default) follows `engine` / `method`.
#' @param engine Approach used to fit CAMS: `"bayesmeta"` (Bayesian, default)
#'   or `"metafor"` (frequentist). Changing `engine` should not silently change
#'   the Bayesian model; use [cam_prior()] for shared priors. `"bayes"` is
#'   accepted as a synonym of `"bayesmeta"`. Omit `engine` when `method` or
#'   `approach` already identifies the fit.
#' @param tau.prior,tau.interaction.prior Heterogeneity priors for \eqn{\tau_m}
#'   (overall) and \eqn{\tau_g} (interaction). Character shorthands or
#'   [cam_prior()] objects. If `tau.interaction.prior = NULL` and the user
#'   supplied `tau.prior`, the interaction prior inherits `tau.prior`.
#' @param fixed.prior Shared normal (or improper-uniform) prior for fixed
#'   effects. See [cam_prior()].
#' @param engine_args List of **computational** options for the selected backend
#'   (e.g. quadrature `delta` / `epsilon`, metafor `test`). Model choices such
#'   as priors do not belong here.
#' @param silent Sampler verbosity; default `2` suppresses messages.
#' @param seed Random seed for reproducible draws (default from option
#'   `"interactmeta.seed"`, currently `42L`). Set to `NA` for a non-reproducible draw.
#'   If `seed` is also passed in `...`, a warning is given and this argument wins.
#' @param ... Further **computational** arguments merged into `engine_args`
#'   when not already set.
#' @return A CAMS fit with class `c("interactmeta", "cam", "cams")` plus attribute
#'   `"interactmeta"` (plot settings and `ifrac_data` for [forestplot()]).
#'
#' @examples
#' \dontrun{
#' library(interactmeta)
#'
#' es <- escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
#' fit <- bayescam(es, reference = "Younger")
#' fit_reml <- cam(es, reference = "Younger", engine = "metafor")
#' forestplot(fit)
#' plot_contribution(fit)
#' plot(fit)
#' }
#'
#' @references
#' Panaro RV, Röver C, Friede T (2025).
#' Consistent Bayesian meta-analysis on subgroup specific effects and
#' interactions. arXiv:2512.18785. \doi{10.48550/arXiv.2512.18785}
#'
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \doi{10.48550/arXiv.2508.15531}
#' @seealso [bayescam()], [cam()], [conventional()], [da()], [ad()], [swada()], [forestplot()], [reproducibility()], [audit()]
#' @aliases cams subgroupmeta
#' @export
cam <- function(x,
                 reference = NULL,
                 modifier = NULL,
                 trial = NULL,
                 subgroup = NULL,
                 trial_var = "trial",
                 subgroup_var = "subgroup",
                 contrast = "all",
                 contrast_sep = "-",
                 big_var = 100^2,
                 include_overall = FALSE,
                 include_interaction = FALSE,
                 method = c("bayes", "REML", "ML"),
                 ecological_slope = TRUE,
                 approach = NULL,
                 engine = c("bayesmeta", "metafor"),
                 tau.prior = NULL,
                 tau.interaction.prior = NULL,
                 fixed.prior = NULL,
                 engine_args = list(),
                 silent = 2,
                 seed = getOption("interactmeta.seed", 42L),
                 ...) {
  mc <- match.call()
  engine <- interactmeta_engine_from_approach(
    engine,
    method,
    approach = approach,
    allowed = c(interactmeta_public_engines(), interactmeta_internal_engines())
  )
  if (is.null(method) || (is.character(method) && length(method) > 1L)) {
    method <- interactmeta_default_method(engine, approach)
  }
  method <- interactmeta_resolve_method(engine, method)
  dots <- list(...)
  if ("ecological" %in% names(dots)) {
    warning(
      "`ecological` is deprecated; use `ecological_slope` instead.",
      call. = FALSE
    )
    ecological_slope <- dots$ecological
    dots$ecological <- NULL
  }
  ecological_slope <- isTRUE(ecological_slope)
  ecological <- ecological_slope
  if (identical(engine, "bayesmeta")) {
    interactmeta_reject_nonpositive_sampling_vi(x)
  }
  modifier_name <- interactmeta_resolve_modifier_column(
    modifier, subgroup, subgroup_var
  )
  data <- prepare_ifrac_input(
    x,
    reference = reference,
    trial = if (is.null(trial)) trial_var else trial,
    subgroup = modifier_name,
    contrast = contrast,
    contrast_sep = contrast_sep,
    big_var = big_var,
    include_overall = include_overall,
    include_interaction = include_interaction
  )
  # ifrac() writes the level column as `subgroup`. Keep the caller's
  # column name for display headings only.
  if ("subgroup" %in% names(data)) {
    subgroup_var <- "subgroup"
  } else if (modifier_name %in% names(data)) {
    subgroup_var <- modifier_name
  }
  subgroup <- subgroup_var

  fit_data <- data[
    is_ifrac_component_row(data[[subgroup_var]]),
    ,
    drop = FALSE
  ]

  design <- resolve_comparison_design(
    fit_data,
    reference = reference,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  if (!ecological_slope) {
    design$formula <- comparison_model_formula(
      design$n_nonref,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      ecological = FALSE
    )
  }

  data_fit <- standardize_comparison_fit_data(fit_data, design)

  if ("seed" %in% names(dots)) {
    warning(
      "Both `seed` and `seed` in `...`; using the `seed` argument.",
      call. = FALSE
    )
    dots$seed <- NULL
  }
  if ("silent" %in% names(dots)) {
    silent <- dots$silent
    dots$silent <- NULL
  }
  if ("engine" %in% names(dots)) {
    warning(
      "`engine` in `...` is ignored; pass `engine` as a named argument to `cam()`.",
      call. = FALSE
    )
    dots$engine <- NULL
  }
  if ("method" %in% names(dots)) {
    warning(
      "`method` in `...` is ignored; pass `method` as a named argument to `cam()`.",
      call. = FALSE
    )
    dots$method <- NULL
  }
  if ("ecological_slope" %in% names(dots) || "ecological" %in% names(dots)) {
    warning(
      "`ecological_slope` in `...` is ignored; pass it as a named argument to `cam()`.",
      call. = FALSE
    )
    dots$ecological_slope <- NULL
    dots$ecological <- NULL
  }
  retired <- intersect(
    names(dots),
    c("formula", "correlated_re", "prior", "use_precompiled")
  )
  if (length(retired)) {
    stop(
      "Argument(s) no longer supported: ",
      paste(retired, collapse = ", "),
      ". Use `tau.prior` / `tau.interaction.prior` / `fixed.prior` for Bayesian priors ",
      "and `engine_args` for computational settings.",
      call. = FALSE
    )
  }

  leftover_prior <- intersect(
    names(dots),
    c(
      "tau.prior.overall", "tau.prior.interaction",
      "mu.prior", "beta.prior.mean", "beta.prior.sd", "beta.prior.cov"
    )
  )
  if (length(leftover_prior)) {
    stop(
      "Argument(s) ",
      paste(leftover_prior, collapse = ", "),
      " are not computational `...` / `engine_args`. ",
      "Use `tau.prior`, `tau.interaction.prior`, and `fixed.prior`.",
      call. = FALSE
    )
  }

  user_tau_prior <- !is.null(tau.prior)
  user_tau_interaction_prior <- !is.null(tau.interaction.prior)
  user_fixed_prior <- !is.null(fixed.prior)
  priors <- cam_priors_resolve(
    tau.prior = tau.prior,
    tau.interaction.prior = tau.interaction.prior,
    fixed.prior = fixed.prior,
    user_tau = user_tau_prior
  )
  cam_priors_check_engine(priors, engine)
  engine_args <- interactmeta_merge_engine_args(engine_args, dots)
  if (identical(engine, "metafor")) {
    .interactmeta_warn_unused_priors(
      user_tau_prior || user_tau_interaction_prior || user_fixed_prior
    )
  }

  if (identical(engine, "bayesmeta")) {
    if (!identical(method, "bayes")) {
      stop(
        "`engine = \"bayesmeta\"` only supports `method = \"bayes\"`.",
        call. = FALSE
      )
    }
    fit <- .interactmeta_fit_bayesmeta(
      data_fit = data_fit,
      design = design,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      ecological = ecological,
      seed = seed,
      silent = silent,
      priors = priors,
      engine_args = engine_args
    )
  } else if (identical(engine, "brms")) {
    if (!identical(method, "bayes")) {
      stop("`engine = \"brms\"` only supports `method = \"bayes\"`.", call. = FALSE)
    }
    fit <- .interactmeta_fit_brms(
      data_fit = data_fit,
      design = design,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      ecological = ecological,
      seed = seed,
      silent = silent,
      priors = priors,
      engine_args = engine_args
    )
  } else if (identical(engine, "metafor")) {
    fit <- .interactmeta_fit_metafor(
      data_fit = data_fit,
      design = design,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      ecological = ecological,
      seed = seed,
      silent = silent,
      method = method,
      engine_args = engine_args
    )
  } else {
    stop(
      "`engine = \"",
      engine,
      "\"` is not a CAMS backend. ",
      "Use `engine = \"bayesmeta\"` or `engine = \"metafor\"`.",
      call. = FALSE
    )
  }

  fit$priors <- priors

  fit <- comparison_restore_interactmeta_coef_names(fit, data, design, contrast_sep = contrast_sep)

  fit <- interactmeta_attach_audit_metadata(
    fit = fit,
    call = mc,
    design = design,
    data = data,
    data_fit = data_fit,
    engine = engine,
    method = method,
    priors = priors,
    engine_args = engine_args,
    seed = seed,
    ecological_slope = ecological_slope,
    big_var = big_var,
    contrast = contrast,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )

  attr(fit, "interactmeta") <- list(
    reference = design$reference,
    ifrac_var = design$ifrac_var,
    ifrac_vars = design$ifrac_vars,
    ifrac2_var = design$ifrac2_var,
    ifrac_vars_std = design$ifrac_vars,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    modifier = modifier_name,
    effects = design$effects,
    coef_map = design$coef_map,
    correlated_re = FALSE,
    data_fit = data_fit,
    ifrac_data = data,
    engine = engine,
    method = method,
    algorithm = fit$algorithm,
    ecological_slope = ecological_slope,
    ecological = ecological_slope,
    model = if (ecological_slope) "comparison" else "specific"
  )

  fit
}

#' @rdname cam
#' @export
cams <- cam

#' @rdname cam
#' @export
subgroupmeta <- cam

resolve_comparison_design <- function(data,
                                reference,
                                trial_var,
                                subgroup_var) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  need <- c("yi", "vi", subgroup_var, trial_var)
  missing <- setdiff(need, names(data))
  if (length(missing)) {
    stop(
      "`data` is missing required column(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  ref <- reference
  if (is.null(ref)) {
    if (!"reference" %in% names(data)) {
      stop(
        "Supply `reference` or include a `reference` column from `ifrac()`.",
        call. = FALSE
      )
    }
    ref <- unique(data$reference)
    ref <- ref[!is.na(ref) & nzchar(ref)]
    if (length(ref) != 1L) {
      stop("`reference` must be a single subgroup label.", call. = FALSE)
    }
    ref <- ref[[1]]
  }

  ifrac_cols <- grep("^ifrac_", names(data), value = TRUE)
  ifrac_cols <- ifrac_cols[!grepl("^ifrac2_", ifrac_cols)]
  ifrac_cols <- sort(ifrac_cols)

  if (length(ifrac_cols) == 0L) {
    stop(
      "No `ifrac_<subgroup>` columns found after running `ifrac()`.",
      call. = FALSE
    )
  }

  n_nonref <- length(ifrac_cols)
  non_ref_labels <- sub("^ifrac_", "", ifrac_cols)
  ifrac2_vars <- paste0("ifrac2_", non_ref_labels)
  missing_ifrac2 <- setdiff(ifrac2_vars, names(data))
  if (length(missing_ifrac2)) {
    data <- add_ifrac2_design(data, subgroup_col = subgroup_var)
  }

  if (any(non_ref_labels == ref)) {
    stop(
      "Non-reference `ifrac_*` columns must differ from `reference` (",
      ref,
      ").",
      call. = FALSE
    )
  }

  subgroup_levels <- sort(unique(as.character(data[[subgroup_var]])))
  if (!(ref %in% subgroup_levels)) {
    stop(
      sprintf("`reference` (%s) must appear in `%s`.", ref, subgroup_var),
      call. = FALSE
    )
  }

  spec <- comparison_model_spec(n_nonref, trial_var, subgroup_var)
  others <- setdiff(subgroup_levels, ref)
  sg_fit <- comparison_subgroup_fit_levels(n_nonref + 1L)
  coef_map <- list(
    ifrac_orig = ifrac_cols,
    ifrac_std = comparison_std_ifrac_names(n_nonref),
    ifrac2_orig = ifrac2_vars,
    ifrac2_std = comparison_std_ifrac2_names(n_nonref),
    subgroup_fit = sg_fit,
    subgroup_orig = c(ref, others)
  )

  list(
    reference = ref,
    ifrac_var = ifrac_cols[[1L]],
    ifrac_vars = ifrac_cols,
    ifrac2_var = ifrac2_vars[[1L]],
    ifrac2_vars = ifrac2_vars,
    non_ref_labels = non_ref_labels,
    subgroup_levels = subgroup_levels,
    subgroup_var = subgroup_var,
    trial_var = trial_var,
    n_nonref = n_nonref,
    n_subgroups = n_nonref + 1L,
    effects = paste(
      paste0(ifrac_cols, ":", subgroup_var),
      collapse = " + "
    ),
    formula = spec$formula,
    prior = spec$prior,
    correlated_re = FALSE,
    coef_map = coef_map
  )
}

#' Documented CAMS priors (brms / mixed-model representation)
#' @noRd
default_comparison_prior <- function(ifrac2_coef, group = "trial") {
  group <- as.character(group)[[1]]
  ifrac2_coef <- as.character(ifrac2_coef)
  rows <- list(
    data.frame(
      prior = "normal(0, 0.5)",
      class = "sd",
      coef = "Intercept",
      group = group,
      stringsAsFactors = FALSE
    )
  )
  for (coef in ifrac2_coef) {
    rows[[length(rows) + 1L]] <- data.frame(
      prior = "normal(0, 0.1)",
      class = "sd",
      coef = coef,
      group = group,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

interactmeta_settings <- function(object) {
  st <- attr(object, "interactmeta", exact = TRUE)
  if (!is.null(st)) {
    return(st)
  }

  if (!is_comparison_fit(object)) {
    stop("`object` must be an `interactmeta` from `cam()`.", call. = FALSE)
  }

  dat <- object$data
  ifrac_var <- detect_ifrac_var(dat)
  ifrac2_var <- paste0("ifrac2_", sub("^ifrac_", "", ifrac_var))
  if (!(ifrac2_var %in% names(dat))) {
    ifrac2_var <- NA_character_
  }

  trial_var <- if ("trial" %in% names(dat)) {
    "trial"
  } else if ("study" %in% names(dat)) {
    "study"
  } else {
    stop("Model data must contain `trial` or `study`.", call. = FALSE)
  }

  subgroup_var <- if ("subgroup" %in% names(dat)) "subgroup" else stop("Model data must contain `subgroup`.")

  list(
    ifrac_var = ifrac_var,
    ifrac2_var = ifrac2_var,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    effects = paste0(ifrac_var, ":", subgroup_var)
  )
}

detect_ifrac_var <- function(dat) {
  cols <- grep("^ifrac_", names(dat), value = TRUE)
  cols <- cols[!grepl("^ifrac2_", cols)]
  if (length(cols) == 1L) {
    return(cols[[1]])
  }
  if ("ifrac" %in% names(dat)) {
    return("ifrac")
  }
  if (length(cols) > 1L) {
    stop(
      "Multiple `ifrac_*` columns found; fit with `cam()` or pass `ifrac_var` to `plot_contribution()`.",
      call. = FALSE
    )
  }
  stop("No `ifrac_*` column found in model data.", call. = FALSE)
}

#' Test whether an object is a CAMS model fit
#' @noRd
is_comparison_fit <- function(x) {
  inherits(x, "interactmeta")
}

#' Whether a CAMS fit used Stan optimize (`method = "mle"`)
#' @noRd
comparison_is_mle <- function(x) {
  if (!is_comparison_fit(x)) {
    return(FALSE)
  }
  st <- attr(x, "interactmeta", exact = TRUE)
  if (!is.null(st$method) && identical(as.character(st$method)[[1L]], "mle")) {
    return(TRUE)
  }
  identical(x$method, "mle") ||
    identical(x$algorithm, "optimize") ||
    (identical(x$engine, "cmdstanr") && identical(x$algorithm, "optimization"))
}

