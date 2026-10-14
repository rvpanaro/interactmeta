#' Contribution-adjusted meta-analysis of treatment-by-subgroup interactions
#'
#' Fits the CAMS model (Panaro, Röver & Friede 2025) for **consistent
#' inference on treatment-by-subgroup interactions** from subgroup-level
#' aggregate data. Subgroup-specific treatment effects are inputs and
#' derived outputs; the methodological target is not a separate meta-analysis
#' within each subgroup.
#'
#' [cams()] is the fitter. [cam()] is another name for the same function.
#'
#' By default, [ifrac()] is run inside the fitter on the input (e.g. a
#' [`metafor::escalc`] object). You can also pass an existing [ifrac()] tibble.
#'
#' The model for binary subgroups is the orthogonal CAMS factorization
#' (Panaro, Röver & Friede 2025, §3.4): an overall information-fraction
#' meta-regression of \eqn{m_j} plus an interaction random-effects meta-analysis
#' of \eqn{g_j}. Bayesian CAMS is implemented using DIRECT quadrature via
#' **bayesmeta**. Frequentist CAMS is implemented using **metafor**. Markov
#' chain Monte Carlo (MCMC) is not required for CAMS inference.
#'
#' @section Sampling assumption (within-trial covariance):
#' Trial blocks use
#' \eqn{g_j=y_{Bj}-y_{Aj}} with \eqn{\mathrm{Var}(g_j)=v_{Aj}+v_{Bj}} and
#' information fraction
#' \eqn{\pi_j=(1/v_{Bj})/(1/v_{Aj}+1/v_{Bj})} (IV weight of the non-reference
#' subgroup). That choice makes \eqn{\mathrm{Cov}(m_j,g_j)=0} when
#' \eqn{\mathrm{Cov}(y_{Aj},y_{Bj})=0}. This package currently assumes that
#' within-trial subgroup sampling covariance is zero. Separate overall and
#' interaction fits (`bmr`/`bayesmeta` or two `rma`s) rely on that
#' orthogonality. If subgroup estimates share a control arm or other
#' within-trial dependence, both \eqn{\mathrm{Var}(g_j)} and block
#' independence are misspecified unless \eqn{\pi_j} and the variances are
#' redefined with the covariance (not implemented here).
#'
#' @section Model variants:
#' One call; switch with `ecological_slope`:
#' * `ecological_slope = TRUE` (default) — **comparison** model with free
#'   ecological / composition slope \eqn{\beta}
#' * `ecological_slope = FALSE` — **subgroup-specific** model with
#'   \eqn{\beta \equiv 0} (hence \eqn{\theta = \gamma}). Public
#'   `bayesmeta` / `metafor` engines implement this as a **constrained**
#'   joint model for \eqn{(\alpha,\gamma)} with free across-trial
#'   heterogeneity \eqn{\tau_m} and \eqn{\tau_g}. Switching the composition
#'   slope off does **not** turn off across-trial patterns: trial-level
#'   random effects and both heterogeneity SDs remain part of the inference.
#'   Population lines become \eqn{\mu_A=\alpha}, \eqn{\mu_B=\alpha+\gamma},
#'   \eqn{M(\pi)=\alpha+\pi\gamma}.
#'
#' @section Inference:
#' CAMS can be fitted using Bayesian or frequentist inference. These are not
#' two different scientific models. The approaches are **bayesmeta** (DIRECT,
#' the default) and **metafor** (REML or ML). See `vignette("introduction-to-interact")`.
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
#' [bayescams()] for Bayesian CAMS and `cams(..., engine = "metafor")` for
#' frequentist CAMS. The fit stores `fit$inference`, `fit$engine`, `fit$method`,
#' and `fit$algorithm` (`"quadrature"` or `"optimization"`).
#'
#' @section Priors:
#' Bayesian CAMS passes priors through to bayesmeta. They are stored on
#' `fit$priors`. **Heterogeneity defaults match [bayesmeta::bayesmeta()]**:
#' \describe{
#'   \item{fixed effects}{independent \eqn{\mathrm{Normal}(0, 10)} on the
#'     overall-block means \eqn{(\alpha,\theta)} and on the interaction mean
#'     \eqn{\gamma}, with \eqn{\theta=\beta+\gamma}. Stored `b_ifrac` is the
#'     reconstructed ecological slope \eqn{\beta=\theta-\gamma}; that prior is
#'     **not** an independent \eqn{N(0,10)} on \eqn{\beta}. (bayesmeta's
#'     improper flat \eqn{\mu} is not the CAMS default: the overall block is
#'     fit with `bmr`, where a proper normal is the reliable companion to
#'     improper tau)}
#'   \item{`tau.prior`}{improper uniform \eqn{p(\tau_m)\propto 1}
#'     (\eqn{\tau_m \ge 0}), bayesmeta's `tau.prior = "uniform"`}
#'   \item{`tau.interaction.prior`}{the same improper uniform on
#'     \eqn{\tau_g} when both heterogeneity priors are omitted. If
#'     `tau.interaction.prior` is `NULL` and the user supplied
#'     `tau.prior`, the interaction prior **inherits** that density}
#' }
#' Character shorthands `"uniform"` and `"sqrt"` are the **improper** bayesmeta
#' densities \eqn{p(\tau)\propto 1} and \eqn{p(\tau)\propto \tau^{-1/2}}
#' (\eqn{\tau \ge 0}). They are **not** \eqn{U(0, U)}.
#' `engine = "metafor"` ignores Bayesian priors (warning if the user supplied
#' them).
#'
#' **Empirical heterogeneity (not the package default).** CDSR calibration
#' suggests interaction heterogeneity is often substantially smaller than
#' overall treatment-effect heterogeneity, so a common weakly informative
#' half-normal on both blocks can be too diffuse for \eqn{\tau_g}
#' (Panaro, Röver & Friede, arXiv:2606.23968). That paper reports
#' measure-specific empirical scales (e.g. \eqn{\mathrm{HN}(0.5)}-style WIP
#' for overall log-ratio / SMD-type endpoints in Röver et al., 2021;
#' tighter interaction scales such as \eqn{\mathrm{HN}(0.1)} in sensitivity
#' analyses and \eqn{\mathrm{HN}(0.07)} empirically for log-ROR). Vignettes
#' set those half-normal densities **explicitly**; they are not silently
#' applied by [cams()].
#'
#' **Fixed-effect parameterization:** bayesmeta/bmr receive independent
#' \eqn{N(0, 10)} on \eqn{(\alpha,\theta)} (`beta.prior.cov = diag(sd^2)`) and
#' `mu.prior` on \eqn{\gamma}. That prior factors across the two orthogonal
#' blocks. Independent normals on \eqn{(\alpha,\beta,\gamma)} would couple
#' \eqn{\theta} and \eqn{\gamma} a priori and are not what this engine uses.
#'
#' @section Intervals:
#' For Bayesian CAMS, public extractors ([subgroups()], [interactions()],
#' [combined()], [predict_cams_by_if()], [shrinkage()], [plot_contribution()])
#' report shortest credible intervals from the DIRECT-node Gaussian mixture
#' (`.weight` and within-node `.cond_*` moments), matching the default
#' bayesmeta interval rule. Alpha and theta are never convolved with each
#' other.
#'
#' @param x A [`metafor::escalc`] object, a data frame with `yi` / `vi`, or output
#'   from [ifrac()] (must include `ifrac_*` design columns). Also accepts a
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
#' @param contrast,contrast_sep,big_var,include_combined,include_interaction Passed
#'   to [ifrac()] when information fractions are computed inside `cam()`.
#' @param method Inferential method: `"bayes"` (default for Bayesian CAMS),
#'   `"REML"` or `"ML"` (frequentist CAMS). Not a backend name. If `engine` is
#'   omitted, `"bayes"` uses `engine = "bayesmeta"` and `"REML"` / `"ML"` use
#'   `engine = "metafor"`. An error is thrown only when `engine` and `method`
#'   are both specified and incompatible.
#' @param ecological_slope Logical; if `TRUE` (default), include the ecological /
#'   composition slope \eqn{\beta} (`b_ifrac`) — the **comparison** model. If
#'   `FALSE`, constrain \eqn{\beta \equiv 0} (\eqn{\theta=\gamma}) — the
#'   **subgroup-specific** model. Across-trial heterogeneity
#'   (\eqn{\tau_m}, \eqn{\tau_g}) and trial-level reconstruction remain;
#'   only the mean composition slope is switched off.
#' @param approach Inferential approach. `"bayes"` selects Bayesian CAMS
#'   (`engine = "bayesmeta"` when `engine` is omitted). `"frequentist"`,
#'   `"REML"`, or `"ML"` select frequentist CAMS (`engine = "metafor"`).
#'   `NULL` (default) follows `engine` / `method`.
#' @param engine Approach used to fit CAMS: `"bayesmeta"` (Bayesian, default)
#'   or `"metafor"` (frequentist). Changing `engine` should not silently change
#'   the Bayesian model. `"bayes"` is
#'   accepted as a synonym of `"bayesmeta"`. Omit `engine` when `method` or
#'   `approach` already identifies the fit.
#' @param tau.prior,tau.interaction.prior Heterogeneity priors for \eqn{\tau_m}
#'   (overall) and \eqn{\tau_g} (interaction), passed to bayesmeta.
#'   Defaults when both are omitted match bayesmeta: improper
#'   `"uniform"` on both blocks (see **Priors**). `"uniform"` and `"sqrt"`
#'   are the improper bayesmeta densities. A half-normal or half-Cauchy is a
#'   density, for example
#'   `function(tau) bayesmeta::dhalfnormal(tau, scale = 0.5)`.
#'   If `tau.interaction.prior = NULL` and the user supplied `tau.prior`,
#'   the interaction prior inherits `tau.prior`.
#' @param fixed.prior Normal prior for the bayesmeta fixed effects, as
#'   `c(mean = 0, sd = 10)`. `NULL` uses that default (needed for a stable
#'   `bmr` overall block with improper tau; not bayesmeta's improper flat
#'   \eqn{\mu}). Applied as independent normals on \eqn{(\alpha,\theta)} and
#'   on \eqn{\gamma}, not as an independent normal on the reconstructed
#'   \eqn{\beta=\theta-\gamma}. `"uniform"` is the improper flat prior.
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
#' library(interact)
#'
#' es <- escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
#' fit <- bayescams(es, reference = "Younger")
#' fit_reml <- cams(es, reference = "Younger", engine = "metafor")
#' forestplot(fit)
#' plot_contribution(fit)
#' plot(fit)
#' }
#'
#' @references
#' Panaro RV, Röver C, Friede T (2025).
#' Consistent Bayesian meta-analysis on subgroup specific effects and
#' interactions. arXiv:2512.18785. \url{https://doi.org/10.48550/arXiv.2512.18785}
#'
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \url{https://doi.org/10.48550/arXiv.2508.15531}
#'
#' Panaro RV, Röver C, Friede T (2026).
#' Empirical prior distributions for treatment-by-subgroup interaction
#' heterogeneity in random-effects meta-analysis.
#' arXiv:2606.23968. \url{https://doi.org/10.48550/arXiv.2606.23968}
#'
#' Röver C, Bender R, Dias S, Schmid CH, Schmidli H, Sturtz S, Weber S,
#' Friede T (2021). On weakly informative prior distributions for the
#' heterogeneity parameter in Bayesian random-effects meta-analysis.
#' \emph{Research Synthesis Methods}, 12(4), 448--474.
#' \url{https://doi.org/10.1002/jrsm.1475}
#' @seealso [bayescams()], [conventional()], [da()], [ad()], [swada()], [forestplot()]
#' @name cams
#' @export
cams <- function(x,
                 reference = NULL,
                 modifier = NULL,
                 trial = NULL,
                 subgroup = NULL,
                 trial_var = "trial",
                 subgroup_var = "subgroup",
                 contrast = "all",
                 contrast_sep = "-",
                 big_var = 100^2,
                 include_combined = FALSE,
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
    include_combined = include_combined,
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

  fit <- interactmeta_finalize_fit(
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

#' @rdname cams
#' @export
cam <- cams

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

#' Proper prior table for brms / mixed-model scaffolding (not bayesmeta defaults)
#'
#' Package Bayesian defaults are bayesmeta's improper uniforms, which cannot
#' be represented in Stan. This helper keeps the article WIP half-normal
#' scales (arXiv:2606.23968) for code paths that need a proper prior table.
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

