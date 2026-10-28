# ---- comparison.R ----
#' Contribution-adjusted meta-analysis for treatment effect modification
#'
#' Fits the CAMS model (Panaro, Röver & Friede 2025) for **consistent
#' inference on treatment effect modification** from modifier-level
#' aggregate data (vignettes often use classic **subgroup** tables).
#' Modifier-level treatment effects are inputs and derived outputs; the
#' methodological target is not a separate meta-analysis within each modifier level.
#'
#' [cams()] is the fitter.
#'
#' By default, [ifrac()] is run inside the fitter on the input (e.g. a
#' [`metafor::escalc`] object). You can also pass an existing [ifrac()] tibble.
#'
#' For two subgroups per trial, CAMS fits an information-fraction
#' meta-regression of combined coordinates \eqn{m_j} and a random-effects
#' meta-analysis of within-trial contrasts \eqn{g_j} (Panaro, Röver &
#' Friede 2025, §3.4). Default inference is **Bayesian** (**bayesmeta**);
#' **`method = "REML"`** (or other **`metafor::rma()`** methods) gives
#' frequentist inference.
#'
#' @section Sampling assumption (within-trial covariance):
#' Within each trial,
#' \eqn{g_j=y_{Bj}-y_{Aj}} with \eqn{\mathrm{Var}(g_j)=v_{Aj}+v_{Bj}} and
#' information fraction
#' \eqn{\pi_j=(1/v_{Bj})/(1/v_{Aj}+1/v_{Bj})} (IV weight of the non-reference
#' subgroup). That choice makes \eqn{\mathrm{Cov}(m_j,g_j)=0} when
#' \eqn{\mathrm{Cov}(y_{Aj},y_{Bj})=0}. This package currently assumes that
#' within-trial subgroup sampling covariance is zero. Separate overall and
#' interaction fits (`bmr`/`bayesmeta` or two `rma`s) rely on that
#' orthogonality. If subgroup estimates share a control arm or other
#' within-trial dependence, both \eqn{\mathrm{Var}(g_j)} and separate
#' fits for \eqn{m_j} and \eqn{g_j} are misspecified unless \eqn{\pi_j} and the variances are
#' redefined with the covariance (not implemented here).
#'
#' @section Model:
#' [cams()] fits the **comparison** CAMS model with free ecological /
#' composition slope \eqn{\beta} (`b_ifrac`): population lines
#' \eqn{\mu_A=\alpha}, \eqn{\mu_B=\alpha+\gamma+\beta}, \eqn{M(\pi)=\alpha+\pi\beta+\gamma}.
#'
#' @section Inference:
#' CAMS can be fitted using Bayesian or frequentist inference. These are not
#' two different scientific models. Set **`method = "Bayesian"`** (default) for
#' **bayesmeta**, or any **`metafor::rma()`** `method` for the two
#' frequentist components (default frequentist choice is **`"REML"`**). See
#' `vignette("cams")`.
#'
#' For Bayesian fits, bayesmeta options may go in `fit_args` or `...`.
#' For metafor fits, further **`rma()`** arguments (e.g. `test`, `control`) may
#' be passed in `...` or `fit_args` (`fit_args` wins on conflicts).
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
#' **`method`** selects inference: `"Bayesian"` uses bayesmeta; values documented
#' for [`metafor::rma()`] (`"REML"`, `"ML"`, `"DL"`, `"PM"`, …) use metafor on
#' both components.
#' The fit stores `fit$method` (e.g. `"Bayesian"`, `"REML"`), `fit$inference`,
#' and `fit$support` (DIRECT quadrature support or one-row metafor estimate).
#' [print.interactmeta()] shows the same wording as
#' the `cams(..., method = )` argument (`Bayesian`, `REML`, …).
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
#' Frequentist `method` values ignore Bayesian priors (warning if supplied).
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
#' For Bayesian CAMS, public extractors ([modifiers()], [interactions()],
#' [combined()], [predict_cams_by_if()], [shrinkage()]) report credible
#' intervals from the DIRECT-node Gaussian mixture (`.weight` and within-node
#' `.cond_*` moments). Default `interval.type` is `"central"` (equal-tailed);
#' pass `interval.type = "shortest"` to [cams()] for the narrowest interval rule.
#' Default [combined()] and [interactions()] at `type = "Overall"` with no
#' `exponentiate` use stored `bmr` / `bayesmeta` intervals when possible;
#' other requests still use the mixture grid.
#' [plot_contribution()] uses central (equal-tailed) mixture intervals for
#' ribbon limits by default (`ribbon_interval.type`). Alpha and theta are never
#' convolved with each other.
#'
#' @param x A [`metafor::escalc`] object, a data frame with `yi` / `vi`, or output
#'   from [ifrac()] (must include `ifrac_*` design columns). Also accepts a
#'   data frame built from trial-level columns (`trial`, `subgroup`, `yi`, `vi`).
#' @param reference Subgroup label used as reference in [ifrac()] and in the
#'   model. Required when `x` is not already an [ifrac()] result.
#' @param modifier Column that holds the two levels of the within-trial
#'   contrast (default `"subgroup"`, so the shipped examples keep working).
#'   Overrides `subgroup` and `subgroup_var` when set. Level labels stay as
#'   stored. Display headings capitalize the first letter of this name
#'   (`dose` is shown as `Dose`).
#' @param trial,subgroup When `x` is `escalc` / a data frame, column names passed
#'   to [ifrac()] (defaults `"trial"` and `"subgroup"`). `subgroup` remains
#'   accepted for the default column name; prefer `modifier` for the argument.
#' @param trial_var,subgroup_var Column names used in the model (defaults
#'   match `trial` and `subgroup`). `subgroup_var` remains accepted.
#' @param contrast,contrast_sep,big_var Passed to [ifrac()] when information
#'   fractions are computed inside `cams()`. To append trial-level combined or
#'   interaction contrast rows to an `ifrac()` table, call [ifrac()] with
#'   `include_combined` / `include_interaction` before fitting.
#' @param method Inferential method. `"Bayesian"` (default) fits with bayesmeta.
#'   Any value accepted by [`metafor::rma()`] for its `method` argument fits
#'   both orthogonal frequentist blocks (e.g. `"REML"`, `"ML"`, `"DL"`,
#'   `"PM"`, `"EB"`, `"SJ"`, `"HE"` / `"CO"`, `"HS"`, `"HSk"`, `"EE"`,
#'   `"GENQ"`, `"PMM"`, `"GENQM"`). Aliases `"bayes"` and `"bayesian"` match
#'   `"Bayesian"`. `"frequentist"` selects `"REML"`.
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
#' @param fit_args List passed to the backend (bayesmeta quadrature or
#'   `metafor::rma()`). Model choices such as priors do not belong here.
#'   For Bayesian fits, optional `n_support` (integer) is forwarded to both
#'   `bayesmeta::bmr()` and `bayesmeta::bayesmeta()` as `n.support` (fewer
#'   `\tau` nodes → smaller DIRECT product grid).
#' @param precompute If `TRUE`, build the default parameter matrix and warm
#'   the effect cache before returning (legacy eager behaviour). Default `FALSE`:
#'   [print.interactmeta()], [summary.interactmeta()], and [tidy()] trigger
#'   reporting work on first use.
#' @param lazy_latent For Bayesian fits, skip per-trial `convolve()` marginals
#'   at fit time (default `TRUE`). Trial-level `$latent` objects are built when
#'   shrinkage or other code needs them. Set `FALSE` for the full latent at fit.
#' @param silent Sampler verbosity; default `2` suppresses messages.
#' @param seed Random seed for reproducible draws (default from option
#'   `"interactmeta.seed"`, currently `42L`). Set to `NA` for a non-reproducible draw.
#'   If `seed` is also passed in `...`, a warning is given and this argument wins.
#' @param ... For **`method = "Bayesian"`**, bayesmeta quadrature options (e.g.
#'   `delta`, `interval.type`). For metafor methods, additional `rma()` arguments
#'   (e.g. `test`, `control`). Merged into `fit_args` when not already set.
#' @return A CAMS fit with class `"interactmeta"` plus attribute
#'   `"interactmeta"` (plot settings and `ifrac_data` for [forestplot()]).
#'
#' @examples
#' \dontrun{
#' library(interact)
#'
#' es <- escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
#' fit <- cams(es, reference = "Lower", method = "Bayesian")
#' fit_reml <- cams(es, reference = "Lower", method = "REML")
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
#' @seealso [conventional()], [da()], [ad()], [swada()], [forestplot()]
#' @name cams
#' @export
cams <- function(x,
                 reference = NULL,
                 modifier = NULL,
                 trial = NULL,
                 subgroup = NULL,
                 trial_var = "study",
                 subgroup_var = "level",
                 contrast = "all",
                 contrast_sep = "-",
                 big_var = 100^2,
                 method = "Bayesian",
                 tau.prior = NULL,
                 tau.interaction.prior = NULL,
                 fixed.prior = NULL,
                 fit_args = list(),
                 precompute = FALSE,
                 lazy_latent = TRUE,
                 silent = 2,
                 seed = getOption("interactmeta.seed", 42L),
                 ...) {
  mc <- match.call()
  resolved <- interactmeta_resolve_cams_method(method)
  backend <- resolved$backend
  method_display <- resolved$method
  method_key <- resolved$method_key
  dots <- list(...)
  if ("engine_args" %in% names(dots)) {
    stop(
      "`engine_args` has been renamed to `fit_args` in `cams()`.",
      call. = FALSE
    )
  }
  if ("engine" %in% names(dots)) {
    stop(
      "`engine` is no longer an argument to `cams()`; use `method = \"Bayesian\"` ",
      "or a `metafor::rma()` method such as `\"REML\"`.",
      call. = FALSE
    )
  }
  if ("approach" %in% names(dots)) {
    stop(
      "`approach` is no longer an argument to `cams()`; use `method` instead.",
      call. = FALSE
    )
  }
  if ("ecological" %in% names(dots) || "ecological_slope" %in% names(dots)) {
    stop(
      "`ecological_slope` is no longer an argument to `cams()`; ",
      "the comparison model (free composition slope) is always used.",
      call. = FALSE
    )
  }
  for (stale in c("include_combined", "include_interaction")) {
    if (stale %in% names(dots)) {
      warning(
        "`", stale, "` is not used by `cams()`; call `ifrac(..., ", stale,
        " = TRUE)` when you need those rows in the design table.",
        call. = FALSE
      )
      dots[[stale]] <- NULL
    }
  }
  if (identical(backend, "bayesmeta")) {
    interactmeta_reject_nonpositive_sampling_vi(x)
  }
  modifier_name <- interactmeta_resolve_modifier_column(
    modifier,
    subgroup,
    fallback = subgroup_var,
    dat = as.data.frame(x)
  )
  trial_col <- if (is.null(trial)) trial_var else trial
  data <- prepare_ifrac_input(
    x,
    reference = reference,
    trial = trial_col,
    subgroup = modifier_name,
    contrast = contrast,
    contrast_sep = contrast_sep,
    big_var = big_var
  )
  subgroup_var <- if (modifier_name %in% names(data)) {
    modifier_name
  } else {
    interactmeta_resolve_level_col(data, fallback = modifier_name)
  }
  trial_var <- interactmeta_resolve_study_col(data, fallback = trial_col)
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
  if ("method" %in% names(dots)) {
    stop(
      "`method` must be passed as a named argument to `cams()`, not in `...`.",
      call. = FALSE
    )
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
      "and `fit_args` for computational settings.",
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
      " are not computational `...` / `fit_args`. ",
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
  cam_priors_check_engine(priors, backend)
  rma_exclude <- c(
    "seed", "silent", "precompute", "lazy_latent",
    "include_combined", "include_interaction", "fit_args"
  )
  if (identical(backend, "metafor")) {
    .interactmeta_warn_unused_priors(
      user_tau_prior || user_tau_interaction_prior || user_fixed_prior
    )
    fit_args <- interactmeta_merge_rma_args(fit_args, dots, exclude = rma_exclude)
  } else {
    fit_args <- interactmeta_merge_engine_args(fit_args, dots)
  }

  if (identical(backend, "bayesmeta")) {
    fit <- .interactmeta_fit_bayesmeta(
      data_fit = data_fit,
      design = design,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      seed = seed,
      silent = silent,
      priors = priors,
      fit_args = fit_args,
      lazy_latent = isTRUE(lazy_latent)
    )
  } else if (identical(backend, "metafor")) {
    fit <- .interactmeta_fit_metafor(
      data_fit = data_fit,
      design = design,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      seed = seed,
      silent = silent,
      method = method_key,
      fit_args = fit_args
    )
  } else {
    stop("Unknown CAMS backend `", backend, "`.", call. = FALSE)
  }

  fit$priors <- priors

  fit <- comparison_restore_interactmeta_coef_names(fit, data, design, contrast_sep = contrast_sep)

  fit <- interactmeta_finalize_fit(
    fit = fit,
    call = mc,
    design = design,
    data = data,
    data_fit = data_fit,
    backend = backend,
    method_display = method_display,
    method_key = method_key,
    priors = priors,
    fit_args = fit_args,
    seed = seed,
    ecological_slope = TRUE,
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
    backend = backend,
    method = method_display,
    algorithm = fit$algorithm,
    ecological_slope = TRUE,
    ecological = TRUE,
    model = "comparison"
  )

  fit <- interactmeta_init_fit_sidecar(fit)

  if (isTRUE(precompute)) {
    fit <- interactmeta_precompute_reporting(fit)
  }

  fit
}

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

#' Proper prior table retained in the model specification
#'
#' Package Bayesian defaults are bayesmeta's improper uniforms, which cannot
#' be represented as a finite table. This helper records the article WIP
#' half-normal scales (arXiv:2606.23968) for model metadata.
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
    stop("`object` must be an `interactmeta` from `cams()`.", call. = FALSE)
  }

  dat <- object$data
  ifrac_vars <- detect_ifrac_vars(dat)
  ifrac_var <- if (length(ifrac_vars) == 1L) {
    ifrac_vars[[1L]]
  } else if ("ifrac" %in% names(dat)) {
    "ifrac"
  } else if (length(ifrac_vars) > 1L) {
    stop(
      "Multiple `ifrac_*` columns found; fit with `cams()` or pass `ifrac_var` to `plot_contribution()`.",
      call. = FALSE
    )
  } else {
    stop("No `ifrac_*` column found in model data.", call. = FALSE)
  }
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

  subgroup_var <- interactmeta_resolve_level_col(dat)
  if (!subgroup_var %in% names(dat)) {
    stop(
      "Model data must contain a modifier column (`level` or legacy `subgroup`).",
      call. = FALSE
    )
  }

  list(
    ifrac_var = ifrac_var,
    ifrac2_var = ifrac2_var,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    effects = paste0(ifrac_var, ":", subgroup_var)
  )
}

#' Test whether an object is a CAMS model fit
#' @noRd
is_comparison_fit <- function(x) {
  inherits(x, "interactmeta")
}



# ---- comparison_model.R ----
comparison_std_ifrac_names <- function(n_nonref) {
  paste0("ifrac_", seq_len(n_nonref))
}

comparison_std_ifrac2_names <- function(n_nonref) {
  paste0("ifrac2_", seq_len(n_nonref))
}

comparison_subgroup_fit_levels <- function(n_subgroups) {
  paste0("interactmeta_sg_", seq_len(n_subgroups))
}

comparison_model_formula <- function(n_nonref,
                               trial_var = "study",
                               subgroup_var = "level",
                               ecological = TRUE) {
  if (n_nonref < 1L) {
    stop("CAMS requires at least two subgroup levels (one reference).", call. = FALSE)
  }
  ifrac2_terms <- paste(comparison_std_ifrac2_names(n_nonref), collapse = " + ")
  fe_terms <- if (isTRUE(ecological)) {
    ifrac_terms <- paste(comparison_std_ifrac_names(n_nonref), collapse = " + ")
    paste0(ifrac_terms, " + ", subgroup_var)
  } else {
    as.character(subgroup_var)
  }
  paste0(
    "yi | se(sqrt(vi)) ~ 1 + ",
    fe_terms,
    " + (1 + ",
    ifrac2_terms,
    " || ",
    trial_var,
    ")"
  )
}

comparison_model_spec <- function(n_nonref,
                            trial_var = "study",
                            subgroup_var = "level") {
  ifrac2 <- comparison_std_ifrac2_names(n_nonref)
  list(
    n_nonref = n_nonref,
    n_subgroups = n_nonref + 1L,
    formula = comparison_model_formula(n_nonref, trial_var, subgroup_var),
    correlated_re = FALSE,
    ifrac_vars_std = comparison_std_ifrac_names(n_nonref),
    ifrac2_vars_std = ifrac2,
    prior = default_comparison_prior(ifrac2, group = trial_var)
  )
}

standardize_comparison_fit_data <- function(data, design) {
  out <- data
  n <- design$n_nonref
  ifrac_std <- comparison_std_ifrac_names(n)
  ifrac2_std <- comparison_std_ifrac2_names(n)

  for (i in seq_len(n)) {
    orig_i <- design$ifrac_vars[[i]]
    out[[ifrac_std[i]]] <- out[[orig_i]]
    out[[ifrac2_std[i]]] <- out[[design$ifrac2_vars[[i]]]]
  }

  sg_chr <- as.character(out[[design$subgroup_var]])
  ref <- design$reference
  other <- setdiff(design$subgroup_levels, ref)
  level_order <- c(ref, other)
  out[[design$subgroup_var]] <- factor(
    paste0("interactmeta_sg_", match(sg_chr, level_order)),
    levels = comparison_subgroup_fit_levels(design$n_subgroups)
  )
  out
}

#' Map canonical fit labels (ifrac_1, interactmeta_sg_2, ...) back to data names
#' Name map from standardized fit names to user-facing labels
#' @noRd
comparison_coef_name_replacements <- function(design) {
  cm <- design$coef_map
  reps <- stats::setNames(cm$ifrac_orig, cm$ifrac_std)
  reps <- c(reps, stats::setNames(cm$ifrac2_orig, cm$ifrac2_std))

  ref <- design$reference
  others <- setdiff(design$subgroup_levels, ref)
  sg_fit <- comparison_subgroup_fit_levels(design$n_subgroups)
  sg_orig <- stats::setNames(c(ref, others), sg_fit)
  c(reps, sg_orig)
}

#' @noRd
comparison_replace_fixed_strings <- function(x, replacements) {
  if (!length(replacements)) {
    return(x)
  }
  ord <- order(-nchar(names(replacements)))
  replacements <- replacements[ord]
  for (i in seq_along(replacements)) {
    x <- gsub(names(replacements)[i], replacements[i], x, fixed = TRUE)
  }
  x
}

#' Subgroup rows only, with reference-first factor levels
#' @noRd
comparison_fit_display_data <- function(data, design, contrast_sep = "-") {
  subgroup_var <- design$subgroup_var
  comp <- is_ifrac_component_row(data[[subgroup_var]], contrast_sep = contrast_sep)
  out <- data[comp, , drop = FALSE]
  ref <- design$reference
  others <- setdiff(design$subgroup_levels, ref)
  out[[subgroup_var]] <- factor(
    as.character(out[[subgroup_var]]),
    levels = c(ref, others)
  )
  out
}

#' Restore user-facing names on an `interactmeta` object
#' @noRd
comparison_restore_interactmeta_coef_names <- function(fit, data, design, contrast_sep = "-") {
  replacements <- comparison_coef_name_replacements(design)
  if (length(replacements) && is.data.frame(fit$support)) {
    names(fit$support) <- comparison_replace_fixed_strings(names(fit$support), replacements)
  }
  if (!is.null(fit$formula)) {
    txt <- paste(deparse(fit$formula), collapse = " ")
    fit$formula <- stats::as.formula(
      comparison_replace_fixed_strings(txt, replacements),
      env = environment(fit$formula)
    )
  }
  fit$data <- comparison_fit_display_data(data, design, contrast_sep = contrast_sep)
  fit
}

# ---- modifier.R ----
#' Column that holds the two levels of a within-trial contrast
#'
#' `modifier` is that column name. The default column is `"level"`; legacy
#' `"subgroup"` is still accepted in user data. A display heading capitalizes
#' only the first letter of the column name, so `dose` is shown as `Dose`.
#'
#' @noRd
NULL

#' Standard study identifier column in long-format data
#' @noRd
interactmeta_study_column <- function() {
  "study"
}

#' Study column in a data frame (`study`, else legacy `trial`)
#' @noRd
interactmeta_resolve_study_col <- function(dat, fallback = NULL) {
  col <- interactmeta_study_column()
  if (is.data.frame(dat) && col %in% names(dat)) {
    return(col)
  }
  if (is.data.frame(dat) && "trial" %in% names(dat)) {
    return("trial")
  }
  fallback %||% col
}

#' Rename legacy `trial` column to `study` when needed
#' @noRd
interactmeta_rename_trial_col_to_study <- function(dat) {
  col <- interactmeta_study_column()
  if ("trial" %in% names(dat) && !col %in% names(dat)) {
    names(dat)[names(dat) == "trial"] <- col
  }
  dat
}

#' brms / metafor trial RE columns always use the `trial` group prefix
#' @noRd
interactmeta_sd_re_group <- function() {
  "trial"
}

#' Standard modifier-level column name in long-format data and effect tables
#' @noRd
interactmeta_level_column <- function() {
  "level"
}

#' Modifier column present in a data frame (`level`, else legacy `subgroup`)
#' @noRd
interactmeta_resolve_level_col <- function(dat, fallback = NULL) {
  col <- interactmeta_level_column()
  if (is.data.frame(dat) && col %in% names(dat)) {
    return(col)
  }
  if (is.data.frame(dat) && "subgroup" %in% names(dat)) {
    return("subgroup")
  }
  fallback %||% col
}

#' Rename legacy `subgroup` column to `level` when needed
#' @noRd
interactmeta_rename_modifier_col_to_level <- function(dat) {
  col <- interactmeta_level_column()
  if ("subgroup" %in% names(dat) && !col %in% names(dat)) {
    names(dat)[names(dat) == "subgroup"] <- col
  }
  dat
}

#' Identifier column in modifier / shrinkage effect summaries
#' @noRd
interactmeta_modifier_summary_col <- function(summary) {
  col <- interactmeta_level_column()
  if (is.data.frame(summary) && col %in% names(summary)) {
    return(col)
  }
  if (is.data.frame(summary) && "subgroup" %in% names(summary)) {
    return("subgroup")
  }
  col
}

#' Resolve the modifier column name (`modifier`, then `subgroup`, then fallback)
#' @noRd
interactmeta_resolve_modifier_column <- function(modifier = NULL,
                                                subgroup = NULL,
                                                fallback = NULL,
                                                dat = NULL) {
  pick <- function(x) {
    if (is.null(x) || !length(x)) {
      return(NULL)
    }
    x <- as.character(x)[[1L]]
    if (!nzchar(x)) NULL else x
  }
  pick_in_data <- function(x) {
    nm <- pick(x)
    if (is.null(nm)) {
      return(NULL)
    }
    if (is.data.frame(dat) && !(nm %in% names(dat))) {
      return(NULL)
    }
    nm
  }
  fb <- pick(fallback)
  if (is.null(fb) && is.data.frame(dat)) {
    fb <- interactmeta_resolve_level_col(dat)
  }
  fb <- fb %||% interactmeta_level_column()
  pick_in_data(modifier) %||% pick_in_data(subgroup) %||% fb
}

#' Display heading for a modifier column
#'
#' Capitalizes the first letter and leaves the rest unchanged.
#' @noRd
interactmeta_modifier_heading <- function(name) {
  name <- as.character(name %||% interactmeta_level_column())[[1L]]
  if (!nzchar(name)) {
    name <- interactmeta_level_column()
  }
  paste0(toupper(substring(name, 1L, 1L)), substring(name, 2L))
}

#' @noRd
interactmeta_modifier_specific_effects_title <- function(name) {
  paste0(interactmeta_modifier_heading(name), "-specific effects")
}

# ---- trial_blocks.R ----
#' Trial-level IV overall and interaction contrast blocks
#'
#' Single source for quoted (observed) trial-level CAMS quantities.
#'
#' CAMS overall rotation (IV pool):
#' `m_j = sum_g w_g y_gj` with `w_g = 1 / v_gj`, `Var(m_j) = 1 / sum_g 1/v_gj`.
#'
#' CAMS interaction contrast:
#' `g_j = y_Bj - y_Aj`.
#' Variance uses `v_Aj + v_Bj` (within-trial subgroup covariance **assumed zero**).
#'
#' Information fraction of B:
#' `pi_j = (1/v_Bj) / (1/v_Aj + 1/v_Bj)`.
#'
#' CAMS orthogonality: when `pi_j` is that IV weight,
#' `Cov(m_j, g_j) = -(1-pi_j) v_Aj + pi_j v_Bj = 0`.
#'
#' @noRd
comparison_trial_effect_blocks <- function(data,
                                     trial_var = "study",
                                     subgroup_var = "level",
                                     reference,
                                     nonref = NULL,
                                     ifrac_var = NULL,
                                     min_vi = 1e-12) {
  dat <- as.data.frame(data)
  dat[[trial_var]] <- as.character(dat[[trial_var]])
  dat[[subgroup_var]] <- as.character(dat[[subgroup_var]])
  ref <- as.character(reference)[[1L]]
  trials <- unique(dat[[trial_var]])

  if (is.null(nonref)) {
    sg_all <- unique(dat[[subgroup_var]])
    nonrefs <- setdiff(sg_all, ref)
    if (!length(nonrefs)) {
      stop("No non-reference subgroup found for trial effect blocks.", call. = FALSE)
    }
    nonref <- nonrefs[[1L]]
  }
  nonref <- as.character(nonref)[[1L]]

  overall_rows <- vector("list", length(trials))
  interaction_rows <- vector("list", length(trials))
  pi <- stats::setNames(rep(NA_real_, length(trials)), trials)

  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    d <- dat[dat[[trial_var]] == tr, , drop = FALSE]
    if (!nrow(d)) {
      next
    }
    yi <- as.numeric(d$yi)
    vi_raw <- pmax(as.numeric(d$vi), 0)
    if (all(vi_raw <= 0)) {
      ov_yi <- mean(yi)
      ov_vi <- 0
    } else {
      # CAMS overall rotation: precision-weighted mean of subgroup estimates.
      vi <- pmax(vi_raw, min_vi)
      w <- 1 / vi
      ov_yi <- sum(w * yi) / sum(w)
      ov_vi <- 1 / sum(w)
      if (all(vi_raw == 0)) {
        ov_vi <- 0
      }
    }
    ov <- data.frame(
      trial = tr,
      yi = ov_yi,
      vi = ov_vi,
      stringsAsFactors = FALSE
    )
    if (!is.null(ifrac_var) && ifrac_var %in% names(d)) {
      pi_tr <- mean(as.numeric(d[[ifrac_var]]), na.rm = TRUE)
      ov$ifrac <- pi_tr
      pi[[tr]] <- pi_tr
    }
    overall_rows[[j]] <- ov

    i_ref <- which(d[[subgroup_var]] == ref)
    i_non <- which(d[[subgroup_var]] == nonref)
    if (length(i_ref) != 1L || length(i_non) != 1L) {
      stop(
        "Each trial must have exactly one row for reference and non-reference ",
        "subgroups; trial `",
        tr,
        "` does not.",
        call. = FALSE
      )
    }
    # CAMS interaction: g_j = y_Bj - y_Aj
    # Var(g_j) = v_Aj + v_Bj when Cov(y_Aj, y_Bj) is assumed zero.
    interaction_rows[[j]] <- data.frame(
      trial = tr,
      yi = as.numeric(d$yi[[i_non]]) - as.numeric(d$yi[[i_ref]]),
      vi = pmax(as.numeric(d$vi[[i_non]]), 0) + pmax(as.numeric(d$vi[[i_ref]]), 0),
      stringsAsFactors = FALSE
    )
  }

  list(
    overall = dplyr::bind_rows(overall_rows),
    interaction = dplyr::bind_rows(interaction_rows),
    pi = pi,
    reference = ref,
    nonref = nonref,
    ifrac_var = ifrac_var
  )
}
