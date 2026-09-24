#' Fit CAMS via orthogonal metafor random-effects blocks
#'
#' Frequentist analogue of the binary CAMS factorization:
#' overall IF meta-regression and interaction RE meta-analysis.
#' Bayesian prior arguments are ignored.
#' @noRd
.interactmeta_fit_metafor <- function(data_fit,
                                      design,
                                      trial_var = "trial",
                                      subgroup_var = "subgroup",
                                      method = "REML",
                                      engine_args = list(),
                                      ecological = TRUE,
                                      silent = 2,
                                      priors = NULL,
                                      priors_user_supplied = FALSE,
                                      ...) {
  if (design$n_nonref != 1L) {
    stop(
      "`engine = \"metafor\"` currently supports binary subgroups only ",
      "(one non-reference `ifrac_*` column).",
      call. = FALSE
    )
  }
  if (!isTRUE(ecological)) {
    stop(
      "`ecological_slope = FALSE` is not implemented for `engine = \"metafor\"`.",
      call. = FALSE
    )
  }
  if (isTRUE(priors_user_supplied)) {
    warning(
      "`engine = \"metafor\"` is frequentist; Bayesian prior arguments are ignored.",
      call. = FALSE
    )
  }

  engine_args <- interactmeta_merge_engine_args(engine_args, list(...))
  engine_args[c(
    "method", "tau.prior", "tau.interaction.prior", "fixed.prior",
    "chains", "iter", "warmup", "thin", "refresh", "cores",
    "convolve_derived"
  )] <- NULL

  blocks <- interactmeta_bayesmeta_trial_blocks(
    data_fit = data_fit,
    design = design,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  ov <- blocks$overall
  int <- blocks$interaction

  overall_fit <- interactmeta_rma(
    yi = ov$yi,
    vi = ov$vi,
    mods = cbind(intercept = 1, ifrac = ov$ifrac),
    intercept = FALSE,
    method = method,
    slab = as.character(ov$trial),
    engine_args = engine_args
  )
  interaction_fit <- interactmeta_rma(
    yi = int$yi,
    vi = int$vi,
    method = method,
    slab = as.character(int$trial),
    engine_args = engine_args
  )

  alpha <- as.numeric(overall_fit$beta[[1L]])
  theta <- as.numeric(overall_fit$beta[[2L]])
  gamma <- as.numeric(interaction_fit$beta[[1L]])
  beta <- theta - gamma
  tau_m <- sqrt(max(as.numeric(overall_fit$tau2), 0))
  tau_g <- sqrt(max(as.numeric(interaction_fit$tau2), 0))

  sg_nonref <- design$coef_map$subgroup_fit[[2L]]
  ifrac2_std <- design$coef_map$ifrac2_std[[1L]]
  draws <- data.frame(
    .chain = 1L,
    .iteration = 1L,
    .draw = 1L,
    b_Intercept = alpha,
    b_ifrac_1 = beta,
    sd_trial__Intercept = tau_m,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  draws[[paste0("b_subgroup", sg_nonref)]] <- gamma
  draws[[paste0("sd_trial__", ifrac2_std)]] <- tau_g

  mu <- interactmeta_bayesmeta_mu_matrix(draws, data_fit, design)

  if (!isTRUE(as.integer(silent)[[1L]] >= 2L)) {
    message(
      "interactmeta metafor approach: orthogonal overall (rma meta-regression) + ",
      "interaction (rma) with method = ",
      method,
      "."
    )
  }

  structure(
    list(
      draws = draws,
      mu = mu,
      data = data_fit,
      formula = stats::as.formula(design$formula),
      metafor = list(
        overall = overall_fit,
        interaction = interaction_fit,
        trial_blocks = blocks
      ),
      engine = "metafor",
      prior = NULL,
      priors = priors,
      algorithm = "optimization",
      method = method
    ),
    class = c("interactmeta", "cam", "cams", "list")
  )
}

#' Call `metafor::rma`, pinning `tau2 = 0` when that variance is not identified.
#'
#' REML/ML count residual \(\tau^2\) as a free parameter. The overall CAMS
#' block is intercept + slope, so two trials saturate the mean model and
#' `rma()` refuses to estimate \(\tau_m^2\). The same point estimates are
#' available with residual heterogeneity fixed at 0.
#' @noRd
interactmeta_rma <- function(yi,
                             vi,
                             mods = NULL,
                             intercept = TRUE,
                             method,
                             slab,
                             engine_args = list()) {
  args <- c(
    list(
      yi = yi,
      vi = vi,
      method = method,
      slab = slab
    ),
    engine_args
  )
  if (!is.null(mods)) {
    args$mods <- mods
    args$intercept <- intercept
  }
  tryCatch(
    do.call(metafor::rma, args),
    error = function(e) {
      if (!grepl(
        "larger than the number of observations",
        conditionMessage(e),
        ignore.case = TRUE
      )) {
        stop(e)
      }
      if (!is.null(args$tau2)) {
        stop(e)
      }
      args$tau2 <- 0
      do.call(metafor::rma, args)
    }
  )
}
