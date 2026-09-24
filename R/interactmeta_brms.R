#' Fit CAMS via brms / Stan MCMC
#'
#' Undocumented CAMS backend kept for internal equivalence tests. Not part of
#' the public engine API (`engine = "bayesmeta"` or `"metafor"`).
#' Priors come from the CAMS-level specification (`prior_to_brms()`), not
#' brms defaults.
#' @noRd
.interactmeta_fit_brms <- function(data_fit,
                                   design,
                                   trial_var = "trial",
                                   subgroup_var = "subgroup",
                                   priors,
                                   engine_args = list(),
                                   ecological = TRUE,
                                   seed = 42L,
                                   silent = 2,
                                   ...) {
  check_brms()
  engine_args <- interactmeta_merge_engine_args(engine_args, list(...))
  engine_args[c(
    "method", "tau.prior", "tau.interaction.prior", "fixed.prior",
    "tau.prior.overall", "ecological", "ecological_slope",
    "formula", "data", "family", "prior", "priors"
  )] <- NULL

  data_brms <- data_fit
  data_brms$sei <- sqrt(pmax(as.numeric(data_brms$vi), 0))
  fml_chr <- if (inherits(design$formula, "formula")) {
    paste(deparse(design$formula), collapse = " ")
  } else {
    paste(as.character(design$formula), collapse = " ")
  }
  fml_txt <- gsub("se(sqrt(vi))", "se(sei)", fml_chr, fixed = TRUE)
  fml <- stats::as.formula(fml_txt)
  brms_prior <- prior_to_brms(
    priors,
    trial_var = trial_var,
    ifrac2_vars = design$coef_map$ifrac2_std
  )

  if (!"open_progress" %in% names(engine_args)) {
    engine_args$open_progress <- FALSE
  }
  if (!"refresh" %in% names(engine_args)) {
    engine_args$refresh <- 0
  }
  if (!"silent" %in% names(engine_args) && isTRUE(as.integer(silent)[[1L]] >= 2L)) {
    engine_args$silent <- 2
  }
  seed_arg <- if (is.na(seed)) NULL else as.integer(seed)[[1L]]
  if (!"seed" %in% names(engine_args) && !is.null(seed_arg)) {
    engine_args$seed <- seed_arg
  }

  brms_fit <- do.call(
    brms::brm,
    c(
      list(
        formula = fml,
        data = data_brms,
        family = stats::gaussian(),
        prior = brms_prior
      ),
      engine_args
    )
  )

  draws <- interactmeta_brms_draws_df(brms_fit, data_brms, design, trial_var)
  if (!isTRUE(ecological)) {
    draws <- comparison_zero_ifrac_slope_draws(draws, K = design$n_nonref)
  }
  mu <- interactmeta_brms_mu_matrix(brms_fit)

  structure(
    list(
      draws = draws,
      mu = mu,
      data = data_fit,
      formula = fml,
      brmsfit = brms_fit,
      engine = "brms",
      prior = prior_to_brms_rows(
        priors,
        trial_var = trial_var,
        ifrac2_vars = design$coef_map$ifrac2_std
      ),
      priors = priors,
      algorithm = "MCMC",
      method = "bayes"
    ),
    class = c("interactmeta", "cam", "cams", "list")
  )
}

interactmeta_brms_draws_df <- function(brms_fit, data_fit, design, trial_var) {
  dr <- as.data.frame(brms::as_draws_df(brms_fit))
  out <- as.data.frame(dr)
  if ("Intercept" %in% names(out) && !"b_Intercept" %in% names(out)) {
    names(out)[names(out) == "Intercept"] <- "b_Intercept"
  }
  drop <- grepl("^(lprior|lp__|sigma$)", names(out))
  out[, !drop, drop = FALSE]
}

interactmeta_brms_mu_matrix <- function(brms_fit) {
  ep <- brms::posterior_epred(brms_fit, re_formula = NULL)
  if (length(dim(ep)) == 3L) {
    ep <- ep[, , 1L]
  }
  matrix(as.numeric(ep), nrow = nrow(ep), ncol = ncol(ep))
}
