# ---- interactmeta_metafor_latent.R ----
#' Frequentist CAMS latent structure (Wald / BLUP; mirrors bayesmeta `$latent` roles)
#' @noRd
interactmeta_metafor_latent <- function(overall_fit,
                                        interaction_fit,
                                        blocks,
                                        design,
                                        prob = 0.95) {
  z <- stats::qnorm(1 - (1 - as.numeric(prob)[[1L]]) / 2)
  wald <- function(est, se) {
    est <- as.numeric(est)[[1L]]
    se <- as.numeric(se)[[1L]]
    list(
      estimate = est,
      se = se,
      lower = if (is.finite(est) && is.finite(se)) est - z * se else NA_real_,
      upper = if (is.finite(est) && is.finite(se)) est + z * se else NA_real_
    )
  }
  alpha <- wald(overall_fit$beta[[1L]], overall_fit$se[[1L]])
  theta <- wald(overall_fit$beta[[2L]], overall_fit$se[[2L]])
  gamma <- wald(interaction_fit$beta[[1L]], interaction_fit$se[[1L]])
  eco <- as.numeric(theta$estimate) - as.numeric(gamma$estimate)
  se_eco <- sqrt(max(as.numeric(theta$se)^2 + as.numeric(gamma$se)^2, 0))
  beta <- wald(eco, se_eco)
  tau_m <- sqrt(max(as.numeric(overall_fit$tau2), 0))
  tau_g <- sqrt(max(as.numeric(interaction_fit$tau2), 0))
  trials <- as.character(blocks$overall$trial)
  pi <- stats::setNames(
    as.numeric(blocks$overall$ifrac),
    trials
  )
  list(
    alpha = alpha,
    theta = theta,
    gamma = gamma,
    beta = beta,
    tau_m = wald(tau_m, overall_fit$se.tau2 / (2 * max(tau_m, .Machine$double.eps))),
    tau_g = wald(tau_g, interaction_fit$se.tau2 / (2 * max(tau_g, .Machine$double.eps))),
    pi = pi,
    trials = trials,
    interval.type = "central",
    cache = new.env(parent = emptyenv())
  )
}
