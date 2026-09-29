#' CAMS study-level effects via the canonical linpred
#'
#' Default **shrinkage** (`evaluation = "trial"`) evaluates at each trial's own
#' information fraction \eqn{\pi_j} with stored random effects (same estimand as
#' the observed black squares):
#' \deqn{\widetilde\mu_{ij}
#'   =\alpha+u_j+\beta\pi_j+(x_i-\pi_j)v_j+\gamma\,1(i=B).}
#' Optional `evaluation = "overall"` evaluates the same stored REs at the pooled
#' IF \eqn{\omega} (trial-specific standardization — not shrinkage).
#' Orthogonal summaries use the same evaluation IF.
#'
#' @param evaluation `"trial"` (default, shrinkage at \eqn{\pi_j}) or
#'   `"overall"` (standardized to \eqn{\omega}).
#' @noRd
comparison_shrinkage_rotation_bundle <- function(fit,
                                           max_draws = NULL,
                                           contrast_sep = "-",
                                           ifrac = c("interaction", "overall"),
                                           overall_ifrac = NULL,
                                           evaluation = c("trial", "overall")) {
  ifrac <- match.arg(ifrac)
  evaluation <- match.arg(evaluation)
  check_comparison_fit(fit)
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st) || is.null(st$ifrac_data)) {
    stop(
      "Fit has no `ifrac_data` in attribute `\"interactmeta\"`; refit with `cam()`.",
      call. = FALSE
    )
  }
  trial_var <- st$trial_var %||% "trial"
  subgroup_var <- st$subgroup_var %||% "subgroup"
  ref <- interactmeta_reference(fit)

  prep <- comparison_prepare_marginal_draws(
    fit,
    overall_ifrac = overall_ifrac,
    ifrac = ifrac,
    max_draws = max_draws,
    contrast_sep = contrast_sep
  )
  post <- prep$post
  omega <- prep$omega
  n_draw <- nrow(post)
  design <- prep$design

  data_eval <- fit$data
  if (is.null(data_eval) || !nrow(data_eval)) {
    stop("Fit has no model data; refit with `cam()`.", call. = FALSE)
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
    # Canonical linpred: Shrinkage at π_j; Standardized at ω.
    # Used when stored trial REs exist (Stan). Missing r_trial columns
    # are treated as 0, which is *not* trial-level shrinkage — native
    # bayesmeta / metafor latents are handled above.
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
        re = "stored",
        fit = fit,
        trial = tr
      )
      y_B[, j] <- comparison_effect_linpred_draws(
        post = post,
        evaluation_if = evaluation_if,
        design = design,
        subgroup_label = nonref,
        subgroup_coef = prep$sg_coef[[nonref]],
        re = "stored",
        fit = fit,
        trial = tr
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
#' Returns `NULL` when the engine has stored `r_trial` columns (Stan) so
#' the linpred path remains authoritative there.
#' @noRd
comparison_shrinkage_native_trial_latents <- function(fit,
                                                post,
                                                trials,
                                                pi,
                                                reference,
                                                nonref) {
  if (any(grepl("^r_trial\\[", names(post)))) {
    return(NULL)
  }
  if (!is.null(fit$latent) &&
      !is.null(fit$latent$m) &&
      !is.null(fit$latent$g)) {
    return(comparison_shrinkage_from_bayesmeta_latent(
      fit = fit,
      post = post,
      trials = trials,
      pi = pi
    ))
  }
  if (identical(fit$engine, "metafor") && !is.null(fit$metafor)) {
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

#' Draws × trial matrices of orthogonal CAMS shrinkage (m, g)
#'
#' Forms the FE+RE linpreds at the reference / non-reference design rows, then
#' the orthogonal CAMS coordinates at the trial information fraction
#' \eqn{\pi_j}:
#' \deqn{m_j=(1-\pi_j)\mu_{A,j}+\pi_j\mu_{B,j},\qquad g_j=\mu_{B,j}-\mu_{A,j}.}
#' With the CAMS design \eqn{ifrac=\pi} and \eqn{ifrac2=x-\pi}, this equals
#' \eqn{m_j=\alpha+u_j+\beta\pi_j+\pi_j\gamma} and \eqn{g_j=\gamma+v_j}.
#' Missing trial REs are treated as 0.
#' @noRd
comparison_orthogonal_shrinkage_matrices <- function(post,
                                               trials,
                                               pi,
                                               nonref,
                                               pi_A,
                                               pi_B,
                                               i2_A,
                                               i2_B,
                                               dat = NULL,
                                               trial_var = "trial",
                                               subgroup_var = "subgroup",
                                               reference = NULL,
                                               ifrac_var = NULL) {
  # Prefer shared linpred when full design rows are available
  if (!is.null(dat) && !is.null(reference) && !is.null(ifrac_var)) {
    ifrac_cols <- if (startsWith(ifrac_var, "ifrac_")) {
      ifrac_var
    } else {
      paste0("ifrac_", ifrac_var)
    }
    mu <- comparison_linpred_matrix(
      post = post,
      dat = dat,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = reference,
      ifrac_cols = ifrac_cols,
      re = "stored"
    )
    n_draw <- nrow(mu)
    m <- matrix(NA_real_, nrow = n_draw, ncol = length(trials))
    g <- matrix(NA_real_, nrow = n_draw, ncol = length(trials))
    colnames(m) <- colnames(g) <- trials
    ref <- as.character(reference)[[1L]]
    for (j in seq_along(trials)) {
      tr <- trials[[j]]
      pij <- pi[[tr]]
      ii <- which(as.character(dat[[trial_var]]) == tr)
      sg <- as.character(dat[[subgroup_var]][ii])
      i_ref <- ii[sg == ref][[1L]]
      i_non <- ii[sg == nonref][[1L]]
      mu_A <- mu[, i_ref]
      mu_B <- mu[, i_non]
      m[, j] <- (1 - pij) * mu_A + pij * mu_B
      g[, j] <- mu_B - mu_A
    }
    return(list(m = m, g = g))
  }

  need <- c("b_Intercept", paste0("b_ifrac_", nonref), paste0("b_subgroup", nonref))
  missing <- setdiff(need, names(post))
  if (length(missing)) {
    stop(
      "Fit draws are missing orthogonal-shrinkage coefficients: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  alpha <- as.numeric(post[["b_Intercept"]])
  beta <- as.numeric(post[[paste0("b_ifrac_", nonref)]])
  gamma <- as.numeric(post[[paste0("b_subgroup", nonref)]])
  n_draw <- length(alpha)

  m <- matrix(NA_real_, nrow = n_draw, ncol = length(trials))
  g <- matrix(NA_real_, nrow = n_draw, ncol = length(trials))
  colnames(m) <- colnames(g) <- trials

  for (j in seq_along(trials)) {
    tr <- trials[[j]]
    pij <- pi[[tr]]
    u_col <- paste0("r_trial[", tr, ",Intercept]")
    v_col <- paste0("r_trial[", tr, ",ifrac2_", nonref, "]")
    u <- if (u_col %in% names(post)) {
      as.numeric(post[[u_col]])
    } else {
      rep(0, n_draw)
    }
    v <- if (v_col %in% names(post)) {
      as.numeric(post[[v_col]])
    } else {
      rep(0, n_draw)
    }
    mu_A <- alpha + u + beta * pi_A[[tr]] + v * i2_A[[tr]]
    mu_B <- alpha + u + beta * pi_B[[tr]] + v * i2_B[[tr]] + gamma
    m[, j] <- (1 - pij) * mu_A + pij * mu_B
    g[, j] <- mu_B - mu_A
  }
  list(m = m, g = g)
}

#' Fitted draws × obs matrix aligned to `post` rows
#' @noRd
comparison_shrinkage_fitted_matrix <- function(fit, post, data_eval) {
  fit_eval <- fit
  fit_eval$data <- data_eval
  fits <- fitted(fit_eval, summary = FALSE)
  post_full <- comparison_posterior_df(fit)
  n_obs <- nrow(data_eval)
  if (ncol(fits) == n_obs && nrow(fits) == nrow(post_full)) {
    fit_mat <- fits
  } else if (nrow(fits) == n_obs && ncol(fits) == nrow(post_full)) {
    fit_mat <- t(fits)
  } else {
    stop(
      "Unexpected dimensions from `fitted()`: expected ",
      nrow(post_full),
      " draws and ",
      n_obs,
      " observations.",
      call. = FALSE
    )
  }
  draw_idx <- match(post$.draw, post_full$.draw)
  if (any(is.na(draw_idx))) {
    stop("Could not align shrinkage posterior draws with `fitted()` output.", call. = FALSE)
  }
  fit_mat[draw_idx, , drop = FALSE]
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
  for (cand in c("subgroup", "group", "sg")) {
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
    # draws × trial matrices on the model scale (π_j + stored REs / native latents)
    draws = list(
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
