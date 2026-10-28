# ---- forestplot_construction.R ----
#' Forest plots for DA, AD, SWADA, and conventional fits
#' @noRd
NULL

#' Resolve a construction fit for forest plotting (conventional → DA)
#' @noRd
construction_forest_resolve_fit <- function(fit) {
  if (inherits(fit, "interactmeta_conventional")) {
    if (is.null(fit$da)) {
      stop("conventional() fit has no `da` component.", call. = FALSE)
    }
    return(fit$da)
  }
  fit
}

#' Long trial × subgroup data for forest panels (CAMS-like `ifrac_data` layout)
#' @noRd
construction_forest_long_data <- function(fit) {
  fit <- construction_forest_resolve_fit(fit)
  dat <- fit$data
  if (is.null(dat) || !nrow(dat)) {
    stop("Fit has no `data` for the forest plot.", call. = FALSE)
  }
  study_col <- interactmeta_resolve_study_col(dat)
  level_col <- fit$modifier %||% interactmeta_resolve_level_col(dat)
  out <- data.frame(
    trial = as.character(dat[[study_col]]),
    subgroup = as.character(dat[[level_col]]),
    yi = as.numeric(dat$yi),
    vi = as.numeric(dat$vi),
    reference = fit$reference,
    stringsAsFactors = FALSE
  )
  out[[study_col]] <- out$trial
  if (!identical(level_col, "subgroup")) {
    out[[level_col]] <- out$subgroup
  }
  if ("n" %in% names(dat)) {
    out$n <- dat$n
  }
  out
}

#' Attach minimal `interactmeta` attribute for shared forest builders
#' @noRd
construction_forest_with_st <- function(fit) {
  fit <- construction_forest_resolve_fit(fit)
  long <- construction_forest_long_data(fit)
  study_col <- interactmeta_resolve_study_col(fit$data)
  mod <- fit$modifier %||% interactmeta_resolve_level_col(fit$data)
  st <- list(
    reference = fit$reference,
    trial_var = study_col,
    subgroup_var = mod,
    subgroup_levels = c(fit$reference, fit$nonref),
    ifrac_data = long,
    coef_map = list(ifrac_orig = fit$nonref)
  )
  out <- fit
  attr(out, "interactmeta") <- st
  out
}

#' Reference line for construction forests
#' @noRd
construction_forest_reference <- function(fit) {
  fit <- construction_forest_resolve_fit(fit)
  as.character(fit$reference)[[1L]]
}

#' RE-MA of within-trial contrasts (AD model; used for DA interaction layers)
#' @noRd
construction_forest_interaction_rma <- function(fit) {
  fit <- construction_forest_resolve_fit(fit)
  if (!is.null(fit$fits$interaction)) {
    return(fit$fits$interaction)
  }
  pairs <- fit$pairs
  if (is.null(pairs) || !nrow(pairs)) {
    return(NULL)
  }
  method <- fit$method %||% "REML"
  if (interactmeta_is_bayesian_fit(fit)) {
    return(NULL)
  }
  suppressWarnings(
    metafor::rma(
      yi = pairs$g,
      vi = pairs$vg,
      slab = pairs$trial,
      method = method
    )
  )
}

#' RE-MA of combined trial effects m_j
#' @noRd
construction_forest_combined_rma <- function(fit) {
  fit <- construction_forest_resolve_fit(fit)
  if (!is.null(fit$fits$combined)) {
    return(fit$fits$combined)
  }
  pairs <- fit$pairs
  if (is.null(pairs) || !nrow(pairs)) {
    return(NULL)
  }
  vi_m <- with(pairs, v_ref * v_nonref / (v_ref + v_nonref))
  method <- fit$method %||% "REML"
  if (interactmeta_is_bayesian_fit(fit)) {
    return(NULL)
  }
  fit_m <- suppressWarnings(
    metafor::rma(
      yi = pairs$m,
      vi = vi_m,
      slab = pairs$trial,
      method = method
    )
  )
  fit$fits$combined <- fit_m
  fit_m
}

#' Merge shrinkage onto forest template rows (unmatched rows stay NA)
#' @noRd
construction_forest_merge_shrinkage_yivi <- function(template,
                                                     summary,
                                                     trial_var,
                                                     label_var,
                                                     summary_label_var = label_var,
                                                     use_posterior_qi = FALSE) {
  if (!nrow(template)) {
    return(template)
  }
  out <- template
  out$yi <- NA_real_
  out$vi <- NA_real_
  out$ci.lb <- NA_real_
  out$ci.ub <- NA_real_
  if (is.null(summary) || !nrow(summary)) {
    return(out)
  }
  lab_col <- summary_label_var
  if (!lab_col %in% names(summary)) {
    lab_col <- label_var
  }
  key_t <- paste(
    as.character(template[[trial_var]]),
    as.character(template[[label_var]]),
    sep = "\r"
  )
  key_s <- paste(
    as.character(summary$trial),
    as.character(summary[[lab_col]]),
    sep = "\r"
  )
  idx <- match(key_t, key_s)
  ok <- !is.na(idx)
  if (!any(ok)) {
    return(out)
  }
  s <- summary[idx[ok], , drop = FALSE]
  if (isTRUE(use_posterior_qi) &&
        all(c("estimate", "lci", "uci") %in% names(s))) {
    out$yi[ok] <- as.numeric(s$estimate)
    out$vi[ok] <- pmax(as.numeric(s$sd)^2, 1e-12, na.rm = TRUE)
    out$ci.lb[ok] <- as.numeric(s$lci)
    out$ci.ub[ok] <- as.numeric(s$uci)
  } else {
    out$yi[ok] <- s$mean
    out$vi[ok] <- pmax(s$sd^2, 1e-12, na.rm = TRUE)
    if (all(c("lci", "uci") %in% names(s))) {
      out$ci.lb[ok] <- as.numeric(s$lci)
      out$ci.ub[ok] <- as.numeric(s$uci)
    }
  }
  out
}

construction_forest_na_shrinkage_summary <- function(sg_dat) {
  data.frame(
    trial = as.character(sg_dat$trial),
    subgroup = as.character(sg_dat$subgroup),
    type = "Shrinkage",
    mean = NA_real_,
    sd = NA_real_,
    estimate = NA_real_,
    lci = NA_real_,
    uci = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Shrinkage summary rows for one metafor `rma` fit
#' @noRd
construction_forest_rma_shrinkage_summary <- function(rma_fit,
                                                      subgroup_label,
                                                      prob = 0.95) {
  if (is.null(rma_fit) || !inherits(rma_fit, "rma")) {
    return(NULL)
  }
  pr <- suppressWarnings(
    predict(rma_fit, level = 100 * prob)
  )
  data.frame(
    trial = as.character(rma_fit$slab),
    subgroup = subgroup_label,
    type = "Shrinkage",
    mean = as.numeric(pr$pred),
    sd = as.numeric(pr$se),
    estimate = as.numeric(pr$pred),
    lci = as.numeric(pr$ci.lb),
    uci = as.numeric(pr$ci.ub),
    stringsAsFactors = FALSE
  )
}

#' Shrinkage summary for one bayesmeta trial block
#' @noRd
construction_forest_bayesmeta_shrinkage_summary <- function(bm_fit,
                                                            subgroup_label,
                                                            prob = 0.95) {
  if (is.null(bm_fit) || !inherits(bm_fit, "bayesmeta")) {
    return(NULL)
  }
  n <- as.integer(bm_fit$k)[[1L]]
  if (!is.finite(n) || n < 1L) {
    return(NULL)
  }
  trials <- as.character(bm_fit$slab)
  if (length(trials) != n) {
    labels <- bm_fit$labels
    if (length(labels) == n) {
      trials <- as.character(labels)
    } else if (!is.null(bm_fit$data$y) && length(bm_fit$data$y) == n) {
      trials <- as.character(seq_len(n))
    } else {
      trials <- as.character(seq_len(n))
    }
  }
  rows <- lapply(seq_len(n), function(i) {
    loc <- interactmeta_native_location(
      bm_fit,
      prob = prob,
      interval.type = "shortest",
      individual = i
    )
    sd_i <- (loc$upper - loc$lower) /
      (2 * stats::qnorm(1 - (1 - prob) / 2))
    data.frame(
      trial = trials[[i]],
      subgroup = subgroup_label,
      type = "Shrinkage",
      mean = loc$median,
      sd = sd_i,
      estimate = loc$median,
      lci = loc$lower,
      uci = loc$upper,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Pooled marginal row (Overall / Prediction) for one subgroup label
#' @noRd
construction_forest_marginal_subgroup_row <- function(label,
                                                      estimate,
                                                      se,
                                                      ci.lb,
                                                      ci.ub,
                                                      row_kind,
                                                      ref) {
  data.frame(
    trial = label,
    subgroup = row_kind,
    colour_label = label,
    row_kind = row_kind,
    yi = estimate,
    vi = se^2,
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    stringsAsFactors = FALSE
  )
}

#' DA interaction prediction: \eqn{\theta_B - \theta_A} from independent subgroup RE-MAs
#' @noRd
construction_forest_da_interaction_prediction <- function(fit, prob) {
  fit <- construction_forest_resolve_fit(fit)
  fit_ref <- fit$fits$reference
  fit_non <- fit$fits$nonref
  if (is.null(fit_ref) || is.null(fit_non)) {
    return(NULL)
  }
  z <- stats::qnorm(1 - (1 - prob) / 2)
  if (interactmeta_is_bayesian_fit(fit)) {
    if (!inherits(fit_ref, "bayesmeta") || !inherits(fit_non, "bayesmeta")) {
      return(NULL)
    }
    loc_ref <- interactmeta_native_location(
      fit_ref,
      prob = prob,
      interval.type = "shortest",
      predict = TRUE
    )
    loc_non <- interactmeta_native_location(
      fit_non,
      prob = prob,
      interval.type = "shortest",
      predict = TRUE
    )
    se_ref <- (loc_ref$upper - loc_ref$lower) / (2 * z)
    se_non <- (loc_non$upper - loc_non$lower) / (2 * z)
    est <- loc_non$median - loc_ref$median
    se <- sqrt(se_non^2 + se_ref^2)
    return(list(
      estimate = est,
      se = se,
      ci.lb = est - z * se,
      ci.ub = est + z * se
    ))
  }
  pr_ref <- construction_forest_rma_marginal(fit_ref, prob, "Prediction")
  pr_non <- construction_forest_rma_marginal(fit_non, prob, "Prediction")
  if (is.null(pr_ref) || is.null(pr_non)) {
    return(NULL)
  }
  est <- pr_non$estimate - pr_ref$estimate
  se <- sqrt(pr_non$se^2 + pr_ref$se^2)
  list(
    estimate = est,
    se = se,
    ci.lb = est - z * se,
    ci.ub = est + z * se
  )
}

#' Overall / prediction intervals from a stored `rma` or summary row
#' @noRd
construction_forest_rma_marginal <- function(rma_fit, prob, kind = c("Overall", "Prediction")) {
  kind <- match.arg(kind)
  if (is.null(rma_fit)) {
    return(NULL)
  }
  if (kind == "Overall") {
    z <- stats::qnorm(1 - (1 - prob) / 2)
    est <- as.numeric(rma_fit$beta)[[1L]]
    se <- as.numeric(rma_fit$se)[[1L]]
    return(list(
      estimate = est,
      se = se,
      ci.lb = est - z * se,
      ci.ub = est + z * se
    ))
  }
  if (!inherits(rma_fit, "rma")) {
    return(NULL)
  }
  pr <- suppressWarnings(predict(rma_fit, level = 100 * prob))
  list(
    estimate = as.numeric(pr$pred),
    se = as.numeric(pr$se),
    ci.lb = as.numeric(pr$ci.lb),
    ci.ub = as.numeric(pr$ci.ub)
  )
}

#' SWADA WLS prediction interval (common weights, interaction tau)
#' @noRd
construction_forest_swada_prediction_row <- function(fit, label, prob) {
  interactmeta_swada_prediction(
    fit,
    target = "subgroup",
    label = label,
    prob = prob
  )
}

#' SWADA interaction prediction (same WLS + \eqn{\tau_g} as subgroup panels)
#' @noRd
construction_forest_swada_interaction_prediction <- function(fit, prob) {
  interactmeta_swada_prediction(fit, target = "interaction", prob = prob)
}

#' Overall / Prediction rows when only interaction marginals exist (e.g. AD)
#'
#' Panel A carries structure/labels; panel B carries pooled interaction estimates.
#' @noRd
construction_forest_interaction_marginals_block <- function(overall_int,
                                                            pred_int,
                                                            include_overall = TRUE,
                                                            include_prediction = TRUE) {
  empty_sg <- data.frame(
    trial = character(0),
    subgroup = character(0),
    yi = numeric(0),
    vi = numeric(0),
    stringsAsFactors = FALSE
  )
  empty_b <- data.frame(
    yi = numeric(0),
    vi = numeric(0),
    ci.lb = numeric(0),
    ci.ub = numeric(0),
    stringsAsFactors = FALSE
  )
  add_pair <- function(int_row, row_kind, slot_label) {
    if (is.null(int_row) || !nrow(int_row)) {
      return(NULL)
    }
    sg <- int_row[1L, , drop = FALSE]
    sg$trial <- "Interaction"
    sg$subgroup <- slot_label
    sg$colour_label <- "Interaction"
    sg$row_kind <- row_kind
    sg$yi <- NA_real_
    sg$vi <- NA_real_
    sg$ci.lb <- NA_real_
    sg$ci.ub <- NA_real_
    b <- int_row[1L, , drop = FALSE]
    b$row_kind <- row_kind
    list(
      subgroup = sg,
      panel_b = data.frame(
        yi = b$yi,
        vi = b$vi,
        ci.lb = b$ci.lb,
        ci.ub = b$ci.ub,
        row_kind = row_kind,
        stringsAsFactors = FALSE
      )
    )
  }
  rows_sg <- list()
  rows_b <- list()
  if (isTRUE(include_overall)) {
    pair <- add_pair(overall_int, "overall", "Overall")
    if (!is.null(pair)) {
      rows_sg <- c(rows_sg, list(pair$subgroup))
      rows_b <- c(rows_b, list(pair$panel_b))
    }
  }
  if (isTRUE(include_prediction)) {
    pair <- add_pair(pred_int, "prediction", "Prediction")
    if (!is.null(pair)) {
      rows_sg <- c(rows_sg, list(pair$subgroup))
      rows_b <- c(rows_b, list(pair$panel_b))
    }
  }
  if (!length(rows_sg)) {
    return(list(subgroup = empty_sg, panel_b = empty_b))
  }
  list(
    subgroup = dplyr::bind_rows(rows_sg),
    panel_b = dplyr::bind_rows(rows_b)
  )
}

#' Build Overall / Prediction blocks for construction forests
#' @noRd
construction_forest_summary_blocks <- function(fit,
                                               prob,
                                               contrast_sep,
                                               ref,
                                               include_overall = TRUE,
                                               include_prediction = TRUE,
                                               include_combined = TRUE,
                                               model_names = NULL) {
  fit <- construction_forest_resolve_fit(fit)
  analysis <- interactmeta_fit_analysis(fit)
  show_overall <- comparison_forest_layer_on(include_overall)
  show_prediction <- comparison_forest_layer_on(include_prediction)
  empty_sg <- data.frame(
    trial = character(0),
    subgroup = character(0),
    yi = numeric(0),
    vi = numeric(0),
    stringsAsFactors = FALSE
  )
  overall_sg <- empty_sg
  overall_int <- empty_sg
  pred_sg <- empty_sg
  pred_int <- empty_sg

  add_sg_rows <- function(kind) {
    if (is.null(fit$subgroups) || !nrow(fit$subgroups)) {
      return(empty_sg)
    }
    sg <- fit$subgroups
    rows <- lapply(seq_len(nrow(sg)), function(i) {
      lab <- as.character(sg$label[[i]])
      if (kind == "Overall") {
        construction_forest_marginal_subgroup_row(
          label = lab,
          estimate = sg$estimate[[i]],
          se = sg$se[[i]],
          ci.lb = sg$ci.lb[[i]],
          ci.ub = sg$ci.ub[[i]],
          row_kind = "overall",
          ref = ref
        )
      } else if (identical(analysis, "swada")) {
        pr <- construction_forest_swada_prediction_row(fit, lab, prob)
        construction_forest_marginal_subgroup_row(
          label = lab,
          estimate = pr$estimate,
          se = pr$se,
          ci.lb = pr$ci.lb,
          ci.ub = pr$ci.ub,
          row_kind = "prediction",
          ref = ref
        )
      } else {
        rma_fit <- if (identical(lab, fit$reference)) {
          fit$fits$reference
        } else {
          fit$fits$nonref
        }
        pr <- if (interactmeta_is_bayesian_fit(fit) &&
            inherits(rma_fit, "bayesmeta")) {
          loc <- interactmeta_native_location(
            rma_fit,
            prob = prob,
            interval.type = "shortest",
            predict = TRUE
          )
          se <- (loc$upper - loc$lower) /
            (2 * stats::qnorm(1 - (1 - prob) / 2))
          list(
            estimate = loc$median,
            se = se,
            ci.lb = loc$lower,
            ci.ub = loc$upper
          )
        } else {
          construction_forest_rma_marginal(rma_fit, prob, "Prediction")
        }
        if (is.null(pr)) {
          return(NULL)
        }
        construction_forest_marginal_subgroup_row(
          label = lab,
          estimate = pr$estimate,
          se = pr$se,
          ci.lb = pr$ci.lb,
          ci.ub = pr$ci.ub,
          row_kind = "prediction",
          ref = ref
        )
      }
    })
    rows <- rows[!vapply(rows, is.null, logical(1L))]
    if (!length(rows)) {
      return(empty_sg)
    }
    dplyr::bind_rows(rows)
  }

  int_rma <- construction_forest_interaction_rma(fit)
  add_int_row <- function(kind) {
    int <- fit$interaction
    if (is.null(int) || !nrow(int)) {
      return(empty_sg)
    }
    if (kind == "Overall") {
      construction_forest_marginal_subgroup_row(
        label = "Interaction",
        estimate = int$estimate[[1L]],
        se = int$se[[1L]],
        ci.lb = int$ci.lb[[1L]],
        ci.ub = int$ci.ub[[1L]],
        row_kind = "overall",
        ref = ref
      )
    } else {
      pr <- if (identical(analysis, "da")) {
        construction_forest_da_interaction_prediction(fit, prob)
      } else if (identical(analysis, "swada")) {
        construction_forest_swada_interaction_prediction(fit, prob)
      } else if (!is.null(int_rma)) {
        if (inherits(int_rma, "bayesmeta")) {
          loc <- interactmeta_native_location(
            int_rma,
            prob = prob,
            interval.type = "shortest",
            predict = TRUE
          )
          se <- (loc$upper - loc$lower) /
            (2 * stats::qnorm(1 - (1 - prob) / 2))
          list(
            estimate = loc$median,
            se = se,
            ci.lb = loc$lower,
            ci.ub = loc$upper
          )
        } else {
          construction_forest_rma_marginal(int_rma, prob, "Prediction")
        }
      } else {
        NULL
      }
      if (is.null(pr) || !length(pr$estimate)) {
        if (interactmeta_is_bayesian_fit(fit) && !is.null(fit$posterior$interaction)) {
          sm <- fit$interaction
          return(data.frame(
            trial = "Interaction",
            subgroup = "Prediction",
            colour_label = "Interaction",
            row_kind = "prediction",
            yi = sm$estimate,
            vi = sm$se^2,
            ci.lb = sm$ci.lb,
            ci.ub = sm$ci.ub,
            stringsAsFactors = FALSE
          ))
        }
        return(empty_sg)
      }
      data.frame(
        trial = "Interaction",
        subgroup = "Prediction",
        colour_label = "Interaction",
        row_kind = "prediction",
        yi = pr$estimate,
        vi = pr$se^2,
        ci.lb = pr$ci.lb,
        ci.ub = pr$ci.ub,
        stringsAsFactors = FALSE
      )
    }
  }

  if (show_overall) {
    overall_sg <- add_sg_rows("Overall")
    overall_int <- add_int_row("Overall")
  }
  if (show_prediction) {
    pred_sg <- add_sg_rows("Prediction")
    pred_int <- add_int_row("Prediction")
  }

  interaction_only_marginals <- identical(analysis, "ad") ||
    (!nrow(overall_sg) && !nrow(pred_sg) &&
      (nrow(overall_int) > 0L || nrow(pred_int) > 0L))
  combined <- if (isTRUE(interaction_only_marginals)) {
    construction_forest_interaction_marginals_block(
      overall_int = overall_int,
      pred_int = pred_int,
      include_overall = show_overall,
      include_prediction = show_prediction
    )
  } else {
    comparison_forest_interleave_overall_prediction(
      overall_sg = overall_sg,
      pred_sg = pred_sg,
      overall_int = overall_int,
      pred_int = pred_int,
      reference = ref,
      include_overall = show_overall,
      include_prediction = show_prediction
    )
  }

  if (comparison_forest_layer_on(include_combined) && !identical(analysis, "ad")) {
    overall_types <- character(0)
    if (show_overall || show_prediction) {
      overall_types <- c("Overall", "Prediction")
    }
    if (length(overall_types)) {
      comb_rma <- construction_forest_combined_rma(fit)
      rows <- lapply(overall_types, function(ty) {
        kind <- if (ty == "Overall") "Overall" else "Prediction"
        if (is.null(comb_rma)) {
          return(data.frame(
            trial = "Combined",
            subgroup = ty,
            colour_label = forestplot_combined_block_label(),
            row_kind = "combined",
            yi = NA_real_,
            vi = NA_real_,
            ci.lb = NA_real_,
            ci.ub = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
        pr <- construction_forest_rma_marginal(comb_rma, prob, kind)
        data.frame(
          trial = "Combined",
          subgroup = ty,
          colour_label = forestplot_combined_block_label(),
          row_kind = "combined",
          yi = pr$estimate,
          vi = pr$se^2,
          ci.lb = pr$ci.lb,
          ci.ub = pr$ci.ub,
          stringsAsFactors = FALSE
        )
      })
      extra_sg <- dplyr::bind_rows(rows)
      extra_b <- data.frame(
        yi = extra_sg$yi,
        vi = extra_sg$vi,
        ci.lb = extra_sg$ci.lb,
        ci.ub = extra_sg$ci.ub,
        row_kind = extra_sg$row_kind,
        stringsAsFactors = FALSE
      )
      combined$subgroup <- dplyr::bind_rows(combined$subgroup, extra_sg)
      combined$panel_b <- dplyr::bind_rows(combined$panel_b, extra_b)
    }
  }

  list(list(
    name = "",
    subgroup = combined$subgroup,
    panel_b = combined$panel_b,
    interaction = NULL
  ))
}

#' Shrinkage data for construction forests
#' @noRd
construction_forest_shrinkage_data <- function(fit,
                                               prob,
                                               contrast_sep,
                                               group,
                                               shrinkage_effects = NULL) {
  fit <- construction_forest_resolve_fit(fit)
  st <- attr(fit, "interactmeta", exact = TRUE)
  ref <- st$reference
  trial_var <- st$trial_var
  subgroup_var <- st$subgroup_var
  sg_dat <- comparison_forest_subgroup_rows(
    st$ifrac_data,
    trial_var,
    subgroup_var,
    contrast_sep
  )
  if (is.null(shrinkage_effects)) {
    stop("Internal error: missing shrinkage effects.", call. = FALSE)
  }
  me_sg <- shrinkage_effects$subgroups
  me_int <- shrinkage_effects$interaction
  me_ov <- shrinkage_effects$overall
  sg <- construction_forest_merge_shrinkage_yivi(
    template = sg_dat,
    summary = me_sg,
    trial_var = trial_var,
    label_var = subgroup_var,
    use_posterior_qi = interactmeta_is_bayesian_fit(fit)
  )
  sg <- comparison_forest_order_to_template(
    sg,
    split_ifrac_for_forest(
      validate_forestplot_data(st$ifrac_data),
      contrast_sep = contrast_sep,
      group = group
    )$subgroup,
    trial_var = trial_var,
    subgroup_var = subgroup_var
  )
  int_template <- build_contrast_rows_from_subgroups(
    sg,
    contrast_sep = contrast_sep,
    reference = ref
  )
  int <- construction_forest_merge_shrinkage_yivi(
    template = int_template,
    summary = me_int,
    trial_var = trial_var,
    label_var = subgroup_var,
    summary_label_var = "interaction",
    use_posterior_qi = interactmeta_is_bayesian_fit(fit)
  )
  ov <- comparison_forest_overall_shrinkage_yivi(
    summary = me_ov,
    trial_var = trial_var
  )
  list(subgroup = sg, interaction = int, overall = ov)
}

#' Shrinkage effect summaries for construction fits
#' @noRd
construction_forest_shrinkage_effects <- function(fit, prob, contrast_sep = "-") {
  fit <- construction_forest_resolve_fit(fit)
  ref <- fit$reference
  nonref <- fit$nonref
  is_bayes <- interactmeta_is_bayesian_fit(fit)
  analysis <- interactmeta_fit_analysis(fit)
  sg_summ <- NULL
  if (identical(analysis, "ad")) {
    st <- attr(construction_forest_with_st(fit), "interactmeta", exact = TRUE)
    sg_dat <- comparison_forest_subgroup_rows(
      st$ifrac_data,
      st$trial_var,
      st$subgroup_var,
      contrast_sep
    )
    if (nrow(sg_dat)) {
      sg_summ <- construction_forest_na_shrinkage_summary(sg_dat)
    }
  } else if (identical(analysis, "swada")) {
    st <- attr(construction_forest_with_st(fit), "interactmeta", exact = TRUE)
    sg_dat <- comparison_forest_subgroup_rows(
      st$ifrac_data,
      st$trial_var,
      st$subgroup_var,
      contrast_sep
    )
    if (nrow(sg_dat)) {
      sg_summ <- construction_forest_na_shrinkage_summary(sg_dat)
    }
  } else if (!is.null(fit$subgroups)) {
    if (is_bayes) {
      s_ref <- construction_forest_bayesmeta_shrinkage_summary(
        fit$fits$reference, ref, prob
      )
      s_non <- construction_forest_bayesmeta_shrinkage_summary(
        fit$fits$nonref, nonref, prob
      )
      sg_summ <- dplyr::bind_rows(s_ref, s_non)
    } else if (!is.null(fit$fits$reference) && !is.null(fit$fits$nonref)) {
      s_ref <- construction_forest_rma_shrinkage_summary(
        fit$fits$reference, ref, prob
      )
      s_non <- construction_forest_rma_shrinkage_summary(
        fit$fits$nonref, nonref, prob
      )
      sg_summ <- dplyr::bind_rows(s_ref, s_non)
    }
  }
  int_fit <- construction_forest_interaction_rma(fit)
  int_lab <- paste0(nonref, contrast_sep, ref)
  int_summ <- if (is_bayes && !is.null(fit$fits$interaction)) {
    tmp <- construction_forest_bayesmeta_shrinkage_summary(
      fit$fits$interaction,
      int_lab,
      prob
    )
    tmp$interaction <- int_lab
    tmp
  } else if (!is.null(int_fit)) {
    tmp <- construction_forest_rma_shrinkage_summary(
      int_fit,
      int_lab,
      prob
    )
    tmp$interaction <- int_lab
    tmp
  } else {
    NULL
  }
  comb <- construction_forest_combined_rma(fit)
  ov_summ <- if (!is.null(comb)) {
    construction_forest_rma_shrinkage_summary(comb, "Combined", prob)
  } else {
    NULL
  }
  list(
    subgroups = sg_summ,
    interaction = int_summ,
    overall = ov_summ
  )
}

#' Heterogeneity footnotes for construction forests (subgroup-specific tau)
#' @noRd
construction_forest_annotations <- function(fit,
                                            prob = 0.95,
                                            digits = 2,
                                            show_heterogeneity = TRUE) {
  if (!isTRUE(show_heterogeneity)) {
    return(list(
      lines = character(),
      ifrac = character(),
      subgroup_ifrac = character(),
      heterogeneity = character()
    ))
  }
  het <- tryCatch(
    interactmeta_extract_heterogeneity(fit, prob = prob),
    error = function(e) NULL
  )
  het_notes <- character()
  if (!is.null(het) && !is.null(het$summary) && nrow(het$summary)) {
    summ <- het$summary
    fmt <- function(v) {
      format(round(as.numeric(v), digits), nsmall = digits, trim = TRUE)
    }
    sg_parts <- character()
    int_part <- character()
    for (i in seq_len(nrow(summ))) {
      par <- as.character(summ$parameter[[i]])
      piece <- paste0(
        par,
        ": ",
        fmt(summ$estimate[[i]])
      )
      if (is.finite(summ$lci[[i]]) && is.finite(summ$uci[[i]])) {
        piece <- paste0(
          piece,
          " (",
          fmt(summ$lci[[i]]),
          ", ",
          fmt(summ$uci[[i]]),
          ")"
        )
      }
      if (grepl("^tau\\[interaction\\]", par)) {
        int_part <- piece
      } else {
        sg_parts <- c(sg_parts, piece)
      }
    }
    if (length(sg_parts)) {
      het_notes[["combined"]] <- paste(sg_parts, collapse = " \u00b7 ")
    }
    if (length(int_part) == 1L && nzchar(int_part)) {
      het_notes[["interaction"]] <- int_part
    }
  }
  list(
    lines = character(),
    ifrac = character(),
    subgroup_ifrac = character(),
    heterogeneity = het_notes
  )
}

#' Effect bundle for construction forests
#' @noRd
construction_forest_effect_bundle <- function(fit,
                                              prob = 0.95,
                                              contrast_sep = "-",
                                              include_shrinkage = TRUE,
                                              include_overall = TRUE,
                                              include_prediction = TRUE,
                                              include_combined = TRUE,
                                              show_annotations = TRUE,
                                              annotation_digits = 2,
                                              model_names = NULL) {
  fit_shell <- construction_forest_with_st(fit)
  ref <- construction_forest_reference(fit)
  need_shrinkage <- comparison_forest_layer_on(include_shrinkage)
  shrinkage_effects <- if (need_shrinkage) {
    construction_forest_shrinkage_effects(
      fit_shell,
      prob = prob,
      contrast_sep = contrast_sep
    )
  } else {
    NULL
  }
  summary_blocks <- construction_forest_summary_blocks(
    fit = fit_shell,
    prob = prob,
    contrast_sep = contrast_sep,
    ref = ref,
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
    model_names = model_names
  )
  annotations <- construction_forest_annotations(
    fit = construction_forest_resolve_fit(fit),
    prob = prob,
    digits = as.integer(annotation_digits)[[1L]],
    show_heterogeneity = show_annotations
  )
  list(
    overall_ifrac = NULL,
    shrinkage = shrinkage_effects,
    summary_blocks = summary_blocks,
    annotations = annotations
  )
}

#' Build construction forest panels (shared layout with CAMS)
#' @noRd
construction_forest_build_fit_panels <- function(fit,
                                                 prob,
                                                 contrast_sep,
                                                 model_names,
                                                 include_shrinkage = TRUE,
                                                 include_overall = TRUE,
                                                 include_prediction = TRUE,
                                                 include_combined = TRUE,
                                                 n_var = FALSE,
                                                 group = "subgroup",
                                                 effect_bundle = NULL) {
  fit <- construction_forest_with_st(fit)
  effect_bundle <- effect_bundle %||% construction_forest_effect_bundle(
    fit,
    prob = prob,
    contrast_sep = contrast_sep,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
    show_annotations = TRUE,
    model_names = model_names
  )
  st <- attr(fit, "interactmeta", exact = TRUE)
  group <- match.arg(group, c("subgroup", "trial"))
  trial_var <- st$trial_var
  subgroup_var <- st$subgroup_var
  model_names <- comparison_forest_resolve_model_names(
    model_names = model_names,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction
  )
  data <- validate_forestplot_data(st$ifrac_data)
  panels <- split_ifrac_for_forest(
    data,
    contrast_sep = contrast_sep,
    group = group
  )
  ref <- construction_forest_reference(fit)
  aligned <- align_interaction_to_subgroups(
    panels$subgroup,
    panels$interaction,
    reference = ref
  )
  ifrac_disp <- forestplot_row_information_fraction(
    aligned$subgroup,
    contrast_sep = contrast_sep,
    digits = 2
  )
  resultsA <- forestplot_panel_df(
    aligned$subgroup,
    prob = prob,
    ifrac_display = ifrac_disp
  )
  resultsA <- forestplot_attach_trial_size(
    resultsA,
    aligned$subgroup,
    trial_col = "trial",
    subgroup_col = "subgroup",
    n_var = n_var,
    contrast_sep = contrast_sep
  )
  resultsB <- forestplot_panel_df(aligned$interaction, prob = prob)
  n_obs <- nrow(resultsA)
  row_labels <- forestplot_study_row_labels(
    trial = aligned$subgroup$trial,
    subgroup = aligned$subgroup$subgroup,
    group = group
  )
  resultsA$fill <- "white"
  resultsA$colour <- comparison_forest_marker_colour()
  resultsA$bold <- FALSE
  resultsB$fill <- "white"
  resultsB$colour <- comparison_forest_marker_colour()
  resultsB$bold <- FALSE
  resultsA$row_kind <- "study"
  resultsB$row_kind <- "study"
  ix_obs <- seq_len(n_obs)
  ix_sg_shr <- integer(0)
  ix_int_shr <- integer(0)
  ix_shr <- integer(0)
  overall_shrinkage <- NULL
  if (comparison_forest_layer_on(include_shrinkage) &&
        !is.null(effect_bundle$shrinkage)) {
    shrinkage <- construction_forest_shrinkage_data(
      fit,
      prob = prob,
      contrast_sep = contrast_sep,
      group = group,
      shrinkage_effects = effect_bundle$shrinkage
    )
    overall_shrinkage <- shrinkage$overall
    inter <- interleave_comparison_shrinkage_rows(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      shrinkage = shrinkage,
      prob = prob,
      ref = ref,
      contrast_sep = contrast_sep,
      ifrac_display = ifrac_disp
    )
    resultsA <- inter$resultsA
    resultsB <- inter$resultsB
    row_labels <- inter$row_labels
    ix_obs <- inter$obs_row_index
    ix_sg_shr <- inter$subgroup_shrinkage_row_index
    ix_int_shr <- inter$interaction_shrinkage_row_index
    ix_shr <- c(ix_sg_shr, ix_int_shr)
    styled_shr <- comparison_forest_apply_shrinkage_colours(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      subgroup_shrinkage_row_index = ix_sg_shr,
      interaction_shrinkage_row_index = ix_int_shr,
      ref = ref
    )
    resultsA <- styled_shr$resultsA
    resultsB <- styled_shr$resultsB
    resultsA$row_kind[ix_shr] <- "shrinkage"
    resultsB$row_kind[ix_shr] <- "shrinkage"
  }
  resultsA$shape <- NA_integer_
  resultsB$shape <- NA_integer_
  resultsA$shape[ix_obs] <- 22L
  resultsB$shape[ix_obs] <- 21L
  if (length(ix_sg_shr)) {
    resultsA$shape[ix_sg_shr] <- 15L
    resultsB$shape[ix_sg_shr] <- NA_integer_
  }
  if (length(ix_int_shr)) {
    resultsA$shape[ix_int_shr] <- NA_integer_
    resultsB$shape[ix_int_shr] <- 16L
  }
  resultsA$point_scale <- 1
  resultsB$point_scale <- 1
  if (length(ix_shr)) {
    shr_scale <- comparison_forest_shrinkage_point_scale()
    resultsA$point_scale[ix_shr] <- shr_scale
    resultsB$point_scale[ix_shr] <- shr_scale
  }
  marginal_names <- setdiff(model_names, "Shrinkage")
  blocks <- effect_bundle$summary_blocks
  if (length(blocks)) {
    if (!"colour_label" %in% names(row_labels)) {
      row_labels$colour_label <- row_labels$label
    }
    extra <- append_comparison_forest_blocks(
      resultsA = resultsA,
      resultsB = resultsB,
      row_labels = row_labels,
      blocks = blocks,
      prob = prob,
      ref = ref,
      contrast_sep = contrast_sep
    )
    resultsA <- extra$resultsA
    resultsB <- extra$resultsB
    row_labels <- extra$row_labels
  }
  styled <- comparison_forest_apply_effect_styles(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    ref = ref,
    contrast_sep = contrast_sep
  )
  resultsA <- styled$resultsA
  resultsB <- styled$resultsB
  if (!is.null(overall_shrinkage)) {
    attr(resultsA, "overall_shrinkage") <- overall_shrinkage
  }
  list(
    resultsA = resultsA,
    resultsB = resultsB,
    row_labels = row_labels,
    n_study = n_obs,
    ref = ref,
    obs_row_index = ix_obs,
    include_shrinkage = comparison_forest_layer_on(include_shrinkage),
    include_prediction = comparison_forest_layer_on(include_prediction),
    include_combined = comparison_forest_layer_on(include_combined),
    overall_shrinkage = overall_shrinkage,
    annotations = effect_bundle$annotations
  )
}

#' Forest plot for DA / AD / SWADA / conventional fits
#' @noRd
construction_forestplot_from_fit <- function(fit,
                                             model_names = NULL,
                                             include_shrinkage = TRUE,
                                             include_overall = TRUE,
                                             include_prediction = TRUE,
                                             include_combined = TRUE,
                                             exponentiate = FALSE,
                                             n_var = FALSE,
                                             group = "subgroup",
                                             show_annotations = TRUE,
                                             annotation_digits = 2,
                                             prob = 0.95,
                                             contrast_sep = "-",
                                             ...) {
  analysis <- interactmeta_fit_analysis(fit)
  if (inherits(fit, "interactmeta_conventional")) {
    analysis <- "conventional"
  }
  dots <- list(...)
  if (identical(analysis, "ad")) {
    type <- "interaction"
    if ("type" %in% names(dots)) {
      type <- forestplot_resolve_type(as.character(dots$type)[[1L]])
      dots$type <- NULL
      if (!identical(type, "interaction")) {
        stop(
          "AD forests support only `type = \"interaction\"`.",
          call. = FALSE
        )
      }
    }
  } else {
    type <- "modifiers"
    if ("type" %in% names(dots)) {
      type <- forestplot_resolve_type(as.character(dots$type)[[1L]])
      dots$type <- NULL
    }
  }
  if ("exponentiate" %in% names(dots)) {
    exponentiate <- dots$exponentiate
    dots$exponentiate <- NULL
  }
  exponentiate <- isTRUE(exponentiate)
  if ("show_heterogeneity_footnote" %in% names(dots)) {
    show_annotations <- isTRUE(dots$show_heterogeneity_footnote)
    dots$show_heterogeneity_footnote <- NULL
  }
  group <- match.arg(group, c("subgroup", "trial"))
  validate_forestplot_args(
    prob = prob,
    which_panel = dots$which_panel,
    xlim = dots$xlim,
    xticks = dots$xticks,
    fill = if (!is.null(dots$fill)) dots$fill else "white",
    n_var = n_var
  )
  model_names <- comparison_forest_resolve_model_names(
    model_names = model_names,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction
  )
  effect_bundle <- construction_forest_effect_bundle(
    fit = fit,
    prob = prob,
    contrast_sep = contrast_sep,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
    show_annotations = show_annotations,
    annotation_digits = annotation_digits,
    model_names = model_names
  )
  render_args <- construction_forest_build_fit_panels(
    fit = fit,
    prob = prob,
    contrast_sep = contrast_sep,
    model_names = model_names,
    include_shrinkage = include_shrinkage,
    include_overall = include_overall,
    include_prediction = include_prediction,
    include_combined = include_combined,
    n_var = n_var,
    group = group,
    effect_bundle = effect_bundle
  )
  annotations <- effect_bundle$annotations
  render_args$annotations <- NULL
  fp <- do.call(
    forestplot_render,
    c(
      render_args,
      list(
        prob = prob,
        contrast_sep = contrast_sep,
        exponentiate = exponentiate,
        annotations = annotations,
        group = group,
        type = if (identical(type, "subgroups")) "modifiers" else type,
        fit_analysis = analysis
      ),
      dots
    )
  )
  attr(fp, "footnote") <- if (length(annotations$lines) ||
    length(annotations$heterogeneity)) {
    paste(
      c(
        annotations$lines,
        unname(annotations$heterogeneity)
      ),
      collapse = "\n"
    )
  } else {
    NULL
  }
  attr(fp, "annotations") <- annotations
  attr(fp, "forestplot_fit_render") <- render_args
  fp
}
