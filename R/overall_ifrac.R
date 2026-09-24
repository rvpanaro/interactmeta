#' @rdname ifrac
#' @export
ifrac.interactmeta <- function(x,
                               ifrac = c("interaction", "overall"),
                               ifrac_var = NULL,
                               summary = FALSE,
                               prob = 0.95,
                               ...) {
  overall_ifrac_from_fit(
    x,
    ifrac_var = ifrac_var,
    ifrac = match.arg(ifrac),
    summary = summary,
    prob = prob
  )
}


#' @rdname ifrac
#' @export
ifrac.interactmeta_overall_ifrac <- function(x, ...) {
  x
}

#' @noRd
overall_ifrac_from_fit <- function(fit,
                                   ifrac_var = NULL,
                                   ifrac = c("interaction", "overall"),
                                   summary = FALSE,
                                   prob = 0.95) {
  ifrac <- match.arg(ifrac)
  check_comparison_fit(fit)
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(st)) {
    stop(
      "Fit has no `\"interactmeta\"` attribute; refit with `cam()`.",
      call. = FALSE
    )
  }

  trial_var <- if (is.null(st$trial_var)) "trial" else st$trial_var
  subgroup_var <- if (is.null(st$subgroup_var)) "subgroup" else st$subgroup_var
  ref <- st$reference
  if (is.null(ref) || !nzchar(ref)) {
    stop("Fit has no `reference` in attribute `\"interactmeta\"`.", call. = FALSE)
  }

  d <- st$data_fit
  if (is.null(d) && !is.null(st$ifrac_data)) {
    d <- st$ifrac_data
    if (is.null(ifrac_var)) {
      ifrac_var <- st$ifrac_var
      if (is.null(ifrac_var) || !ifrac_var %in% names(d)) {
        vars <- detect_ifrac_vars(d)
        ifrac_var <- vars[[1]]
      }
    }
  } else if (is.null(ifrac_var)) {
    ifrac_var <- st$ifrac_vars_std[[1]]
    if (is.null(ifrac_var) && !is.null(st$ifrac_var)) {
      ifrac_var <- st$ifrac_var
    }
    if (is.null(ifrac_var) || !length(ifrac_var) || !ifrac_var %in% names(d)) {
      if (!is.null(d)) {
        vars <- detect_ifrac_vars(d)
        ifrac_var <- vars[[1]]
      }
    }
  }

  if (is.null(d)) {
    stop("Fit has no `data_fit` or `ifrac_data`; refit with `cam()`.", call. = FALSE)
  }

  draws <- if (identical(ifrac, "interaction")) {
    tau_gamma <- comparison_tau_gamma_samples(fit, ifrac_var = ifrac_var, trial_var = trial_var)
    interaction_ifrac_compute(
      data = d,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = ref,
      ifrac_var = ifrac_var,
      tau_gamma = tau_gamma,
      contrast_sep = "-"
    )
  } else {
    tau_hyp <- sprintf("%s__Intercept = 0", trial_var)
    tau <- comparison_sd_samples(fit, tau_hyp)
    overall_ifrac_compute(
      data = d,
      trial_var = trial_var,
      subgroup_var = subgroup_var,
      reference = ref,
      ifrac_var = ifrac_var,
      tau = tau,
      contrast_sep = "-"
    )
  }

  trial_pi <- overall_ifrac_trial_pi(
    d,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = ref,
    ifrac_var = ifrac_var
  )
  bundle <- list(
    draws = draws,
    trial_pi = trial_pi,
    reference = ref,
    ifrac_var = ifrac_var,
    subgroup_label = overall_ifrac_subgroup_label(ifrac_var, interactmeta = st),
    n_trials = nrow(trial_pi),
    summary = overall_ifrac_summary_df(draws, prob = prob, weight = interactmeta_ifrac_draw_weights(fit, draws)),
    prob = prob,
    ifrac = ifrac,
    fit = fit
  )
  obj <- overall_ifrac_result(draws, summary = summary, prob = prob, bundle = bundle)
  if (inherits(obj, "interactmeta_overall_ifrac")) {
    w <- interactmeta_ifrac_draw_weights(fit, draws)
    if (!is.null(w)) {
      attr(obj, ".weight") <- w
    }
    obj <- interactmeta_attach_origin(obj, fit, quantity = "overall_ifrac")
  }
  obj
}

#' @noRd
overall_ifrac_from_data <- function(data,
                                    fit = NULL,
                                    tau = NULL,
                                    tau_gamma = NULL,
                                    ifrac = c("interaction", "overall"),
                                    reference = NULL,
                                    trial = "trial",
                                    subgroup = "subgroup",
                                    ifrac_var = NULL,
                                    contrast_sep = "-",
                                    summary = FALSE,
                                    prob = 0.95) {
  ifrac <- match.arg(ifrac)
  if (!is.null(fit)) {
    check_comparison_fit(fit)
    if (!is_comparison_fit(fit)) {
      stop("`fit` must be an `interactmeta` from `cam()`.", call. = FALSE)
    }
    st <- attr(fit, "interactmeta", exact = TRUE)
    trial_var <- if (!is.null(st) && !is.null(st$trial_var)) st$trial_var else trial
    if (identical(ifrac, "interaction")) {
      if (is.null(ifrac_var)) {
        ifrac_var <- detect_ifrac_vars(data)[[1]]
      }
      tau_gamma <- comparison_tau_gamma_samples(fit, ifrac_var = ifrac_var, trial_var = trial_var)
    } else {
      tau_hyp <- sprintf("%s__Intercept = 0", trial_var)
      tau <- comparison_sd_samples(fit, tau_hyp)
    }
  } else if (identical(ifrac, "interaction")) {
    if (is.null(tau_gamma)) {
      stop(
        "Supply `fit` or `tau_gamma` when `object` is a data frame and `ifrac = \"interaction\"`.",
        call. = FALSE
      )
    }
  } else if (is.null(tau)) {
    stop(
      "Supply `fit` or `tau` when `object` is a data frame and `ifrac = \"overall\"`.",
      call. = FALSE
    )
  }

  if (identical(ifrac, "interaction")) {
    tau_gamma <- as.numeric(tau_gamma)
    if (!length(tau_gamma) || !all(is.finite(tau_gamma))) {
      stop("`tau_gamma` must be one or more finite numeric values.", call. = FALSE)
    }
  } else {
    tau <- as.numeric(tau)
    if (!length(tau) || !all(is.finite(tau))) {
      stop("`tau` must be one or more finite numeric values.", call. = FALSE)
    }
  }

  ref <- reference
  if (is.null(ref)) {
    if ("reference" %in% names(data)) {
      ref <- unique(data$reference)
      ref <- ref[!is.na(ref) & nzchar(ref)]
    }
    if (length(ref) != 1L) {
      stop(
        "Supply `reference` or include a single `reference` column in `object`.",
        call. = FALSE
      )
    }
    ref <- ref[[1]]
  }

  if (is.null(ifrac_var)) {
    ifrac_var <- detect_ifrac_vars(data)[[1]]
  }

  draws <- if (identical(ifrac, "interaction")) {
    interaction_ifrac_compute(
      data = data,
      trial_var = trial,
      subgroup_var = subgroup,
      reference = ref,
      ifrac_var = ifrac_var,
      tau_gamma = tau_gamma,
      contrast_sep = contrast_sep
    )
  } else {
    overall_ifrac_compute(
      data = data,
      trial_var = trial,
      subgroup_var = subgroup,
      reference = ref,
      ifrac_var = ifrac_var,
      tau = tau,
      contrast_sep = contrast_sep
    )
  }

  trial_pi <- overall_ifrac_trial_pi(
    data,
    trial_var = trial,
    subgroup_var = subgroup,
    reference = ref,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  bundle <- list(
    draws = draws,
    trial_pi = trial_pi,
    reference = ref,
    ifrac_var = ifrac_var,
    subgroup_label = overall_ifrac_subgroup_label(ifrac_var),
    n_trials = nrow(trial_pi),
    summary = overall_ifrac_summary_df(draws, prob = prob),
    prob = prob,
    ifrac = ifrac
  )
  overall_ifrac_result(draws, summary = summary, prob = prob, bundle = bundle)
}

#' @noRd
overall_ifrac_prepare_trials <- function(data,
                                         trial_var,
                                         subgroup_var,
                                         reference,
                                         ifrac_var,
                                         contrast_sep = "-") {
  if (!ifrac_var %in% names(data)) {
    stop(sprintf("Column `%s` not found in data.", ifrac_var), call. = FALSE)
  }

  d <- data[
    is_ifrac_component_row(data[[subgroup_var]], contrast_sep = contrast_sep),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) {
    stop("No subgroup rows found for overall IF computation.", call. = FALSE)
  }

  sg <- as.character(d[[subgroup_var]])
  reference <- overall_ifrac_resolve_reference(sg, reference)
  d <- d[order(d[[trial_var]], sg != reference, sg), , drop = FALSE]
  d$se <- sqrt(d$vi)

  trials <- unique(as.character(d[[trial_var]]))
  n_trials <- length(trials)
  n_per <- nrow(d) / n_trials
  if (abs(n_per - round(n_per)) > 1e-8 || round(n_per) != n_per) {
    stop("Each trial must have the same number of subgroup rows.", call. = FALSE)
  }

  pi <- vapply(
    trials,
    function(tr) {
      d[[ifrac_var]][d[[trial_var]] == tr][[1]]
    },
    numeric(1)
  )

  list(data = d, trials = trials, pi = pi, reference = reference)
}

#' Match `reference` to subgroup labels in the IF table
#'
#' Fitted `data_fit` recodes subgroups to `interactmeta_sg_1` (reference)
#' / `interactmeta_sg_2` while `st$reference` stays the original label.
#' @noRd
overall_ifrac_resolve_reference <- function(sg, reference) {
  sg <- unique(as.character(sg))
  ref <- as.character(reference)[[1L]]
  if (length(ref) && ref %in% sg) {
    return(ref)
  }
  canon <- "interactmeta_sg_1"
  if (canon %in% sg) {
    return(canon)
  }
  ref
}

#' @noRd
overall_ifrac_weighted_pool <- function(pi, den, hetero_sd) {
  hetero_sd <- as.numeric(hetero_sd)
  wi <- vapply(hetero_sd, function(h) 1 / (h^2 + den), numeric(length(pi)))
  if (!is.matrix(wi)) {
    wi <- matrix(wi, nrow = length(pi), ncol = length(hetero_sd))
  }
  vapply(
    seq_len(ncol(wi)),
    function(j) {
      w <- wi[, j]
      sum(pi * w) / sum(w)
    },
    numeric(1)
  )
}

#' Trial-level IF coordinates and the two pooling denominators
#'
#' `pi` is the trial information fraction (same under both weighting schemes).
#' Overall denominator is \(\mathrm{Var}(m_j)=1/\sum_g 1/\sigma_{gj}^2\),
#' the overall-block sampling variance
#' \((1-\pi_j)^2\sigma_{Aj}^2+\pi_j^2\sigma_{Bj}^2\).
#' Interaction denominator is \(\sigma_{Aj}^2+\sigma_{Bj}^2=\mathrm{Var}(g_j)\).
#' Weights are \(w_j=1/(\tau^2+v_j)\) at a supplied heterogeneity SD.
#' @noRd
overall_ifrac_trial_components <- function(data,
                                           trial_var,
                                           subgroup_var,
                                           reference,
                                           ifrac_var,
                                           contrast_sep = "-") {
  prep <- overall_ifrac_prepare_trials(
    data = data,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = reference,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  d <- prep$data
  trials <- prep$trials
  pi <- unname(prep$pi)
  reference <- prep$reference

  den_overall <- vapply(
    trials,
    function(tr) {
      ix <- d[[trial_var]] == tr
      vi <- pmax(as.numeric(d$vi[ix]), 1e-12)
      1 / sum(1 / vi)
    },
    numeric(1)
  )
  den_interaction <- vapply(
    trials,
    function(tr) {
      ix <- d[[trial_var]] == tr
      sum(d$vi[ix])
    },
    numeric(1)
  )
  max_vi <- vapply(
    trials,
    function(tr) {
      max(d$vi[d[[trial_var]] == tr], na.rm = TRUE)
    },
    numeric(1)
  )

  data.frame(
    trial = as.character(trials),
    pi = as.numeric(pi),
    den_overall = unname(as.numeric(den_overall)),
    den_interaction = unname(as.numeric(den_interaction)),
    max_vi = unname(as.numeric(max_vi)),
    stringsAsFactors = FALSE
  )
}

#' Inverse-variance IF weights at a heterogeneity SD
#' @noRd
overall_ifrac_weights_at <- function(den, hetero_sd) {
  h <- as.numeric(hetero_sd)[[1L]]
  1 / (h^2 + as.numeric(den))
}

#' @noRd
overall_ifrac_compute <- function(data,
                                  trial_var,
                                  subgroup_var,
                                  reference,
                                  ifrac_var,
                                  tau,
                                  contrast_sep = "-") {
  comp <- overall_ifrac_trial_components(
    data = data,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = reference,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  overall_ifrac_weighted_pool(
    pi = comp$pi,
    den = comp$den_overall,
    hetero_sd = tau
  )
}

#' Interaction-IF weights: inverse of tau_gamma^2 + subgroup sampling variances
#' @noRd
interaction_ifrac_compute <- function(data,
                                      trial_var,
                                      subgroup_var,
                                      reference,
                                      ifrac_var,
                                      tau_gamma,
                                      contrast_sep = "-") {
  comp <- overall_ifrac_trial_components(
    data = data,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = reference,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  overall_ifrac_weighted_pool(
    pi = comp$pi,
    den = comp$den_interaction,
    hetero_sd = tau_gamma
  )
}

#' Posterior draws of interaction heterogeneity (trial ifrac2 SD)
#' @noRd
comparison_tau_gamma_samples <- function(fit, ifrac_var, trial_var = NULL) {
  st <- attr(fit, "interactmeta", exact = TRUE)
  if (is.null(trial_var) || !nzchar(trial_var)) {
    trial_var <- if (!is.null(st) && !is.null(st$trial_var)) st$trial_var else "trial"
  }
  subgroup <- sub("^ifrac_", "", ifrac_var)
  interaction_param <- paste0("Interaction_", subgroup)
  legacy_ifrac2_param <- paste0("ifrac2_", subgroup)
  specs <- comparison_heterogeneity_specs(fit, trial_var = trial_var, component = "slope")
  if (!nrow(specs)) {
    stop(
      "No interaction heterogeneity (`Interaction_*` / `ifrac2_*` slope SD) found in the fit.",
      call. = FALSE
    )
  }
  hit <- specs$parameter %in% c(interaction_param, legacy_ifrac2_param)
  if (!any(hit)) {
    hit <- specs$type == "slope"
    if (!any(hit)) {
      stop(
        sprintf(
          "No slope heterogeneity parameter for `%s`; available: %s.",
          interaction_param,
          paste(specs$parameter, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    ix <- which(hit)[[1L]]
  } else {
    ix <- which(hit)[[1L]]
  }
  comparison_sd_samples(fit, specs$hypothesis[[ix]])
}

#' @noRd
overall_ifrac_weight_label <- function(ifrac = c("interaction", "overall")) {
  ifrac <- match.arg(ifrac)
  if (identical(ifrac, "interaction")) {
    "Interaction-weighted information fraction"
  } else {
    "Overall-weighted information fraction"
  }
}

#' Trial-level information fractions used in the overall IF weights
#' @noRd
overall_ifrac_trial_pi <- function(data,
                                   trial_var,
                                   subgroup_var,
                                   reference,
                                   ifrac_var,
                                   contrast_sep = "-") {
  d <- data[
    is_ifrac_component_row(data[[subgroup_var]], contrast_sep = contrast_sep),
    ,
    drop = FALSE
  ]
  trials <- unique(as.character(d[[trial_var]]))
  data.frame(
    trial = trials,
    pi = vapply(
      trials,
      function(tr) {
        d[[ifrac_var]][d[[trial_var]] == tr][[1]]
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

#' Draws plus metadata for plotting
#' @noRd
overall_ifrac_bundle <- function(object,
                                 fit = NULL,
                                 tau = NULL,
                                 tau_gamma = NULL,
                                 ifrac = c("interaction", "overall"),
                                 reference = NULL,
                                 trial = "trial",
                                 subgroup = "subgroup",
                                 ifrac_var = NULL,
                                 contrast_sep = "-",
                                 prob = 0.95,
                                 ...) {
  ifrac <- match.arg(ifrac)
  if (inherits(object, "interactmeta_overall_ifrac")) {
    return(list(
      draws = as.numeric(object),
      trial_pi = attr(object, "trial_pi", exact = TRUE),
      reference = attr(object, "reference", exact = TRUE),
      ifrac_var = attr(object, "ifrac_var", exact = TRUE),
      subgroup_label = attr(object, "subgroup_label", exact = TRUE),
      n_trials = attr(object, "n_trials", exact = TRUE),
      summary = attr(object, "summary", exact = TRUE),
      prob = attr(object, "prob", exact = TRUE) %||% prob,
      ifrac = attr(object, "ifrac", exact = TRUE) %||% "interaction",
      fit = attr(object, "interactmeta", exact = TRUE)
    ))
  }

  if (is.numeric(object) && !is.data.frame(object) && !inherits(object, "interactmeta_overall_ifrac")) {
    draws <- as.numeric(object)
    if (!length(draws) || !all(is.finite(draws))) {
      stop("`object` must be a finite numeric vector of overall IF draws.", call. = FALSE)
    }
    summ <- overall_ifrac_summary_df(draws, prob = prob)
    return(list(
      draws = draws,
      trial_pi = NULL,
      reference = reference,
      ifrac_var = ifrac_var,
      subgroup_label = overall_ifrac_subgroup_label(ifrac_var),
      n_trials = NA_integer_,
      summary = summ,
      prob = prob,
      ifrac = ifrac
    ))
  }

  if (is_comparison_fit(object)) {
    check_comparison_fit(fit)
    st <- attr(object, "interactmeta", exact = TRUE)
    trial_var <- if (is.null(st$trial_var)) "trial" else st$trial_var
    subgroup_var <- if (is.null(st$subgroup_var)) "subgroup" else st$subgroup_var
    ref <- st$reference
    d <- st$data_fit
    if (is.null(d) && !is.null(st$ifrac_data)) {
      d <- st$ifrac_data
      if (is.null(ifrac_var)) {
        ifrac_var <- st$ifrac_var
        if (is.null(ifrac_var) || !ifrac_var %in% names(d)) {
          ifrac_var <- detect_ifrac_vars(d)[[1]]
        }
      }
    } else if (is.null(ifrac_var)) {
      ifrac_var <- st$ifrac_vars_std[[1]]
      if (is.null(ifrac_var) && !is.null(st$ifrac_var)) {
        ifrac_var <- st$ifrac_var
      }
      if (is.null(ifrac_var) || !ifrac_var %in% names(d)) {
        ifrac_var <- detect_ifrac_vars(d)[[1]]
      }
    }
    draws <- overall_ifrac_from_fit(
      object,
      ifrac_var = ifrac_var,
      ifrac = ifrac,
      summary = FALSE,
      prob = prob
    )
    # Prefer attributes already attached by overall_ifrac_from_fit
    trial_pi <- attr(draws, "trial_pi", exact = TRUE)
    if (is.null(trial_pi)) {
      trial_pi <- overall_ifrac_trial_pi(
        d,
        trial_var = trial_var,
        subgroup_var = subgroup_var,
        reference = ref,
        ifrac_var = ifrac_var
      )
    }
    summ <- attr(draws, "summary", exact = TRUE)
    if (is.null(summ)) {
      summ <- overall_ifrac_summary_df(draws, prob = prob)
    }
    return(list(
      draws = draws,
      trial_pi = trial_pi,
      reference = ref,
      ifrac_var = ifrac_var,
      subgroup_label = overall_ifrac_subgroup_label(ifrac_var, interactmeta = st),
      n_trials = if (is.data.frame(trial_pi)) nrow(trial_pi) else length(trial_pi),
      summary = summ,
      prob = prob,
      ifrac = ifrac,
      fit = object
    ))
  }

  ref <- reference
  if (is.null(ref) && "reference" %in% names(object)) {
    ref <- unique(object$reference)
    ref <- ref[!is.na(ref) & nzchar(ref)]
    if (length(ref) == 1L) {
      ref <- ref[[1]]
    }
  }
  if (is.null(ifrac_var)) {
    ifrac_var <- detect_ifrac_vars(object)[[1]]
  }
  draws <- overall_ifrac_from_data(
    object,
    fit = fit,
    tau = tau,
    tau_gamma = tau_gamma,
    ifrac = ifrac,
    reference = ref,
    trial = trial,
    subgroup = subgroup,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep,
    summary = FALSE,
    prob = prob
  )
  trial_pi <- overall_ifrac_trial_pi(
    object,
    trial_var = trial,
    subgroup_var = subgroup,
    reference = ref,
    ifrac_var = ifrac_var,
    contrast_sep = contrast_sep
  )
  list(
    draws = draws,
    trial_pi = trial_pi,
    reference = ref,
    ifrac_var = ifrac_var,
    subgroup_label = overall_ifrac_subgroup_label(ifrac_var),
    n_trials = nrow(trial_pi),
    summary = overall_ifrac_summary_df(draws, prob = prob),
    prob = prob,
    ifrac = ifrac
  )
}

#' @noRd
overall_ifrac_summary_df <- function(draws, prob = 0.95, weight = NULL) {
  if (is.null(weight)) {
    weight <- attr(draws, ".weight", exact = TRUE)
  }
  mom <- comparison_draw_moments(as.numeric(draws), prob = prob, weight = weight)
  data.frame(
    estimate = mom$estimate,
    lci = mom$lci,
    uci = mom$uci,
    prob = mom$prob,
    stringsAsFactors = FALSE
  )
}

#' Align quadrature node weights with IF draws from the same fit
#' @noRd
interactmeta_ifrac_draw_weights <- function(fit, draws) {
  w <- as.data.frame(fit$draws)$.weight
  if (is.null(w)) {
    return(NULL)
  }
  if (length(w) == length(draws)) {
    return(as.numeric(w))
  }
  NULL
}

#' One-line summary for overall-IF plots (estimate, interval, studies, reference)
#' @noRd
overall_ifrac_footnote <- function(summary,
                                   prob = NULL,
                                   n_trials = NULL,
                                   reference = NULL,
                                   digits = 2,
                                   plugin = FALSE) {
  if (is.null(prob)) {
    prob <- if (!is.null(summary$prob)) summary$prob[[1L]] else 0.95
  }
  fmt <- function(x) {
    format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE)
  }
  est <- fmt(summary$estimate)
  has_int <- !isTRUE(plugin) &&
    "lci" %in% names(summary) &&
    "uci" %in% names(summary) &&
    is.finite(as.numeric(summary$lci)[[1L]]) &&
    is.finite(as.numeric(summary$uci)[[1L]])
  out <- if (has_int) {
    paste0(
      "IF ",
      est,
      " (95% interval: ",
      fmt(summary$lci),
      ", ",
      fmt(summary$uci),
      ")"
    )
  } else {
    paste0("IF ", est, " (plug-in at fitted tau)")
  }
  if (!is.null(n_trials) && !is.na(n_trials)) {
    out <- paste0(out, " \u00b7 k = ", n_trials, " studies")
  }
  if (!is.null(reference) && nzchar(as.character(reference))) {
    out <- paste0(out, " \u00b7 reference: ", reference)
  }
  out
}

#' @noRd
overall_ifrac_footnote_bundle <- function(bundle, digits = 2) {
  plugin <- FALSE
  if (!is.null(bundle$fit)) {
    plugin <- isTRUE(interactmeta_ifrac_is_plugin(bundle$fit))
  }
  overall_ifrac_footnote(
    summary = bundle$summary,
    prob = bundle$prob %||% 0.95,
    n_trials = bundle$n_trials,
    reference = bundle$reference,
    digits = digits,
    plugin = plugin
  )
}

#' @noRd
new_interactmeta_overall_ifrac <- function(draws, bundle, prob = 0.95) {
  obj <- draws
  class(obj) <- c("interactmeta_overall_ifrac", "numeric")
  attr(obj, "trial_pi") <- bundle$trial_pi
  attr(obj, "summary") <- bundle$summary
  attr(obj, "reference") <- bundle$reference
  attr(obj, "ifrac_var") <- bundle$ifrac_var
  attr(obj, "subgroup_label") <- bundle$subgroup_label
  attr(obj, "n_trials") <- bundle$n_trials
  attr(obj, "prob") <- prob
  attr(obj, "ifrac") <- bundle$ifrac %||% "interaction"
  attr(obj, "interactmeta") <- bundle$fit
  obj
}

#' @noRd
overall_ifrac_result <- function(draws, summary = FALSE, prob = 0.95, bundle = NULL) {
  if (summary) {
    return(overall_ifrac_summary_df(draws, prob = prob))
  }
  if (is.null(bundle)) {
    bundle <- list(
      draws = draws,
      trial_pi = NULL,
      reference = NULL,
      ifrac_var = NULL,
      subgroup_label = NULL,
      n_trials = NA_integer_,
      summary = overall_ifrac_summary_df(draws, prob = prob),
      prob = prob
    )
  }
  new_interactmeta_overall_ifrac(draws, bundle = bundle, prob = prob)
}

#' Print overall information fraction summary
#'
#' @param x A `interactmeta_overall_ifrac` object.
#' @param digits Decimal places (default `2`).
#' @param ... Ignored.
#'
#' @export
print.interactmeta_overall_ifrac <- function(x, digits = 2, ...) {
  summ <- attr(x, "summary", exact = TRUE)
  prob <- attr(x, "prob", exact = TRUE) %||% 0.95
  if (is.null(summ)) {
    summ <- overall_ifrac_summary_df(x, prob = prob)
  }
  scheme <- attr(x, "ifrac", exact = TRUE) %||% "interaction"
  style <- interactmeta_interval_style(x)
  plugin_fit <- identical(style, "frequentist")
  cat(overall_ifrac_weight_label(scheme), "\n", sep = "")
  comparison_print_effects_table(
    summ,
    digits = digits,
    prob = prob,
    id_cols = character(),
    title = NULL,
    footnote = !plugin_fit,
    style = style
  )
  if (plugin_fit) {
    tau_lab <- if (identical(scheme, "interaction")) {
      "plug-in at fitted interaction-block tau_gamma"
    } else {
      "plug-in at fitted overall-block tau"
    }
    cat("(", tau_lab, ")\n", sep = "")
  }
  invisible(x)
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
