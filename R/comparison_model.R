comparison_user_cache_dir <- function() {
  if (getRversion() >= "4.0") {
    tools::R_user_dir("interactmeta", which = "cache")
  } else {
    file.path(path.expand("~"), ".cache", "interactmeta")
  }
}

# Alias kept for any lingering call sites
comparison_user_precompile_dir <- comparison_user_cache_dir

comparison_nonref_from_cols <- function(ifrac_cols) {
  sort(sub("^ifrac_", "", ifrac_cols))
}

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
                               trial_var = "trial",
                               subgroup_var = "subgroup",
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
                            trial_var = "trial",
                            subgroup_var = "subgroup") {
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

comparison_dummy_data <- function(n_nonref, n_trials = 2L) {
  n_sg <- n_nonref + 1L
  sg_levels <- comparison_subgroup_fit_levels(n_sg)
  ifrac_std <- comparison_std_ifrac_names(n_nonref)
  ifrac2_std <- comparison_std_ifrac2_names(n_nonref)

  rows <- lapply(seq_len(n_trials), function(i) {
    trial <- paste0("trial_", i)
    out <- data.frame(
      trial = trial,
      subgroup = sg_levels,
      yi = seq_along(sg_levels) * 0.1,
      vi = rep(0.05, n_sg),
      stringsAsFactors = FALSE
    )
    for (j in seq_len(n_nonref)) {
      frac <- j / (n_sg + 0.5)
      out[[ifrac_std[j]]] <- frac
      out[[ifrac2_std[j]]] <- as.numeric(out$subgroup == sg_levels[j + 1L]) - frac
    }
    out$subgroup <- factor(out$subgroup, levels = sg_levels)
    out
  })
  dplyr::bind_rows(rows)
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

ifrac_var_std <- function(ifrac_var, coef_map) {
  if (is.null(coef_map) || !length(coef_map)) {
    return(ifrac_var)
  }
  if (ifrac_var %in% coef_map$ifrac_std) {
    return(ifrac_var)
  }
  hit <- coef_map$ifrac_std[coef_map$ifrac_orig == ifrac_var]
  if (length(hit) == 1L) {
    return(unname(hit))
  }
  ifrac_var
}

ifrac_var_orig <- function(ifrac_var, coef_map) {
  if (is.null(coef_map) || !length(coef_map)) {
    return(ifrac_var)
  }
  if (ifrac_var %in% coef_map$ifrac_orig) {
    return(ifrac_var)
  }
  hit <- coef_map$ifrac_orig[coef_map$ifrac_std == ifrac_var]
  if (length(hit) == 1L) {
    return(unname(hit))
  }
  ifrac_var
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
  if (length(replacements) && is.data.frame(fit$draws)) {
    names(fit$draws) <- comparison_replace_fixed_strings(names(fit$draws), replacements)
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
