#' Observation-level CAMS linear predictor (draws × observations)
#'
#' Single implementation of the Stan mean structure:
#' \deqn{\mu = \alpha + u_j + \beta\cdot\mathrm{ifrac} + v_j\cdot\mathrm{ifrac2}
#'   + \gamma\,1(\mathrm{nonref}).}
#'
#' @param post Draws data frame with `b_Intercept`, `b_ifrac_*`, `b_subgroup*`,
#'   and optional `r_trial[...]` / `sd_trial__*` columns.
#' @param dat Model / prediction data with trial, subgroup, and IF columns.
#' @param trial_var,subgroup_var,reference Column / label metadata.
#' @param ifrac_cols Character vector of `ifrac_*` column names on `dat` /
#'   coefficient names (`b_<ifrac_col>`).
#' @param re One of `"stored"` (use `r_trial` when present), `"zero"` (FE only),
#'   or `"population"` (draw new trial REs from `sd_trial__*`).
#' @param seed,stream_prefix RNG controls for `re = "population"`.
#' @return Numeric matrix (draws × observations).
#' @noRd
comparison_linpred_matrix <- function(post,
                                dat,
                                trial_var = "trial",
                                subgroup_var = "subgroup",
                                reference,
                                ifrac_cols,
                                re = c("stored", "zero", "population"),
                                seed = NULL,
                                stream_prefix = "linpred") {
  re <- match.arg(re)
  ref <- as.character(reference)[[1L]]
  n_draw <- nrow(post)
  n_obs <- nrow(dat)
  if (!n_obs) {
    return(matrix(numeric(0), nrow = n_draw, ncol = 0L))
  }
  if (!"b_Intercept" %in% names(post)) {
    stop("Fit draws are missing `b_Intercept`.", call. = FALSE)
  }
  ifrac_cols <- as.character(ifrac_cols)
  trials <- as.character(dat[[trial_var]])
  subgroups <- as.character(dat[[subgroup_var]])
  unique_trials <- unique(trials)

  re_int <- matrix(0, nrow = n_draw, ncol = length(unique_trials))
  colnames(re_int) <- unique_trials
  re_slope <- list()
  for (ifc in ifrac_cols) {
    lbl <- sub("^ifrac_", "", ifc)
    re_slope[[lbl]] <- matrix(0, nrow = n_draw, ncol = length(unique_trials))
    colnames(re_slope[[lbl]]) <- unique_trials
  }

  if (identical(re, "population")) {
    tau <- if ("sd_trial__Intercept" %in% names(post)) {
      pmax(as.numeric(post$sd_trial__Intercept), 0)
    } else {
      rep(0, n_draw)
    }
    re_int[] <- comparison_predict_rnorm_matrix(
      n_draw,
      length(unique_trials),
      sd = tau,
      seed = seed,
      stream = paste0(stream_prefix, "_re_intercept")
    )
    for (ifc in ifrac_cols) {
      lbl <- sub("^ifrac_", "", ifc)
      sd_col <- paste0("sd_trial__ifrac2_", lbl)
      tau_g <- if (sd_col %in% names(post)) {
        pmax(as.numeric(post[[sd_col]]), 0)
      } else {
        rep(0, n_draw)
      }
      re_slope[[lbl]][] <- comparison_predict_rnorm_matrix(
        n_draw,
        length(unique_trials),
        sd = tau_g,
        seed = seed,
        stream = paste0(stream_prefix, "_re_ifrac_", lbl)
      )
    }
  } else if (identical(re, "stored")) {
    for (tr in unique_trials) {
      col <- paste0("r_trial[", tr, ",Intercept]")
      if (col %in% names(post)) {
        re_int[, tr] <- as.numeric(post[[col]])
      }
    }
    for (ifc in ifrac_cols) {
      lbl <- sub("^ifrac_", "", ifc)
      for (tr in unique_trials) {
        col <- paste0("r_trial[", tr, ",ifrac2_", lbl, "]")
        if (col %in% names(post)) {
          re_slope[[lbl]][, tr] <- as.numeric(post[[col]])
        }
      }
    }
  }

  out <- matrix(0, nrow = n_draw, ncol = n_obs)
  intercept <- as.numeric(post$b_Intercept)
  for (j in seq_len(n_obs)) {
    tr <- trials[[j]]
    sg <- subgroups[[j]]
    lp <- intercept + re_int[, tr]
    if (!identical(sg, ref)) {
      sg_col <- paste0("b_subgroup", sg)
      if (!sg_col %in% names(post)) {
        stop("Missing subgroup FE column `", sg_col, "`.", call. = FALSE)
      }
      lp <- lp + as.numeric(post[[sg_col]])
    }
    for (ifc in ifrac_cols) {
      lbl <- sub("^ifrac_", "", ifc)
      xcol <- if (ifc %in% names(dat)) {
        ifc
      } else if (paste0("ifrac_", lbl) %in% names(dat)) {
        paste0("ifrac_", lbl)
      } else {
        stop("Missing IF column `", ifc, "` in prediction data.", call. = FALSE)
      }
      bcol <- paste0("b_", ifc)
      if (!bcol %in% names(post)) {
        # bayesmeta standardized: b_ifrac_1 while column may be ifrac_Older
        bcol_alt <- paste0("b_ifrac_", lbl)
        if (bcol_alt %in% names(post)) {
          bcol <- bcol_alt
        } else if ("b_ifrac_1" %in% names(post) && length(ifrac_cols) == 1L) {
          bcol <- "b_ifrac_1"
        } else {
          stop("Missing IF slope `", paste0("b_", ifc), "`.", call. = FALSE)
        }
      }
      x <- as.numeric(dat[[xcol]][[j]])
      lp <- lp + as.numeric(post[[bcol]]) * x
      x2_name <- paste0("ifrac2_", lbl)
      if (x2_name %in% names(dat) && !is.null(re_slope[[lbl]])) {
        x2 <- as.numeric(dat[[x2_name]][[j]])
        lp <- lp + re_slope[[lbl]][, tr] * x2
      }
    }
    out[, j] <- lp
  }
  out
}
