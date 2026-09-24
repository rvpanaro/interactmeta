#' Observed study-level treatment-effect estimates and variances
#'
#' Stable extractor for the observed estimates used in CAMS post-fit analysis.
#' Prefer this over reading undocumented fields on `object$data`.
#'
#' @param object A [`interactmeta`] from [cam()].
#' @param type Estimand:
#'   - `"subgroups"` (default; `"subgroup"` accepted): one entry per
#'     trial × subgroup row
#'   - `"overall"` or `"combined"`: trial-specific combined-population
#'     effect \eqn{m_j} (inverse-variance combination within each trial)
#'   - `"interaction"`: trial-specific treatment-by-subgroup interaction
#'     \eqn{g_j} (non-reference minus reference within each trial)
#' @param contrast_sep Separator for interaction labels (default `"-"`).
#'
#' @return A list with:
#'   - `estimate`: named numeric vector
#'   - `variance`: named numeric vector (sampling variances)
#'   - `type`: resolved type
#'   - `labels`: character vector of names (same order)
#'
#' For `type = "subgroups"`, names are `"<trial> / <subgroup>"`. For
#' `"overall"` and `"interaction"`, names are trial IDs.
#'
#' @seealso [shrinkage()], [fitted.interactmeta()]
#' @examples
#' \dontrun{
#' observed <- observed_effects(fit, type = "overall")
#' observed$estimate["RESPIRE1"]
#' observed$variance["RESPIRE1"]
#' }
#' @export
observed_effects <- function(object,
                             type = c("subgroups", "overall", "interaction"),
                             contrast_sep = "-") {
  check_comparison_fit(object)
  type <- comparison_resolve_effect_type(type)
  dat <- comparison_predict_data(object, newdata = NULL)
  st <- attr(object, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% "trial"
  subgroup_var <- st$subgroup_var %||% "subgroup"
  ref <- as.character(st$reference)[[1L]]

  if (identical(type, "subgroups")) {
    labs <- paste(
      as.character(dat[[trial_var]]),
      as.character(dat[[subgroup_var]]),
      sep = " / "
    )
    if (anyDuplicated(labs)) {
      labs <- make.unique(labs, sep = " / ")
    }
    est <- stats::setNames(as.numeric(dat$yi), labs)
    var <- stats::setNames(pmax(as.numeric(dat$vi), 0), labs)
    return(list(
      estimate = est,
      variance = var,
      type = type,
      labels = labs
    ))
  }

  trials <- unique(as.character(dat[[trial_var]]))
  comparison_assert_unique_names(trials, what = "trial")

  nonrefs <- setdiff(unique(as.character(dat[[subgroup_var]])), ref)
  blocks <- comparison_trial_effect_blocks(
    data = dat,
    trial_var = trial_var,
    subgroup_var = subgroup_var,
    reference = ref,
    nonref = if (length(nonrefs)) nonrefs[[1L]] else NULL
  )

  if (identical(type, "overall")) {
    est <- stats::setNames(blocks$overall$yi, blocks$overall$trial)
    var <- stats::setNames(blocks$overall$vi, blocks$overall$trial)
    return(list(
      estimate = est[trials],
      variance = var[trials],
      type = type,
      labels = trials
    ))
  }

  # interaction
  est <- stats::setNames(blocks$interaction$yi, blocks$interaction$trial)
  var <- stats::setNames(blocks$interaction$vi, blocks$interaction$trial)
  list(
    estimate = est[trials],
    variance = var[trials],
    type = type,
    labels = trials
  )
}

#' Resolve estimand `type`, accepting `combined` as an alias of `overall`
#'
#' Internals keep `"overall"` for trial-level combined-population quantities
#' (`m_j`) and forest/shrinkage blocks. `"combined"` is the public synonym.
#' @noRd
comparison_resolve_effect_type <- function(type) {
  type <- as.character(type)[[1L]]
  if (identical(type, "subgroup")) {
    type <- "subgroups"
  }
  if (type %in% c("combined", "combined_population")) {
    type <- "overall"
  }
  match.arg(type, c("subgroups", "overall", "interaction"))
}

#' @noRd
comparison_predict_data <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    dat <- object$data
  } else {
    dat <- newdata
  }
  if (is.null(dat) || !nrow(dat)) {
    stop("No data available for prediction.", call. = FALSE)
  }
  need <- c("yi", "vi")
  miss <- setdiff(need, names(dat))
  if (length(miss)) {
    stop(
      "Prediction data must include columns: ",
      paste(need, collapse = ", "),
      ". Missing: ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  st <- attr(object, "interactmeta", exact = TRUE)
  trial_var <- st$trial_var %||% "trial"
  subgroup_var <- st$subgroup_var %||% "subgroup"
  if (!trial_var %in% names(dat) || !subgroup_var %in% names(dat)) {
    stop(
      "Prediction data must include `",
      trial_var,
      "` and `",
      subgroup_var,
      "`.",
      call. = FALSE
    )
  }
  dat
}

#' @noRd
comparison_assert_unique_names <- function(x, what = "name") {
  x <- as.character(x)
  if (!length(x)) {
    stop("No ", what, " names found.", call. = FALSE)
  }
  if (anyNA(x) || any(!nzchar(x))) {
    stop(what, " names must be non-missing and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop(what, " names must be unique.", call. = FALSE)
  }
  invisible(x)
}

#' Draw N(0, sd_i) matrices with optional seed + stream offset
#' @noRd
comparison_predict_rnorm_matrix <- function(n_draw,
                                      n_col,
                                      sd,
                                      seed = NULL,
                                      stream = "default") {
  sd <- as.numeric(sd)
  if (length(sd) == 1L) {
    sd <- rep(sd, n_draw)
  }
  if (length(sd) != n_draw) {
    stop("Length of `sd` must equal number of draws.", call. = FALSE)
  }
  gen <- function() {
    z <- matrix(stats::rnorm(n_draw * n_col), nrow = n_draw, ncol = n_col)
    z * sd
  }
  if (is.null(seed)) {
    return(gen())
  }
  # Distinct streams for different noise sources under the same user seed
  stream_id <- sum(utf8ToInt(as.character(stream))) %% 10000L
  comparison_with_seed(as.integer(seed) + as.integer(stream_id), gen())
}
