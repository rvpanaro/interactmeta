#' Assemble a data frame from yi / vi / subgroup / trial
#' @noRd
interactmeta_build_input <- function(data,
                                     yi_expr,
                                     vi_expr,
                                     subgroup_expr,
                                     trial_expr) {
  eval_one <- function(expr, default_name = NULL) {
    if (is.null(expr)) {
      if (!is.null(default_name) && default_name %in% names(data)) {
        return(data[[default_name]])
      }
      return(NULL)
    }
    if (is.character(expr) && length(expr) == 1L) {
      if (expr %in% names(data)) {
        return(data[[expr]])
      }
      stop("Column `", expr, "` not found in `data`.", call. = FALSE)
    }
    if (is.symbol(expr)) {
      nm <- as.character(expr)
      if (!is.null(data) && nm %in% names(data)) {
        return(data[[nm]])
      }
    }
    if (is.null(data)) {
      return(eval(expr, envir = parent.frame(2L)))
    }
    eval(expr, envir = data, enclos = parent.frame(2L))
  }

  if (missing(data) || is.null(data)) {
    yi_v <- eval_one(yi_expr)
    vi_v <- eval_one(vi_expr)
    sg_v <- eval_one(subgroup_expr)
    tr_v <- eval_one(trial_expr)
    if (is.null(tr_v) || (is.character(tr_v) && length(tr_v) == 1L)) {
      stop("`trial` must be a vector when `data` is omitted.", call. = FALSE)
    }
    return(data.frame(
      trial = tr_v,
      subgroup = sg_v,
      yi = as.numeric(yi_v),
      vi = as.numeric(vi_v),
      stringsAsFactors = FALSE
    ))
  }

  data <- as.data.frame(data)
  yi_v <- eval_one(yi_expr, "yi")
  vi_v <- eval_one(vi_expr, "vi")
  sg_v <- eval_one(subgroup_expr, "subgroup")
  tr_v <- eval_one(trial_expr, "trial")
  if (is.null(yi_v) || is.null(vi_v) || is.null(sg_v) || is.null(tr_v)) {
    stop(
      "`cam()` needs `yi`, `vi`, `subgroup`, and `trial` ",
      "(vectors or column names in `data`).",
      call. = FALSE
    )
  }
  data.frame(
    trial = tr_v,
    subgroup = sg_v,
    yi = as.numeric(yi_v),
    vi = as.numeric(vi_v),
    stringsAsFactors = FALSE
  )
}

#' Bayesian CAMS fit
#'
#' Convenience wrapper for [cams()] with Bayesian inference
#' (`method = "bayes"`), implemented using bayesmeta/DIRECT.
#' Frequentist CAMS uses `cams(..., engine = "metafor")`.
#'
#' @inheritParams cams
#' @param engine Approach for Bayesian CAMS. The public implementation is
#'   `"bayesmeta"`.
#' @param ... Further arguments passed to [cams()] (e.g. `tau.prior`,
#'   `engine_args`, `ecological_slope`, `silent`, `seed`).
#'   `method` is forced to `"bayes"`.
#' @return A fit of class `"interactmeta"`.
#'
#' @examples
#' \dontrun{
#' library(interactmeta)
#' es <- escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
#' fit <- bayescams(es, reference = "Younger")
#' }
#'
#' @seealso [cams()], [forestplot()]
#' @name bayescams
#' @export
bayescams <- function(x,
                     reference = NULL,
                     modifier = NULL,
                     trial = NULL,
                     subgroup = NULL,
                     trial_var = "trial",
                     subgroup_var = "subgroup",
                     contrast = "all",
                     contrast_sep = "-",
                     big_var = 100^2,
                     include_overall = FALSE,
                     include_interaction = FALSE,
                     ecological_slope = TRUE,
                     engine = "bayesmeta",
                     silent = 2,
                     seed = getOption("interactmeta.seed", 42L),
                     ...) {
  dots <- list(...)
  if ("method" %in% names(dots)) {
    warning(
      "`method` is ignored by `bayescams()`; Bayesian CAMS is always used. ",
      "For frequentist CAMS, call `cams(..., engine = \"metafor\")`.",
      call. = FALSE
    )
    dots$method <- NULL
  }
  do.call(
    cam,
    c(
      list(
        x = x,
        reference = reference,
        modifier = modifier,
        trial = trial,
        subgroup = subgroup,
        trial_var = trial_var,
        subgroup_var = subgroup_var,
        contrast = contrast,
        contrast_sep = contrast_sep,
        big_var = big_var,
        include_overall = include_overall,
        include_interaction = include_interaction,
        method = "bayes",
        ecological_slope = ecological_slope,
        engine = engine,
        silent = silent,
        seed = seed
      ),
      dots
    )
  )
}

#' @rdname bayescams
#' @export
bayescam <- bayescams
