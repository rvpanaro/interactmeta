# ---- sensitivity.R ----
#' Fit multiple meta-analysis constructions on the same data
#'
#' Runs **DA**, **AD**, **SWADA**, and/or **CAMS** on one prepared dataset
#' (typically after [ifrac()]). Use [compare()] on the result for side-by-side
#' tables.
#'
#' @param data_or_fit A data frame with `yi`, `vi`, trial and subgroup columns,
#'   or an existing [`interactmeta`] from [cams()] (data and reference are read
#'   from the fit).
#' @param reference Reference subgroup level. Required for a data frame; inferred
#'   from `data_or_fit` when it is a CAMS fit.
#' @param constructions Character vector of constructions to fit (`"da"`,
#'   `"ad"`, `"swada"`, `"cams"`). Aliased as `methods` in some docs; do not
#'   confuse with `method = "REML"` passed via `...` to fitters.
#' @param trial Trial id column (default `"study"`).
#' @param subgroup Modifier level column (default `"level"`).
#' @param ... Passed to each construction fitter ([da()], [ad()], [swada()],
#'   [cams()]).
#' @return A list with class `interactmeta_sensitivity` containing `$fits`.
#' @export
#' @seealso [compare()], [compare_interactions()], [cams()]
sensitivity <- function(data_or_fit,
                        reference = NULL,
                        constructions = c("da", "ad", "swada", "cams"),
                        trial = "study",
                        subgroup = "level",
                        ...) {
  constructions <- match.arg(
    constructions,
    c("da", "ad", "swada", "cams"),
    several.ok = TRUE
  )
  fit_args <- list(...)
  method_arg <- fit_args$method
  if (!is.null(method_arg)) {
    fit_args$method <- NULL
  }
  silent_arg <- fit_args$silent
  if (!is.null(silent_arg)) {
    fit_args$silent <- NULL
  }
  if (inherits(data_or_fit, "interactmeta")) {
    st <- attr(data_or_fit, "interactmeta", exact = TRUE)
    if (!is.list(st)) {
      stop("CAMS fit is missing the `interactmeta` attribute.", call. = FALSE)
    }
    data <- st$ifrac_data
    reference <- reference %||% st$reference
    trial <- st$trial_var %||% trial
    subgroup <- st$subgroup_var %||% subgroup
    if (is.null(reference)) {
      stop("Could not infer `reference` from the fit.", call. = FALSE)
    }
  } else {
    data <- data_or_fit
    if (is.null(reference)) {
      stop("`reference` is required when `data_or_fit` is a data frame.", call. = FALSE)
    }
  }
  if (!is.data.frame(data)) {
    stop("`data_or_fit` must be a data frame or CAMS fit.", call. = FALSE)
  }
  da_args <- c(
    list(
      data = data,
      reference = reference,
      trial = trial,
      subgroup = subgroup,
      method = method_arg %||% "Bayesian"
    ),
    fit_args
  )
  fits <- list()
  if ("da" %in% constructions) {
    fits$da <- do.call(da, da_args)
  }
  if ("ad" %in% constructions) {
    fits$ad <- do.call(ad, da_args)
  }
  if ("swada" %in% constructions) {
    fits$swada <- do.call(swada, da_args)
  }
  if ("cams" %in% constructions) {
    cam_args <- c(
      list(
        x = data,
        reference = reference,
        trial_var = trial,
        subgroup_var = subgroup,
        method = method_arg %||% "Bayesian",
        silent = silent_arg %||% 2
      ),
      fit_args
    )
    fits$cams <- do.call(cams, cam_args)
  }
  structure(
    list(
      data = data,
      reference = reference,
      trial = trial,
      subgroup = subgroup,
      constructions = constructions,
      fits = fits,
      fit_args = fit_args
    ),
    class = c("interactmeta_sensitivity", "list")
  )
}

#' Compare estimands across a sensitivity analysis or multiple fits
#'
#' @param x An `interactmeta_sensitivity` object, or two or more fitted objects
#'   passed via `...` (delegates to legacy compare helpers).
#' @param estimand `"interaction"`, `"modifiers"`, or `"combined"`.
#' @param ... For `interactmeta_sensitivity`, passed to the underlying compare
#'   builder (`exponentiate`, `type`, …). For multiple fits, forwarded to the
#'   internal interaction compare builder.
#' @export
compare <- function(x, ...) {
  UseMethod("compare")
}

#' @rdname compare
#' @export
compare.default <- function(x, y = NULL, ..., exponentiate = FALSE, type = "Overall") {
  dots <- list(...)
  extra_fits <- Filter(function(o) is.object(o), c(list(y), dots))
  extra_fits <- extra_fits[!vapply(extra_fits, is.null, logical(1L))]
  fits <- do.call(interactmeta_compact_fits, c(list(x), extra_fits))
  if (length(fits) < 2L) {
    stop(
      "`compare()` needs at least two fits or a `sensitivity()` object.",
      call. = FALSE
    )
  }
  if (interactmeta_all_cams(fits)) {
    return(interactmeta_compare_interactions_fits(
      fits,
      exponentiate = exponentiate,
      type = type
    ))
  }
  do.call(
    compare_interactions,
    c(fits, list(exponentiate = exponentiate, type = type))
  )
}

#' @rdname compare
#' @param estimand Estimand to tabulate when `x` is an `interactmeta_sensitivity`
#'   object.
#' @export
compare.interactmeta_sensitivity <- function(x,
                                             estimand = c(
                                               "interaction",
                                               "modifiers",
                                               "combined"
                                             ),
                                             ...) {
  estimand <- match.arg(estimand)
  f <- x$fits
  switch(
    estimand,
    interaction = {
      need <- c("da", "ad", "swada", "cams")
      miss <- setdiff(need, names(f))
      if (length(miss)) {
        stop(
          "Interaction sensitivity needs DA, AD, SWADA, and CAMS fits; missing: ",
          paste(miss, collapse = ", "),
          call. = FALSE
        )
      }
      compare_interactions(f$da, f$ad, f$swada, f$cams, ...)
    },
    modifiers = {
      need <- c("da", "swada", "cams")
      miss <- setdiff(need, names(f))
      if (length(miss)) {
        stop(
          "Modifier sensitivity needs DA, SWADA, and CAMS fits; missing: ",
          paste(miss, collapse = ", "),
          call. = FALSE
        )
      }
      compare_modifiers(f$da, f$swada, cam = f$cams, ...)
    },
    combined = {
      cams_fits <- Filter(function(x) inherits(x, "interactmeta"), f)
      if (length(cams_fits) < 2L) {
        stop(
          "`estimand = \"combined\"` needs two or more CAMS fits; ",
          "use `tidy(combined(fits$cams))` for a single CAMS export row.",
          call. = FALSE
        )
      }
      do.call(compare_combined, c(cams_fits, list(...)))
    }
  )
}

#' @rdname compare
#' @export
compare.interactmeta <- function(x,
                                 y = NULL,
                                 ...,
                                 exponentiate = FALSE,
                                 type = "Overall",
                                 estimand = "interaction") {
  dots <- list(...)
  extra_fits <- Filter(function(o) is.object(o), c(list(y), dots))
  extra_fits <- extra_fits[!vapply(extra_fits, is.null, logical(1L))]
  if (length(extra_fits)) {
    fits <- do.call(interactmeta_compact_fits, c(list(x), extra_fits))
    if (length(fits) < 2L) {
      stop(
        "`compare()` needs at least two fits or a `sensitivity()` object.",
        call. = FALSE
      )
    }
    if (interactmeta_all_cams(fits)) {
      return(interactmeta_compare_interactions_fits(
        fits,
        exponentiate = exponentiate,
        type = type
      ))
    }
    return(do.call(
      compare_interactions,
      c(fits, list(exponentiate = exponentiate, type = type))
    ))
  }
  stop(
    "Pass a `sensitivity()` object to `compare()`, e.g. ",
    "`compare(sensitivity(data, reference = ...), estimand = \"interaction\")`.",
    call. = FALSE
  )
}
