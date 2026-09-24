#' Re-exported helpers from dependencies
#'
#' [metafor::escalc()] builds effect-size / sampling-variance data for CAMS.
#' [bayesmeta::traceplot()] draws Rubin-style shrinkage-vs-τ diagnostics
#' (including the [traceplot.interactmeta] method).
#' [generics::tidy()] turns CAMS effect / heterogeneity summaries into tibbles
#' (see [tidy.interactmeta_subgroup_effects] and related methods). Tibbles print with
#' the usual tibble formatting; classed `interactmeta_*` objects use a compact
#' biostatistical `print()` instead.
#'
#' @name interactmeta-reexports
#' @aliases escalc tidy traceplot
#' @importFrom metafor escalc
#' @export escalc
#' @importFrom generics tidy
#' @export tidy
#' @importFrom bayesmeta traceplot
#' @export traceplot
NULL
