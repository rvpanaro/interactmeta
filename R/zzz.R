#' @importFrom stats setNames
#' @importFrom paletteer paletteer_d
#' @importFrom scales alpha
NULL

.onLoad <- function(libname, pkgname) {
  # bayesmeta registers S3 methods forestplot.escalc / forestplot.bayesmeta /
  # forestplot.bmr when its namespace is loaded. Keep them available for the
  # session once interactmeta is loaded (Imports already ensures the package is
  # installed).
  loadNamespace("bayesmeta")
  loadNamespace("forestplot")
  loadNamespace("metafor")
  invisible()
}

.onAttach <- function(libname, pkgname) {
  invisible()
}
