# Cached fits for the slow coverage sweep (one object per engine / pipeline).

fixture_coverage_suite <- function() {
  key <- "coverage_suite_v1"
  if (exists(key, envir = .fixture_cache, inherits = FALSE)) {
    return(get(key, envir = .fixture_cache))
  }
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("forestplot")

  sens <- sensitivity(
    respire_age,
    reference = "<65",
    method = "REML",
    silent = 2
  )
  bayes_tiny <- fixture_cam_tiny_fit()
  reml_respire <- suppressWarnings(
    cams(respire_age, reference = "<65", method = "REML", silent = 2)
  )
  reml_ml <- suppressWarnings(
    cams(respire_age, reference = "<65", method = "ML", silent = 2)
  )

  suite <- list(
    bayes_tiny = bayes_tiny,
    reml_respire = reml_respire,
    reml_ml = reml_ml,
    sens = sens,
    da = sens$fits$da,
    ad = sens$fits$ad,
    swada = sens$fits$swada,
    cams = sens$fits$cams
  )
  assign(key, suite, envir = .fixture_cache)
  suite
}