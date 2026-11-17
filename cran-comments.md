## interactmeta 0.1.0

### CRAN submission

* First submission.

### Test environments

* Local: macOS, R 4.4.x — `R CMD check` on source tarball (0 errors, 0 warnings).
* GitHub Actions: `R-CMD-check` on ubuntu-latest (error-on: warning).
* Coverage: `covr::codecov()` with `INTERACT_SLOW_TESTS=true` (graphics draw tests).

### Notes

* Default [`cams()`] fit is Bayesian (Suggests/Imports **bayesmeta**); frequentist path uses **metafor**.
* Vignettes knit four articles; forest-plot chunks can be slow; full check run exceeded
  typical notebook cache on maintainer machine — allow vignette build time on Linux builders.

### Downstream dependencies

* None known at first release.