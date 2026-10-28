
<!-- README.md is generated from README.Rmd. Please edit that file -->

# interactmeta

<img src="man/figures/logo.png" align="right" width="139" alt="interactmeta logo" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/rvpanaro/interactmeta/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rvpanaro/interactmeta/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/rvpanaro/interactmeta/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/rvpanaro/interactmeta/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

**interactmeta** implements **contribution-adjusted meta-analysis
(CAMS)** for treatment effect modification from modifier-level aggregate
data (Panaro, Röver & Friede 2025): [*Consistent Bayesian meta-analysis
on subgroup specific effects and
interactions*](https://arxiv.org/abs/2512.18785).

`cams()` fits contribution-adjusted meta-analysis (CAMS) from
modifier-level aggregate data. `interactions()`, `modifiers()`, and
`combined()` report the main estimands; `forestplot()` and `plot()` and
`tidy()` visualize and export a fitted model.

| Start here | More |
|----|----|
| `vignette("cams")` — **Getting started** (estimands, assumptions, Bayesian & REML on `steroids`) | `vignette("export-and-qc")` — tables & sensitivity QC; `vignette("worked-examples")` — `age`, `edss`, and other datasets |

## Installation

``` r
# install.packages("remotes")
remotes::install_github("rvpanaro/interactmeta")
```

Documentation: <https://rvpanaro.github.io/interactmeta/>

CAMS supports **Bayesian** (default) and **frequentist** inference
(`method = "REML"`). See the getting-started vignette for the full
workflow.

## Example

``` r
library(interactmeta)

st <- steroids[steroids$n > 0, ]
es <- escalc(
  measure = "OR",
  ai = treat.events,
  bi = treat.total - treat.events,
  ci = cont.events,
  di = cont.total - cont.events,
  data = st
)
both <- c("CAPE COVID", "COVID STEROID", "REMAP-CAP", "Steroids-SARI")
fit <- cams(es[es$study %in% both, ], reference = "NIV", method = "REML")
interactions(fit, exponentiate = TRUE)
modifiers(fit, exponentiate = TRUE)
forestplot(fit, exponentiate = TRUE)
```

## Citation

``` bibtex
@article{panaro2025cams,
  title   = {Consistent Bayesian meta-analysis on subgroup specific effects and interactions},
  author  = {Panaro, Renato and R{\"o}ver, Christian and Friede, Tim},
  journal = {arXiv preprint arXiv:2512.18785},
  year    = {2025},
  doi     = {10.48550/arXiv.2512.18785}
}
```
