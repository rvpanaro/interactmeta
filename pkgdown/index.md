# interactmeta

<p class="badges">
<a href="https://github.com/rvpanaro/interactmeta/actions/workflows/R-CMD-check.yaml"><img src="https://github.com/rvpanaro/interactmeta/actions/workflows/R-CMD-check.yaml/badge.svg" alt="R-CMD-check"></a>
<a href="https://github.com/rvpanaro/interactmeta/actions/workflows/pkgdown.yaml"><img src="https://github.com/rvpanaro/interactmeta/actions/workflows/pkgdown.yaml/badge.svg" alt="pkgdown"></a>
<a href="https://codecov.io/gh/rvpanaro/interactmeta"><img src="https://codecov.io/gh/rvpanaro/interactmeta/branch/master/graph/badge.svg" alt="Codecov"></a>
</p>

**Contribution-adjusted meta-analysis (CAMS)** for treatment effect modification
from modifier-level aggregate data ([Panaro, Röver & Friede 2025](https://arxiv.org/abs/2512.18785)).

`cams()` fits CAMS from paired subgroup summaries; `interactions()`, `modifiers()`,
and `combined()` report estimands; `forestplot()` and `plot()` visualize results.

## Articles

| | |
|---|---|
| [Getting started with CAMS](articles/cams.html) | Paper walk-through on `age` (Bayesian) and `steroids` (REML) |
| [Heterogeneity after CAMS](articles/heterogeneity-after-cams.html) | $\tau$ / $\tau_\gamma$, heterogeneity plots, traceplots |
| [DA, AD, and SWADA vs CAMS](articles/comparing-da-ad-swada.html) | Constructions after the main CAMS workflow |
| [Example datasets](articles/worked-examples.html) | Package data tables (`steroids`, `age`, RESPIRE, …) |

## Install

```r
# install.packages("remotes")
remotes::install_github("rvpanaro/interactmeta")
```

See the [Reference](reference/index.html) for all functions and [NEWS](news/index.html) for changes.