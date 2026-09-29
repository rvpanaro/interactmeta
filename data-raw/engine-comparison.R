# Optional script: CAMS engine × method comparison on the bundled
# `edss` data. The vignette `vignettes/engine-comparison.Rmd` fits the same
# objects live when the backends are installed. This script is the
# documented procedure for regenerating those fits outside a vignette
# build (especially the optional internal brms backend).
#
# From the package root, after installing interactmeta:
#
#   source("data-raw/engine-comparison.R")
#
# Seed and engine_args are recorded so a third party can regenerate
# summaries rather than relying on hard-coded numbers.

seed <- 42L
mcmc_args <- list(
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 1,
  refresh = 0,
  seed = seed
)

library(interactmeta)
library(metafor)

data(edss, package = "interactmeta")
edss_es <- escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)

fit_da <- da(edss_es, reference = "Lower", method = "REML")
fit_ad <- ad(edss_es, reference = "Lower", method = "REML")
fit_sw <- swada(edss_es, reference = "Lower", method = "REML")

fit_bm <- cam(
  edss_es,
  reference = "Lower",
  engine = "bayesmeta",
  silent = 2,
  seed = seed
)
fit_reml <- cam(
  edss_es,
  reference = "Lower",
  engine = "metafor",
  method = "REML",
  silent = 2
)
fit_ml <- cam(
  edss_es,
  reference = "Lower",
  engine = "metafor",
  method = "ML",
  silent = 2
)

print(tidy(interactions(fit_bm, type = "Pooled")))
print(summary(fit_bm)$parameters)
print(tidy(fit_ad))
print(tidy(fit_sw))

if (requireNamespace("brms", quietly = TRUE)) {
  fit_brms <- bayescam(
    edss_es,
    reference = "Lower",
    engine = "brms",
    silent = 2,
    seed = seed,
    engine_args = list(
      chains = mcmc_args$chains,
      iter = mcmc_args$iter,
      warmup = mcmc_args$warmup,
      cores = mcmc_args$cores,
      refresh = mcmc_args$refresh,
      open_progress = FALSE
    )
  )
  print(tidy(interactions(fit_brms, type = "Pooled")))
  print(summary(fit_brms)$parameters)
}
