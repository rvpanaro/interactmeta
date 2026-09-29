# Optional script: Bayesian CAMS engine comparison on RESPIRE age.
# Not run during package checks or vignette building (brms MCMC is slow).
#
# From the package root, after installing interactmeta:
#
#   source("data-raw/respire-brms.R")
#
# Seed and engine_args are recorded so a third party can regenerate the
# posterior summaries rather than relying on hard-coded numbers.

seed <- 42L
engine_args_brms <- list(
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 1,
  refresh = 0,
  seed = seed
)

library(interactmeta)

data(respire_age, package = "interactmeta")

fit_bm <- cam(
  respire_age,
  reference = "<65",
  engine = "bayesmeta",
  silent = 2
)

fit_brms <- cam(
  respire_age,
  reference = "<65",
  engine = "brms",
  seed = seed,
  silent = 2,
  engine_args = engine_args_brms
)

print(compare_subgroups(
  da(respire_age, reference = "<65"),
  swada(respire_age, reference = "<65"),
  fit_bm
))

print(tidy(subgroups(fit_bm, type = "Pooled", exponentiate = TRUE)))
print(tidy(subgroups(fit_brms, type = "Pooled", exponentiate = TRUE)))
print(tidy(interactions(fit_bm, type = "Pooled")))
print(tidy(interactions(fit_brms, type = "Pooled")))
