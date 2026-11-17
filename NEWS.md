# interactmeta (development)

## interactmeta 0.1.0

### New features

- **`tidy()`** method for **`interactmeta`** fits: one long export table across
  interaction, modifier, combined, and heterogeneity estimands.
- **`plot(fit, type = ...)`** routes contribution (default), heterogeneity, and
  forest plots.
- **`sensitivity()`** and **`compare()`** for construction sensitivity analyses.

### Performance

- **`cams()`** defaults to **`precompute = FALSE`** and **`lazy_latent = TRUE`**:
  fit returns after the bayesmeta blocks and DIRECT support, without warming
  marginal caches or per-trial latent convolves. **`print()`**, **`summary()`**,
  and **`tidy()`** build reporting on first use; pass **`precompute = TRUE`** for
  the previous eager behaviour.
- **`fit_args = list(n_support = …)`** is forwarded to both **`bmr`** and
  **`bayesmeta`** as **`n.support`** to reduce quadrature grid size.
- On-demand trial latents for **`shrinkage()`** and Bayesian **`forestplot()`**
  shrinkage rows use the **independent-block mixture product** only
  (`trial_marginals = "trials"`). **`lazy_latent = FALSE`** still builds trial
  blocks at fit time with the same mixture product as on-demand shrinkage (no
  per-trial **`bayesmeta::convolve()`** at fit).

### Documentation

- **`vignette("comparing-da-ad-swada")`** — DA, AD, SWADA vs CAMS after
  **`vignette("cams")`** (`respire_age`, `sensitivity()`, `compare()`).
- Split **`vignette("diagnostics")`** into four articles (pkgdown **After CAMS**):
  **`heterogeneity-after-cams`**, **`trial-specific-effects-and-shrinkage`**,
  **`traceplots`**, **`single-subgroup-trials`** (each refits `steroids` like
  **`vignette("cams")`**). Old `diagnostics.html` redirects to `cams.html`.

### Breaking changes

- **`traceplot(fit, type = "modifier")`** removed (use `type = "interaction"` or
  `"combined"`).
- **Unexported:** `plot_contribution()`, `plot_heterogeneity()` — use
  **`plot(fit)`** or **`plot(fit, type = "heterogeneity")`**.
- **Unexported:** `compare_interactions()`, `compare_modifiers()`,
  `compare_combined()`, `compare_ifrac()`, `compare_weights()` — use
  **`compare(sensitivity(...), estimand = ...)`** or **`compare(fit1, fit2, ...)`**
  for two CAMS approaches.
- Marginal extractors **`modifiers()`**, **`interactions()`**, **`combined()`**
  inherit **`prob`**, **`ifrac`**, **`contrast_sep`**, and **`exponentiate`**
  from **`fit$reporting`** when those arguments are omitted (`NULL`).

### Internal

- **`interactmeta_marginal()`** canonical marginal dispatch for reporting code paths.
