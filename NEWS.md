# interactmeta 0.1.0

## New features

- **`tidy()`** method for **`interactmeta`** fits: one long export table across
  interaction, modifier, combined, and heterogeneity estimands.
- **`plot(fit, type = ...)`** routes contribution (default), heterogeneity, and
  forest plots.
- **`sensitivity()`** and **`compare()`** for construction sensitivity analyses.

## Performance

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

## Documentation

- **`vignette("comparing-da-ad-swada")`** — DA, AD, SWADA vs CAMS after
  **`vignette("cams")`** (`respire_age`, `sensitivity()`, `compare()`).
- **`vignette("heterogeneity-after-cams")`** — $\tau$ / $\tau_\gamma$
  summaries, heterogeneity plots, and Bayesian `traceplot()` on the `age` fit
  (after **`vignette("cams")`**). Old `diagnostics.html` and `traceplots.html`
  redirect on the site.

## Breaking changes

- **`subgroupmeta_palette()`** and **`subgroupmeta_theme_report()`** are no longer
  exported; use [interactmeta_palette()] and [interactmeta_theme_report()].
- **`traceplot(fit, type = "modifier")`** removed (use `type = "interaction"` or
  `"combined"`).
- **Unexported:** `plot_contribution()`, `plot_heterogeneity()` — use
  **`plot(fit)`** or **`plot(fit, type = "heterogeneity")`**.
- **Unexported:** `compare_interactions()`, `compare_modifiers()`,
  `compare_combined()`, and `compare_ifrac()` — use **`compare(sensitivity(...),
  estimand = ...)`** or **`compare(fit1, fit2, ...)`** for two CAMS approaches.
  **`compare_weights()`** and **`compare_heterogeneity()`** remain exported for
  SWADA-weight and τ side-by-side tables.
- Marginal extractors **`modifiers()`**, **`interactions()`**, **`combined()`**
  inherit **`prob`**, **`ifrac`**, **`contrast_sep`**, and **`exponentiate`**
  from **`fit$reporting`** when those arguments are omitted (`NULL`).

## Internal

- **`interactmeta_marginal()`** canonical marginal dispatch for reporting code paths.