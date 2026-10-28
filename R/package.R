# ---- interactmeta-package.R ----
#' Meta-analysis of treatment effect modification from modifier-level aggregate data
#'
#' **interact** implements **contribution-adjusted meta-analysis
#' (CAMS)** for treatment effect modification from
#' modifier-level aggregate data (Panaro, Röver & Friede 2025). Each trial
#' contributes two modifier-level estimates; CAMS derives an information fraction
#' \eqn{\pi_j}, fits a random-effects meta-analysis of within-trial
#' contrasts and a random-effects meta-regression of precision-weighted
#' combined effects, then reports the interaction \eqn{\gamma}, modifier-level
#' means \eqn{\mu_A(\pi)} and \eqn{\mu_B(\pi)}, and the combined-population
#' mean \eqn{M(\pi)} at a chosen reporting fraction.
#'
#' Fitter: [cams()]. Bayesian inference (default) uses **bayesmeta**;
#' frequentist inference uses **metafor** (REML or ML). The statistical
#' model is the same; intervals are credible or confidence intervals
#' respectively. See `vignette("cams")`.
#'
#' [da()] (Difference of Averages) and [ad()] (Average Difference) are
#' **different constructions** (Panaro, Röver & Friede 2026), not two names
#' for one estimand and not alternative CAMS models. [swada()] applies the
#' Same Weighting Across Different Analyses (SWADA) principle so those
#' constructions share study contributions.
#' [conventional()] is a compatibility bundle of DA and AD with
#' analysis-specific weights. A Bayesian posterior on DA or AD does not add
#' information-fraction estimands.
#'
#' @section CAMS workflow:
#'
#' 1. [ifrac()] (optional before fit) builds trial information fractions
#'    \eqn{\pi_j=\sigma_{Bj}^{-2}/(\sigma_{Aj}^{-2}+\sigma_{Bj}^{-2})};
#'    [cams()] computes them from `vi` by default.
#' 2. Within each trial, form contrast \eqn{g_j=y_{Bj}-y_{Aj}} and combined
#'    coordinate \eqn{m_j=(1-\pi_j)y_{Aj}+\pi_j y_{Bj}}.
#' 3. Fit a random-effects meta-analysis of \eqn{g_j} and a random-effects
#'    meta-regression of \eqn{m_j} on \eqn{\pi_j}.
#' 4. [interactions()], [modifiers()], and [combined()] summarize
#'    \eqn{\gamma} and effects at the reporting information fraction.
#' 5. [plot()] and [forestplot()] visualize effects across \eqn{\pi} and by
#'    trial (see `vignette("cams")` and `vignette("worked-examples")`).
#'
#' @section Documentation:
#' \describe{
#'   \item{Getting started}{`vignette("cams")` — estimands, assumptions, Bayesian and REML on the IV/NIV `steroids` example.}
#'   \item{Export and QC}{`vignette("export-and-qc")` — `tidy()` exports and DA/AD/SWADA/CAMS sensitivity tables.}
#'   \item{Other datasets}{`vignette("worked-examples")` — brief forests for `age`, `edss`, `il6`, and RESPIRE tables (including incomplete-subgroup patterns).}
#' }
#'
#' @section Inference:
#'
#' CAMS can be fitted using Bayesian or frequentist inference. These are
#' not two different scientific models. Default Bayesian fitting uses
#' **bayesmeta** (Röver 2020); pass `method = "REML"` (or `"ML"`) for
#' **metafor** (Viechtbauer 2010). See `vignette("cams")`.
#'
#' Default Bayesian heterogeneity priors match [bayesmeta::bayesmeta()]:
#' improper uniform on \eqn{\tau_m} and \eqn{\tau_g}. Fixed effects use
#' independent \eqn{\mathrm{Normal}(0,10)} on \eqn{(\alpha,\theta)} and on
#' \eqn{\gamma}.
#' Empirical work indicates interaction heterogeneity is often smaller than
#' overall heterogeneity (Panaro, Röver & Friede, arXiv:2606.23968); vignettes
#' set half-normal scales from that literature explicitly. See [cams()].
#'
#' @section Reporting:
#'
#' [plot()] (default contribution view) plots modifier-level and combined effects
#' across the information fraction. [forestplot()] shows observed estimates and
#' trial-level CAMS shrinkage (grey markers; `include_shrinkage = TRUE` by
#' default). Typography helpers
#' ([interactmeta_theme_report()], [interactmeta_palette()]) are documented
#' on their help pages.
#'
#' @section Reference:
#' The package implements CAMS:
#'
#' Panaro RV, Röver C, Friede T (2025).
#' *Consistent Bayesian meta-analysis on subgroup specific effects and interactions.*
#' arXiv:2512.18785.
#' \url{https://doi.org/10.48550/arXiv.2512.18785}
#'
#' Inference backends:
#'
#' Viechtbauer W (2010).
#' Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software}, 36(3), 1--48.
#' \url{https://doi.org/10.18637/jss.v036.i03}
#'
#' Röver C (2020).
#' Bayesian random-effects meta-analysis using the bayesmeta R package.
#' \emph{Journal of Statistical Software}, 93(6), 1--51.
#' \url{https://doi.org/10.18637/jss.v093.i06}
#'
#' Related methods (not the CAMS model). Difference of Averages (DA) and
#' Average Difference (AD) are different constructions:
#'
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \url{https://doi.org/10.48550/arXiv.2508.15531}
#'
#' Empirical heterogeneity priors (optional; not package defaults):
#'
#' Panaro RV, Röver C, Friede T (2026).
#' Empirical prior distributions for treatment-by-subgroup interaction
#' heterogeneity in random-effects meta-analysis.
#' arXiv:2606.23968. \url{https://doi.org/10.48550/arXiv.2606.23968}
#'
#' Fisher DJ, Carpenter JR, Morris TP, Freeman SC, Tierney JF (2017).
#' Meta-analytical methods to identify who benefits most from treatments:
#' daft, deluded, or deft approach? \emph{BMJ}, 356, j573.
#' \url{https://doi.org/10.1136/bmj.j573}
#'
#' van Houwelingen HC, Arends LR, Stijnen T (2002).
#' Advanced methods in meta-analysis: multivariate approach and meta-regression.
#' \emph{Statistics in Medicine}, 21(4), 589--624. \url{https://doi.org/10.1002/sim.1040}
#'
#' van Houwelingen HC, Zwinderman KH, Stijnen T (1993).
#' A bivariate approach to meta-analysis. \emph{Statistics in Medicine},
#' 12(24), 2273--2284. \url{https://doi.org/10.1002/sim.4780122405}
#'
#' @seealso [cams()], [ifrac()], [interactions()], [modifiers()],
#'   [combined()], [tidy()], [sensitivity()], [compare()], [predict_cams_by_if()],
#'   [da()], [ad()],
#'   [swada()], `vignette("cams", package = "interact")`,
#'   `vignette("export-and-qc", package = "interact")`,
#'   `vignette("worked-examples", package = "interact")`
#'
#' @aliases interact interact-package interactmeta
"_PACKAGE"

# ---- reexports.R ----
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

# ---- zzz.R ----
#' @importFrom stats setNames
#' @importFrom paletteer paletteer_d
#' @importFrom scales alpha
NULL

utils::globalVariables(c(".data", ".prec"))

.onLoad <- function(libname, pkgname) {
  # bayesmeta registers S3 methods forestplot.escalc / forestplot.bayesmeta /
  # forestplot.bmr when its namespace is loaded. Keep them available for the
  # session once interact is loaded (Imports already ensures the package is
  # installed).
  loadNamespace("bayesmeta")
  loadNamespace("forestplot")
  loadNamespace("metafor")
  invisible()
}

.onAttach <- function(libname, pkgname) {
  invisible()
}
