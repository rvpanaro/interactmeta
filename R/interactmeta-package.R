#' Meta-analysis of treatment-by-subgroup interactions using subgroup-level aggregate data
#'
#' **interactmeta** implements **contribution-adjusted meta-analysis
#' (CAMS)** for treatment-by-subgroup interactions from
#' subgroup-level aggregate data (Panaro, Röver & Friede 2025). The
#' methodological map is
#' \deqn{(y_{Aj},y_{Bj}) \to (m_j,g_j,\pi_j) \to
#' \{\gamma,\,\mu_A(\pi),\,\mu_B(\pi),\,M(\pi)\}.}
#' Subgroup-specific effects are inputs and derived outputs. The target is
#' not a separate meta-analysis within each subgroup.
#'
#' Preferred fitter: [cam()]. [cams()] and [subgroupmeta()] are aliases.
#' CAMS can be fitted using Bayesian or frequentist inference. Bayesian CAMS
#' is implemented using **bayesmeta**/DIRECT quadrature on the two
#' orthogonal blocks; frequentist CAMS is implemented using **metafor**.
#' Those packages fit random-effects series. interactmeta constructs the
#' information fraction \(\pi_j\), rotates the observations, reconstructs
#' contribution-adjusted estimands, and preserves
#' \(\mu_B(\pi)-\mu_A(\pi)=\gamma\) and
#' \(M(\pi)=(1-\pi)\mu_A(\pi)+\pi\mu_B(\pi)\).
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
#' 1. [ifrac()] builds trial information fractions
#'    \(\pi_j=\sigma_{Bj}^{-2}/(\sigma_{Aj}^{-2}+\sigma_{Bj}^{-2})\).
#' 2. The observations rotate to \(g_j=y_{Bj}-y_{Aj}\) and
#'    \(m_j=(1-\pi_j)y_{Aj}+\pi_j y_{Bj}\).
#' 3. The overall block is a meta-regression of \(m_j\) on \(\pi_j\); the
#'    interaction block is a random-effects analysis of \(g_j\).
#' 4. Extractors return the interaction [interactions()], subgroup effects
#'    [subgroups()], and combined-population effect [combined()]
#'    ([overall()] is an alias) at a reporting \(\pi\).
#' 5. [plot_contribution()] / [predict_cam_by_if()] evaluate
#'    \(\mu_A(\pi)\), \(\mu_B(\pi)\), and \(M(\pi)\) across \(\pi\).
#'    Trial-level reconstruction \(\mu_{Aj}=m_j-\pi_j g_j\) is stored on
#'    `fit$latent` (see `vignette("convolution")`).
#'
#' @section Inference:
#'
#' CAMS can be fitted using Bayesian or frequentist inference. These are
#' not two different scientific models. The approaches are **bayesmeta**
#' (DIRECT, the default) and **metafor** (REML or ML). See
#' `vignette("engine-comparison")`.
#'
#' Default Bayesian priors are \(\mathrm{Normal}(0,10)\) on the fixed
#' effects, \(\mathrm{Normal}^+(0,0.5)\) on \(\tau_m\), and
#' \(\mathrm{Normal}^+(0,0.1)\) on \(\tau_g\). See [cam_prior()].
#'
#' @section Reporting:
#'
#' [plot_contribution()] is the contribution-adjusted meta-regression.
#' [forestplot()] overlays trial-level CAMS reconstruction (grey shrinkage;
#' `include_shrinkage = TRUE` by default). Typography helpers
#' ([interactmeta_theme_report()], [interactmeta_palette()]) are documented
#' on their help pages.
#'
#' Fitted objects store provenance for [audit()]; see
#' [interactmeta-auditability].
#'
#' @section Reference:
#' The package implements CAMS:
#'
#' Panaro RV, Röver C, Friede T (2025).
#' *Consistent Bayesian meta-analysis on subgroup specific effects and interactions.*
#' arXiv:2512.18785.
#' \doi{10.48550/arXiv.2512.18785}
#'
#' Related methods (not the CAMS model). Difference of Averages (DA) and
#' Average Difference (AD) are different constructions:
#'
#' Panaro RV, Röver C, Friede T (2026).
#' Subgroup comparisons within and across studies in meta-analysis.
#' \emph{Research Synthesis Methods} (in press).
#' arXiv:2508.15531. \doi{10.48550/arXiv.2508.15531}
#'
#' Fisher DJ, Carpenter JR, Morris TP, Freeman SC, Tierney JF (2017).
#' Meta-analytical methods to identify who benefits most from treatments:
#' daft, deluded, or deft approach? \emph{BMJ}, 356, j573.
#' \doi{10.1136/bmj.j573}
#'
#' van Houwelingen HC, Arends LR, Stijnen T (2002).
#' Advanced methods in meta-analysis: multivariate approach and meta-regression.
#' \emph{Statistics in Medicine}, 21(4), 589--624. \doi{10.1002/sim.1040}
#'
#' van Houwelingen HC, Zwinderman KH, Stijnen T (1993).
#' A bivariate approach to meta-analysis. \emph{Statistics in Medicine},
#' 12(24), 2273--2284. \doi{10.1002/sim.4780122405}
#'
#' @seealso [cam()], [bayescam()], [ifrac()], [interactions()], [subgroups()],
#'   [combined()], [plot_contribution()], [predict_cam_by_if()], [da()], [ad()],
#'   [swada()], `vignette("cam", package = "interactmeta")`,
#'   `vignette("consistent-interactions", package = "interactmeta")`,
#'   `vignette("single-subgroup-trials", package = "interactmeta")`,
#'   `vignette("further", package = "interactmeta")`,
#'   `vignette("convolution", package = "interactmeta")`,
#'   `vignette("engine-comparison", package = "interactmeta")`
#'
#' @aliases interactmeta
"_PACKAGE"
