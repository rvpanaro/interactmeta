#' Auditability of CAMS fits
#'
#' A [cam()] fit is intended to be independently inspectable: a
#' third party should be able to trace reported estimates from the input data
#' through transformations, model assumptions, the chosen inferential
#' approach, and derived quantities.
#'
#' @section Provenance chain:
#'
#' \preformatted{
#' user input
#'   -> ifrac() design (information fractions)
#'   -> standardized fit data
#'   -> trial-effect blocks (binary CAMS: m_j, g_j, pi_j)
#'   -> model_spec + priors
#'   -> inference (Bayesian: bayesmeta/DIRECT; frequentist: metafor)
#'   -> fit (native backend + common $draws)
#'   -> derived summaries / plots
#' }
#'
#' Inspect a fit with [audit()]:
#'
#' \preformatted{
#' a <- audit(fit)
#' a$model
#' a$priors
#' a$computation
#' a$assumptions
#' a$transformed_data
#' a$draws_semantics
#' a$software
#' }
#'
#' @section Stored fields:
#'
#' \describe{
#'   \item{`$call`}{The [cam()] call.}
#'   \item{`$model_spec`}{Structured model choices (subgroups, IF definition,
#'     ecological slope, RE covariance, parameterization, method).}
#'   \item{`$inference` / `$engine` / `$method` / `$algorithm`}{Inferential
#'     approach, computational backend, method, and algorithm (not
#'     interchangeable). `$engine` records the R implementation
#'     (`"bayesmeta"` or `"metafor"`).}
#'   \item{`$priors`}{Normalized CAMS priors (`$common`) plus backend translation
#'     (`$engine`). Top-level `$tau` / `$tau_interaction` / `$fixed` remain for
#'     compatibility.}
#'   \item{`$computation`}{Seed, stochastic vs deterministic, and `engine_args`
#'     actually used.}
#'   \item{`$assumptions`}{Including `subgroup_covariance = "assumed_zero"`.}
#'   \item{`$transformed_data`}{Observed subgroup estimates and, for binary CAMS,
#'     trial blocks `m_j`, `g_j`, `pi_j` with variances.}
#'   \item{`$draws_semantics`}{What a `$draws` row is: MCMC draw, quadrature
#'     node, or point estimate. Bayesmeta rows are **not** IID posterior samples.}
#'   \item{`$provenance`}{R and package versions, seed, UTC timestamp.}
#' }
#'
#' Native backend objects are retained: `$bayesmeta`, `$metafor`.
#'
#' @section Inferential implementations:
#'
#' **Bayesian CAMS** (implemented using bayesmeta/DIRECT). Overall
#' `bmr(m_j ~ pi_j)` and interaction `bayesmeta(g_j)`, combined by independent
#' product quadrature. Posterior inference uses the DIRECT Gaussian-mixture
#' support, including heterogeneity uncertainty for \(\tau_m\) and
#' \(\tau_g\). Intervals are **credible intervals**. Reconstruction of CAMS
#' quantities keeps the required dependence: the joint \((\alpha,\theta)\)
#' support of the overall block with the independent \(\gamma\) block.
#' \eqn{\alpha} = overall intercept, \eqn{\theta} = IF slope, \eqn{\gamma} =
#' interaction mean, stored `b_ifrac` is \eqn{\beta = \theta - \gamma}. `$mu`
#' is the fixed-effect linear predictor. `fit$latent$mu_A[[j]]` / `$mu_B[[j]]`
#' are **trial-specific** \eqn{\mu_{Aj}=m_j-\pi_j g_j},
#' \eqn{\mu_{Bj}=m_j+(1-\pi_j)g_j}. Population effects at a reporting IF are
#' [subgroups()]:
#' \eqn{\mu_A(\pi)=\alpha+\pi b_{\mathrm{ifrac}}},
#' \eqn{\mu_B(\pi)=\alpha+\pi b_{\mathrm{ifrac}}+\gamma}.
#'
#' **Frequentist CAMS** (implemented using metafor). Same factorization and
#' the same CAMS estimands, with `rma` (REML/ML). `$draws` is one
#' point-estimate row. Intervals are **confidence intervals** from the stored
#' `rma` covariance.
#'
#' Bayesian CAMS places independent \eqn{N} on \eqn{(\alpha,\theta)} and
#' \eqn{\gamma} (bmr/bayesmeta).
#'
#' @section Derived inference:
#'
#' After one CAMS fit, [interactions()], [subgroups()], [combined()],
#' [heterogeneity()], [ifrac()], [summary()], [tidy()], and [forestplot()]
#' are extractors or deterministic transformations of that fit. They do
#' not switch inferential approaches. Uncertainty uses the representation
#' of the fit (DIRECT mixture for Bayesian CAMS, `rma` covariance for
#' frequentist CAMS). [audit()] and [interactmeta_origin()] record that
#' chain.
#'
#' @section Statistically relevant defaults:
#'
#' Defined once (not per backend), unless noted:
#'
#' \tabular{lll}{
#' Default \tab Value \tab Where \cr
#' `engine` \tab `"metafor"` \tab [cam()] \cr
#' `method` \tab `"REML"` (or `"bayes"` if bayesmeta) \tab [cam()] \cr
#' `ecological_slope` \tab `TRUE` \tab [cam()] \cr
#' `tau.prior` \tab Half-Normal(0, 0.5) \tab [cam_prior()] defaults \cr
#' `tau.interaction.prior` \tab Half-Normal(0, 0.1) if `tau.prior` omitted;
#' inherit `tau.prior` if the user supplied it \tab `cam_priors_resolve()` \cr
#' `fixed.prior` \tab Normal(0, 10) \tab [cam_prior()] defaults \cr
#' RE covariance \tab uncorrelated (`||`) \tab design / block factorization \cr
#' subgroup covariance \tab assumed zero (`v_A + v_B`) \tab trial blocks \cr
#' `big_var` \tab `100^2` \tab missing-cell imputation in [ifrac()] \cr
#' `"uniform"` \tab improper \eqn{p(\tau)\propto 1}, \eqn{\tau\ge 0} \tab not \eqn{U(0,U)} \cr
#' }
#'
#' @section Analysis, inference, and implementation:
#'
#' These are three separate axes, recorded on `$analysis`, `$inference`,
#' `$engine`, `$method`, `$algorithm`, and `$capabilities`. DA and AD are
#' not inherently frequentist; CAMS is not inherently Bayesian. IF
#' functionality (`$capabilities$ifrac`) is `TRUE` only for CAMS. See
#' [interactmeta_capability_matrix()].
#'
#' One selected analysis + inference + backend produces one fit. Downstream
#' extractors use that fit only: Bayesian DA does not borrow from CAMS or
#' AD; Bayesian AD does not borrow from DA or CAMS.
#'
#' @section Reporting quantities:
#'
#' [interactions()], [subgroups()], [combined()], [heterogeneity()],
#' [ifrac()], and [forestplot()] summarize the **selected CAMS fit**.
#' For Bayesian CAMS that is a weighted Gaussian mixture from the stored
#' DIRECT support (`.weight` and within-tau moments), aligned with
#' `$latent`. For frequentist CAMS it is the stored `rma` coefficients and
#' covariance (delta method for derived estimands). Intervals from
#' quadrature are **credible intervals**.
#'
#' @seealso [audit()], [native_fit()], [interactmeta_origin()], [cam()],
#'   [cam_prior()]
#' @name interactmeta-auditability
NULL
