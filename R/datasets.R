# ---- age.R ----
#' Age subgroup data for contribution-aware meta-analysis
#'
#' Trial-level subgroup estimates from a systematic review of immunotherapies
#' for relapsing-remitting multiple sclerosis (Heesen et al. 2024). Each
#' trial reports treatment effects in **Younger** and **Older** age strata
#' (log relapse **rate ratio** and variance—not mean differences in counts).
#' Trial **`CONFIRM-B`** / **`CONFIRM-C`** in this table correspond to a single
#' **`CONFIRM`** row in [edss]. Subgroup labels do not include numeric age cuts.
#' See `vignette("worked-examples")`.
#'
#' @format A data frame with 14 rows and 6 variables:
#' \describe{
#'   \item{study}{Study identifier.}
#'   \item{subgroup}{Age subgroup label.}
#'   \item{n}{Subgroup sample size.}
#'   \item{treatment}{Treatment label.}
#'   \item{yi}{Effect estimate.}
#'   \item{vi}{Variance of \code{yi} (\code{se^2} from the source data).}
#' }
#' @source
#' Heesen C, Röver C, Salem S, et al. (2024). Treatment effect modifiers of
#' immunotherapies for relapsing-remitting multiple sclerosis: a systematic review
#' and meta-analysis. \emph{Mult Scler J Exp Transl Clin}, 10(4), 20552173241274618.
#' \url{https://doi.org/10.1177/20552173241274618}
#'
#' \url{https://pubmed.ncbi.nlm.nih.gov/39493424/}
#' @usage data(age)
"age"

# ---- edss.R ----
#' Subgroup relapse data for contribution-aware meta-analysis
#'
#' Same MS immunotherapy review as [age]: trials split by **Lower** vs
#' **Higher EDSS** (disability) strata. Log relapse rate ratios and
#' variances (not mean differences). Numeric EDSS thresholds are not stored.
#' CONFIRM BG-12 and copaxone arms are separate trial ids (**`CONFIRM-B`**,
#' **`CONFIRM-C`**) so CAMS does not pool them under one comparison (shared
#' control is not modelled). Includes the **COP1 pilot** copaxone comparison
#' (**`COP1 PILOT`**) used in the EDSS CAMS analyses in Panaro, Röver & Friede.
#' See `vignette("worked-examples")`.
#'
#' @format A data frame with 16 rows and 6 variables:
#' \describe{
#'   \item{study}{Study identifier.}
#'   \item{subgroup}{EDSS subgroup factor with levels \code{Lower} (default
#'     reference) then \code{Higher}.}
#'   \item{n}{Subgroup sample size.}
#'   \item{treatment}{Treatment label.}
#'   \item{yi}{Log effect estimate.}
#'   \item{vi}{Variance of \code{yi}.}
#' }
#' @source
#' Heesen C, Röver C, Salem S, et al. (2024). Treatment effect modifiers of
#' immunotherapies for relapsing-remitting multiple sclerosis: a systematic review
#' and meta-analysis. \emph{Mult Scler J Exp Transl Clin}, 10(4), 20552173241274618.
#' \url{https://doi.org/10.1177/20552173241274618}
#'
#' \url{https://pubmed.ncbi.nlm.nih.gov/39493424/}
#' @usage data(edss)
"edss"

# ---- il6.R ----
#' IL-6 antagonist COVID-19 mortality data (REACT meta-analysis)
#'
#' Event counts by trial and corticosteroid subgroup from the WHO REACT
#' prospective meta-analysis of IL-6 antagonists in patients hospitalized
#' for COVID-19. Some trials report only one corticosteroid subgroup
#' (`n = 0` patients in the other subgroup); drop those empty rows before [escalc()] and
#' let [ifrac()] / [cams()] complete the missing cells. See
#' `vignette("worked-examples")`.
#'
#' @format A data frame with 30 rows and 7 variables:
#' \describe{
#'   \item{study}{Study identifier.}
#'   \item{subgroup}{Corticosteroid subgroup (\code{yes} / \code{no}).}
#'   \item{n}{Total patients (`treat.total` + `cont.total`).}
#'   \item{treat.events}{Events in the treatment arm.}
#'   \item{treat.total}{Total subjects in the treatment arm.}
#'   \item{cont.events}{Events in the control arm.}
#'   \item{cont.total}{Total subjects in the control arm.}
#' }
#' @source
#' WHO Rapid Evidence Appraisal for COVID-19 Therapies (REACT) Working Group,
#' Sterne JAC, Murthy S, et al. (2021). Association between administration of
#' IL-6 antagonists and mortality among patients hospitalized for COVID-19:
#' a meta-analysis. \emph{JAMA}, 326(6), 499--518.
#' \url{https://doi.org/10.1001/jama.2021.11330}
#'
#' \url{https://pubmed.ncbi.nlm.nih.gov/34228774/}
#' @usage data(il6)
"il6"

# ---- respire.R ----
#' Log-hazard ratio and SE from a reported HR confidence interval
#'
#' For a confidence level \eqn{c},
#' \eqn{\alpha = 1-c}, \eqn{z = \Phi^{-1}(1-\alpha/2)}, and
#' \deqn{y = \log(\mathrm{HR}),\qquad
#' \mathrm{SE} = \{\log(U)-\log(L)\}/(2z).}
#'
#' Do **not** substitute \eqn{z=1.96} when \eqn{c \ne 0.95}. RESPIRE 1
#' subgroup intervals use \eqn{c=0.975}; RESPIRE 2 uses \eqn{c=0.951}.
#'
#' @param hr Hazard ratio.
#' @param ci_lower,ci_upper Confidence limits on the HR scale.
#' @param ci_level Confidence level in \eqn{(0,1)}, recycled to `hr`.
#' @return A data frame with `yi`, `sei`, `vi`, and `z`.
#' @export
#' @examples
#' log_hr_from_ci(hr = 0.61, ci_lower = 0.32, ci_upper = 1.16, ci_level = 0.975)
log_hr_from_ci <- function(hr, ci_lower, ci_upper, ci_level) {
  hr <- as.numeric(hr)
  ci_lower <- as.numeric(ci_lower)
  ci_upper <- as.numeric(ci_upper)
  ci_level <- as.numeric(ci_level)
  n <- length(hr)
  if (length(ci_level) == 1L) {
    ci_level <- rep(ci_level, n)
  }
  if (length(ci_lower) != n || length(ci_upper) != n || length(ci_level) != n) {
    stop("`hr`, `ci_lower`, `ci_upper`, and `ci_level` must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(hr) | hr <= 0)) {
    stop("`hr` must be positive and finite.", call. = FALSE)
  }
  if (any(!is.finite(ci_lower) | ci_lower <= 0 | !is.finite(ci_upper) | ci_upper <= 0)) {
    stop("CI limits must be positive and finite.", call. = FALSE)
  }
  if (any(ci_upper <= ci_lower)) {
    stop("`ci_upper` must exceed `ci_lower`.", call. = FALSE)
  }
  if (any(!is.finite(ci_level) | ci_level <= 0 | ci_level >= 1)) {
    stop("`ci_level` must lie in (0, 1).", call. = FALSE)
  }
  alpha <- 1 - ci_level
  z <- stats::qnorm(1 - alpha / 2)
  yi <- log(hr)
  sei <- (log(ci_upper) - log(ci_lower)) / (2 * z)
  data.frame(
    yi = yi,
    sei = sei,
    vi = sei^2,
    z = z,
    stringsAsFactors = FALSE
  )
}

#' Reconstruct an HR confidence interval from `yi` and `sei`
#'
#' Inverse of [log_hr_from_ci()] on the HR scale.
#'
#' @param yi Log hazard ratio.
#' @param sei Standard error of `yi`.
#' @param ci_level Confidence level used to form `sei`.
#' @return A data frame with `hr`, `ci_lower`, and `ci_upper`.
#' @export
hr_ci_from_log <- function(yi, sei, ci_level) {
  yi <- as.numeric(yi)
  sei <- as.numeric(sei)
  ci_level <- as.numeric(ci_level)
  n <- length(yi)
  if (length(ci_level) == 1L) {
    ci_level <- rep(ci_level, n)
  }
  alpha <- 1 - ci_level
  z <- stats::qnorm(1 - alpha / 2)
  data.frame(
    hr = exp(yi),
    ci_lower = exp(yi - z * sei),
    ci_upper = exp(yi + z * sei),
    stringsAsFactors = FALSE
  )
}

#' RESPIRE age-subgroup hazard ratios (ciprofloxacin DPI)
#'
#' Published subgroup results for **time to first exacerbation**,
#' ciprofloxacin DPI 14-day on/off versus pooled placebo, in the replicate
#' Phase III RESPIRE 1 and RESPIRE 2 trials in non-CF bronchiectasis.
#'
#' Columns `hr`, `ci_lower`, `ci_upper`, and `ci_level` are the original
#' published values. `yi`, `sei`, and `vi` are derived by [log_hr_from_ci()];
#' RESPIRE 1 uses 97.5% intervals and RESPIRE 2 uses 95.1% intervals.
#' See `vignette("worked-examples")`.
#'
#' @format A data frame with 4 rows and 10 variables:
#' \describe{
#'   \item{study}{`RESPIRE 1` or `RESPIRE 2`.}
#'   \item{subgroup}{Age stratum (`<65`, `>=65`).}
#'   \item{n}{Subgroup sample size.}
#'   \item{hr}{Published hazard ratio.}
#'   \item{ci_lower,ci_upper}{Published HR confidence limits.}
#'   \item{ci_level}{Confidence level of those limits (`0.975` or `0.951`).}
#'   \item{yi}{\eqn{\log(\mathrm{HR})}.}
#'   \item{sei}{Standard error of `yi` from the published CI.}
#'   \item{vi}{`sei^2`.}
#' }
#' @source
#' Subgroup estimates were extracted from FDA review / briefing materials for
#' the ciprofloxacin dry-powder inhalation (DPI) RESPIRE program (14-day on/off
#' regimen versus pooled placebo). Interval levels match the primary-endpoint
#' multiplicity adjustments in the publications:
#'
#' De Soyza A, Aksamit T, Bandel T-J, et al. (2018). RESPIRE 1: a phase III
#' placebo-controlled randomised trial of ciprofloxacin dry powder for inhalation
#' in non-cystic fibrosis bronchiectasis. \emph{Eur Respir J}, 51(1), 1702052.
#' \url{https://doi.org/10.1183/13993003.02052-2017}
#'
#' Aksamit T, De Soyza A, Bandel T-J, et al. (2018). RESPIRE 2: a phase III
#' placebo-controlled randomised trial of ciprofloxacin dry powder for inhalation
#' in non-cystic fibrosis bronchiectasis. \emph{Eur Respir J}, 51(1), 1702053.
#' \url{https://doi.org/10.1183/13993003.02053-2017}
#'
#' \url{https://www.accessdata.fda.gov/}
#' @usage data(respire_age)
#' @seealso [respire_sex], [log_hr_from_ci()]
"respire_age"

#' RESPIRE sex-subgroup hazard ratios (ciprofloxacin DPI)
#'
#' Same endpoint, treatment comparison, and source as [respire_age()], for
#' sex subgroups (Female / Male).
#'
#' @format A data frame with 4 rows and 10 variables (see [respire_age()]).
#' @source See [respire_age()].
#' @usage data(respire_sex)
#' @seealso [respire_age], [log_hr_from_ci()]
"respire_sex"

# ---- steroids.R ----
#' Corticosteroid COVID-19 mortality data (REACT meta-analysis)
#'
#' Event counts by trial, ventilation subgroup, and treatment from the WHO REACT
#' prospective meta-analysis of systemic corticosteroids in critically ill
#' patients with COVID-19. Some trials report only one ventilation subgroup
#' (`n = 0` patients in the other subgroup); drop those empty rows before [escalc()] and
#' let [ifrac()] / [cams()] complete the missing cells. See
#' `vignette("cams")`.
#'
#' @format A data frame with 14 rows and 8 variables:
#' \describe{
#'   \item{study}{Study identifier.}
#'   \item{subgroup}{Ventilation support category (\code{IV}, \code{NIV}, etc.).}
#'   \item{treatment}{Treatment label (drug).}
#'   \item{n}{Total patients (`treat.total` + `cont.total`).}
#'   \item{treat.events}{Events in the treatment arm.}
#'   \item{treat.total}{Total subjects in the treatment arm.}
#'   \item{cont.events}{Events in the control arm.}
#'   \item{cont.total}{Total subjects in the control arm.}
#' }
#' @source
#' WHO Rapid Evidence Appraisal for COVID-19 Therapies (REACT) Working Group,
#' Sterne JAC, Murthy S, et al. (2020). Association between administration of
#' systemic corticosteroids and mortality among critically ill patients with
#' COVID-19: a meta-analysis. \emph{JAMA}, 324(13), 1330--1341.
#' \url{https://doi.org/10.1001/jama.2020.17023}
#'
#' \url{https://pubmed.ncbi.nlm.nih.gov/32876694/}
#' @usage data(steroids)
"steroids"
