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

respire_build_dataset <- function(trial,
                                  subgroup,
                                  n,
                                  hr,
                                  ci_lower,
                                  ci_upper,
                                  ci_level,
                                  subgroup_levels) {
  trans <- log_hr_from_ci(hr, ci_lower, ci_upper, ci_level)
  data.frame(
    trial = as.character(trial),
    subgroup = factor(as.character(subgroup), levels = subgroup_levels),
    n = as.integer(n),
    hr = as.numeric(hr),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    ci_level = as.numeric(ci_level),
    yi = trans$yi,
    sei = trans$sei,
    vi = trans$vi,
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
#' See `vignette("consistent-interactions")`.
#'
#' @format A data frame with 4 rows and 10 variables:
#' \describe{
#'   \item{trial}{`RESPIRE 1` or `RESPIRE 2`.}
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
