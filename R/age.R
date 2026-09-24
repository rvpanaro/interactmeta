#' Age subgroup data for contribution-aware meta-analysis
#'
#' Trial-level subgroup estimates by age stratum with contrast-level summaries.
#'
#' @format A data frame with 14 rows and 6 variables:
#' \describe{
#'   \item{trial}{Trial identifier.}
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
#' \doi{10.1177/20552173241274618}
#'
#' \url{https://pubmed.ncbi.nlm.nih.gov/39493424/}
#' @usage data(age)
"age"
