#' P-values from abstracts from articles in 5 biomedical journals (American Journal of Epidemiology, BMJ, JAMA, Lancet, New England Journal of Medicine), over 11 years (2000-2010).
#'
#' A dataset containing 15,653 p-values.
#'
#' @docType data
#'
#' @usage data(journals_pVals)
#'
#' @format A data frame with 15,653 rows and 7 variables:
#' \describe{
#'   \item{pvalue}{P-value}
#'   \item{pvalueTruncated}{Equals to 1 if the p-value is truncated, 0 otherwise}
#'   \item{pubmedID}{Pubmed ID of the article}
#'   \item{year}{Year of publication}
#'   \item{abstract}{Abstract}
#'   \item{title}{Title}
#'   \item{journal}{Journal}
#' }
#' 
#' @keywords datasets
#' 
#' @source Code for extracting p-values at: \url{https://github.com/jtleek/swfdr/blob/master/getPvalues.R}
#'
"journals_pVals"
