#' Subset of SNPs from meta-analysis of BMI GWAS study.
#'
#' A dataset containing 50,000 SNPs and results for their associations with BMI.
#' 
#'
#' @format A data frame with 50,000 rows and 9 variables:
#' \describe{
#'   \item{SNP}{ID for SNP (single nucleotide polymorphism)}
#'   \item{A1}{Allele 1 for SNP}
#'   \item{A2}{Allele 2 for SNP}
#'   \item{Freq_MAF_Hapmap}{Frequency of minor allele (MAF) in Hapmap project}
#'   \item{b}{Estimated beta for association between SNP and BMI}
#'   \item{se}{Estimated standard error (se) for association between SNP and BMI}
#'   \item{p}{P-value for association between SNP and BMI}
#'   \item{N}{Total sample size considered for association of SNP and BMI}
#' }
#' @source \url{https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2015_BMI}
"BMI_GIANT_GWAS_sample"

#' P-values from abstracts from articles in 5 biomedical journals (American Journal of Epidemiology, BMJ, Jama, Lancet, New England Journal of Medicine), over 11 years (2000-2010).
#'
#' A dataset containing 15,653 p-values.
#' 
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
#' @source Code for extracting p-values at: \url{https://github.com/jtleek/swfdr/blob/master/getPvalues.R}

