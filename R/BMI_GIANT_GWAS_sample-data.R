#' Subset of SNPs from meta-analysis of BMI GWAS study.
#'
#' A dataset containing 50,000 SNPs and results for their associations with BMI.
#' 
#' @docType data
#'
#' @usage data(BMI_GIANT_GWAS_sample)
#'
#' @return Object of class tbl_df, tbl, data.frame.
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
#'   \item{Freq_MAF_Int_Hapmap}{Three approximately equal intervals for the Hapmap MAFs}
#' }
#'
#' @keywords datasets
#' 
#' @source \url{https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2015_BMI}
#'
"BMI_GIANT_GWAS_sample"