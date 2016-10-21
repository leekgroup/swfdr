## ------------------------------------------------------------------------
library(swfdr)

## ------------------------------------------------------------------------
data(BMI_GIANT_GWAS_sample)
head(BMI_GIANT_GWAS_sample)
dim(BMI_GIANT_GWAS_sample)

## ------------------------------------------------------------------------
table(BMI_GIANT_GWAS_sample$Freq_MAF_Int_Hapmap)

## ------------------------------------------------------------------------
X <- model.matrix(~ splines::ns(N,5) + Freq_MAF_Int_Hapmap, data = BMI_GIANT_GWAS_sample)[,-1]
head(X)

