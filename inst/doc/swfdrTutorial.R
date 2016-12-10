## ------------------------------------------------------------------------
library(swfdr)

## ------------------------------------------------------------------------
data(journals_pVals)
colnames(journals_pVals)

## ------------------------------------------------------------------------
table(journals_pVals$year)
table(journals_pVals$journal)

## ------------------------------------------------------------------------
journals_pVals1 <- journals_pVals[journals_pVals$year==2005 & 
                                    journals_pVals$journal == "American Journal of Epidemiology" &
                                    journals_pVals$pvalue < 0.05,]
dim(journals_pVals1)

## ------------------------------------------------------------------------
tt <- journals_pVals1[,2]
rr <- rep(0,length(tt))
rr[tt == 0] <- (journals_pVals1[tt==0,1] == round(journals_pVals1[tt==0,1],2))
pVals <- journals_pVals1[,1]
resSwfdr <- calculateSwfdr(pValues = pVals, truncated = tt, rounded = rr, numEmIterations=100)
names(resSwfdr)

## ------------------------------------------------------------------------
resSwfdr

## ------------------------------------------------------------------------
data(BMI_GIANT_GWAS_sample)
head(BMI_GIANT_GWAS_sample)
dim(BMI_GIANT_GWAS_sample)

## ------------------------------------------------------------------------
table(BMI_GIANT_GWAS_sample$Freq_MAF_Int_Hapmap)

## ------------------------------------------------------------------------
X <- model.matrix(~ splines::ns(N,5) + Freq_MAF_Int_Hapmap, data = BMI_GIANT_GWAS_sample)[,-1]
head(X)

## ------------------------------------------------------------------------
pi0x <- lm_pi0(pValues=BMI_GIANT_GWAS_sample$p, 
               X=X, smooth.df=3)
names(pi0x)

## ------------------------------------------------------------------------
BMI_GIANT_GWAS_sample$fitted0.8 <- pi0x$pi0.lambda[,round(pi0x$lambda,2)==0.8]
BMI_GIANT_GWAS_sample$fitted0.9 <- pi0x$pi0.lambda[,round(pi0x$lambda,2)==0.9]
BMI_GIANT_GWAS_sample$fitted.final.smooth <- pi0x$pi0

## ------------------------------------------------------------------------
ldf <- reshape2::melt(BMI_GIANT_GWAS_sample,
                      id.vars=colnames(BMI_GIANT_GWAS_sample)[-grep("fitted",
                                                                    colnames(BMI_GIANT_GWAS_sample))],
                      value.name = "pi0",variable.name = "lambda")
ldf$lambda <- as.character(ldf$lambda)
ldf$lambda[ldf$lambda=="fitted0.8"] <- "lambda=0.8"
ldf$lambda[ldf$lambda=="fitted0.9"] <- "lambda=0.9"
ldf$lambda[ldf$lambda=="fitted.final.smooth"] <- "final smoothed pi0(x)"

head(ldf)

## ---- BMI_GWAS_plot------------------------------------------------------
library(ggplot2)
ggplot(ldf, aes(x=N, y=pi0))+
  geom_line(aes(col=Freq_MAF_Int_Hapmap, linetype=lambda)) +
  ylab("Estimated proportion of nulls") +
  guides(color=guide_legend(title="MAF in HapMap CEU population"),
         linetype=guide_legend(title="Estimate"))


