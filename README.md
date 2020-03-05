# swfdr
Science-wide false discovery rate

`swfdr` is an R package for dealing with multiple hypothesis testing.



## Installation

The package is available through Bioconductor.

```
library(BiocManager)
BiocManager::install("swfdr")
```

It is also possible to install in-development versions from github.

```
library(devtools)
install_github("leekgroup/swfdr")
```




## Example

A minimal working example for using the `lm_pi0` and `lm_qvalue` functions
is as follows

```
library(swfdr)
x <- c(rbeta(1000, 0.2, 1), runif(1000, 0, 1))
x.covariate = rep(c(1,2), each=1000)
# compute covariate-adjusted pi0 - probability of null hypothesis to be true
x.pi0 = lm_pi0(x, X=x.covariate)
x.pi0
# compute covariate-adjusted qvalues
x.qvalues <- lm_qvalue(x, X=x.covariate)
x.qvalues
```

