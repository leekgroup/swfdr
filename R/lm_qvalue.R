# Esimation of qvalues from pvalues and matrix of covariates


#' Compute qvalues taking into account a matrix of covariates
#'
#' The recipe for turning pvalues into qvalues is adapted from package
#' 'qvalue' and articles by Storey, Tibshirani, Taylor, Siegmund.
#' 
#'
#' @param p numeric vector of p-values
#' @param X matrix of covariates (can be missing if pi0 is specified instead)
#' @param pfdr logical
#' @param pi0 list with pi0 estimates from lm_pi0
#' @param ... other parameters (passed on to lm_pi0 if pi0 is not provided)
#' 
#' @return list
lm_qvalue <- function(p, X, pfdr=FALSE, pi0=NULL, ...) {
  
  # check inputs
  prange <- check_p(p)
  
  # pi0 is required - compute it if not provided
  if (is.null(pi0)) {
    pi0 <- lm_pi0(p, X=X, ...)
  }
  
  # block to compute qvalues (adapted from package qvalue)
  # This is almost the same as pi0s*p.adjust(p, method="BH")
  # However, this implementation allows setting pfdr=TRUE similarly as in qvalue()
  n <- length(p)
  i <- n:1
  o <- order(p, decreasing=TRUE)
  ro = order(o, decreasing=FALSE)
  if (pfdr) {
    q <- pmin(1, cummin( p[o]*n/ (i*(1- (1-p[o])^n))))
  } else {
    q <- pmin(1, cummin( p[o]*n/ i ))
  }
  q <- (q* pi0$pi0[o])[ro]
  
  # create output 
  result <- c(list(call=match.call()), pi0, list(pvalues = p, qvalues = q))
  class(result) <- c("lm_qvalue")
  result
}

