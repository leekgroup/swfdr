# Computation of local false discovery rate
#


#' Estimation of local false discovery rate
#'
#' Local FDR is a comparison of a p-value distribution with a flat/uniform distribution.
#' This first half of this implementation is similar to the code for lfdr in package qvalue.
#' The second half deviates from qvalue::lfdr in that it avoids comparing the right tail
#' of the p-values with dnorm. The overall effect is to avoid producing pathologically low
#' estimates of lfdr for p-values close to unity. This pathological behavior is avoided
#' in qvalue::lfdr by using argument monotone=TRUE; here, the monotone transformation is
#' omitted to allow adjustment with a vector of pi0, and thus the pathological behavior
#' is handled through the block containing dnorm_flat.
#'
#' @param p numeric vector of p-values
#' @param pi0 numeric vector of pi0 estimates (same length as p)
#' @param eps numeric, regularization of p-values, same effect as in qvalue::lfdr
#' @param adj numeric, adjustment of smoothing bandwidth in density estimte, same effect
#' as in qvalue::lfdr
#'
#' @return numeric vector of local fdr values, truncated into interval [0, 1]
lfdr <- function(p, pi0, eps=1e-8, adj=1.5) {

  # estimation of density of empirical p-values in units of the normal distribution
  p <- pmin(pmax(p, eps), 1-eps)
  x <- qnorm(p)
  xd <- density(x, adjust = adj)
  xs <- smooth.spline(x = xd$x, y = xd$y)
  y <- predict(xs, x)$y
  
  # The calculation in qvalue::lfdr compares y to dnorm(x) 
  # Here, use an adjusted function that does not decrease when x>0
  dnorm_flat = function(x) {
    result = dnorm(x)
    result[x>0] = dnorm(0)
    result
  }
  
  result <- dnorm_flat(x) / y
  pmin(pi0*result, 1)
}

