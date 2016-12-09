#' Estimate pi0(x)
#' 
#' @param pValues Numerical vector of p-values
#' @param lambda Numerical vector of thresholds. Must be in [0,1).
#' @param X Design matrix (one test per row, one variable per column). Do not include the intercept.
#' @param smooth.df Number of degrees of freedom when estimating pi0(x) with a smoother.
#' @param threshold If TRUE (default), all estimates are thresholded at 0 and 1, if FALSE, none of them are.
#' 
#' @return pi0 Numerical vector of smoothed estimate of pi0(x). The length is the number of rows in X.
#' @return pi0.lambda Numerical matrix of estimated pi0(x) for each value of lambda. The number of columns is the number of tests, the number of rows is the length of lambda.
#' @return lambda Vector of the values of lambda used in calculating pi0.lambda
#' @return pi0.smooth Matrix of fitted values from the smoother fit to the pi0(x) estimates at each value of lambda (same number of rows and columns as pi0.lambda)
#'
#' @export
lm_pi0 <- function(pValues, lambda = seq(0.05, 0.95, 0.05), X, smooth.df=3, threshold=TRUE)
{
  ##if X is a vector, change it into a matrix
  if(is.null(dim(X)))
  {
    X <- matrix(X, ncol=1)
  }
  
  ##number of tests
  n <- nrow(X)
  ##number of lambdas
  nLambda <- length(lambda)
  
  ##sort lambdas from smallest to largest and take only unique values
  lambda <- sort(unique(lambda))
  
  ##make a design matrix with the intercept
  Xint <- cbind(1, X)
  
  ##get the estimate for each value of lambda 
  pi0.lambda <- matrix(NA, nrow=n, ncol=nLambda)
  for(i in 1:nLambda)
  {
    lambda.i <- lambda[i]
    y <- pValues > lambda.i
    
    ##fit regression
    regFit <- lsfit(X, y)
    
    ##get the estimated values of pi0
    pi0.lambda[,i] <- (Xint %*% matrix(regFit$coefficients, ncol=1))[,1]/(1-lambda.i)
    
    if(threshold){
      pi0.lambda[,i] <- ifelse(pi0.lambda[,i] > 1, 1, pi0.lambda[,i])
      pi0.lambda[,i] <- ifelse(pi0.lambda[,i] < 0, 0, pi0.lambda[,i])
    }
  }
  
  ##smooth over values of lambda (do this for each test in part)
  pi0.smooth <- matrix(NA, nrow=n, ncol=nLambda)
  ##also save final estimate (maximum of 0 and minimum of 1 and smoothed value at largest lambda)
  pi0 <- rep(NA, length=n)
  for(i in 1:n)
  {
    if(i %% 10000==0)
    {
      message(paste("At test #:",i))
    }
    spi0 <- smooth.spline(lambda, pi0.lambda[i,], df=smooth.df)
    pi0.smooth[i, ] <- spi0$y
    pi0[i] <- pi0.smooth[i,nLambda]
  }
  
  if(threshold){ 
    pi0 <- ifelse(pi0 > 1, 1, pi0)
    pi0 <- ifelse(pi0 < 0, 0, pi0)
  }
  
  return(list(pi0=pi0, pi0.lambda=pi0.lambda, lambda=lambda, pi0.smooth=pi0.smooth))
}

