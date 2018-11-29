# Computation of pi0 from pvalues and matrix of covariates
#
# Notes:
#  - a previous version (v1.2.1) of this function is available in a sister file
#
# Changes in lm_pi0_1.3 w.r.t version 1.2.1:
#  - added argument value matching for "type"
#  - added checks on arguments lambda, X, smooth.df
#  - moved fitting pi0.lambda to separate functions
#  - removed some import decorators from documentation 
#  - added some importFrom decorators
#  - allowed X to be missing (uses a constant/uniformative covariate)
#
# Additional changes in lm_pi0
#  - removed argument smooth.df and uses df=3 always
#  - removed argument threshold (always force probabilities into [0,1])
#  - pi0, pi0.lambda, as well as pi0.smooth are now truncated into [0,1]
#  - removed pi0.smooth from output list
#  - fitted splines using a customized function fast.spline instead of smooth.spline


#' Estimate pi0(x)
#' 
#' @param p Numerical vector of p-values
#' @param lambda Numerical vector of thresholds. Must be in [0,1).
#' @param X Design matrix (one test per row, one variable per column). Do not include the intercept.
#' @param type Type of regression, "logistic" or "linear." Default is logistic.
#' @param smooth.df Number of degrees of freedom when estimating pi0(x) with a smoother.
#' 
#' @return pi0 Numerical vector of smoothed estimate of pi0(x). The length is the number of rows in X.
#' @return pi0.lambda Numerical matrix of estimated pi0(x) for each value of lambda. The number of columns is the number of tests, the number of rows is the length of lambda.
#' @return lambda Vector of the values of lambda used in calculating pi0.lambda
#' @return pi0.smooth Matrix of fitted values from the smoother fit to the pi0(x) estimates at each value of lambda (same number of rows and columns as pi0.lambda)
#'
#' @importFrom stats binomial glm
#'
#' @examples
#' X <- seq(-1,2,length=1000) ##covariate
#' pi0 <- 1/4*X + 1/2 ##probability of being null
#' nullI <- rbinom(1000,prob=pi0,size=1)> 0 ##generate null/alternative p-values
#' pValues <- rep(NA,1000) ##vector of p-values
#' pValues[nullI] <- runif(sum(nullI)) ##null from U(0,1)
#' pValues[!nullI] <- rbeta(sum(!nullI),1,2) ##alternative from Beta
#' pi0x <- lm_pi0(pValues, X=X, smooth.df=3)
#'
#' @export
lm_pi0_1.3 <- function(p, lambda = seq(0.05, 0.95, 0.05), X,
                   type=c("logistic", "linear"), smooth.df=3, threshold=TRUE) {
  
  # check validity of inputs
  type <- match.arg(type)
  prange <- check_p(p)
  lambda <- check_lambda(lambda, prange[2])
  n.lambda <- length(lambda)
  smooth.df <- check_df(smooth.df, n.lambda)
  X <- check_X(X=X, p=p)
  
  # pick a modeling function, fit odels for each lambda
  available.functions <- list(logistic=fit_logistic, linear=fit_linear)
  fit.function <- available.functions[[type]]  
  pi0.lambda <- matrix(NA, nrow=nrow(X), ncol=n.lambda)
  for (i in 1:n.lambda) {
    y <- (p > lambda[i])
    pi0.lambda[, i] <- fit.function(y, X)/(1-lambda[i])
  }
  if (threshold) {
    pi0.lambda <- force.unit.interval(pi0.lambda)
  }
  
  # smooth over values of lambda (for each p-value/ row in X)
  smooth.function = function(y) {
    spline <- smooth.spline(lambda, y, df=smooth.df,
                           tol=1e-7, keep.data=FALSE)
    spline$y
  }
  pi0.smooth = t(apply(pi0.lambda, 1, smooth.function))
  if (threshold) {
    pi0.smooth = force.unit.interval(pi0.smooth)
  }
  
  # get final estimate of pi0
  # (instead of taking limit lambda->1, use largest available lambda)
  pi0 <- force.unit.interval(pi0.smooth[, n.lambda])
  
  list(pi0=pi0, pi0.lambda=pi0.lambda, lambda=lambda, pi0.smooth=pi0.smooth)
}



#' Estimation of pi0, proportion of p-values consistent with a null hypothesis
#'
#' @param p numeric vector, p-values
#' @param lambda numeric vector, thresholds used to bin pvalues, must be in [0,1).
#' @param X numeric matrix, covariates that might be related to p values
#' (one test per row, one variable per column). 
#' @param type character, type of regression used to fit features to pvalues
#'
#' @return pi0 numerical vector of smoothed estimate of pi0(x).
#' The length is the number of rows in X.
#' @return pi0.lambda numeric matrix of estimated pi0(x) for each value of lambda.
#' The number of columns is the number of tests, the number of rows is the length of lambda.
#' @return lambda numeric vector of the thresholds used in calculating pi0.lambda
#' 
#' @importFrom stats binomial glm
#'
#' @examples
#' X <- seq(-1,2,length=1000) ##covariate
#' pi0 <- 1/4*X + 1/2 ##probability of being null
#' nullI <- rbinom(1000,prob=pi0,size=1)> 0 ##generate null/alternative p-values
#' pValues <- rep(NA,1000) ##vector of p-values
#' pValues[nullI] <- runif(sum(nullI)) ##null from U(0,1)
#' pValues[!nullI] <- rbeta(sum(!nullI),1,2) ##alternative from Beta
#' pi0x <- lm_pi0(pValues, X=X, smooth.df=3)
#'
#' 
#' @return list with pi0, pi0.lambda, lambda
#'
#' @export
lm_pi0 <- function(p, lambda = seq(0.05, 0.95, 0.05), X,
                   type=c("logistic", "linear")) {
  
  # check validity of inputs
  type <- match.arg(type)
  prange <- check_p(p)
  lambda <- check_lambda(lambda, prange[2])
  n.lambda <- length(lambda)
  X <- check_X(X=X, p=p)
  
  # pick a modeling function, fit odels for each lambda
  available.functions <- list(logistic=fit_logistic, linear=fit_linear)
  fit.function <- available.functions[[type]]  
  pi0.lambda <- matrix(NA, nrow=nrow(X), ncol=n.lambda)
  for (i in 1:n.lambda) {
    y <- (p > lambda[i])
    pi0.lambda[, i] <- fit.function(y, X)/(1-lambda[i])
  }
  pi0.lambda <- regularize.interval(pi0.lambda)
  
  # smooth over values of lambda (for each p-value/ row in X)
  # (instead of taking limit lambda->1, use largest available lambda)
  smooth.large.lambda = function(y) {
    fast.spline(lambda, y)[n.lambda]
  }
  pi0 <- regularize.interval(apply(pi0.lambda, 1, smooth.large.lambda))
  
  list(pi0=pi0, pi0.lambda=pi0.lambda, lambda=lambda)
}





# #############################################################################
# modeling of a response vector with covariates
#
# all functions take as input:
#
# @param y vector of response values
# @param X matrix of covariates with nrow(X) = length(y)
#
# @return numeric vector of length(y)


# fit response values using a binomial/logit model
fit_logistic <- function(y, X) {
  glm(y ~ X, family=binomial)$fitted.values
}


# fit response values using a linear regression
# 
# (This could be implemnted via glm(... family=gaussian) but the
# implementation with lsfit is much faster
fit_linear <- function(y, X) {
  regFit <- lsfit(X, y)$coefficients
  regFit[1] + (X %*% matrix(regFit[-1], ncol=1))
}




# #############################################################################
# other helper functions


# helper to force a set of values in a regular interval
#
# (This is a simpler implementation than ifelse)
#
# @param x numeric vector or matrix
#
# @return same object like x, with values truncated by [0,1]
regularize.interval <- function(x, interval = c(0, 1)) {
  x[x < interval[1]] <- interval[1]
  x[x > interval[2]] <- interval[2]
  x
}



