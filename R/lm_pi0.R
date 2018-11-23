# Computation of pi0 from pvalues and matrix of covariates
#
# Notes:
#  - a previous version (v1.2.1) of this function is available in a sister file
#
# Changes in lm_pi0 w.r.t version 1.2.1:
#  - added argument value matching for "type"
#  - added checks on arguments lambda, X, smooth.df
#  - moved fitting pi0.lambda to separate functions
#  - removed some import decorators from documentation 
#  - added some importFrom decorators
#
# Additional changes in lm_pi0_2
#  - removed argument smooth.df and uses df=3 always
#  - removed argument threshold (always force probabilities into [0,1])
#  - pi0, pi0.lambda, as well as pi0.smooth are now truncated into [0,1]
#  - shifted final estimate to using limit lambda=1 (instead of largest lambda)
#  - removed pi0.smooth from output list


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
lm_pi0 <- function(p, lambda = seq(0.05, 0.95, 0.05), X,
                   type=c("logistic", "linear"), smooth.df=3, threshold=TRUE) {
  
  # check validity of inputs
  type <- match.arg(type)
  prange <- check_p(p)
  lambda <- check_lambda(lambda, prange[2])
  n.lambda <- length(lambda)
  smooth.df <- check_df(smooth.df, n.lambda)
  X <- check_X(X, p)
  
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



#' Another version of lm_pi0
#'
#' @param p numeric vector
#' @param lambda numeric vector
#' @param X matrix of covariates
#' @param type character
#' @param smooth character
#'
#' @return list with pi0, pi0.lambda, lambda
lm_pi0_2 <- function(p, lambda = seq(0.05, 0.95, 0.05), X,
                     type=c("logistic", "linear")) {
  
  # check validity of inputs
  type <- match.arg(type)
  prange <- check_p(p)
  lambda <- check_lambda(lambda, prange[2])
  n.lambda <- length(lambda)
  X <- check_X(X, p)
    
  # pick a modeling function, fit odels for each lambda
  available.functions <- list(logistic=fit_logistic, linear=fit_linear)
  fit.function <- available.functions[[type]]  
  pi0.lambda <- matrix(NA, nrow=nrow(X), ncol=n.lambda)
  for (i in 1:n.lambda) {
    y <- (p > lambda[i])
    pi0.lambda[, i] <- fit.function(y, X)/(1-lambda[i])
  }
  pi0.lambda <- force.unit.interval(pi0.lambda)
  
  # smooth over values of lambda (for each p-value/ row in X)
  # (instead of taking limit lambda->1, use largest available lambda)
  smooth.large.lambda = function(y) {
    spline <- smooth.spline(lambda, y,
                            df=3, keep.data=FALSE, tol=1e-7)
    spline$y[n.lambda]
  }
  pi0 <- force.unit.interval(apply(pi0.lambda, 1, smooth.large.lambda))
  
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
# Input validation checks


# verify quality of pvalues
#
# @param p numeric vector of p-values
#
# @return numeric vector of length 2 with (min(p), max(p))
check_p <- function(p) {
  if (class(p) != "numeric") {
    stop("p must be a numeric vector\n", call. = FALSE)
  }
  prange <- range(p, na.rm=TRUE)
  if (prange[1] < 0 | prange[2] > 1) {
    stop("p must be numeric in range [0, 1]\n", call. = FALSE)
  }
  if (any(is.na(p))) {
    stop("p must not have missinge values\n", call. = FALSE)
  }
  prange
}


# verify the suitability of lambda
#
# @param x vector of lambda values
# @param pmax numeric, maximal pvalue
#
# @return numeric vector of sorted unique x, or stop if x does not satisfy criteria
check_lambda <- function(x, pmax) {
  if (class(x)!="numeric") {
    stop("lambda must be a numeric vector \n", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("lambda must not contain NAs, NULL, or non-finite elements\n", call. = FALSE)
  }
  x <- sort(unique(x))
  if (length(x)<4) {
    stop("lambda must be of length >=4\n", call. = FALSE)
  }
  if (min(x)<=0 | max(x)>=1) {
    stop("lambda values must all be in range (0, 1)\n", call. = FALSE)
  }
  if (pmax < max(x)) {
    warning("maximal p is smaller than maximal lambda", call. = FALSE)
  }
  x
}


# verify the suitability of a value as a number of degrees of freedrom
#
# @param x expect a single number
# @param max.value numeric, maximal value allowed for x
#
# @return integer derived from x
check_df <- function(x, max.value) {
  if (class(x) != "numeric" & class(x) != "integer") {
    stop("df must be a number")
  }
  if (length(x) != 1 | any(!is.finite(x))) {
    stop("df must be a single finite number", call. = FALSE)
  }
  x <- round(x)
  if (x <= 1 | x > max.value) {
    stop("df must be in range 1 < df < length(lambda)\n", call. = FALSE)
  }
  x
}


# verify that X is a matrix of covariates compatible with a vector of pvalues
#
# @param X vector or matrix of covariates
# @param pvalues vector of values
#
# @return matrix
check_X <- function(X, pvalues) {
  # allow for null input (no covariates)
  if (is.null(X)) {
    X <- cbind(rep(1, length(pvalues)))
    rownames(X) <- names(pvalues)
  }
  # allow for a single covariates specified as a vector
  if (is.null(dim(X))) {
    X <- cbind(X)
  }
  # ensure that X and pvalues are compatible
  if (length(pvalues)!=nrow(X)) {
    stop("incompatible X and p - different lengths\n", call. = FALSE)
  }
  if (class(X) != "matrix") {
    warning(paste0("coercing X info a matrix from a ", class(X)), call.=FALSE)
    X <- as.matrix(X)
  }
  if (!identical(rownames(X), names(pvalues))) {
    warning("X and p have different names", call. = FALSE)
  }
  # ensure that all columns in X are numeric
  if (!all(apply(X, 2, class) %in% c("numeric", "integer", "factor"))) {
    stop("X must be a numeric vector or numeric matrix\n", call. = FALSE)
  }
  X
}


# helper to force a set of values in a unit interval
#
# (This is a simpler implementation than ifelse)
#
# @param x numeric vector or matrix
#
# @return same object like x, with values truncated by [0,1]
force.unit.interval <- function(x) {
  x[x<0] <- 0
  x[x>1] <- 1
  x
}



