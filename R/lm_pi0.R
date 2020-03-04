# Computation of pi0 from pvalues and matrix of covariates
#


#' Estimation of pi0, proportion of p-values consistent with a null hypothesis
#'
#' @param p numeric vector, p-values
#' @param lambda numeric vector, thresholds used to bin pvalues, must be in [0,1).
#' @param X numeric matrix, covariates that might be related to p values
#' (one test per row, one variable per column). 
#' @param type character, type of regression used to fit features to pvalues
#' @param smooth.df Number of degrees of freedom when estimating pi0(x) with a smoother.
#' @param threshold logical, if TRUE, all estimates are thresholded into unit interval;
#' if FALSE, all estimates are left as they are computed
#' @param smoothing character, type of smoothing used to fit pi0
#'
#' @return pi0 numerical vector of smoothed estimate of pi0(x).
#' The length is the number of rows in X.
#' @return pi0.lambda numeric matrix of estimated pi0(x) for each value of lambda.
#' The number of columns is the number of tests, the number of rows is the length of lambda.
#' @return lambda numeric vector of the thresholds used in calculating pi0.lambda
#' @return pi0.smooth (only output with smoothing="smooth.spline") Matrix of fitted values
#' from the smoother fit to the pi0(x) estimates at each value of lambda (same number of
#' rows and columns as pi0.lambda)
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
#' pi0x <- lm_pi0(pValues, X=X)
#'
#' @export
lm_pi0 <- function(p, lambda = seq(0.05, 0.95, 0.05), X,
                   type=c("logistic", "linear"), smooth.df=3, threshold=TRUE,
                   smoothing=c("unit.spline", "smooth.spline")) {
  
  # check validity of inputs
  type <- match.arg(type)
  smoothing <- match.arg(smoothing)
  prange <- check_p(p)
  lambda <- check_lambda(lambda, prange[2])
  n.lambda <- length(lambda)
  smooth.df <- check_df(smooth.df, n.lambda)
  X <- check_X(X=X, p=p)
  
  # pick a modeling function, fit odels for each lambda
  available.regressions <- list(logistic=fit_logistic, linear=fit_linear)
  fit.function <- available.regressions[[type]]
  pi0.lambda <- matrix(NA, nrow=nrow(X), ncol=n.lambda)
  for (i in 1:n.lambda) {
    y <- (p >= lambda[i])
    pi0.lambda[, i] <- fit.function(y, X)/(1-lambda[i])
  }
  if (threshold) {
    pi0.lambda <- regularize.interval(pi0.lambda)
  }
  
  # smooth over values of lambda (for each p-value/ row in X)
  # (instead of taking limit lambda->1, use largest available lambda)
  available.smoothings <- list(smooth.spline=smooth.spline.pi0,
                               unit.spline=unit.spline.pi0)
  smoothing.function <- available.smoothings[[smoothing]]
  pi0 <- smoothing.function(lambda, pi0.lambda, smooth.df)
  if (threshold) {
    pi0 <- regularize.interval(pi0)
  }
  
  result <- c(list(call=match.call(), lambda=lambda, X.names = colnames(X),
                 pi0.lambda=pi0.lambda), pi0)
  
  class(result) <- "lm_pi0"
  result
}





# #############################################################################
# modeling of a response vector with covariates


#' Fit response values using a binomial/logit model
#'
#' @keywords internal
#' @noRd
#' @param y numeric vector
#' @param X numeric matrix (covariates)
#'
#' @return numeric vector
#'
#' @importFrom stats binomial glm
fit_logistic <- function(y, X) {
  glm(y ~ X, family=binomial)$fitted.values
}


#' Fit response values using a linear regression
#' 
#' (This could be implemnted via glm(... family=gaussian) but the
#' implementation with lsfit is much faster
#'
#' @keywords internal
#' @noRd
#' @param y numeric vector
#' @param X numeric matrix (covariates)
#'
#' @return numeric vector
#'
#' @importFrom stats lsfit
fit_linear <- function(y, X) {
  regFit <- lsfit(X, y)$coefficients
  regFit[1] + (X %*% matrix(regFit[-1], ncol=1))
}




# #############################################################################
# other helper functions


#' Force a set of values in a regular interval
#'
#' (This is a simpler/faster implementation than with ifelse)
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector or matrix
#' @param interval numeric vector of length 2 with min/max values
#'
#' @return same object like x, with values truncated by [0,1]
regularize.interval <- function(x, interval = c(0, 1)) {
  if (is(x, "list")) {
    return (lapply(x, regularize.interval, interval=interval))
  }
  x[x < interval[1]] <- interval[1]
  x[x > interval[2]] <- interval[2]
  x
}

