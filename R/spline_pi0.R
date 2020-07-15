# Estimation of a pi0 curves using splines 
# This is used in lm_pi0 to estimate pi0 based on estimates at various lambda


#' Fit smoothing splines to several curves 
#'
#' This implementation uses smooth.spline. It gives a slight optimization
#' by not reporting original data (keep.data=FALSE) and by setting numeric
#' tolerance to a fixed number (tol=1e-7). In the latter, the precise number
#' does not matter, but it is important that it is not computed from x
#' at each iteration in the loop.
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector, position of knots
#' @param ymat numeric matrix, ncol(ymat) should match length(x)
#' @param df integer, degrees of freedom for the spline
#'
#' @return list with components pi0.smooth and pi0
smooth.spline.pi0 <- function(x, ymat, df=3) {
  nx <- length(x)
  pi0.smooth = matrix(NA, nrow=nrow(ymat), ncol=nx)
  pi0 <- rep(NA, nrow(ymat))
  for (i in seq_len(nrow(ymat))) {
    yfit <- smooth.spline(x, ymat[i,], df=df, tol=1e-7, keep.data=FALSE)$y
    pi0.smooth[i,] = yfit
    pi0[i] <- yfit[nx]
  }
  list(pi0.smooth=pi0.smooth, pi0=pi0)
}


#' Fit smoothing unit splines to several curves
#'
#' In contrast to smooth.spline.pi0, this function uses the bs() bases
#' and uses boundary knots at x=(0, 1)
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector, position of knots
#' @param ymat numeric matrix, ncol(ymat) should match length(x)
#' @param df integer, degrees of freedom for the spline
#' 
#' @return list with component pi0 (numeric vector)
unit.spline.pi0 <- function(x, ymat, df=3) {
  transform.matrix <- unit.spline.matrix(x, df=df)
  ## make prediction only from the last row of the transformation matrix
  transform.last <- transform.matrix[nrow(transform.matrix), ]
  pi0 <- apply(ymat, 1, function(z) { sum(transform.last * z) })
  list(pi0=pi0)
}


#' Construct a transformation matrix for a unit spline
#'
#' This function uses library::bs to construct a b-spline basis
#' from knots, using [0,1] as the boundary knots.
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector, position of internal knots within [0,1]
#' @param df integer, degrees of freedom for the spline
#' 
#' @importFrom splines bs
#'
#' @return matrix of dimension c(length(x), length(x))
unit.spline.matrix <- function(x, df=3) {
  # construct a splines::bs basis from x
  x.bs = bs(x, knots=x, Boundary.knots=c(0, 1))
  x.bsTbs = t(x.bs) %*% x.bs
  n.bs = ncol(x.bs)
  n.x = length(x)
  
  # construct a smoothing component
  omega.sqrt = matrix(0, ncol=n.bs, nrow=n.bs)
  for (i in seq(3, n.bs)) {
    omega.sqrt[i, c(i-2, i-1, i)] = c(1,-2,1)
  }
  omega = t(omega.sqrt) %*% omega.sqrt
  
  # optimize smoothing weight to get number of degrees of freedom right
  multiplier = optimize.multiplier(x.bsTbs, omega, target=2+df)
  
  # output a transformation matrix
  x.bs %*% solve(x.bsTbs + multiplier*omega) %*% t(x.bs)
}


#' Fit smoothing unit spline to one curve
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector
#' @param y numeric vector
#'
#' @return vector of same length as y
unit.spline <- function(x, y, df=3) {
  transform.matrix = unit.spline.matrix(x, df=df)
  as.numeric(transform.matrix %*% cbind(y))
}


#' recursive internal function to find an optimal value of lambda/multiplier
#' so that trace(X+lambda*omega) = target
#'
#' This implementation relies on the knowledge that a larger lagrange
#' multiplier leads to lower trace(X+lambda*omega).
#'
#' The implementation starts in an interval (0, 1, 2), which are lower,
#' middle, and upper bounds, and then expands the interval or zooms in to
#' find a reasonable multiplier
#'
#' @keywords internal
#' @noRd
#' @param X matrix
#' @param omega matrix
#' @param target numeric, target for trace(X+lambda*omega)
#' @param interval numeric vector of length 2, current range for lambda
#' (internal use)
#' @param values numeric vector of length 2, current target estimates
#' (internal use)
#' @param tol numeric, numerical tolerance, does not need to be very small
#'
#' @return numeric, lagrange multiplier that brings criterion close to its
#' target
optimize.multiplier <- function(X, omega, target=5,
                            interval=c(1e-6, 1, 2), values=c(NA, NA, NA),
                            tol=1e-4) {
  # exit early if the interval is very narrow
  interval.width <- interval[3]-interval[1]
  if (interval.width/interval[2] < tol) {
    return(interval[2])
  }
  # estimate values for the target based on the test interval range
  for (i in c(1,2,3)) {
    if (is.na(values[i])) {
      values[i] <- sum(diag(solve(X + interval[i]*omega)))
    }
  }
  # exit if satisfactory solution is found
  if (abs(values[1] - values[3]) < tol) {
    return(mean(interval))
  }
  # adjust the interval range and look further
  if (target < values[3]) {
    interval.new <- interval[3] + c(0, 2*interval.width, 4*interval.width)
    values.new <- c(values[2], NA, NA)
  } else {
    if (target < values[2]) {
      interval.new <- c(interval[2], mean(interval[2:3]), interval[3])
      values.new <- c(values[2], NA, values[3])      
    } else {
      interval.new <- c(interval[1], mean(interval[1:2]), interval[2])
      values.new <- c(values[1], NA, values[2])
    }
  }
  optimize.multiplier(X, omega, target, interval.new, values.new)
}

