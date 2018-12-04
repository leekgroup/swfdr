#  Fast calculation of smoothing curves using splines
# 
#  This is a simplified version of function 'smooth.spline' from the
#  stats package (file smspline.R).
#
#  It is simplified by removing many tuning parameters, validation
#  checks on the input, and components from the output. Because the
#  validation checks in the original function included calculations
#  O(N) and O(NlogN) steps, this stripped-down version can run
#  substantially faster in some cases. However, this function now
#  relies on external validations of the inputs and can fail abruptly.
#  For this reason, this function should be kept internal to the package
#  where it can be called with safe inputs.
#
#  Copied and modified from file src/library/stats/R/smspline.R
#  Part of the R package, https://www.R-project.org
#
#
#
#  License notice from original file:
#
#  Copyright (C) 1995-2016 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/



#' Fast smoothing splines
#'
#' Smoothing is computed using the same code as smooth.spline from the stats
#' package. This version can run substantially faster because it does not
#' perform validation checks on the inputs. The inputs must be prepared
#' outside of the function. This setup can be a benefit if smoothing
#' is performed many times with the same values of x.
#'
#' This version also supports a very limited number of arguments. The
#' behavior of this function follows the default settings for smooth.spline.
#'
#' Note: this function uses a function internal to the stats package to
#' perform its work. In principle, this means fast.spline can break
#' if that function from stats changes. However, that function from stats
#' has remained fixed for many years (now 2018) and is thus not likely to
#' change. If that happens, the results from fast.spline can be reproduced
#' via 'smooth.spline(x, y, df=3, tol=4.5e-7)$y'
#'
#' @keywords internal
#' @param x numeric vector, must have sorted, unique, finite values
#' @param y numeric vector, must have finite values
#' @param df integer, number of degrees of freedom
#' @param tol numeric, tolerance that determines when smoothing stops.
#'
#' @return numeric vector with smoothed y values
fast.spline <- function(x, y = NULL, df = 3, tol = 4.5e-7) {
  
  ctrl.Num <- c(low=-1.5, high=1.5, tol=1e-4, eps=2e-8, ratio=-1.0)
  nx <- length(x)
  wbar <- rep(1, nx)
  r.ux <- x[nx] - x[1L]
  xbar <- (x - x[1L])/r.ux
  nknots <- nx
  knot <- c(rep(xbar[1], 3), xbar, rep(xbar[nx], 3))
  nk <- nknots + 2L   
  ispar <- 0L
  spar <- 0  
  icrit <- 3L
  dofoff <- df
  iparms <- c(icrit=icrit, ispar=ispar, iter = as.integer(500), FALSE)

  ## warnings might appear because this part calls C_rbart, which is
  ## an internal function within the stats package that is not usually
  ## exported, i.e. not made available for use outside the stats package. 
  suppressWarnings(.Fortran(stats:::C_rbart,
                            as.double(1.0),
                            as.double(dofoff),
                            x = as.double(xbar),
                            y = as.double(y),
                            w = as.double(wbar), 
                            ssw = as.double(0),
                            as.integer(nx),
                            as.double(knot),
                            as.integer(nk),
                            coef = double(nk),
                            ty = double(nx),
                            lev = double(nx),
                            crit = double(1),
                            iparms = iparms,
                            spar = spar,
                            parms = ctrl.Num, 
                            scratch = double((17L+1L) * nk + 1L),
                            ld4  = 4L,
                            ldnk = 1L,
                            ier = integer(1L)
                            )[["ty"]])
}

