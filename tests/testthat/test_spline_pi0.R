## Tests for estimating pi0 using smoothing splines

cat("\ntest_spline_pi0.R\n")


# generate a noisy curve with y values around ~1 (ballpark scale)
set.seed(12345)
xx <- seq(0.05, 0.95, by=0.05)
xx.len <- length(xx)
yy <- ((xx-0.4)^2) + runif(length(xx), 0, 0.1)
yy <- (yy - mean(yy)) / sd(yy)

# dataset with multiple curves:
#  yy defined in header
#  yy with flipped sign
#  flat curve at 0
#  flat curve at 1
#  non-noisy curve of slope 1
yymat <- rbind(yy,
               -yy,
               rep(0, xx.len),
               rep(1, xx.len),
               seq(1, xx.len))
rownames(yymat) <- NULL


# fit the curve using R's smooth.spline
yy.smooth <- smooth.spline(xx, yy, df=3)$y




# #############################################################################
# Last-point estimate using smooth.spline


test_that("smooth.spline.last with matrix (1 row)", {
  result <- smooth.spline.last(xx, matrix(yy, nrow=1))
  expect_equal(result, yy.smooth[length(xx)], tol=1e-6)
})


test_that("smooth.spline.last with matrix (multiple rows)", {
  result <- smooth.spline.last(xx, yymat)
  smooth.last = yy.smooth[xx.len]
  expect_equal(result, c(smooth.last, -smooth.last, 0, 1, xx.len), tol=0.1)
})



# #############################################################################
# Last-point estimate using b-splines


test_that("entire unit.spline is similar to smooth.spline", {
  # difference between smooth spline and data
  yy.smooth.maxerr <- max(abs(yy.smooth-yy))
  # difference between unit spline and data
  yy.bs <- unit.spline(xx, yy)
  yy.bs.maxerr <- max(abs(yy.bs-yy))
  # the errors should be comparable between two approaches
  expect_lt(abs(yy.bs.maxerr-yy.smooth.maxerr)/yy.smooth.maxerr, 0.1)
})


test_that("can estimate last point on unit.spline", {
  result <- unit.spline.last(xx, matrix(yy, nrow=1))
  # just check ballpark estimate
  expect_equal(result, yy.smooth[length(xx)], tol=0.1)
})


test_that("can estimate last point on unit.spline", {
  result <- unit.spline.last(xx, yymat)
  smooth.last = yy.smooth[xx.len]
  ## first two items in result are based on a noisy curve ->  tolerant comparison
  expect_equal(result[1:2], smooth.last * c(1, -1), tol=1e-1)
  # last three items in result are based on non-noisy curve -> stricter
  expect_equal(result[3:5], c(0, 1, xx.len), tol=1e-2)
})

