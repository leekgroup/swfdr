## Tests for R/fastspline.R

cat("\ntest_fastspline.R\n")


# generate a noisy curve
xx = 1:19
yy = ((xx-10)^2)/20 + runif(length(xx), 0, 4)


test_that("fast.spline produces similar results to smooth.spline", {
  s1 = smooth.spline(xx, yy, df=3)$y
  s2 = fast.spline(xx, yy)
  expect_equal(s1, s2, tol=1e-6)
})

