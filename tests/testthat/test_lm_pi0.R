# Tests for pi0  lm_pi0

# This outputs a label for devtools::test()
# (useful if test suite has multiple files)
cat("\ntest_lm_pi0.R\n")




# #############################################################################
# Data objects used within the tests

# Create some artificial p-values for tests
p.len <- 20 

# "Uniform" dataset has p-values from a uniform distribution
p.uniform <- seq(0.01, 0.99, length=p.len)
X.flat <- rep(0, p.len)
X.informative <- rep(c(0, 1), each=p.len/2)
X.alt <- rep(c(1,2), p.len/2)
X.matrix <- cbind(flat=X.flat, alt=X.alt, informative=X.informative)

# A shorter-than-default vector for lambda
# (many tests below can use this shorter vector)
lambda.5 <- seq(0.05, 0.95, length=5)




# #############################################################################
# validation of input data


test_that("lambda must be numeric, without NAs, and within (0,1)", {
  temp <- seq(0.1, 0.9, length=9)
  # non-numeric
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=c("0.05", temp)), "numeric")
  # with NA or non-finite
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=c(NA, temp)), "finite")
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=c(temp, Inf)), "finite")
  # well outside range
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=c(-0.05, temp)),"range")
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=c(temp, 1.2)), "range")
  # exactly boundary values 0, 1
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=c(-0.1, 0, temp)),"range")
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=c(temp, 1)), "range")
})


test_that("lambda must be at least a vector of length 4", {
  # length 0
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=NULL))
  # length 1, 2, 3
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=0.5), "4")
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=seq(0.05, 0.95, length=2)), "4")
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=seq(0.05, 0.95, length=3)), "4")
  # repeated values are not ok
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=rep(0.5, 10)), "4")
  # length 4 is ok
  expect_silent(lm_pi0(p.uniform, X=X.flat, lambda=seq(0.05, 0.95, length=4)))  
})


test_that("p must be a numeric vector", {
  expect_error(lm_pi0(), "required")
  temp.char <- letters[1:length(X.flat)]
  expect_error(lm_pi0(temp.char, X=X.flat), "numeric")
  temp.range <- seq(-0.2, 1.2, length=length(X.flat))
  expect_error(lm_pi0(temp.range, X=X.flat), "range")
  temp.NA <- c(NA, seq(0,1, length=length(X.flat)-1))
  expect_error(lm_pi0(temp.NA, X=X.flat), "missing")
})


test_that("some p should be higher than lambda range", {
  p.narrow <- seq(0, 0.5, length=length(p.uniform))
  expect_warning(lm_pi0(p.narrow, X=X.flat, lambda=lambda.5))
})


test_that("df must be a single finite number", {
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df="3"), "number")
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=NA), "number")
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=NULL), "number")
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=rep(4, length(p.uniform))), "number")
})


test_that("df must be in proper range", {
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=1), "range")
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=0), "range")
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=-2), "range")
  expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=100), "lambda")
})


test_that("argument 'type' must convey a supported method name", {
  expect_error(lm_pi0(p.uniform, X=X.flat, type="unsupported"))
}) 


test_that("pvalues and covariates must match", {
  expect_error(lm_pi0(p.uniform, X=c(X.flat, 0.5)), "length")
  expect_error(lm_pi0(p.uniform, X=rbind(X.matrix, 0)), "length")
})


test_that("warnings if covariates are NULL or missing", {
  expect_warning(lm_pi0(p.uniform, lambda=lambda.5), "covariate")
})


test_that("warnings if pvalues and covariates have different names", {
  temp.p <- setNames(seq(0, 1, length=10), letters[1:10])
  temp.X <- setNames(rep(0, length=10), letters[1:10])
  expect_silent(lm_pi0(temp.p, X=temp.X, lambda=lambda.5))
  temp.X2 <- temp.X
  names(temp.X2)[2] = "z"
  expect_error(lm_pi0(temp.p, X=temp.X2, lambda=lambda.5), "names")
})


test_that("covariates must be a numeric vector or matrix", {
  temp.p <- seq(0, 1, length=100)
  expect_error(lm_pi0(temp.p, X=rep(letters[1:10], 10)), "numeric")
  temp.mat <- cbind(A=1:10, B=1)
  temp.df <- data.frame(A=rep(1:10, 10), B=1:5, C=letters[1:10], stringsAsFactors=F)
  ## data frame with numeric columns is OK
  expect_silent(lm_pi0(temp.p, X=temp.df[, c("A", "B")], lambda=lambda.5))
  ## data frame with non-numeric columns is not ok
  expect_error(suppressWarnings(lm_pi0(temp.p, X=temp.df[, c("A", "C")], lambda=lambda.5)),
               "numeric matrix")
})


test_that("covariates must not contain bad values", {
  X.na = X.matrix
  X.na[1,2] = NA
  expect_error(lm_pi0(p.uniform, X=X.na), "missing")
  X.inf = X.matrix
  X.inf[3,1] = Inf
  expect_error(lm_pi0(p.uniform, X=X.inf), "finite")
})




# #############################################################################
# structure of output


test_that("lm_pi0 gives an object of class lm_pi0", {
  result <- lm_pi0(p.uniform, X=X.flat, lambda=lambda.5)
  expect_equal(class(result), "lm_pi0")
})


test_that("all components of lm_qvalue must have unique names", {
  result <- lm_pi0(p.uniform, X=X.flat, lambda=lambda.5)
  expect_equal(names(result), unique(names(result)))
})




# #############################################################################
# toggle between smoothing methods


test_that("smoothing with smooth.spline outputs pi0 for each lambda", {
  result <- lm_pi0(p.uniform, X=X.flat, lambda=lambda.5, smoothing="smooth")
  expect_true("pi0.smooth" %in% names(result))
  expect_equal(dim(result$pi0.smooth), dim(result$pi0.lambda))
})


test_that("smoothing with unit.spline gives pi0 but omits pi0.smooth", {
  result <- lm_pi0(p.uniform, X=X.flat, lambda=lambda.5, smoothing="unit")
  expect_true("pi0" %in% names(result))
  expect_equal(length(result$pi0), length(p.uniform))
  expect_false("pi0.smooth" %in% names(result))
})





# #############################################################################
# uniform p values and uninformative covariates


test_that("uniform pvalues with uninformative covariate vector yield 1", {
  result <- lm_pi0(p.uniform, X=X.flat, lambda=lambda.5)
  expect_equal(result$pi0, rep(1, length(p.uniform)), tol=1e-3)
})


test_that("uniform pvalues with missing covariates yield 1", {
  result <- suppressWarnings(lm_pi0(p.uniform, lambda=lambda.5))
  expect_equal(result$pi0, rep(1, length(p.uniform)), tol=1e-3)
})


test_that("uniform pvalues with null covariates yield 1", {
  result <- suppressWarnings(lm_pi0(p.uniform, X=NULL, lambda=lambda.5))
  expect_equal(result$pi0, rep(1, length(p.uniform)), tol=1e-3)
})


test_that("ordering of lambda does not matter", {
  result.fwd <- lm_pi0(p.uniform, X=X.flat, lambda=lambda.5)
  result.rev <- lm_pi0(p.uniform, X=X.flat, lambda=rev(lambda.5))
  expect_equal(result.fwd$pi0.lambda, result.rev$pi0.lambda)
})


test_that("uniform pvals with uninformative covariate matrix yield 1 (dim 2)", {
  # a two-variable covariate matrix
  temp.X <- cbind(A=rep(0, length(p.uniform)),
                  B=rep(2, length(p.uniform)))
  result <- lm_pi0(p.uniform, X=temp.X, lambda=lambda.5)
  expect_equal(result$pi0, rep(1, length(p.uniform)), tol=1e-3)
})


test_that("uniform pvals with uninformative covariate matrix yield 1 (dim 1)", {
  # a two-variable covariate matrix
  temp.X <- cbind(A=X.flat, B=X.flat+2)        
  result <- lm_pi0(p.uniform, X=temp.X, lambda=lambda.5)
  expect_equal(result$pi0, rep(1, length(p.uniform)), tol=1e-3)
})




# #############################################################################
# uniform p values and informative covariates


test_that("uniform pvals with informative covariate matrix yields non-1s", {
  result <- lm_pi0(p.uniform, X=X.informative, lambda=lambda.5)
  expect_lt(mean(head(result$pi0, p.len/2)), 0.5)
  expect_gt(mean(tail(result$pi0, p.len/2)), 0.5)
})


test_that("adding uninformative covariate does not impact on pi0", {
  result.0 <- lm_pi0(p.uniform, X=X.matrix[, c("informative")],
                     lambda=lambda.5)
  result.1 <- lm_pi0(p.uniform, X=X.matrix[, c("informative", "flat")],
                     lambda=lambda.5)
  expect_equal(result.0$pi0, result.1$pi0)
})


test_that("uniform pvals, informative covariate, all pi0 values are in range (0, 1)", {
  result <- lm_pi0(p.uniform, X=X.informative, lambda=lambda.5)
  expect_lte(max(result$pi0), 1)
  expect_lte(max(result$pi0.lambda), 1)
  expect_gte(min(result$pi0), 0)
  expect_gte(min(result$pi0.lambda), 0)
})




# #############################################################################
# use logistic or linear modeling


test_that("non-informative covariate lead to eqiv pi0 under logistic/linear models", {
  result.logistic <- suppressWarnings(lm_pi0(p.uniform, X=X.flat, type="logistic"))
  result.linear <- suppressWarnings(lm_pi0(p.uniform, X=X.flat, type="linear"))
  expect_equal(result.logistic$pi0, result.linear$pi0)
})


test_that("informative covariate lead to non-eqiv pi0 w. logistic/linear models", {
  result.logistic <- suppressWarnings(lm_pi0(p.uniform, X=p.uniform, type="logistic"))
  result.linear <- suppressWarnings(lm_pi0(p.uniform, X=p.uniform, type="linear"))
  expect_gt(sum(abs(result.logistic$pi0-result.linear$pi0)), 0)
})




# #############################################################################
# use unusual arguments


test_that("smooth.df can accept fractional values for smooth.df", {
  set.seed(1234)
  # create a setup that yields non-unit pi0 estimates
  X.random = p.uniform + runif(length(p.uniform))
  lambda.8 = seq(0.5, 0.95, length=8)
  result.30 <- lm_pi0(p.uniform, X=X.random, smooth.df=3, lambda=lambda.8)
  result.32 <- lm_pi0(p.uniform, X=X.random, smooth.df=3.2, lambda=lambda.8)
  result.39 <- lm_pi0(p.uniform, X=X.random, smooth.df=3.9, lambda=lambda.8)
  # 3.0 and 3.2 should round to the same value, hence exactly equal result
  expect_equal(result.30$pi0, result.32$pi0)
  # 3.9 should round to 4, hence slightly different results
  expect_gt(sum(abs(result.30$pi0-result.39$pi0)), 0)
})





# #############################################################################
# types of p-value inputs

test_that("p value can be a single value", {
  # here suppress warnings in case p<lambda warning
  result.0 <- suppressWarnings(lm_pi0(0.5, X=0, lambda=lambda.5))
  expect_equal(result.0$pi0, 0)
  # here suppress warnings in case p<lambda warning
  result.1 <- suppressWarnings(lm_pi0(1, X=0, lambda=lambda.5))
  expect_equal(result.1$pi0, 1)
})


test_that("p value can be a short vector", {
  result <- lm_pi0(c(0.5, 1), X=c(0,0), lambda=lambda.5)
  expect_equal(result$pi0, c(1, 1))
})

