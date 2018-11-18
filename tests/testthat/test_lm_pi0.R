# Testing lm_pi0

# This outputs a label for devtools::test()
# (useful if test suite has multiple files)
cat("\ntest_lm_pi0.R\n")




# #############################################################################
# Data objects used within the tests

# Create some artificial p-values for tests
p.len <- 20 

# "Uniform" dataset has p-values from a uniform distribution
p.uniform <- seq(0, 1, length=p.len)
X.flat <- rep(0, p.len)
X.informative <- rep(c(0, 1), each=p.len/2)
X.alt <- rep(c(1,2), p.len/2)
X.matrix <- cbind(flat=X.flat, alt=X.alt, informative=X.informative)

# A shorter-than-default vector for lambda
# (many tests below can use this shorter vector)
lambda.5 = seq(0.05, 0.95, length=5)




# #############################################################################
# compatibility of input data


test_that("pvalues and covariates must match", {
  expect_error(lm_pi0(p.uniform, X=c(X.flat, 0.5)))
  expect_error(lm_pi0(p.uniform, X=rbind(X.matrix, 0)))
})


test_that("lambda must be at least a vector of length 4", {
  # length 0
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=NULL))
  # length 1
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=0.5))
  # length 2,3
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=seq(0.05, 0.95, length=2)), "four")
  expect_error(lm_pi0(p.uniform, X=X.flat, lambda=seq(0.05, 0.95, length=3)), "four")
  # length 4 is ok
  expect_silent(lm_pi0(p.uniform, X=X.flat, lambda=seq(0.05, 0.95, length=4)))  
})


if (FALSE) {
  # TO DO
  test_that("df must be greater than one", {
    expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=1))
    expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=0))
    expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=-2))
  })

  test_that("df must be smaller than length of lambda", {
    expect_error(lm_pi0(p.uniform, X=X.flat, smooth.df=6))
  })
}


if (FALSE) {
  # TO DO 
  test_that("argument 'type' must convey a supported method name", {
    expect_error(lm_pi0(p.uniform, X=X.flat, type="unsupported"))
  }) 
}




# #############################################################################
# uniform p values and uninformative covariates


test_that("uniform pvalues with uninformative covariate vector yield 1", {
  result <- lm_pi0(p.uniform, X=X.flat, lambda=lambda.5)
  expect_equal(result$pi0, rep(1, length(p.uniform)))
})

test_that("uniform pvals with uninformative covariate matrix yield 1 (dim 2)", {
  # a two-variable covariate matrix
  temp.X <- cbind(A=rep(0, length(p.uniform)),
                 B=rep(2, length(p.uniform)))
  result <- lm_pi0(p.uniform, X=temp.X, lambda=lambda.5)
  expect_equal(result$pi0, rep(1, length(p.uniform)))
})

test_that("uniform pvals with uninformative covariate matrix yield 1 (dim 1)", {
  # a two-variable covariate matrix
  temp.X <- cbind(A=X.flat, B=X.flat+2)        
  result <- lm_pi0(p.uniform, X=temp.X, lambda=lambda.5)
  expect_equal(result$pi0, rep(1, length(p.uniform)))
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
# types of p-value inputs


test_that("p value can be a single value", {
  result <- lm_pi0(0.5, X=0, lambda=lambda.5)
  expect_equal(result$pi0, 0)
})


if (FALSE) {
  ## TO-DO
  test_that("p values can include NA", {
    p.temp = c(0.1, NA, 0.5, 0.9)
    expect_silent(lm_pi0(p.temp, X=rep(0, length(p.temp)), lambda=lambda.5))
  })
}

