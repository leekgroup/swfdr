# Tests for lm_qvalue

## This outputs a label during devtools::test()
## (useful if test suite has multiple files)
cat("\ntest_lm_qvalue.R\n")




# #############################################################################
# Data objects used within the tests

# A shorter-than-default vector for lambda
lambda.5 <- seq(0.05, 0.95, length=5)
lambda.10 <- seq(0.05, 0.95, by=0.1)


# short examples
len.short <- 20
p.uniform <- seq(1/len.short, 1, length=len.short)
X.flat = rep(0, len.short)
X.informative <- rep(c(0, 1), each=len.short/2)


# longer examples (divisible by 2, 3, 4, 5, 6)
len.long <- 900
i.long <- seq(1, len.long)
p.uniform.long <- seq(1/len.long, 1, length=len.long)
X.flat.long <- rep(0, len.long)


# helper function to count odd and even numbers in a vector
count.even <- function(x) {
  if (length(x)==0) return (0)
  sum(x %% 2 == 0)
}
count.odd <- function(x) {
  if (length(x)==0) return (0)
  sum(x %% 2 == 1)
}
count.in.set <- function(x, s) {
  if (length(x)==0) return (0)
  sum(x %in% s)
}


## This script uses randomization
## For reproducibility and testing flexibility, each test uses its own seed
## To try out different seeds, use a base seed here
base.seed = 555




# #############################################################################
# validation of input data
# (extensive tests are not required - they are carried out in lm_pi0)


test_that("pvalues must be numeric", {
  expect_error(lm_qvalue(), "required")
  expect_error(lm_qvalue(c(p.uniform, NA)), "missing")
})




# #############################################################################
# structure of output


test_that("lm_qvalue gives an object of class lm_qvalue", {
  result <- lm_qvalue(p.uniform, X=X.flat, lambda=lambda.5)
  expect_is(result, "lm_qvalue")
})


test_that("all components of lm_qvalue must have unique names", {
  result <- lm_qvalue(p.uniform, X=X.flat, lambda=lambda.5)
  expect_equal(names(result), unique(names(result)))
})




# #############################################################################
# uniform p values and uninformative covariates


test_that("uniform pvalues with uninformative covariate give large q", {
  result <- lm_qvalue(p.uniform, X=X.flat)
  expect_gt(min(result$qvalue), 0.1)
})


test_that("uniform pvalues and no covariatesgive large q", {
  result.none <- suppressWarnings(lm_qvalue(p.uniform))
  result.flat <- lm_qvalue(p.uniform, X=X.flat)
  expect_equal(result.none$qvalue, result.flat$qvalue)
})


test_that("uniform pvalues give large q (fpdr)", {
  result <- lm_qvalue(p.uniform, X=X.flat, pfdr=TRUE)
  expect_gt(min(result$qvalue), 0.1)
})




# #############################################################################
# uniform p values and stratifying but uninformative covariates


test_that("uniform pvalues with stratifying or random covariate", {
  set.seed(base.seed + 38382)
  p.uniform.medium <- seq(2/len.long, 1, length=len.long/2)
  X.strat <- rep(c(1,2), len.long/4)
  X.rand <- runif(len.long/2)
  result.strat <- lm_qvalue(p.uniform.medium, X=X.strat, lambda=lambda.10)$qvalues
  result.rand <- lm_qvalue(p.uniform.medium, X=X.rand, lambda=lambda.10)$qvalues
  # since p values were uniform, expect most qvalues to be about the same
  expect_lt(sd(result.strat), 0.1)
  expect_lt(sd(result.rand), 0.2)
  # random covariate might have some variability in q
  expect_lt(sd(result.strat), sd(result.rand))
})




# #############################################################################
# non-unifrom p values without covariates


test_that("non-uniform pvalues with uninformative covariate give reasonable q", {
  p.hits <- c(1e-12, 1e-11, p.uniform)
  result <- lm_qvalue(p.hits, X=rep(0, length(p.hits)), lambda=lambda.5, pfdr=FALSE)
  result2 <- lm_qvalue(p.hits, X=rep(0, length(p.hits)), lambda=lambda.5, pfdr=TRUE)
  expect_lt(min(result$qvalues), 0.05)
})


test_that("non-uniform pvalues with uninformative covariate give reasonable q", {
  # a large fraction of hit with low p valuess
  p.hits <- c(p.uniform/1e2, p.uniform)
  result <- suppressWarnings(lm_qvalue(p.hits, lambda=lambda.5))
  expect_gte(sum(result$qvalues<0.05), length(p.uniform))
})




# #############################################################################
# non-unifrom p values with informative covariates


test_that("non-uniform pvalues with easy covariate", {
  set.seed(base.seed + 92761)
  n.hits <- 50
  p.hits <- runif(n.hits, 0, 1e-4)
  p <- c(runif(len.long, 0, 1), p.hits)
  X <- c(runif(len.long, 0, 1), runif(n.hits, 2, 3))
  i.hits <- len.long + 1:n.hits
  result.0 <- suppressWarnings(lm_qvalue(p, lambda=lambda.10)$qvalue)
  result.X <- lm_qvalue(p, X=X, lambda=lambda.10)$qvalue
  # qvalues using covariates should be lower for the hits
  q.hits.ratio <- result.X[i.hits]/result.0[i.hits]
  expect_lt(max(q.hits.ratio), 1)  
})


test_that("non-uniform pvalues with very easy covariate", {
  set.seed(base.seed + 81331)
  n.hits <- 50
  p.hits <- runif(n.hits, 0, 1e-4)
  p <- c(runif(len.long, 0, 1), p.hits)
  i.hits <- len.long + 1:n.hits  
  # here introduce a p-dependent component into the X covariate
  X <- c(runif(len.long, 0, 1), runif(n.hits, 1.8, 2.2)*log(p[i.hits]))
  result.0 <- suppressWarnings(lm_qvalue(p, lambda=lambda.10)$qvalues)
  result.X <- lm_qvalue(p, X=X, lambda=lambda.10)$qvalues
  # qvalues using covariates should be much lower for the hits
  q.hits.ratio <- result.X[i.hits]/result.0[i.hits]
  expect_lt(max(q.hits.ratio), 1)
  expect_lt(median(q.hits.ratio), 0.1)
  # qvalues among the non-hits should be about the same
  q.nonhits.ratio <- result.X[i.long]/result.0[i.long]
  expect_gt(median(q.nonhits.ratio), 0.8)
  expect_lt(median(q.nonhits.ratio), 1.2)
})


test_that("pvalues from overlapping distributions with informative covariate", {
  set.seed(base.seed + 30001)
  p <- rep(NA, length=len.long)
  # let odd indexes be noise and even indexes be hits
  i.A <- seq(1, len.long, by=2)
  i.B <- seq(2, len.long, by=2)
  # make the p-values uniform + skewed
  p[i.A] <- runif(len.long/2, 0, 1)
  p[i.B] <- rbeta(len.long/2, 1, 10)
  # here introduce a p-dependent component into the X covariate
  X <- rep(c(1,2), len.long/2)
  result.0 <- suppressWarnings(lm_qvalue(p, lambda=lambda.10)$qvalues)
  result.X <- lm_qvalue(p, X=X, lambda=lambda.10)$qvalues
  # qvalues using covariates should be much lower for the hits
  # using covariates should detect many more of the B distribution
  expect_gt(sum(result.X<0.01)/sum(result.0<0.01), 10)
  expect_gt(count.even(which(result.X<0.01)), count.even(which(result.0<0.01)))
  ## should correspond to qvalue interpretation
  hits.X.10 = which(result.X<0.10)
  hits.X.25 = which(result.X<0.25)
  hits.X.50 = which(result.X<0.50)
  # in each test, set a larger number on RHS to allow some noise
  # the test is to catch gross errors, not slight tuning
  expect_lt(count.odd(hits.X.10)/length(hits.X.10), 0.20)
  expect_lt(count.odd(hits.X.25)/length(hits.X.25), 0.35)
  expect_lt(count.odd(hits.X.50)/length(hits.X.50), 0.6)
})


test_that("pvalues from overlapping distributions with noisy informative covariate", {
  set.seed(base.seed + 96732)
  p <- rep(NA, length=len.long)
  # let odd indexes be noise and even indexes be hits
  i.A <- seq(1, len.long, by=2)
  i.B <- seq(2, len.long, by=2)
  p[i.A] <- runif(len.long/2, 0, 1)
  p[i.B] <- rbeta(len.long/2, 0.5, 20)
  # here introduce a p-dependent component into the X covariate
  X <- rep(c(1,2), len.long/2) + rnorm(len.long)
  result.0 <- suppressWarnings(lm_qvalue(p, lambda=lambda.10)$qvalues)
  result.X <- lm_qvalue(p, X=X, lambda=lambda.10)$qvalues
  # there should be more true hits than false hits 
  expect_gt(count.even(which(result.X<0.01)), count.odd(which(result.X<0.01)))
  # should detect more hits with covariate than without
  expect_gt(sum(result.X<0.01), sum(result.0<0.01))
  expect_gt(count.even(which(result.X<0.01)), count.even(which(result.0<0.01)))
  # should correspond to qvalue interpretation
  hits.X.10 <- which(result.X<0.1)
  hits.X.25 <- which(result.X<0.25)
  hits.X.50 <- which(result.X<0.5)
  # in each test, set a larger number on RHS to allow some noise
  # the test is to catch gross errors, not slight tuning
  expect_lt(count.odd(hits.X.10)/length(hits.X.10), 0.20)
  expect_lt(count.odd(hits.X.25)/length(hits.X.25), 0.35)
  expect_lt(count.odd(hits.X.50)/length(hits.X.50), 0.6)
})


test_that("pvalues with a few hits with noisy informative covariate", {
  set.seed(base.seed + 71921)
  len.hits = len.long/6
  p <- rep(NA, length=len.long)
  i.hits = 1:len.hits
  i.nonhits = (len.hits+1):len.long
  p[i.long] <- runif(len.long, 0, 1)
  p[i.hits] <- rbeta(len.hits, 0.5, 10)
  # here introduce a p-dependent component into the X covariate
  X <- rnorm(len.long)
  X[i.hits] <- log(p[i.hits])*rnorm(len.hits, 1, 1)
  result.0 <- suppressWarnings(lm_qvalue(p, lambda=lambda.10)$qvalues)
  result.X <- lm_qvalue(p, X=X, lambda=lambda.10)$qvalues
  # make sure starting with some naive p values below threshold
  FP.naive = count.in.set(which(p<0.05), i.nonhits)
  expect_gt(FP.naive, 0)
  # number of false low-qs should be lower than number of false low-ps
  expect_lt(count.in.set(which(result.0<0.05), i.nonhits), FP.naive)
  # FP rates should correspond to qvalue interpretation
  hits.X.10 <- which(result.X<0.1)
  hits.X.25 <- which(result.X<0.25)
  hits.X.50 <- which(result.X<0.5)
  expect_lt(count.in.set(hits.X.10, i.nonhits)/length(hits.X.10), 0.20)
  expect_lt(count.in.set(hits.X.25, i.nonhits)/length(hits.X.25), 0.35)
  expect_lt(count.in.set(hits.X.50, i.nonhits)/length(hits.X.50), 0.6)
  # TP should be higher with the covariate
  hits.0.10 <- which(result.0<0.1)
  hits.0.25 <- which(result.0<0.25)
  hits.0.50 <- which(result.0<0.5)
  expect_gt(count.in.set(hits.X.10, i.hits), count.in.set(hits.0.10, i.hits))
  expect_gt(count.in.set(hits.X.25, i.hits), count.in.set(hits.0.25, i.hits))
  expect_gt(count.in.set(hits.X.50, i.hits), count.in.set(hits.0.50, i.hits))
})

