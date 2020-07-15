# Tests for pretty-printing package objects


# This outputs a label for devtools::test()
# (useful if test suite has multiple files)
cat("\ntest_print.R\n")




# #############################################################################
# Data objects used within the tests

# Create some artificial p-values for tests
p.len <- 20 


# "Uniform" dataset has p-values from a uniform distribution
p.uniform <- seq(0.01, 0.99, length=p.len)
X.flat = rep(0, p.len)
X.binary = rep(c(1,2), each=p.len/2)
X.matrix = cbind(flat=X.flat, binary=X.binary)
lambda.5 = seq(0.1, 0.95, length=5)


# instead of recreating lm_pi0 and lm_qvalue objects, have some ready
pi0.flat = lm_pi0(p.uniform, X=X.flat, lambda=lambda.5)
pi0.binary = lm_pi0(p.uniform, X=X.binary, lambda=lambda.5)
pi0.matrix = lm_pi0(p.uniform, X=X.binary, lambda=lambda.5)

q.flat = lm_qvalue(p.uniform, X=X.flat, lambda=lambda.5)
q.binary = lm_qvalue(p.uniform, X=X.binary, lambda=lambda.5)
q.matrix = lm_qvalue(p.uniform, X=X.matrix, lambda=lambda.5)




# #############################################################################
# validation of input data


test_that("print.lm_pi0 emits error on input that is not lm_pi0", {
  # do not accept primitive types
  expect_error(print.lm_pi0(NULL), "lm_pi0")
  expect_error(print.lm_pi0(c(1,2)), "lm_pi0")
  # does accept other package objects
  expect_message(print.lm_pi0(q.flat))
})


test_that("print.lm_qvalue emits error on input that is not lm_pi0", {
  # do not accept primitive types
  expect_error(print.lm_qvalue(NULL), "lm_qvalue")
  expect_error(print.lm_qvalue(NULL), "lm_qvalue")
  # do not accept other package objects
  expect_error(print.lm_qvalue(pi0.flat), "lm_qvalue")
})



# #############################################################################
# print and summary are the same


test_that("summary.lm_pi0 produces output", {
  expect_message(summary.lm_pi0(q.flat), "pi0")
  expect_message(summary.lm_pi0(pi0.flat), "pi0")
})

test_that("summary.lm_qvalue produces output", {
  expect_error(summary.lm_qvalue(pi0.flat), "lm_qvalue")
  expect_message(summary.lm_qvalue(q.flat), "qvalue")
})



# #############################################################################
# content of printed output


test_that("vec.to.string preoduces equally spaced items", {
  aa = c(1, 100, 23)
  expected = "    1   100    23"
  expect_equal(v2s(aa, 5), expected)
})


test_that("print pi0 gives call, lambda, and pi0", {
  expect_message(print(pi0.flat), "lm_pi0")
  expect_message(print(pi0.flat), "lambda")
  expect_message(print(pi0.flat), "pi0")
})


test_that("print qvalues gives call, lambda, pi0, qvalue hits", {
  expect_message(print(q.flat), "lm_qvalue")
  expect_message(print(q.flat), "lambda")
  expect_message(print(q.matrix), "pi0")
  expect_message(print(q.matrix), "q-value")
  expect_message(print(q.matrix), "0.05")
})


test_that("print can rpdouce partial output", {
  expect_message(print(q.flat, components="lambda"), "lambda")
  expect_message(print(q.flat, components="X"), "covariate")
})

