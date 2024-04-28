context("Fisher-Scoring")
set.seed(2912)

n <- 500
x1 <- runif(n)
x2 <- runif(n)
z1 <- runif(n)
z2 <- runif(n)
locationvec <- cbind(x1, x2)
scalevec <- cbind(z1)

y <- rnorm(n, mean = 3 * x1 - 1.5 * x2, sd = exp(2*z1))

modelWLS <- LocationScaleRegressionWLS$new(y ~ x1 + x2, ~ z1)
modelreg <- LocationScaleRegression$new(y ~ x1 + x2, ~ z1)

gradient_descent(modelreg, stepsize = 0.001, maxit = 1000000)
newton_raphson_WLS(modelWLS, maxit = 10000)

test_that("outputs are of right dimensions", {
  expect_equivalent(length(modelWLS$beta), (ncol(locationvec) + 1))
  expect_equivalent(length(modelWLS$gamma), (ncol(scalevec) + 1))
})

test_that("outputs same/similar as GD in master", {
  expect_equivalent(round(modelWLS$beta, digits = 2), round(modelreg$beta, digits = 2))
  expect_equivalent(round(modelWLS$gamma, digits = 3), round(modelreg$gamma, digits = 3))
})
