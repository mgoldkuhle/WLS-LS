context("Gradient-Descent")
library(numDeriv)

set.seed(2912)

n <- 500
x1 <- runif(n)
x2 <- runif(n)
z1 <- runif(n)
z2 <- runif(n)
locationvec <- cbind(x1, x2)
scalevec <- cbind(z1, z2)

y <- rnorm(n, mean = 3 * x1 - 1.5 * x2, sd = exp(z1 + 2 * z2))

modelWLS <- LocationScaleRegressionWLS$new(y ~ x1 + x2, ~ z1 + z2)
modelreg <- LocationScaleRegression$new(y ~ x1 + x2, ~ z1 + z2)

modelWLS_GD <- gradient_descent_WLS(modelWLS)
modelreg_GD <- gradient_descent(modelreg)

test_that("outputs are of right dimensions", {
  expect_equivalent(length(modelWLS_GD$beta), (ncol(locationvec) + 1))
  expect_equivalent(length(modelWLS_GD$gamma), (ncol(scalevec) + 1))
})

test_that("outputs same/similar as GD in master", {  expect_equivalent(round(modelWLS_GD$beta, digits = 1), round(modelWLS$beta, digits = 1))
  expect_equivalent(round(modelWLS_GD$gamma, digits = 1), round(modelWLS$gamma, digits = 1))
})
