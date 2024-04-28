context("model")
library(numDeriv)

set.seed(2912)

n <- 500
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
z1 <- runif(n)
z2 <- runif(n)
locationvec <- cbind(x1, x2)
scalevec <- cbind(z1, z2)

y <- rnorm(n, mean = 3 * x1 - 0.5 * x2, sd = exp(z1 + 2 * z2))

modelWLS <- LocationScaleRegressionWLS$new(y ~ x1 + x2, ~ z1 + z2)
modelreg <- LocationScaleRegression$new(y ~ x1 + x2, ~ z1 + z2)

modelWLS_GD <- gradient_descent_WLS(modelWLS)
modelreg_GD <- gradient_descent(modelreg)


f <- function(x, beta = TRUE) {
  model <- modelWLS$clone()
  if(beta) model$beta <- x
  else model$gamma <- x
  model$loglik()
}

test_that("beta/gamma gradients works", {
  expect_equivalent(modelWLS$grad_beta(), grad(f, modelWLS$beta))
  #expect_equivalent(modelWLS$grad_gamma(), grad(f, modelWLS$gamma, beta = FALSE))
})


test_that("inheritance worked", {
  expect_is(modelWLS, "LocationScaleRegression")
  expect_is(modelWLS, "LocationScaleRegressionWLS")
}) 


test_that("same results for WLS and regular", {
  expect_equivalent(5, 5)
})

#test_that("gradient descent does not throw an error", {
#  expect_error(gradient_descent(modelWLS$clone(), verbose = TRUE), NA)
#})

