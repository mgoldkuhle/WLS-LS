library(numDeriv)

set.seed(2912)

n <- 500
x1 <- runif(n)
x2 <- runif(n)
z1 <- runif(n)
z2 <- runif(n)
locationvec <- cbind(x1, x2)
scalevec <- cbind(z1, z2)

y <- rnorm(n, mean = 3 * x1 - 0.5 * x2, sd = exp(z1 + 2 * z2))

modelWLS <- LocationScaleRegressionWLS$new(y ~ x1 + x2, ~ z1 + z2)
modelreg <- LocationScaleRegression$new(y ~ x1 + x2, ~ z1 + z2)

#newton_raphson_WLS(modelWLS, maxit = 10000)

#test_that("plot functions work as intended", {
#  expect_error(plot.LocationScaleRegressionWLS(modelWLS), NA)
#})