library(rworldmap)
data(countryExData)
countries <- countryExData[, 2]
n <- length(countries)
x <- runif(n)
y <- rnorm(n, 5*x, sd = exp(1.2 * x))

model <- LocationScaleRegressionWLS$new(y~x, ~x,
                                        labels = countries)
optim_params(model, method = "ensemble")
pdf.plot(model, which.obs = c(10,100,136,131,22), 
         single.plot = TRUE, legend = TRUE,
         legend_inset = c(0.01, 0.01))
