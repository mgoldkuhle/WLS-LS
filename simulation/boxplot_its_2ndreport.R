#note: for this simulation, the optimization algorithms were temporarily 
#changed st they returned i (the number of iterations they needed to converge) 
#rather than the optimized parameter values 

library(ggplot2)
library(tidyr)

reg <- c()
GD <- c()
NR <- c()

n <- 500
for(i in 1:100){
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- runif(n)
  z2 <- runif(n)
  coef <- runif(4, -5, 5)
  
  y <- rnorm(n, mean = coef[1] * x1 - coef[2] * x2, sd = exp(coef[3]*z1 + coef[4] * z2))
  
  modelGD <- LocationScaleRegressionWLS$new(y~x1 + x2, ~z1 + z2)
  modelNR <- LocationScaleRegressionWLS$new(y~x1 + x2, ~z1 + z2)
  modelreg <- LocationScaleRegression$new(y~x1 + x2, ~z1 + z2)
  
  GD[i] <- gradient_descent_WLS(modelGD, maxit = 10000)
  NR[i] <- newton_raphson_WLS(modelNR)
  reg[i] <- gradient_descent(modelreg, maxit = 20000)
}

GD1 <- GD[1:81]
perform <- data.frame(Fisher_Scoring = NR, GD_WLS = GD1, GD_reg = reg)

perform_long <- gather(perform, method, iterations)
perform_long$method <- factor(perform_long$method,
                              levels = c("Fisher_Scoring", "GD_WLS", "GD_reg"), 
                              ordered = TRUE)

ggplot(perform_long, aes(x = method, y = iterations, fill = method)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_brewer(palette="Blues") +
  ylab("number of iterations (log-transformed)") + 
  xlab("optimization method used") +
  guides(fill=FALSE) +
  theme_minimal()