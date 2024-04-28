source("simulation/simdat.R")
source("simulation/sepdat.R")
library(tidyverse)
library(gridExtra)


nits <- 200

###Run with 
n <- 1000
lock <- 5
sclk <- 5
range_loc <- c(0,10)
range_scl <- c(-5,5)
range_noise <- c(-1,1)
random_range <- TRUE
same_preds <- FALSE
noisek <- 5
beta_range <- c(1,5)
gamma_range <- c(-1,1)
seeds <- Sys.time()


its <- rep(nits, times = 1)
res1 <- apply(as.matrix(its), MARGIN = 1, simdata, lock = lock, sclk = sclk, 
             range_loc = range_loc, range_scl = range_scl, random_range = random_range, 
             n = n, same_preds = same_preds, noisek = noisek,
             beta_range = beta_range, gamma_range = gamma_range, smallvar = FALSE, seeds = seeds)

#res1 <- res
#allres <-   sepdata(res[,1], nits = 10, lock = 5, sclk = 5, range_loc = c(0,10), range_scl = c(-5,5), random_range = TRUE, 
#                    n, same_preds = FALSE, noisek = 2, range_noise = c(-1,1),
#                    beta_range = c(-10, 10), gamma_range = c(-1,1))
allres <- data.frame()
for(i in 1:length(its)){
   temp <- sepdata(res1[,i], nits = 10, lock = lock, sclk = sclk, 
                   range_loc = range_loc, range_scl = range_scl, random_range = random_range, 
                   n = n, same_preds = same_preds, noisek = noisek, range_noise = c(-1,1),
                   beta_range = beta_range, gamma_range = gamma_range)
   allres <- rbind(temp, allres)
}
saveRDS(allres, paste0("simulation/results/allres", 12, ".RDS"))

##########Simulation Boxplots#########
#for the boxplots
nits <- 500
sets <- 4 #number of settings

n <- c(1000, 600, 1000, 1000)
lock <- 5
sclk <- 5
range_loc <- c(-10,10)
range_scl <- c(-1,1)
random_range <- TRUE
same_preds <- FALSE
noi <- c(0,0,0,5)

#draw true parameter values
set.seed(4600)
beta <- runif(lock+1, 0, 10)
gamma <- runif(sclk+1, -0.7, 0.7)
gamma_large <- runif(sclk+1, -1.2, 1.2)
gamma_large[5] <- gamma[5]
gamma_list <- list(gamma, gamma, gamma_large, gamma)


resbox <- matrix(NA, nrow = nits*(lock+sclk+2), ncol = sets) #result matrix
for (i in 1:sets){
  temp <- simdata(nits = nits, lock = lock, sclk = sclk, 
                  range_loc = range_loc, range_scl = range_scl, random_range = TRUE, 
                  n = n[i], same_preds = same_preds, noisek = noi[i],
                  beta = beta, gamma = gamma_list[[i]], careful = TRUE, ret = TRUE)
  
  print(paste("done with round", i))
  temp <- temp[,c(1:(lock+1), (lock+noi[i]+2):(lock+noi[i]+sclk+2))]
  resbox[,i] <- as.vector(temp)
  
  #resbox1 <- apply(as.matrix(its), MARGIN = 1, simdata, lock = lock, sclk = sclk, 
  #                range_loc = range_loc, range_scl = range_scl, random_range = random_range, 
  #                n = n, same_preds = same_preds, noisek = noisek,
  #                beta = beta, gamma = gamma, careful = TRUE, ret = TRUE)  
}
#boxplotsimresults <- list(n = n, range_loc = range_loc, range_scl = range_scl, noisek = noi, 
#                          same_preds = FALSE, beta = beta, gamma = gamma_list, nits = nits, results = resbox)
#saveRDS(boxplotsimresults, "simulation/boxplotsimresults3.RDS")


######make boxplots##########
boxdat <- data.frame()
for (i in 1:sets){
  #resbox1 <- readRDS("simulation/boxplotsimresults2.RDS")$results
  tempdat <- data.frame(setting = rep(i, nits),
                        result = resbox[(2*nits+1):(3*nits),i])
  boxdat <- rbind(tempdat, boxdat)
}

plot1 <- ggplot(data = boxdat, aes(x = factor(setting), y = result, 
                       fill = factor(setting))) +
  geom_boxplot(outlier.shape = 1) + 
  scale_fill_brewer(palette="Blues") +
  geom_hline(yintercept = beta[3])+
  ylab("Estimate beta2") +
  xlab("Setting") +
  guides(fill=FALSE) +
  theme_minimal()

boxdat <- data.frame()
for (i in 1:sets){
  #resbox1 <- readRDS("simulation/boxplotsimresults2.RDS")$results
  tempdat <- data.frame(setting = rep(i, nits),
                        result = resbox[(10*nits+1):(11*nits),i])
  boxdat <- rbind(tempdat, boxdat)
}

plot2 <- ggplot(data = boxdat, aes(x = factor(setting), y = result, 
                                   fill = factor(setting))) +
  geom_boxplot(outlier.shape = 1) + 
  scale_fill_brewer(palette="Blues") +
  geom_hline(yintercept = gamma[5])+
  ylab("Estimate gamma4") +
  xlab("Setting") +
  guides(fill=FALSE) +
  theme_minimal()

grid.arrange(plot1, plot2, ncol = 2)






############Test different algorithms ############
#for the boxplots
nits <- 150

n <- 750
lock <- 3
sclk <- 1
range_loc <- c(0,1)
range_scl <- c(0,1)
random_range <- TRUE
same_preds <- FALSE
noi <- 0

#draw true parameter values
set.seed(500)
beta <- runif(lock+1, 0, 5)
gamma <- runif(sclk+1, -0.3, 0.3)

temp <- simdata_diffalgs(nits = nits, lock = lock, sclk = sclk, 
                         range_loc = range_loc, range_scl = range_scl, random_range = TRUE, 
                         n = n, same_preds = same_preds, noisek = noi,
                         beta = beta, gamma = gamma, careful = TRUE, ret = TRUE)