simdata <- function(nits, lock, sclk, range_loc, range_scl, random_range, 
                    n, same_preds = FALSE,
                    beta = NULL, gamma = NULL,
                    beta_range = NULL, gamma_range = NULL,
                    noisek = 0, range_noise = c(-1,1),
                    smallvar = FALSE,
                    ret = FALSE,
                    seeds = 1001){
  #set.seed(202)
  if(same_preds & lock!=sclk){
    warning("same_preds = TRUE, but number of parameters for scale and loc not the same. Setting sclk <- lock")
    sclk <- lock
  }
  noise <- noisek!=0
  #noise <- (!is.null(noisek) | noisek!=0)
  
  if(!noise) {
    noisek <- 0
    range_noise <- NULL
  }

  #predictor names
  names_loc <- sapply(1:lock, function(k1) paste0("x", k1))
  names_scl <- sapply(1:sclk, function(k2) paste0("z", k2))
  if(noise) {
    noise_names <- sapply(1:noisek, function(k3) paste0("v", k3))
  }

  #initialize results list
  results <- list()
  results$est <- matrix(NA, ncol = lock + sclk + 2 * noisek + 2, nrow = nits) #matrix that will contain the estimates
  results$CIcounter <- rep(0, length = lock + sclk + 2 * noisek + 2)
  
  ###make ranges to simulate parameters from####
  if(!is.null(beta_range)){ 
    beta <- runif(lock+1, beta_range[1], beta_range[2])
  }
  
  if(!is.null(gamma_range)){
    gamma <- runif(sclk+1, gamma_range[1], gamma_range[2])
  }
  if(noise){
    beta <- c(beta, rep(0, length = noisek))
    gamma <- c(gamma, rep(0, length = noisek))
  }

  ####make ranges to simulate predictors from###
  #expand ranges into matrices
  range_loc <- matrix(range_loc, nrow = 2, ncol = lock)
  range_scl <- matrix(range_scl, nrow = 2, ncol = sclk)
  if(noise) range_noise <- matrix(range_noise, nrow = 2, ncol = noisek)
  
  if(random_range){ #replace ranges by random values (within the given range)
    #loc_mat <- 1
    set.seed(seeds)
    range_loc <- apply(range_loc, MARGIN = 2, 
                       function(x) sort(runif(2, x[1], x[2])))
    range_scl <- apply(range_scl, MARGIN = 2, 
                       function(x) sort(runif(2, x[1], x[2])))
    
    if(noise) range_noise <- apply(range_noise, MARGIN = 2,
                                   function(x) sort(runif(2, x[1], x[2])))
  }
  
  set.seed(Sys.time())
  ###initialize dataframes for location and scale predictors###
  locdf <- data.frame(rep(NA, length = n))
  scldf <- data.frame(rep(NA, length = n))
  if(noise) noisedf <- data.frame(rep(NA, length = n))
  offs <- rep(1, length = n) #intercept
  
  ###simulation runs###
  for(i in 1:nits){
    
    for (k1 in 1:lock){ #fill location predictors
      locdf[,k1] <- runif(n,range_loc[1, k1], range_loc[2, k1]) 
    }
    names(locdf) <- names_loc
    
    if (same_preds){
      scldf <- locdf
      
    } else {
      for (k2 in 1:sclk){ #fill scale predictors
        scldf[,k2] <- runif(n, range_scl[1, k2], range_scl[2, k2])
      }
      names(scldf) <- names_scl
    }
    
    fulldat <- cbind(locdf, scldf)
    
    #simulate noise preds
    if(noise){
      for (k3 in 1:noisek){
        noisedf[,k3] <- runif(n, range_noise[1, k3], range_noise[2, k3])
      }
      names(noisedf) <- noise_names
    } else {
      noisedf <- NULL
    }
    
    ##simulate dependent variable##
    y_mean = drop(as.matrix(cbind(offs, locdf)) %*% beta[1:(lock+1)]) #subset beta because of noise
    y_sd = exp(drop(as.matrix(cbind(offs, scldf)) %*% gamma[1:(sclk+1)]))
    y <- y_mean + rnorm(n, 0, sd = y_sd)
    
    if(noise){
      fulldat <- cbind(y, fulldat, noisedf)
    } else {
      fulldat <- cbind(y, fulldat)
    }

    #make location and scale formulas
    if(!noise){
      loc_form <- as.formula(paste("y", paste(names(locdf), collapse = " + "), 
                                   sep = " ~ "))
      scl_form <- as.formula(paste(" ", paste(names(scldf), collapse = " + "), 
                                   sep = " ~ "))
    } else {
      loc_form <- as.formula(paste("y", paste(paste(names(locdf), collapse = " + "), 
                                              paste(names(noisedf), collapse = " + "), sep = " + "), 
                                   sep = " ~ "))
      scl_form <- as.formula(paste(" ", paste(paste(names(scldf), collapse = " + "), 
                                              paste(names(noisedf), collapse = " + "), sep = " + "), 
                                   sep = " ~ "))
    }
    
    
    #generate model and optimize parameters
    model <- LocationScaleRegressionWLS$new(loc_form, scl_form, data = fulldat)
    model <- newton_raphson_WLS(model, smallvar = smallvar)
    #model <- optim_params(model, method = "ensemble")
    
    results$est[i, 1:(lock+noisek+1)] <- model$beta
    results$est[i, (lock+noisek+2):(lock+sclk+2*noisek+2)] <- model$gamma
    #results$data <- fulldat
    
    #make CIs matrices (transpose because otherwise mapply doesn't work properly)
    CIs_beta <- t(model$CI_beta())
    CIs_gamma <- t(model$CI_gamma())
                   
    cifun <- function(x, trueparam){ 
      as.numeric(x[1]) < as.numeric(trueparam) & as.numeric(trueparam) < as.numeric(x[2])
    }
    #make binary vectors (is true parameter in CI)
    CIbetatrue <- as.integer(mapply(FUN = cifun, 
                                    x = as.data.frame(CIs_beta), trueparam = beta))
    CIgammatrue <- as.integer(mapply(FUN = cifun, 
                                     x = as.data.frame(CIs_gamma), trueparam = gamma))
    cic <- c(CIbetatrue, CIgammatrue)
    results$CIcounter <- results$CIcounter + cic
  }
  #also return ranges, "true" parameters and (avg) bias
  results$params <- c(beta, gamma)
  results$ranges_preds <- c(range_loc, range_scl, range_noise)
  results$ranges_params <- c(beta_range, gamma_range)
  results$biases <- sweep(results$est, 2, results$params)
  results$avgbias <- apply(results$biases, 2, mean)
  results$relavgbias <- results$avgbias/results$params
  #for the noise parameters, replace relative average bias with just average bias
  results$relavgbias[is.infinite(results$relavgbias)] <- results$avgbias[is.infinite(results$relavgbias)]
  if(ret){
    return(results$est)
  } else  {
    return(c(results$avgbias, results$relavgbias, results$CIcounter))
  }
}
