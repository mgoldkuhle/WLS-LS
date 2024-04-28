sepdata <- function(results, nits, lock, sclk, range_loc, range_scl, random_range, 
                              n, same_preds,
                              beta_range, gamma_range,
                              noisek = NULL, range_noise = c(-1,1)){
  sepres <- length(results)/3
  if(is.null(noisek)){
    noisek <- 0
    range_noise <- c(NA, NA)
  }
  
  tabres <- data.frame(n = n, lock = lock, sclk = sclk, 
                       loc_rangel = range_loc[1],
                       loc_rangeu = range_loc[2],
                       scl_rangel = range_scl[1],
                       scl_rangeu = range_scl[2],
                       beta_rangel = beta_range[1],
                       beta_rangeu = beta_range[2],
                       gamma_rangel = gamma_range[1],
                       gamma_rangeu = gamma_range[2],
                       noisek = noisek,
                       noise_rangel = range_noise[1],
                       noise_rangeu = range_noise[2])
  tabres$same_preds <- factor(same_preds, levels = c(FALSE,TRUE),
                              labels = c("no", "yes"))
  
  tabres$avgbias_intbeta <- results[1]
  tabres$avgbias_betas <- mean(results[2:(1+lock)])
  tabres$relbias_intbeta <- results[sepres + 1]
  tabres$relbias_betas <- mean(results[(sepres + 2):(sepres + lock + 2)])
  tabres$CIcbetas <- mean(results[(2*sepres + 1):(2*sepres + 1 + lock)])
  
  tabres$avgbias_intgamma <- results[lock+noisek+2]
  tabres$avgbias_gammas <- mean(results[(lock+noisek+3):(lock+noisek+sclk+2)])
  tabres$relbias_intgamma <- results[sepres + lock + noisek + 2]
  tabres$relbias_gammas <- mean(results[(sepres + lock + noisek + 3):(sepres+lock+noisek+sclk + 2)])
  tabres$CIcgammas <- mean(results[(2*sepres + lock + noisek + 2):(2*sepres + lock + noisek + sclk + 2)])

  if(noisek != 0){
    tabres$avgbias_noisebeta <- mean(results[(lock+2):(lock + noisek + 1)])
    tabres$avgbias_noisegamma <- mean(results[(lock + noisek + sclk + 3):(lock + 2*noisek + sclk + 2)])
    tabres$CIcnoisebeta <- mean(results[(2*sepres + lock + 2):(2*sepres + lock + noisek + 1)])
    tabres$CIcnoisegamma <- mean(results[(2*sepres + lock + noisek + sclk + 3):(2*sepres + lock + 2*noisek + sclk + 2)])
  } else {
    tabres$avgbias_noisebeta <- NA
    tabres$avgbias_noisegamma <- NA
    tabres$CIcnoisebeta <- NA
    tabres$CIcnoisegamma <- NA
  }
  return(tabres)
}
