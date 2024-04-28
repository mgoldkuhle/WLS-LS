#' @title asp20plot
#' @description The R6 class \code{asp20plot} provides parameter estimation based on 
#' weighted/generalized least squares (WLS/GLS) and corresponding plots functions for model diagnostics 
#' 
#' @field W A numeric matrix with the weight parameters.
#' 
#' @import R6
#' @import asp20model
#' 
#' @export

#library(R6)
#library(asp20model)

LocationScaleRegressionWLS <- R6Class(classname = "LocationScaleRegressionWLS",
  inherit = LocationScaleRegression,
  public = list(
    W = numeric(),
    preds = numeric(),
    inv = numeric(),
    df = numeric(),
    n = numeric(),
    k_loc = numeric(),
    k_scale = numeric(),
    labels = character(),
    
    
    #' @details
    #' Create a new `LocationScaleRegressionWLS` object.
    #' @param location A two-sided formula with the response variable on the
    #'                 LHS and the predictor for the location (i.e. the mean)
    #'                 on the RHS.
    #' @param scale A one-sided formula with the predictor for the scale
    #'              (i.e. the standard deviation) on the RHS.
    #' @param data A data frame (or list or environment) in which to evaluate
    #'             the `location` and `scale` formulas.
    #' @param labels optional, an extra column of labels for the observations
    #'             (e.g. country names or patient numbers). Can be used in plots
    #'             later on.
    #'  
    #' @param ... Passed on to [stats::model.matrix()].
    #'
    #' @return
    #' A `LocationScaleRegressionWLS` object.
    #'
    #' @examples
    #' y <- rnorm(30)
    #' LocationScaleRegressionWLS$new(y ~ 1)
    #'
    #' @importFrom stats model.matrix
    #' @importFrom stats formula
    
    
    
    initialize = function(location,
                          scale = ~1,
                          data = environment(location),
                          labels = NULL,
                          ...) {
      
      if (!is.null(labels) & !is.character(labels)){
        warning("coercing labels to class character")
        labels <- as.character(labels)
      }
      
      #scale <- update(scale, paste(location[[2]], "~ ."))
      #private$y <- eval(location[[2]], data, environment(location))
      #private$X <- model.matrix(location, data, ...)
      #private$Z <- model.matrix(scale, data, ...)
      
      super$initialize(location, scale, data)
      
      private$loc_form <- location
      private$scale_form <- scale
      
      #remove first column in X as it contains y
      private$Xdf <- model.frame(location, data, ...)[,-1]
      private$Zdf <- model.frame(scale, data, ...)
      
      self$beta <- rep.int(0, ncol(private$X))
      self$gamma <- rep.int(0, ncol(private$Z))
      
      invisible(self)
      self$n <- nrow(private$X)
      self$k_loc <- ncol(private$X) - 1
      self$k_scale <- ncol(private$Z) - 1
      self$df <- sum(self$hatvalues)
      
      self$W <- diag(nrow(private$X))
      self$preds <- matrix(private$X[,-1])
      self$inv <- diag(self$k_loc + 1)
      
      if (!is.null(labels) & !is.character(labels)) stop("labels must be character")
      self$labels <- labels
      
    },
    
    #' @details
    #' Returns different types of residuals of a
    #' `LocationScaleRegression` object at the current parameter values.
    #' The working and response residuals are the observed values
    #' minus the fitted values for the location.
    #' The deviance and Pearson residuals are weighted residuals,
    #' divided by the fitted values for the scale.
    #' The standardized residuals are weighted residuals, divided
    #' by the fitted values for the scale and the square root of
    #' (1-h_i), where h_i is the i-th hat value.
    #' The studentized residuals are weighted residuals, divided
    #' by (1), the fitted values for the scale after removing the 
    #' i-th observation and (2), the square root of (1-h_i), where 
    #' h_i is the i-th hat value  
    #'
    #' @param type The type of residuals which should be returned.
    #'             Can be abbreviated.
    #'
    #' @return
    #' A numeric vector.
    #'
    #' @examples
    #' y <- rnorm(30)
    #' model <- LocationScaleRegression$new(y ~ 1)
    #' model$resid()
    #' model$resid("deviance")
    
    resid = function(type = c("working", "response", "deviance", 
                              "pearson", "standard", "student")){
      type <- match.arg(type)
      
      if (type == "working" || type == "response"){
        resid <- private$.resid
        
      } else if (type == "deviance" || type == "pearson"){
        resid <- private$.resid / self$fitted_scale
        
      } else if (type == "standard"){
        denom <- sqrt(1 - self$hatvalues)
        resid <- (1/self$fitted_scale) * (private$.resid/denom)
        
      } else if (type == "student"){
        str <- self$resid(type = "standard")
        fac <- (1/self$fitted_scale)
        rerr <- sqrt(((self$SSEw) - str^2)/(self$n - self$df - 1))
        resid <- fac * (private$.resid/(rerr * sqrt(1 - self$hatvalues)))
      }
      resid
    },
  
  #' @details
  #' Returns confidence intervals for the \eqn{\beta} parameters of a
  #' fitted `LocationScaleRegression` object at the current parameter 
  #' values.
  #' 
  #' @param alpha the desired confidence level, e.g. alpha = 0.05 (the
  #' default)
  #'
  #' @return
  #' A matrix with each row corresponding to one \eqn{\beta} (in order).
  #' The first (second) column corresponds to the left (right) boundary 
  #' of the confidence interval.
  #'
  #' @examples
  #' y <- rnorm(30)
  #' model <- LocationScaleRegression$new(y ~ 1)
  #' optim_params(model)
  #' model$CI_beta(alpha = 0.01)
  
    CI_beta = function(alpha = 0.05){
      if (any(self$grad_beta() > 0.01)){
        return("model does not seem to be fitted. Run optim_params() first.")
      }
      qalph <- 1 - (alpha/2)
      edf <- self$n - sum(self$hatvalues)
      confis <- matrix(NA, nrow = self$k_loc + 1, ncol = 2)
      beta_vars <- diag(solve(self$hess_beta()))
      beta_sds <- sqrt(beta_vars)
      
      confis[,1] <- self$beta - qt(qalph, edf) * beta_sds
      confis[,2] <- self$beta + qt(qalph, edf) * beta_sds
      return(confis)
    },
  
  #' @details
  #' Returns confidence intervals for the \eqn{\gamma} parameters of a
  #' fitted `LocationScaleRegression` object at the current parameter 
  #' values.
  #' 
  #' @param alpha the desired confidence level, e.g. alpha = 0.05 (the
  #' default)
  #'
  #' @return
  #' A matrix with each row corresponding to one \eqn{\gamma} (in order).
  #' The first (second) column corresponds to the left (right) boundary 
  #' of the confidence interval.
  #'
  #' @examples
  #' y <- rnorm(30)
  #' model <- LocationScaleRegression$new(y ~ 1)
  #' optim_params(model)
  #' model$CI_gamma(alpha = 0.01)
  
    CI_gamma = function(alpha = 0.05){
      if (any(self$grad_gamma() > 0.01)){
        return("model does not seem to be fitted. Run optim_params() first.")
      }
      qalph <- 1 - (alpha/2)
      edf <- self$n - sum(self$hatvalues)
      confis <- matrix(NA, nrow = self$k_scale + 1, ncol = 2)
      gamma_vars <- diag(solve(self$hess_gamma()))
      gamma_sds <- sqrt(gamma_vars)
      
      confis[,1] <- self$gamma - qt(qalph, edf) * gamma_sds
      confis[,2] <- self$gamma + qt(qalph, edf) * gamma_sds
      return(confis)
    }
  ),
  
  private = list(
    Xdf = numeric(),
    Zdf = numeric(),
    loc_form = formula(),
    scale_form = formula(),
    
    update_gamma = function(value){
      private$.gamma <- value
      self$fitted_scale <- exp(drop(private$Z %*% self$gamma))
      self$W <- diag(1/(self$fitted_scale)^2)
      self$inv <- chol2inv(chol((t(private$X) %*% self$W) %*% private$X))
      self$beta <- drop(self$inv %*% (t(private$X) %*% self$W %*% private$y))
      invisible(self)
    }
  ),
  
  active = list(
    #' @field hatmatrix The model's hat matrix wrt to the `location` predictors
    hatmatrix = function(value){
      if (missing(value)){
        private$X %*% self$inv %*% (t(private$X) %*% self$W)
      } else {
        stop("can't set `$hatmatrix`", call. = FALSE)
      }
    },
    
    #' @field hatvalues The model's hat (or leverage) values wrt to 
    #' the `location` predictors 
    hatvalues = function(value){
      if (missing(value)){
        diag(self$hatmatrix)
      } else {
        stop("can't set `$hatvalues`", call. = FALSE)
      }
    },
    
    #' @field zhatvalues (informal measure) The model's hat (or leverage) 
    #' values wrt to the `scale` predictors
    zhatvalues = function(value){
      if (missing(value)){
        zinv <- solve(t(private$Z) %*% private$Z)
        diag(private$Z %*% zinv %*% t(private$Z)) 
      } else {
        stop("can't set `$zhatvalues`", call. = FALSE)
      } 
    },
    
    #' @field cooks_dist The model's Cook's Distance values wrt to the 
    #' `location` predictors
    cooks_dist = function(value){
      str <- self$resid(type = "standard")
      hi <- self$hatvalues
      (str^2 / (self$k_loc + 1)) * (hi / (1 - hi))
    },
    
    #' @field residuals The model's response residuals
    residuals = function(value){
      if (missing(value)){
        self$resid()
      } else {
        stop("can't set `$residuals`", call. = FALSE)
      }
    },
    
    #' @field edf The model's effective degrees of freedom, the trace 
    #' of the hat matrix
    edf = function(value){
      if (missing(value)){
        self$df
      } else {
        stop("can't set `$edf`", call. = FALSE)
      }
    },
    
    #' @field AIC The model's Akaike information criterion
    AIC = function(value){
      if (missing(value)){
        -2 * self$loglik() + 2 * (self$df)
      } else {
        stop("can't set `$AIC`", call. = FALSE)
      }
    },
    
    #' @field R2 The model's R squared
    R2 = function(value){
      if (missing(value)){
        1 - (self$SSE/self$TSS)
      } else {
        stop("can't set `$R2`", call. = FALSE)
      }
    },
    
    #' @field SSE The model's sum of squared errors
    SSE = function(value){
      if (missing(value)){
        sum(self$resid()^2)
      } else {
        stop("can't set `$SSE`", call. = FALSE)
      }
    },
    
    #' @field SSEw The model's weighted sum of squared errors. Equal to `n`
    SSEw = function(value){
      if (missing(value)){
        drop(self$resid() %*% self$W %*% self$resid())
      } else {
        stop("can't set `$SSE`", call. = FALSE)
      }
    },
    
    #' @field MSE The model's mean squared error
    MSE = function(value){
      if (missing(value)){
        self$SSE / (self$n - self$df)
      } else {
        stop("can't set `$MSE`", call. = FALSE)
      }
    },
    
    #' @field RMSE The model RMSE (sqrt of the MSE)
    RMSE = function(value){
      if (missing(value)){
        sqrt(self$MSE)
      } else {
        stop("can't set `$RMSE`", call. = FALSE)
      }
    },
    
    #' @field TSS The model's total sum of squares
    TSS = function(value){
      if (missing(value)){
        yy <- private$y
        sum((yy - mean(yy))^2)
      } else {
        stop("can't set `$TSS`", call. = FALSE)
      }
    },
    
    #' @field Xm The design matrix used for estimating the `location`
    Xm = function(value){
      if (missing(value)){
        private$X
      } else {
        stop("can't set `$X`", call. = FALSE)
      }
    },
    
    #' @field Zm The design matrix used for estimating the `scale`
    Zm = function(value){
      if (missing(value)){
        private$Z
      } else {
        stop("can't set `$Z`", call. = FALSE)
      }
    },
    
    #' @field pred_data A data.frame containing all `location` and 
    #' `scale` predictors
    pred_data = function(value){
      if (missing(value)){
        cbind(private$Xdf, private$Zdf)
      } else {
        stop("can't set `$pred_data`", call. = FALSE)
      }
    },
    
    #' @field formula_location The formula used to predict the `location`
    formula_location = function(value){
      if (missing(value)){
        private$loc_form
      } else {
        stop("can't set `$formula_location`", call. = FALSE)
      }
    }, 
    
    #' @field formula_scale The formula used to predict the `scale`
    formula_scale = function(value){
      if (missing(value)){
        private$scale_form
      } else {
        stop("can't set `$formula_scale`", call. = FALSE)
      }
    }
  )
)

tryupdate <- function(mod, step, type){
  mod0 <- mod$clone() #get old values in case of error
  gamma0 <- mod0$gamma
  #gradgam <- mod$grad_gamma()
  
  tryCatch(
    {
      if (type=="GD"){
        mod$gamma <- mod$gamma + step * mod$grad_gamma()
      } else {
        fisher_gamma <- solve(mod$hess_gamma())
        mod$gamma <- mod$gamma + step * drop(fisher_gamma %*% mod$grad_gamma()) 
      }
    },
    
    error = function(e){
      mod$gamma <- gamma0
      step = step/10
      tryupdate(mod = mod, step = step, type = type)
    }
  )
  return(mod)
}

newton_raphson_WLS <- function(model,
                               abstol = 0.0001,
                               maxit = 1000,
                               verbose = FALSE,
                               smallvar = FALSE) {
  
  stepsi = c(rep(10^-8, 2), rep(10^-7,2), rep(10^-6,2), rep(10^-4, 2), 10^-3, 10^-2, 10^-1)
  for (i in seq_len(maxit)){
    if(!smallvar & i<=11){
      tryupdate(model, step = stepsi[i], type = "NR")
    } else {
      tryupdate(model, step = 1, type = "NR")
    }
    
    grad_gamma <- model$grad_gamma()
    
    if (verbose) {
      par_msg <- c(model$beta, model$gamma)
      par_msg <- format(par_msg, trim = TRUE, digits = 3)
      par_msg <- paste(par_msg, collapse = " ")
      
      grad_msg <- c(grad_gamma)
      grad_msg <- format(grad_msg, trim = TRUE, digits = 3)
      grad_msg <- paste(grad_msg, collapse = " ")
      
      loglik_msg <- format(model$loglik(), digits = 3)
      
      message(
        "Iteration:      ", i, "\n",
        "Parameters:     ", par_msg, "\n",
        "Gradient:       ", grad_msg, "\n",
        "Log-likelihood: ", loglik_msg, "\n",
        "==============="
      )
    }
    
    if (any(is.nan(c(model$beta, model$gamma)))){
      stop("Sorry, seems like we've diverged. Try again, setting smallvar = FALSE")
    }
    
    if (all(abs(c(grad_gamma)) <= abstol)) break
  }
  
  message("Finishing after ", i, " iterations")
  invisible(model)
}

gradient_descent_WLS <- function(model,
                                 stepsize = 10^-3,
                                 maxit = 1000,
                                 abstol = 0.001,
                                 verbose = FALSE) {
  
  for (i in seq_len(maxit)) {
    
    tryupdate(model, stepsize, type = "GD")
    
    grad_gamma <- model$grad_gamma()
    
    if (verbose) {
      par_msg <- c(model$beta, model$gamma)
      par_msg <- format(par_msg, trim = TRUE, digits = 3)
      par_msg <- paste(par_msg, collapse = " ")
      
      grad_msg <- c(grad_gamma)
      grad_msg <- format(grad_msg, trim = TRUE, digits = 3)
      grad_msg <- paste(grad_msg, collapse = " ")
      
      loglik_msg <- format(model$loglik(), digits = 3)
      
      message(
        "Iteration:      ", i, "\n",
        "Parameters:     ", par_msg, "\n",
        "Gradient:       ", grad_msg, "\n",
        "Log-likelihood: ", loglik_msg, "\n",
        "==============="
      )
    }
    if (all(abs(c(grad_gamma)) <= abstol)) break
  }
  message("Finishing after ", i, " iterations")
  invisible(model)
}

#' Optimizer for the location and scale parameters in the 
#' `LocationScaleRegressionWLS` model class
#'
#' This function optimizes the log-likelihood of the given location-scale
#' regression model by either Fisher Scoring (the default) or Gradient Descent.
#' The Fisher Information (derivative of the gradient) is approximated by the 
#' differential quotient.  
#' The function has a side effect on the `model` object.
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' @param maxit The maximum number of iterations.
#' @param abstol The absolute convergence tolerance. The algorithm stops if the
#'               absolute value of the gradient drops below this value.
#' @param smallvar logical, set to FALSE by default. If TRUE, will override the 
#'               restriction on updating size that is placed on the Fisher-Scoring
#'               algorithm for the first iterations. Only advisable in model settings
#'               where the variance is suspected or known to be small
#'
#' @return
#' The updated model, invisibly.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#'
#' @export
#' 

optim_params <- function(model, maxit = 1000, abstol = 0.001, 
                         method = c("FS", "GD", "ensemble"),
                         smallvar = FALSE, stepsize = 10^-3){
  
  if (!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  } 
  
  if ((all(abs(c(model$grad_gamma())) <= abstol))){
    stop("model is already fitted given the supplied value for abstol. If a finer fit is desired, reduce abstol")
  } 
  
  if(length(method) > 1) {
    method <- "FS"
  }

  if (method == "FS"){
    tryCatch(
      {
        mod0 <- model$clone()
        newton_raphson_WLS(model, maxit = maxit, abstol = abstol, smallvar = smallvar)
      },
      error = function(e){
        message("Fisher Scoring failed due to initial overshooting. \n
                consider using method = \"ensemble\" or method = \"GD\" instead.")
      }
    )
    
  } else if (method == "GD"){
    warning("Usage of the gradient descent algorithm is in general not advised, 
            as it is very slow. \n Please only use in cases where Fisher Scoring faile,")
    gradient_descent_WLS(model, maxit = maxit, abstol = abstol, stepsize = stepsize)
    
  } else if (method == "ensemble"){
    #model$gamma <- model$gamma + 10^-10 * model$grad_gamma()
    while (all(abs(c(model$grad_gamma())) > abstol)){
      tryCatch(
        {
          mod0 <- model$clone()
          newton_raphson_WLS(model, maxit = maxit, abstol = abstol)
        },
        error = function(e){
          model$gamma <- mod0$gamma
          gradient_descent_WLS(model, maxit = 10, abstol = abstol, stepsize = 10^(-10))
        }
      )
    } 
  } else {
    stop("did not provide a valid convergence algorithm. The options are 
         Fisher Scoring (type = \"FS\") or Gradient Descent (type = \"GD\")")
  }
  invisible(model)
}

#' prediction function for location and scale for new data supplied to 
#' a model in the `LocationScaleRegressionWLS` model class
#'
#' This function makes predictions for the location and scale for a 
#' data frame of new data, given a fitted `LocationScaleRegressionWLS
#' model
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' @param newdata a data frame of new data points to make predictions on.
#'
#' @return
#' A dataframe with the predictions for location and scale.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' newdata <- 1
#' predict(model, newdata)
#'
#' @export
#' 

predict.LocationScaleRegressionWLS <- function(model, newdata){
  
  if (!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  } 
  
  if (missing(newdata) || is.null(newdata)){
    message("since no new data was supplied, returning the original predictions")
    return(cbind(model$fitted_location, model$fitted_scale))
  }
  
  if (nrow(newdata)==0){
    stop("newdata does not contain any observations")
  }
  
  #check if everything in newdata is correct class 
  if (!is.data.frame(newdata)) {
    stop("new data must be of class data.frame")
  }
  if (any ((sapply(newdata, class))=="character")) {
    stop("class character in newdata not allowed \n newdata must be of type numeric")
  }
  if (any ((sapply(newdata, class))=="logical")) {
    stop("class logical in newdata not allowed \n newdata must be of type numeric")
  }
  if (any ((sapply(newdata, class))=="factor")) {
    stop("class factor in newdata not allowed \n newdata must be of type numeric")
  }
  
  ##### prep data #####
  newdata <- newdata[unique(names(newdata))] #get rid of duplicate columns in newdata
  
  #get variable names from model
  loc_names <- all.vars(model$formula_location)[-1]
  scale_names <- all.vars(model$formula_scale)
  modnames <- unique(c(loc_names, scale_names))
  
  #check if newdata has the correct variables
  #do this twice (loc and scale) for more meaningful error messages
  tryCatch(
    {
      locdat <- newdata[loc_names]
    }, 
    error = function(e){
      stop("location variables do not match those in function call")
    }
  )
  tryCatch(
    {
      scaledat <- newdata[scale_names]
    }, 
    error = function(e){
      stop("scale variables do not match those in function call")
    }
  )
  
  #check if the range of parameters is roughly similar
  rg_x <- apply(model$Xm, MARGIN = 2, range)[,-1]
  rg_z <- apply(model$Zm, MARGIN = 2, range)[,-1]
  for (i in 1:ncol(rg_x)){
    if(!all(locdat[,i] >= rg_x[1,i] & locdat[,i] <= rg_x[2,i])){
      which_one <- loc_names[i]
      warning(paste("some values of location predictor", which_one, "fall outside 
                    of the range of the \n", 
                    "predictor's values in model fit. Prediction may not be reliable"))
    } 
  }
  for (i in 1:ncol(rg_z)){
    if(!all(scaledat[,i] >= rg_z[1,i] & scaledat[,i] <= rg_z[2,i])){
      which_one <- scale_names[i]
      warning(paste("some values of scale predictor", which_one, "fall outside of the
                    range of the \n", 
                    "predictor's values in model fit. Prediction may not be reliable"))
    } 
  }
  
  offs <- rep(1, nrow(newdata))
  locmat <- as.matrix(cbind(offs, locdat))
  scalemat <- as.matrix(cbind(offs, scaledat))

  ###make predictions###
  results <- data.frame(
    loc_pred = drop(locmat %*% model$beta),
    scale_pred = exp(drop(scalemat %*% model$gamma))
  )
  return(results)
}



#' Plot function for the `LocationScaleRegressionWLS` model class
#'
#' This function provides plots for checking the results of a fitted
#' LocationScaleRegressionWLS object.
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' @param type type of plots that are created. If nothing is specified,
#' fitted vs actual values, residuals' normal Q-Q, residuals kernel density,
#' residuals vs. fitted are plotted. Options:
#' `influence`: fitted vs studentized residuals with representation of cook's distance
#' to observe the influence of single data points
#' `spread`: fitted values vs absolute studentized residuals with regression line to
#' observe the spread-level and assess constant error variance
#' @param single optional parameter that allows to only plot one plot of a set
#' specified with the `type` parameter. 
#' If no `type` was specified, valid values are:
#' `fitactual`: scatter plot that shows fitted vs actual values
#' `qq`: Normal QQ-Plot of the residuals
#' `density`: kernel density estimation of the ...
#' `fitres`: scatter plot fitted values vs residuals including a loess smoothing line
#' If `type`="influence" was specified:
#' `cooks`: Cook's distance for each observation with labeled outliers
#' `hatvalues`: Hat-Values for each observation with labeled outliers
#' `combined`: Hat-Values vs. Studentized Residuals with the circle size indicating
#' Cook's distance, threshold lines at -2, 0, 2 for studentized residuals and
#' at 2p/n and 3p/n for hat-values and marked and labeled observations that exceed the
#' Cook's distance threshold of 4/(n-k-1)
#' If `type`="spread" was specified, there are no further choices, since this only
#' produces one plot by default.
#'
#' @return
#' Returns the specified plots
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' plot(model)
#'
#' @export

plot.LocationScaleRegressionWLS <- function(model, type = "", single = ""){
  if(any(model$grad_gamma() > 0.01)){
    warning("It seems that the model was not fitted. Run optim_params(model) first.")
  }
  
  resid_reg <- model$resid() #residuals
  resid_dev <- model$resid("deviance") #residuals divided by the scale
  fitted_vals <- model$fitted_location 
  y_vals <- model$fitted_location + resid_reg
  
  #get labels or replace with obs. numbers if no labels were passed in the model spec.
  if(!is.null(model$labels)){
    obs.names <- model$labels
  } else {
    obs.names <- sapply(1:model$n, 
                        FUN = function(x) as.character(x))
  }
  
  #influence plots
  if(type == "influence"){
    hatvalues <- model$hatvalues
    rstudent <- model$resid(type = "student")
    cooks <- model$cooks_dist
    scale <- 8 / max(cooks, na.rm=TRUE)
    p <- ncol(model$preds) + 1
    n <- nrow(model$preds)
    threshold_cooks <- 4/(n-p-2)
    hcooks <- which(cooks > threshold_cooks)
    threshold_hat <- 2*p/n
    hhat <- which(hatvalues > threshold_hat)
    
    if(single == "cooks"){
      #cook's distance per obs.
      plot(cooks, main="Cook's Distance",
           ylab="Cook's Distance", xlab="Obs. Number", 
           sub=paste("Cut-off at", toString(round(threshold_cooks, digits=3))),
           ylim=c(0, max(cooks)+0.05))
      abline(h=c(threshold_cooks, 1), lty=2, col="red")
      if(length(hcooks) > 0){
        text(hcooks, cooks[hcooks], obs.names[hcooks], pos=3, cex=0.7)
      }
    } else if(single == "hatvalues") {
      #hatvalue per obs.
      plot(hatvalues, main="Hat-Values",
           xlab="Obs. Number", ylab="Hat-Value", ylim=c(0, max(hatvalues)+0.05))
      abline(h=c(2, 3)*p/n, col="red", lty=2)
      if(length(hhat>0)){
        text(hhat, hatvalues[hhat], obs.names[hhat], pos=3, cex=0.7)
      }
    } else if(single == "combined"){
      #hatvalues vs stud. residuals with cook's distance
      plot(hatvalues, rstudent, 
           main="Influence",
           xlab="Hat-Values", sub="Circle Size: Cook's Distance", 
           ylab="Studentized Residuals", type="n",
           xlim=c(0, max(hatvalues)+0.05), ylim=c(min(rstudent)*1.2, max(rstudent)*1.2))
      points(hatvalues, rstudent, cex=scale*cooks)
      if(length(hcooks)>0){
        points(hatvalues[hcooks], rstudent[hcooks], cex=scale*cooks[hcooks], col="red")
        text(hatvalues[hcooks], rstudent[hcooks], obs.names[hcooks], col="red", cex=0.7)
      }
      abline(h=c(-2, 0, 2), lty=2)
      abline(v=c(2, 3)*p/n, lty=2)
    } else {
    #combined plots
    #restore plotting device par on exit
    par_orig <- par(ask=T) #click enter to show next plot
    on.exit(par(par_orig))
    
    plot(cooks, main="Cook's Distance",
         ylab="Cook's Distance", xlab="Obs. Number", 
         sub=paste("Cut-off at", toString(round(threshold_cooks, digits=3))),
         ylim=c(0, max(cooks)+0.05))
    abline(h=c(threshold_cooks, 1), lty=2, col="red")
    if(length(hcooks>0)){
      text(hcooks, cooks[hcooks], obs.names[hcooks], pos=3, cex=0.7)
    }
      
    plot(hatvalues, main="Hat-Values",
         xlab="Obs. Number", ylab="Hat-Value", ylim=c(0, max(hatvalues)+0.05))
    abline(h=c(2, 3)*p/n, col="red", lty=2)
    if(length(hhat>0)){
    text(hhat, hatvalues[hhat], obs.names[hhat], pos=3, cex=0.7)
    }
      
    plot(hatvalues, rstudent, 
         main="Hat-Values vs. Stud. Residuals",
         xlab="Hat-Values", sub="Circle Size: Cook's Distance", 
         ylab="Studentized Residuals", type="n",
         xlim=c(0, max(hatvalues)+0.05), ylim=c(min(rstudent)*1.2, max(rstudent)*1.2))
    points(hatvalues, rstudent, cex=scale*cooks)
    if(length(hcooks>0)){
      points(hatvalues[hcooks], rstudent[hcooks], cex=scale*cooks[hcooks], col="red")
      text(hatvalues[hcooks], rstudent[hcooks], obs.names[hcooks], col="red", cex=0.7)
    }
    abline(h=c(-2, 0, 2), lty=2)
    abline(v=c(2, 3)*p/n, lty=2)
    mtext("Influence", outer=T, cex=1.5)
    }
  } 
  
  else if(type == "scale"){
    
    rstand <- model$resid(type = "standard")
    rstudent <- model$resid(type = "student")
    rstudent_abs <- abs(rstudent)
    p <- ncol(model$Zm) + 1
    n <- nrow(model$Zm)
    
    if(single == "scllinear"){
      #restore plotting device par on exit
      par_orig <- par(mfrow = c(1, 1), oma=c(0,0,2,0))
      on.exit(par(par_orig))
      plot(model$fitted_scale, rstand, xlab = "fitted scale", ylab = "standardized residuals",
           main = "fitted scale vs. standardized residuals")
    } else if (single == "zhatvals"){
      #restore plotting device par on exit
      par_orig <- par(mfrow = c(1, 1), oma=c(0,0,2,0))
      on.exit(par(par_orig))
      plot(model$zhatvalues, type = "p", xlab = "observations", ylab = "hatvalues",
           main = "Hatvalues scale predictor")
      abline(h=c(2, 3)*p/n, col="red", lty=2)
    } else if (single == "spread"){
      #restore plotting device par on exit
      par_orig <- par(mfrow = c(1, 1))
      on.exit(par(par_orig))
    
      plot(fitted_vals, rstudent_abs, xlab="Fitted Values", 
           ylab="Absolute Studentized Residuals", main="Spread-Level")
      regline <- lm(rstudent_abs~fitted_vals)
      abline(regline, col="red")
    } else {
      #restore plotting device par on exit
      par_orig <- par(mfrow = c(2, 2), oma=c(0,0,2,0))
      on.exit(par(par_orig))
      plot(model$fitted_scale, rstand, xlab = "fitted scale", ylab = "standardized residuals",
           main = "fitted scale vs. standardized residuals")
      plot(fitted_vals, rstudent_abs, xlab="Fitted Values", 
           ylab="Absolute Studentized Residuals", main="Spread-Level")
      regline <- lm(rstudent_abs~fitted_vals)
      abline(regline, col="red")
      plot(model$zhatvalues, type = "p", xlab = "observations", ylab = "hatvalues",
           main = "Hatvalues scale predictor")
      abline(h=c(2, 3)*p/n, col="red", lty=2)
    }
    
  } else {
    #basic diagnostic plots
    rstand <- model$resid(type = "standard")
    Rsq <- 1 - (model$SSE/model$TSS)
    
    if(single=="fitactual"){
      plot(fitted_vals, y_vals, main = "Fitted vs. actual Values",
           xlab="Fitted Values", ylab="Actual Values")
      legend("topleft", legend=paste("R^2 =", round(Rsq, digits=2)))
      abline(0, 1, col="red")
    } else if(single=="fitres"){
      plot(fitted_vals, rstand, main="Fitted vs. Residuals", xlab = "Fitted Values",
           ylab = "Stand. Residuals")
      lo <- loess(rstand ~ fitted_vals)
      xl <- seq(min(fitted_vals),max(fitted_vals), (max(fitted_vals) - min(fitted_vals))/1000)
      lines(xl, predict(lo,xl), col='red', lwd=2)
    } else if(single=="qq"){
      #calculate theoretical standard distribution quantiles
      ord <- order(rstand)
      rstand_ord <- rstand[ord]
      n <- length(rstand)
      P <- ppoints(n)
      z <- qnorm(P)
      plot(z, rstand_ord, main = "Normal Q-Q Plot of Standardized Residuals", xlab = "Theoretical Quantiles", 
           ylab = "Sample Quantiles", frame.plot=T)
      #diagonal of perfect fit
      coef <- coef(lm(rstand_ord~z))
      a <- coef[1]
      b <- coef[2]
      abline(a, b, col="red")
      #confidence interval bands
      alpha <- 0.95
      Z <- qnorm(1-(1-alpha)/2)
      se <- (b/dnorm(z))*sqrt(P*(1-P)/n)
      vals <- a+b*z
      upper <- vals + Z*se
      lower <- vals - Z*se
      lines(z,upper,lty=2,col="red")
      lines(z,lower,lty=2,col="red")
      legend("topleft", legend=c("0.95 Confidence Interval"),col="red", lty=2, cex=0.7)
    } else if(single=="density"){
      xgrid <- seq(-max(fitted_vals)*1.5, max(fitted_vals)*1.5, by=0.01)
      fitnorm <- dnorm(xgrid, 0, sd(resid_reg))
      plot(density(resid_reg), main = "Kernel Density Estimation", ylab="Density", 
           frame.plot=T, ylim=c(0, max(c(fitnorm, density(resid_reg)$y))))
      rug(resid_reg, col = 'red')
      lines(xgrid, fitnorm, lty=2, col="red")
      legend("topleft", legend=c("Normal Curve"), col="red", lty=2, cex=0.7)
    } else {
    
      #restore plotting device parameters on exit
      par_orig <- par(mfrow = c(2, 2))
      on.exit(par(par_orig))
      
      #Fitted against y
      plot(fitted_vals, y_vals, main = "Fitted vs. Actual Values", 
           xlab=paste0("Fitted Values ", "(R^2 = ", round(Rsq, digits=2), ")"), 
           ylab="Actual Values")
      abline(0, 1, col="red")
      
      #Fitted against Residuals
      plot(fitted_vals, rstand, main="Fitted vs. Residuals", xlab = "Fitted Values",
           ylab = "Stand. Residuals")
      lo <- loess(rstand ~ fitted_vals)
      xl <- seq(min(fitted_vals),max(fitted_vals), (max(fitted_vals) - min(fitted_vals))/1000)
      lines(xl, predict(lo,xl), col='red', lwd=2)
      
      #QQ normal
      #calculation of theoretical normal distribution quantiles
      ord <- order(rstand)
      rstand_ord <- rstand[ord]
      n <- length(rstand)
      P <- ppoints(n)
      z <- qnorm(P)
      plot(z, rstand_ord, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", 
           ylab = "Sample Quantiles", frame.plot=T)
      #diagonal line of perfect fit
      coef <- coef(lm(rstand_ord~z))
      a <- coef[1]
      b <- coef[2]
      abline(a, b, col="red")
      #confidence interval bands
      alpha <- 0.95
      Z <- qnorm(1-(1-alpha)/2)
      se <- (b/dnorm(z))*sqrt(P*(1-P)/n)
      vals <- a+b*z
      upper <- vals + Z*se
      lower <- vals - Z*se
      lines(z,upper,lty=2,col="red")
      lines(z,lower,lty=2,col="red")
      
      #Kernel density
      #fitted normal distribution as comparison
      xgrid <- seq(-max(fitted_vals)*1.5, max(fitted_vals)*1.5, by=0.01)
      fitnorm <- dnorm(xgrid, 0, sd(resid_reg))
      plot(density(resid_reg), main = "Kernel Density Estimation", ylab="Density", 
           frame.plot=T, ylim=c(0, max(c(fitnorm, density(resid_reg)$y))))
      rug(resid_reg, col = 'red')
      lines(xgrid, fitnorm, lty=2, col="red")
    }
  }
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides a plot with fitted vs. actual values. This
#' allows to observe the goodness of fit of the model. As an additional
#' measure, R-squared is displayed.
#' 
#' #' This is a wrapper for `plot(model, single="fitactual")`, where `model` is an 
#' object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the fit plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' fitactual.plot(model)
#'
#' @export

fitactual.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, single="fitactual")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides a QQ-normal plot of the residuals.
#' 
#' #' This is a wrapper for `plot(model, single="qq")`, where `model` is an 
#' object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the QQ plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' qq.plot(model)
#'
#' @export

qq.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, single="qq")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides a kernel density estimation of the residuals.
#' 
#' #' This is a wrapper for `plot(model, single="density")`, where `model` is an 
#' object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the Kernel Density Estimation plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' density.plot(model)
#'
#' @export

density.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, single="density")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides a plot showing Fitted Values vs. Residuals.
#' 
#' #' This is a wrapper for `plot(model, single="fitres")`, where `model` is an 
#' object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the Fitted vs. Residuals plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' fitres.plot(model)
#'
#' @export

fitres.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, single="fitres")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides plots for checking fitted vs studentized residuals 
#' with a representation of cook's distance to observe the influence of 
#' single data points.
#' 
#' #' This is a wrapper for `plot(model, type="influence")`, where `model` is an 
#' object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the influence plots.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' influence.plot(model)
#'
#' @export

influence.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, type="influence")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides a plot showing the Cook's Distance Value 
#' for each observation.
#' 
#' This is a wrapper for `plot(model, type="influence", single="cooks")`, 
#' where `model` is an object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the Cook's Distance plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' cooks.plot(model)
#'
#' @export

cooks.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, type="influence", single="cooks")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides a plot showing the Hat-Value 
#' for each observation.
#' 
#' This is a wrapper for `plot(model, type="influence", single="hatvalues")`, 
#' where `model` is an object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the Hat-Values level plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' hats.plot(model)
#'
#' @export

hats.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, type="influence", single="hatvalues")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides a plot for checking fitted vs studentized residuals 
#' with a representation of Cook's Distance to observe the influence of 
#' single data points.
#' 
#' This is a wrapper for `plot(model, type="influence", single="combined")`, 
#' where `model` is an object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the outliers plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' outliers.plot(model)
#'
#' @export

influence.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, type="influence", single="combined")
}

#' Plot function for the `LocationScaleRegressionWLS` model class
#' 
#' This function provides plots for checking fitted values vs absolute 
#' studentized residuals with regression line to observe the spread-level 
#' and assess constant error variance.
#' 
#' This is a wrapper for `plot(model, type="spread")`, where `model` is an 
#' object of class `LocationScaleRegressionWLS`
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' 
#' @return
#' Shows the spread level plot.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' influence.plot(model)
#'
#' @export

spread.plot <- function(model){
  #only allow correctly specified models as input
  if(!inherits(model, "LocationScaleRegressionWLS")){
    stop("model must be of class LocationScaleRegressionWLS")
  }
  plot(model, type="spread")
}


#' pdf plot function for the `LocationScaleRegressionWLS` model class
#'
#' This function provides plots of the probability distribution function
#' (pdf). This function allows either plotting of the fitted distributions 
#' for up to eight observations or predictive distribution for new data
#' with the fitted parameters in the model object
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#' @param newdata optional, new data for which the predictive pdfs are 
#' to be plotted
#' @param which.obs a vector of up to lentgth eight indicating the case 
#' numbers of the observations for which fitted distributions a (either in
#' the original data used to fit the model or in the new data)
#' @param sdrange vector with two entries, indicating which range of x-values 
#' to the left and right of the smallest and highest predicted mean should be 
#' considered. Formulated in terms of standard deviations of the distribution
#' of the respective mean. The default lower(higher) end is 4(4) standard 
#' deviations to the left(right) of the smallest(largest) mean. 
#' @param single logical (FALSE by default), indicating whether the pdf's 
#' should all be displayed in the same plot. The respecitive pdf's are then
#' color-coded
#' @param obs.names optional, a vector of names for the observations. Will
#' be displayed in the title for single.plot = F and in the legend for
#' single.plot = T
#' @param grid logical, indicating whether the plots should be arranged in 
#' a grid. The default is TRUE for when length(which.obs) > 3. If set to FALSE, 
#' the plots are displayed in one column.
#' @param legend logical, indicating whether a legend should be included 
#' when single.plot = T
#' @param legend_position character indicating where legend should be placed
#' (if single.plot = T). Should be one of “bottomright”, “bottom”, “bottomleft”, 
#' “left”, “topleft”, “top”, “topright”, “right”, “center”. The default is
#' "topright", but this might sometimes overlap with the plot, with the need
#' to change the position interactively
#' @param legend_inset 2-entry vector, indicating how war the legend should 
#' be placed away from the plot's borders (default is c(0.05, 0.05))
#' @param colorblind logical, indicating whether or not the color palette 
#' used for single.plot = T should be colorblind-friendly. The default is 
#' TRUE and thus a sequential palette. If desired, one can use the
#' qualitative palette by setting colorblind = FALSE
#'
#' @return
#' Returns the specified plot
#'
#' @examples
#' library(rworldmap)
#' countries <- countryExData[, 2]
#' n <- length(countries)
#' x <- runif(n)
#' y <- 5 * x + rnorm(n, sd = exp(1.2 * x))
#' model <- LocationScaleRegressionWLS$new(y ~ x, ~ x,
#' +           labels = countries)
#' optim_params(model, method = "ensemble")
#' pdf.plot(model, which.obs = c(10, 100, 136, 131, 22), 
#' +           single.plot = TRUE, legend = TRUE)
#'
#' @export

pdf.plot <- function(model, newdata = NULL, which.obs = NULL, sdrange = c(4,4),
                     single.plot = FALSE, obs.names = NULL,
                     grid = TRUE, legend = FALSE, colorblind = TRUE,
                     legend_position = "topright", legend_inset = c(0.05, 0.05)){
  
  if(is.null(which.obs) & is.null(newdata)){
    stop("please supply a vector which.obs")
  }
  
  if(is.null(which.obs) & !is.null(newdata)){
    if(nrow(newdata) > 8){
      stop("can't produce pdf plots for more than 8 observations. Please supply which.obs")
    } else {
      which.obs <- 1:nrow(newdata)
    }
  }
  
  if(!is.null(model$labels) & !is.null(obs.names) & is.null(newdata)){
    warning("you gave obs.names, but the model already has labels. Using supplied obs.names.")
  }
    
  nop <-  length(which.obs)
  
  #tests for correct input
  if(!is.null(obs.names) & nop != length(obs.names)){
    stop("length of obs.names does not match length of which.obs")
  }
  if(nop > 8) stop("can't produce pdf plots for more than 8 observations")
  if(!single.plot & legend) warning("single.plot = F does not have a legend option")
  
  
  #make obs names in case none were given
  if(is.null(obs.names)){
    if(!is.null(model$labels) & is.null(newdata)){
      obs.names <- model$labels[which.obs]
    } else {
      obs.names <- sapply(which.obs, 
                        FUN = function(x) paste0("obs", x))
    }
  }
  
  #get data and parameters for plotting
  if (!is.null(newdata)){
    params <- as.matrix(predict(model, newdata))
    if(nop == 1){
      mus <- params[1]
      sigmas <- params[2]
    } else {
      params <- params[which.obs, ]
      mus <- params[,1]
      sigmas <- params[,2]
    }
    
  } else {
    mus <- model$fitted_location[which.obs]
    sigmas <- model$fitted_scale[which.obs]
    params <- cbind(mus, sigmas)
  }
  

  ####find x range for the plots
  imin <- which(mus == min(mus))
  imax <- which(mus == max(mus))
  
  lower <- mus[imin] - sdrange[1] * sigmas[imin]
  upper <- mus[imax] + sdrange[2] * sigmas[imax]

  range <- seq(lower, upper, length = 10000)
  
  dfun <- function(x) dnorm(range, x[1], x[2])
  dens <- apply(params, MARGIN = 1,
                FUN = dfun)
  
  if(grid & nop > 3){
    #plots: get grid parameters
    paru <- ceiling(sqrt(nop))
    parl <- floor(sqrt(nop))
    if ((paru * parl) < nop) parl <- paru
    
    #restore plotting device par on exit
    par_orig <- par(mfrow = c(paru, parl))
    on.exit(par(par_orig))
    
  } else {
    par_orig <- par(mfrow = c(nop, 1))
    on.exit(par(par_orig))
  }
  
  if (!single.plot){
    #produce plots
    for (j in 1:ncol(dens)){
      plot(range, dens[,j], type = "l",
           ylab = "pdf, f(y)", xlab = "y",
           main = paste0("\"", obs.names[j], "\"\U003A",
                         " mu = ", round(mus[j], 3), 
                         ", sigma = ", round(sigmas[j], 3))) 
    }
    
  } else {
    par_orig <- par(mfrow = c(1,1))
    on.exit(par(par_orig))
    
    if (colorblind){
      cols <- c(rgb(34/255,94/255,168/255), rgb(127/255,205/255,187/255), 
                rgb(237/255,248/255,177/255), rgb(65/255,182/255,196/255), 
                rgb(12/255,44/255,132/255), rgb(199/255,233/255,180/255),
                rgb(29/255,145/255,192/255), rgb(255/255,255/255,217/255))
    } else {
      cols <- c(rgb(27/255,158/255,119/255), rgb(217/255,95/255,2/255),
                rgb(117/255,112/255,179/255), rgb(231/255,41/255,138/255),
                rgb(102/255,166/255,30/255), rgb(230/255,171/255,2/255),
                rgb(166/255,118/255,29/255), rgb(102/255,102/255,102/255))
    }
    
    maxden <- max(dens)
    plot(range, dens[,1], type = "l",
         ylim = c(0, maxden), col = cols[1], 
         lwd = 2.5, ylab = "pdf, f(y)", xlab = "y")
    if(nop>1){
      for (j in 2:ncol(dens)){
        lines(range, dens[,j], col = cols[j], lwd = 2.5)
      }
    }
    
    if(legend){
      legend(legend_position, 
             legend = paste0(obs.names, ": mu = ", round(mus, digits = 1), 
                             ", sigma = ", round(sigmas, digits = 1)), 
             col = cols, 
             pch = "-", 
             #bty = "n", 
             pt.cex = 2.5, 
             cex = 1.2, 
             text.col = "black", 
             horiz = F, 
             inset = legend_inset)
    }
  }
}

#' Summary function for the `LocationScaleRegressionWLS` model class
#'
#' This function provides a summary for checking the results of a fitted
#' LocationScaleRegressionWLS object.
#'
#' @param model A [`LocationScaleRegressionWLS`] object.
#'
#' @return
#' prints a statistical summary of the input LocationScaleRegressionWLS object
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionWLS$new(y ~ 1)
#' optim_params(model)
#' summary(model)
#'
#' @export

summary.LocationScaleRegressionWLS <- function(model){
  if(any(model$grad_gamma() > 0.01)) stop("can't return summary as model was 
                                          not fitted. Run optim_params(model) first.")
  
  dep_name <- all.vars(model$formula_location)[1]
  loc_names <- c("(I)", all.vars(model$formula_location)[-1])
  scale_names <- c("(I)", all.vars(model$formula_scale))
  
  #get summary values for betas
  betas <- round(model$beta, digits = 4)
  beta_sds <- round(sqrt(diag(solve(model$hess_beta()))), digits = 4)
  beta_ts <- round(betas/beta_sds, digits = 4)
  beta_ps <- as.character(round(1 - round(pt(beta_ts, 
                                       df = (model$n - sum(model$hatvalues))), 
                                       digits = 4)), digits = 4)
  beta_ps[beta_ps == 0] <- "<2e-16"
  
  gammas <- round(model$gamma, digits = 4)
  gamma_sds <- round(sqrt(diag(solve(model$hess_gamma()))), digits = 4)
  gamma_ts <- round(gammas/gamma_sds, digits = 4)
  gamma_ps <- as.character(round(1 - round(pt(gamma_ts, 
                                        df = (model$n - sum(model$hatvalues))), 
                                        digits = 4)), digits = 4)
  gamma_ps[gamma_ps == 0] <- "<2e-16"
  
  weights <- diag(model$W)
  weights_qu <- round(summary(weights), digits = 4)
  
  cat("Call: \n")
  cat("LocationScaleRegressionWLS$new(location = ", paste(dep_name, paste(loc_names[-1], collapse = " + "), sep = " ~ "),",\n")
  cat("                               scale    =  ", paste("", paste(scale_names[-1], collapse = " + "), sep = " ~ "), ")\n")
  cat("\n")
  cat("Family: gaussian\n")
  cat("Link functions: location - identity; scale - log \n")
  cat("\n")
  cat("weights:  \n")
  cat("  Min      1Q     Median    3Q      Max\n")
  cat(weights_qu[1], " ", weights_qu[2], " ",  weights_qu[3], " ",  weights_qu[5], " ",  weights_qu[6], "\n")
  cat("\n")
  cat("beta coefficients:\n")
  cat("                 Estimate     Std. Error      t value      Pr(>|t|)\n")
  for(i in 1:length(betas)){
    cat(loc_names[i], "            ", betas[i],"      ",beta_sds[i], "     ", beta_ts[i], "    ", beta_ps[i], "\n")
  }
  cat("\n")
  cat("gamma coefficients:\n")
  cat("                 Estimate     Std. Error      t value      Pr(>|t|)\n")
  for(i in 1:length(gammas)){
    cat(scale_names[i], "            ", gammas[i],"      ",gamma_sds[i], "     ", gamma_ts[i], "    ", gamma_ps[i], "\n")
  }
  cat("\n")
  #cat("Residual standard error: ", model$SSE/(model$n - sum(model$hatvalues)), "\n")
  cat("AIC: ", model$AIC, "\n")
  cat("R-squared: ", 1 - (model$SSE/model$TSS))
}
