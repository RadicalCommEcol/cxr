
#' Join parameters in a 1d vector
#' 
#' Generate a 1d vector from a series of parameters in a certain order. It also returns the same vector
#' for lower and upper bounds. This function is intended to work with
#' parameters for a single species (i.e. a single lambda value, etc). Note that lambda_cov and alpha_cov must be consistent
#' with num_covariates.
#'
#' @param init_lambda numeric, lambda
#' @param init_sigma numeric, sigma
#' @param init_effect 1d vector, 
#' @param init_response 1d vector, 

#' @param init_lambda_cov 1d vector, initial values for lambda_cov
#' @param init_effect_cov 1d vector
#' @param init_response_cov 1d vector

#' @param lower_lambda lower bound for lambda
#' @param upper_lambda upper bound for lambda
#' @param lower_sigma lower bound for sigma
#' @param upper_sigma upper bound for sigma
#' @param lower_effect 
#' @param upper_effect 
#' @param lower_response 
#' @param upper_response
#' @param lower_lambda_cov lower bound for lambda_cov
#' @param upper_lambda_cov upper bound for lambda_cov
#' @param lower_effect_cov 
#' @param upper_effect_cov 
#' @param lower_response_cov 
#' @param upper_response_cov
#'
#' @return list with three 1d vectors, ready for passing to the optim methods, consistent with the functions model_BH1-5 
#' @export
cxr_init_er_params <- function(init_lambda = NULL,
                       init_sigma = 0,
                       init_effect = NULL,
                       init_response = NULL,
                       init_lambda_cov = NULL,
                       init_effect_cov = NULL,
                       init_response_cov = NULL,
                       lower_lambda = 1,
                       upper_lambda = 1e5,
                       lower_sigma = 1e-5,
                       upper_sigma = 1e5,
                       lower_effect = 1e-5,
                       upper_effect = 1e5,
                       lower_response = 1e-5,
                       upper_response = 1e5,
                       lower_lambda_cov = 1e-5,
                       upper_lambda_cov = 1e5,
                       lower_effect_cov = 1e-5,
                       upper_effect_cov = 1e5,
                       lower_response_cov = 1e-5,
                       upper_response_cov = 1e5
                       ){
  init_par <- NULL
  lower_bounds <- NULL
  upper_bounds <- NULL
  
  # if lambda is not null, it goes first

  if(!is.null(init_lambda)){
    init_par <- init_lambda
    lower_bounds <- rep(lower_lambda,length(init_lambda))
    upper_bounds <- rep(upper_lambda,length(init_lambda))
  }
  
  # effect of covariates on lambda
  if(!is.null(init_lambda_cov)){
    init_par <- c(init_par,init_lambda_cov)
    lower_bounds <- c(lower_bounds,rep(lower_lambda_cov,length(init_lambda_cov)))
    upper_bounds <- c(upper_bounds,rep(upper_lambda_cov,length(init_lambda_cov)))
  }
  
  # competitive effect values
  if(!is.null(init_effect)){
    init_par <- c(init_par,init_effect)
    lower_bounds <- c(lower_bounds,rep(lower_effect,length(init_effect)))
    upper_bounds <- c(upper_bounds,rep(upper_effect,length(init_effect)))                  
  }
  
  # effect of covariates on effect
  if(!is.null(init_effect_cov)){
    init_par <- c(init_par,init_effect_cov)
    lower_bounds <- c(lower_bounds,rep(lower_effect_cov,length(init_effect_cov)))
    upper_bounds <- c(upper_bounds,rep(upper_effect_cov,length(init_effect_cov)))
  }
  
  # competitive response values
  if(!is.null(init_response)){
    init_par <- c(init_par,init_response)
    lower_bounds <- c(lower_bounds,rep(lower_response,length(init_response)))
    upper_bounds <- c(upper_bounds,rep(upper_response,length(init_response)))                  
  }
  
  # effect of covariates on response
  if(!is.null(init_response_cov)){
    init_par <- c(init_par,init_response_cov)
    lower_bounds <- c(lower_bounds,rep(lower_response_cov,length(init_response_cov)))
    upper_bounds <- c(upper_bounds,rep(upper_response_cov,length(init_response_cov)))
  }
  
  # sigma goes last
  init_par <- c(init_par,init_sigma)
  lower_bounds <- c(lower_bounds,lower_sigma)
  upper_bounds <- c(upper_bounds,upper_sigma)
  
  return(list(init_par = init_par, lower_bounds = lower_bounds, upper_bounds = upper_bounds))
}



