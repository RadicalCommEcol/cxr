
#' Join parameters in a 1d vector
#' 
#' Generate a 1d vector from a series of parameters in a certain order. It also returns the same vector
#' for lower and upper bounds. This function is intended to work with
#' parameters for a single species (i.e. a single lambda value, etc). Note that lambda_cov and alpha_cov must be consistent
#' with num_covariates.
#'
#' @param init_lambda numeric, lambda
#' @param init_sigma numeric, sigma
#' @param init_alpha 1d vector, interaction coefficients over the species
#' @param init_lambda_cov 1d vector, initial values for lambda_cov
#' @param init_alpha_cov 1d vector, initial values for alpha_cov
#' @param lower_lambda lower bound for lambda
#' @param upper_lambda upper bound for lambda
#' @param lower_sigma lower bound for sigma
#' @param upper_sigma upper bound for sigma
#' @param lower_alpha lower bound for alpha
#' @param upper_alpha upper bound for alpha
#' @param lower_lambda_cov lower bound for lambda_cov
#' @param upper_lambda_cov upper bound for lambda_cov
#' @param lower_alpha_cov lower bound for alpha_cov
#' @param upper_alpha_cov upper bound for alpha_cov
#'
#' @return list with three 1d vectors, ready for passing to the optim methods, consistent with the functions model_BH1-5 
#' @export
cxr_init_params <- function(init_lambda = NULL,
                       init_sigma = 0,
                       init_alpha = NULL,
                       init_lambda_cov = NULL,
                       init_alpha_cov = NULL,
                       lower_lambda = 1,
                       upper_lambda = 1e5,
                       lower_sigma = 1e-5,
                       upper_sigma = 1e5,
                       lower_alpha = 1e-5,
                       upper_alpha = 1e5,
                       lower_lambda_cov = 1e-5,
                       upper_lambda_cov = 1e5,
                       lower_alpha_cov = 1e-5,
                       upper_alpha_cov = 1e5
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
  
  # alpha value/matrix
  if(!is.null(init_alpha)){
    init_par <- c(init_par,init_alpha)
    lower_bounds <- c(lower_bounds,rep(lower_alpha,length(init_alpha)))
    upper_bounds <- c(upper_bounds,rep(upper_alpha,length(init_alpha)))                  
  }
  
  # effect of covariates on alpha
  if(!is.null(init_alpha_cov)){
    init_par <- c(init_par,init_alpha_cov)
    lower_bounds <- c(lower_bounds,rep(lower_alpha_cov,length(init_alpha_cov)))
    upper_bounds <- c(upper_bounds,rep(upper_alpha_cov,length(init_alpha_cov)))
  }
  
  # sigma goes at the end
  init_par <- c(init_par,init_sigma)
  lower_bounds <- c(lower_bounds,lower_sigma)
  upper_bounds <- c(upper_bounds,upper_sigma)
  
  # lower_bounds <- ifelse(lower_bounds > 0, lower_bounds, 1e-10)
  
  return(list(init_par = init_par, lower_bounds = lower_bounds, upper_bounds = upper_bounds))
}



