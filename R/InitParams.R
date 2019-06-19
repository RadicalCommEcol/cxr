
#' Join parameters in a 1d vector
#' 
#' Generate a 1d vector from a series of parameters in a certain order. It also returns the same vector
#' for lower and upper bounds. This function is intended to work with
#' parameters for a single species (i.e. a single lambda value, etc). Note that lambda.cov and alpha.cov must be consistent
#' with num.covariates.
#'
#' @param init.lambda numeric, lambda
#' @param init.sigma numeric, sigma
#' @param init.alpha 1d vector, interaction coefficients over the species
#' @param init.lambda.cov 1d vector, initial values for lambda.cov
#' @param init.alpha.cov 1d vector, initial values for alpha.cov
#' @param lower.lambda lower bound for lambda
#' @param upper.lambda upper bound for lambda
#' @param lower.sigma lower bound for sigma
#' @param upper.sigma upper bound for sigma
#' @param lower.alpha lower bound for alpha
#' @param upper.alpha upper bound for alpha
#' @param lower.lambda.cov lower bound for lambda.cov
#' @param upper.lambda.cov upper bound for lambda.cov
#' @param lower.alpha.cov lower bound for alpha.cov
#' @param upper.alpha.cov upper bound for alpha.cov
#' @param num.competitors number of competitors
#' @param num.covariates number of covariates
#'
#' @return list with three 1d vectors, ready for passing to the optim methods, consistent with the functions BH_1-5 
#' @export
InitParams <- function(init.lambda = NULL,
                       init.sigma = 0,
                       init.alpha = NULL,
                       init.lambda.cov = NULL,
                       init.alpha.cov = NULL,
                       lower.lambda = 1,
                       upper.lambda = 1e5,
                       lower.sigma = 1e-5,
                       upper.sigma = 1e5,
                       lower.alpha = 1e-5,
                       upper.alpha = 1e5,
                       lower.lambda.cov = 1e-5,
                       upper.lambda.cov = 1e5,
                       lower.alpha.cov = 1e-5,
                       upper.alpha.cov = 1e5,
                       num.competitors,
                       num.covariates
){
  init.par <- NULL
  lower.bounds <- NULL
  upper.bounds <- NULL
  
  # if lambda is not null, it goes first
  if(!is.null(init.lambda)){
    init.par <- init.lambda
    lower.bounds <- rep(lower.lambda,length(init.lambda))
    upper.bounds <- rep(upper.lambda,length(init.lambda))
  }
  
  # effect of covariates on lambda
  if(!is.null(init.lambda.cov)){
    init.par <- c(init.par,init.lambda.cov)
    lower.bounds <- c(lower.bounds,rep(lower.lambda.cov,length(init.lambda.cov)))
    upper.bounds <- c(upper.bounds,rep(upper.lambda.cov,length(init.lambda.cov)))
  }
  
  # alpha value/matrix
  if(!is.null(init.alpha)){
    init.par <- c(init.par,init.alpha)
    lower.bounds <- c(lower.bounds,rep(lower.alpha,length(init.alpha)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha,length(init.alpha)))                  
  }
  
  # effect of covariates on alpha
  if(!is.null(init.alpha.cov)){
    init.par <- c(init.par,init.alpha.cov)
    lower.bounds <- c(lower.bounds,rep(lower.alpha.cov,length(init.alpha.cov)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha.cov,length(init.alpha.cov)))
  }
  
  # sigma goes at the end
  init.par <- c(init.par,init.sigma)
  lower.bounds <- c(lower.bounds,lower.sigma)
  upper.bounds <- c(upper.bounds,upper.sigma)
  
  lower.bounds <- ifelse(lower.bounds > 0, lower.bounds, 1e-10)
  
  return(list(init.par = init.par, lower.bounds = lower.bounds, upper.bounds = upper.bounds))
}



