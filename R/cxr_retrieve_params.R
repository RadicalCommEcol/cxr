
#' Retrieve parameters from the vector returned by the optimization procedures
#'
#' @param optim_params 1d vector, the result of an optimization method
#' @param lambda_length either 0 (lambda not fit) or 1
#' @param alpha_length either 0 (alpha not fit) or a positive number
#' @param lambda_cov_length either 0 (lambda_cov not fit) or a positive number
#' @param alpha_cov_length either 0 (alpha_cov not fit) or a positive number
#'
#' @return list with elements "lambda", "alpha", "lambda_cov", "alpha_cov", "sigma". If one of these elements is not present, returns NULL.
#' @export
cxr_retrieve_params <- function(optim_params, lambda_length, alpha_length, lambda_cov_length, alpha_cov_length){
  
  lambda <- NULL
  alpha <- NULL
  lambda_cov <- NULL
  alpha_cov <- NULL
  sigma <- optim_params[length(optim_params)]
  
  pos <- 1
  
  if(lambda_length == 1){
    lambda <- optim_params[pos]
    pos <- pos + 1
  }
  
  if(lambda_cov_length > 0){
    lambda_cov <- optim_params[pos:(pos+lambda_cov_length-1)]
    pos <- pos + lambda_cov_length
  }
  
  if(alpha_length > 0){
    alpha <- optim_params[pos:(pos+alpha_length-1)]
    pos <- pos + alpha_length
  }
  
  if(alpha_cov_length > 0){
    alpha_cov <- optim_params[pos:(pos+alpha_cov_length-1)]
  }
  
  return(list(lambda = lambda, alpha = alpha, lambda_cov = lambda_cov, alpha_cov = alpha_cov, sigma = sigma))
  
}

