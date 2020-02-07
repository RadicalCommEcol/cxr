#' Internal, retrieve parameters from the vector returned by the optimization procedures
#'
#' @param optim_par 1d vector, the result of an optimization method
#' @param lambda_length either 0 (lambda not fit) or 1
#' @param alpha_length either 0 (alpha not fit) or a positive number
#' @param lambda_cov_length either 0 (lambda_cov not fit) or a positive number
#' @param alpha_cov_length either 0 (alpha_cov not fit) or a positive number
#'
#' @return list with elements "lambda", "alpha", "lambda_cov", "alpha_cov", "sigma". 
#' If one of these elements is not present, returns NULL.
#' @noRd
cxr_retrieve_params <- function(optim_par, 
                                lambda_length = 0, 
                                alpha_length = 0, 
                                lambda_cov_length = 0, 
                                alpha_cov_length = 0){
  
  lambda <- NULL
  alpha <- NULL
  lambda_cov <- NULL
  alpha_cov <- NULL
  sigma <- optim_par[length(optim_par)]
  
  pos <- 1
  
  if(lambda_length == 1){
    lambda <- optim_par[pos]
    pos <- pos + 1
  }
  
  if(lambda_cov_length > 0){
    lambda_cov <- optim_par[pos:(pos+lambda_cov_length-1)]
    pos <- pos + lambda_cov_length
  }
  
  if(alpha_length > 0){
    alpha <- optim_par[pos:(pos+alpha_length-1)]
    pos <- pos + alpha_length
  }
  
  if(alpha_cov_length > 0){
    alpha_cov <- optim_par[pos:(pos+alpha_cov_length-1)]
  }
  
  return(list(lambda = lambda, alpha = alpha, lambda_cov = lambda_cov, alpha_cov = alpha_cov, sigma = sigma))
  
}

