
#' Retrieve parameters from the vector returned by the optimization procedures
#'
#' @param optim_params 1d vector, the result of an optimization method
#' @param lambda_length either 0 (lambda not fit) or a positive number
#' @param effect_length either 0 (effect not fit) or a positive number
#' @param response_length either 0 (response not fit) or a positive number
#' @param lambda_cov_length either 0 (lambda_cov not fit) or a positive number
#' @param effect_cov_length either 0 (effect_cov not fit) or a positive number
#' @param response_cov_length either 0 (response_cov not fit) or a positive number
#'
#' @return list with elements "lambda", "effect", "response", "lambda_cov",
#' "effect_cov", "response_cov", "sigma". 
#' If one of these elements is not present, returns NULL.
#' @export
cxr_retrieve_er_params <- function(optim_params, 
                                lambda_length = 0,
                                effect_length = 0,
                                response_length = 0,
                                lambda_cov_length = 0,
                                effect_cov_length = 0,
                                response_cov_length = 0
                                ){
  
  lambda <- NULL
  effect <- NULL
  response <- NULL
  lambda_cov <- NULL
  effect_cov <- NULL
  response_cov <- NULL  
  sigma <- optim_params[length(optim_params)]
  
  pos <- 1
  
  if(lambda_length > 0){
    lambda <- optim_params[pos:(pos+lambda_length-1)]
    pos <- pos + lambda_length
  }
  
  if(lambda_cov_length > 0){
    lambda_cov <- optim_params[pos:(pos+lambda_cov_length-1)]
    pos <- pos + lambda_cov_length
  }
  
  if(effect_length > 0){
    effect <- optim_params[pos:(pos+effect_length-1)]
    pos <- pos + effect_length
  }
  
  if(effect_cov_length > 0){
    effect_cov <- optim_params[pos:(pos+effect_cov_length-1)]
    pos <- pos + effect_cov_length
  }
  
  if(response_length > 0){
    response <- optim_params[pos:(pos+response_length-1)]
    pos <- pos + response_length
  }
  
  if(response_cov_length > 0){
    response_cov <- optim_params[pos:(pos+response_cov_length-1)]
    pos <- pos + response_cov_length
  }
  
  return(list(lambda = lambda, 
              effect = effect,
              response = response, 
              lambda_cov = lambda_cov,
              effect_cov = effect_cov,
              response_cov = response_cov, 
              sigma = sigma))
  
}

