
#' Retrieve parameters from the vector returned by the optimization procedures
#'
#' @param optim.params 1d vector, the result of an optimization method
#' @param param.list character vector, which parameters are present. Possible elements are "lambda", "lambda.cov", "alpha", "alpha.cov".
#' @param alpha.length if alpha is to be retrieved, its length
#' @param alpha.cov.length if alpha.cov is to be retrieved, its length
#' @param num.competitors how many competitor species
#' @param num.covariates how many covariates
#'
#' @return list with elements "lambda", "alpha", "lambda.cov", "alpha.cov", "sigma". If one of these elements is not present, returns NULL.
#' @export
RetrieveParams <- function(optim.params, param.list, alpha.length, alpha.cov.length, num.competitors, num.covariates){
  
  lambda <- NULL
  alpha <- NULL
  lambda.cov <- NULL
  alpha.cov <- NULL
  sigma <- optim.params[length(optim.params)]
  
  pos <- 1
  
  if("lambda" %in% param.list){
    lambda <- optim.params[pos]
    pos <- pos + 1
  }
  
  if("lambda.cov" %in% param.list){
    lambda.cov <- optim.params[pos:(pos+num.covariates-1)]
    pos <- pos + num.covariates
  }
  
  if("alpha" %in% param.list){
    alpha <- optim.params[pos:(pos+alpha.length-1)]
    pos <- pos + alpha.length
  }
  
  if("alpha.cov" %in% param.list){
    alpha.cov <- optim.params[pos:(pos+alpha.cov.length-1)]
  }
  
  return(list(lambda = lambda, alpha = alpha, lambda.cov = lambda.cov, alpha.cov = alpha.cov, sigma = sigma))
  
}

