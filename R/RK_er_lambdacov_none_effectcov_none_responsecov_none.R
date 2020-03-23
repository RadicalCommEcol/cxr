#' Effect response Ricker model without covariate effects
#'
#' Note that, as e and r are not pair-specific, all species parameters are fit in the same function.
#'
#' @param par 1d vector with initial parameters in the order: 
#' lambda,effect,response,sigma.
#' @param fitness 1d vector with fitness observations.
#' @param target matrix with species in rows, observations in columns. Value is 1 if
#' a species is focal for a given observation, 0 otherwise.
#' @param density matrix with species in rows, observations in columns. Value is 
#' density of each sp as neighbour for each observation.
#' @param covariates included for compatibility, not used in this model.
#' @param fixed_parameters optional list specifying values of fixed parameters, 
#' with components "lambda","effect","response".
#'
#' @return log-likelihood value
#' @export
RK_er_lambdacov_none_effectcov_none_responsecov_none <- function(par,
                                                                 fitness,
                                                                 target,
                                                                 density,
                                                                 covariates,
                                                                 fixed_parameters){
  
  num.sp <- nrow(target)
  
  # parameters to fit are all in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,effect,response,sigma}
  
  # comment or uncomment sections for the different parameters
  # depending on whether your model includes them
  pos <- 1
  
  # if a parameter is passed within the "par" vector,
  # it should be NULL in the "fixed_parameters" list
  
  # lambda
  if(is.null(fixed_parameters[["lambda"]])){
    lambda <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    lambda <- fixed_parameters[["lambda"]]
  }
  
  # lambda_cov
  # if(is.null(fixed_parameters$lambda_cov)){
  #   # the covariate effects are more efficient in a matrix form
  #   # with species in rows (hence byrow = T, because by default
  #   # the vector is sorted first by covariates)
  #   lambda_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],nrow = num.sp,byrow = TRUE)
  #   pos <- pos + ncol(covariates)*num.sp
  # }else{
  #   lambda_cov <- fixed_parameters[["lambda_cov"]]
  # }
  
  # effect
  if(is.null(fixed_parameters[["effect"]])){
    effect <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    effect <- fixed_parameters[["effect"]]
  }
  
  # effect_cov
  # if(is.null(fixed_parameters$effect_cov)){
  #   effect_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],nrow = num.sp,byrow = TRUE)
  #   pos <- pos + ncol(covariates)*num.sp
  # }else{
  #   effect_cov <- fixed_parameters[["effect_cov"]]
  # }
  
  # response
  if(is.null(fixed_parameters[["response"]])){
    response <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    response <- fixed_parameters[["response"]]
  }
  
  # response_cov
  # if(is.null(fixed_parameters[["response_cov"]])){
  #   response_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],nrow = num.sp,byrow = TRUE)
  #   pos <- pos + ncol(covariates)*num.sp
  # }else{
  #   response_cov <- fixed_parameters[["response_cov"]]
  # }
  
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  lambda.part <- colSums(lambda*target)
  r.part <- colSums(response*target)
  e.part <- colSums(effect*density)
  
  pred <- lambda.part * exp(r.part * - e.part)
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
  
}