#' Effect response Beverton-Holt model with covariate effects on lambda, effect, and response
#'
#' The function for calculating fecundity given 
#' effect and response values is taken from Godoy et al. (2014). 
#' Note that, as e and r are not pair-specific, all species parameters are fit in the same function.
#'
#' @param par 1d vector with initial parameters in the order: 
#' lambda,lambda_cov,effect,effect_cov,response,response_cov,sigma
#' @param fitness 1d vector with fitness observations
#' @param target matrix with species in rows, observations in columns. Value is 1 if
#' a species is focal for a given observation, 0 otherwise.
#' @param density matrix with species in rows, observations in columns. Value is 
#' density of each sp as neighbour for each observation.
#' @param covariates numeric dataframe or matrix with observations 
#' in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation
#' @param fixed_parameters optional list specifying values of fixed parameters, 
#' with components "lambda","lambda_cov","effect","effect_cov",
#' "response","response_cov".
#'
#' @return log-likelihood value
#' @export
BH_er_lambdacov_global_effectcov_global_responsecov_global <- function(par,
                                                                 fitness,
                                                                 target,
                                                                 density,
                                                                 covariates,
                                                                 fixed_parameters){
  
  num.sp <- nrow(target)
  
  # parameters to fit are all in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,lambda_cov,effect,effect_cov,response,response_cov,sigma}
  
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
  if(is.null(fixed_parameters$lambda_cov)){
    # the covariate effects are more efficient in a matrix form
    # with species in rows (hence byrow = T, because by default
    # the vector is sorted first by covariates)
    lambda_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],nrow = num.sp,byrow = TRUE)
    pos <- pos + ncol(covariates)*num.sp
  }else{
    lambda_cov <- fixed_parameters[["lambda_cov"]]
  }
  
  # effect
  if(is.null(fixed_parameters[["effect"]])){
    effect <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    effect <- fixed_parameters[["effect"]]
  }
  
  # effect_cov
  if(is.null(fixed_parameters$effect_cov)){
    effect_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],nrow = num.sp,byrow = TRUE)
    pos <- pos + ncol(covariates)*num.sp
  }else{
    effect_cov <- fixed_parameters[["effect_cov"]]
  }
  
  # response
  if(is.null(fixed_parameters[["response"]])){
    response <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    response <- fixed_parameters[["response"]]
  }
  
  # response_cov
  if(is.null(fixed_parameters[["response_cov"]])){
    response_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],nrow = num.sp,byrow = TRUE)
    pos <- pos + ncol(covariates)*num.sp
  }else{
    response_cov <- fixed_parameters[["response_cov"]]
  }
  
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  lambda_cov_all <- matrix(0,nrow = num.sp,ncol = length(fitness))
  response_cov_all <- matrix(0,nrow = num.sp,ncol = length(fitness))
  effect_cov_all <- matrix(0,nrow = num.sp,ncol = length(fitness))
  for(i.sp in 1:num.sp){
    for(i.obs in 1:length(fitness)){
      lambda_cov_all[i.sp,i.obs] <- sum(lambda_cov[i.sp,]*covariates[i.obs,])
      response_cov_all[i.sp,i.obs] <- sum(response_cov[i.sp,]*covariates[i.obs,])
      effect_cov_all[i.sp,i.obs] <- sum(effect_cov[i.sp,]*covariates[i.obs,])
      
    }# for i_cov
  }# for i.sp
  
  lambda.part <- colSums(lambda*(1+lambda_cov_all)*target)
  r.part <- colSums(response*(1+response_cov_all)*target)
  e.part <- colSums(effect*(1+effect_cov_all)*density)
  
  pred <- lambda.part/ (1+ e.part*r.part )
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
  
}