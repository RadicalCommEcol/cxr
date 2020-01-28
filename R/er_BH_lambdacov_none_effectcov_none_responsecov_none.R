er_BH_lambdacov_none_effectcov_none_responsecov_none <- function(par,
                                                                 fitness,
                                                                 target,
                                                                 density,
                                                                 covariates,
                                                                 fixed_parameters){
  
  num.sp <- nrow(target)
  
  # TODO check which fixed
  
  # lambda.vector <- par[1:num.sp]
  # r.vector <- par[(num.sp+1):(num.sp+num.sp)]
  # e.vector <- par[(num.sp+1+num.sp):(length(par)-1)]
  # sigma <- par[length(par)]
  
  # parameters to fit are all in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,lambda_cov,alpha,alpha_cov,sigma}
  
  # comment or uncomment sections for the different parameters
  # depending on whether your model includes them
  pos <- 1
  
  # if a parameter is passed within the "par" vector,
  # it should be NULL in the "fixed_parameters" list
  if(is.null(fixed_parameters[["lambda"]])){
    lambda <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    lambda <- fixed_parameters[["lambda"]]
  }
  
  # if(is.null(fixed_parameters$lambda_cov)){
  #   lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
  #   pos <- pos + ncol(covariates)
  # }else{
  #   lambda_cov <- fixed_parameters[["lambda_cov"]]
  # }
  
  if(is.null(fixed_parameters[["effect"]])){
    effect <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    effect <- fixed_parameters[["effect"]]
  }
  
  # if(is.null(fixed_parameters[["effect_cov"]])){
  #   effect_cov <- par[pos]
  #   pos <- pos + num.sp - 1
  # }else{
  #   effect_cov <- fixed_parameters[["effect_cov"]]
  # }
  
  if(is.null(fixed_parameters[["response"]])){
    response <- par[pos:(pos + num.sp - 1)]
    pos <- pos + num.sp
  }else{
    response <- fixed_parameters[["response"]]
  }
  
  # if(is.null(fixed_parameters[["response_cov"]])){
  #   response_cov <- par[pos]
  #   pos <- pos + num.sp - 1
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
  
  pred <- lambda.part/ (1+ e.part*r.part )
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
  
}