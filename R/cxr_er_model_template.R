
# 1 - MODEL FAMILY
# Model names are of the form shown below. "family" is the acronym of general family of your model,
# as e.g. BH for Beverton-Holt, LW for Law-Wilkinson, etc.
# these acronyms may also be useful to name non-linear versions of the models, e.g. "BH-NL"

# 2 - PARAMETER FORMS
# All effect-response models have, necessarily, three parameters:
# lambda, e, and r. In addition, the effect of covariates on these parameters
# can be included.

# The template for naming your model 
# can be seen below. Note how a model name includes information about the model family,
# and also about how parameters are modelled.
# In particular, the effects of covariates over each parameter can be included by
# naming the parameter as "global". 

# 3 - FIXED TERMS
# if you have independent estimates of a subset of parameters, you may want
# to keep these as fixed, and only fit those you need to. This can be achieved
# through the argument "fixed_parameters", which is a list with as many components 
# as model parameters (lambda, effect, response, lambda_cov, effect_cov, response_cov),
# where each component should be either NULL if the parameter is to be fit, 
# or the fixed value of the parameter in question.

# 4 - Saving your model
# name the R file with the same name as your model, and for using it within cxr, 
# make it available in the global environment. This is as simple as sourcing the file.

# 5 - adding your model to cxr
# document your model, file a pull_request etc

er_BH_lambdacov_form_effectcov_form_responsecov_form <- function(par,
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
  
  # the model should return a "pred" value
  # a function of lambda, effect, response, lambda_cov, effect_cov, response_cov
  pred <- 0
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
  
}