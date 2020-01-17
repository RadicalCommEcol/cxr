
# 1 - MODEL FAMILY
# Model names are of the form shown below. "family" is the acronym of general family of your model,
# as e.g. BH for Beverton-Holt, LW for Law-Wilkinson, etc.

# 2 - PARAMETER FORMS
# every model has, at least, a lambda parameter. Other potential parameters are 
# pairwise interactions (alpha), and the effects of covariates over lambda (lambda_cov)
# and over alpha (alpha_cov).

# How these parameters are modelled is important to differentiate models,
# and so should be reflected in the function name. The template for naming your model 
# can be seen below. Note how a model name includes information about the model family,
# and also about how parameters are modelled.

# alpha_form can be one of "alpha_none", "alpha_global", or "alpha_pairwise" depending
# on whether and how you include pairwise interactions in your model.
# "global" means a single alpha value for every interaction, and "pairwise"
# means specific values for every pairwise interaction.

# lambdacov_form can be one of "lambdacov_none" or "lambdacov_global", where 
# "none" indicates that no effect of covariates on lambda is modelled,
# and "global" indicates that each covariate affects lambda

# alphacov_form can be one of "alphacov_none", "alphacov_global", or "alphacov_pairwise",
# where "none" indicates no effect of covariates on alpha values,
# "global" indicates a single effect term of each covariate on alpha values,
# and "pairwise" indicates a specific effect term of each covariate on each pairwise alpha

# 3 - FIXED TERMS
# if you have independent estimates of a subset of parameters, you may want
# to keep these as fixed, and only fit those you need to. This can be achieved
# through the argument "fixed_parameters", which is a list with four components (lambda, alpha, lambda_cov, alpha_cov),
# where each component should be either NULL if the parameter is to be fit, 
# or the fixed value of the parameter in question.

# 4 - Saving your model
# name the R file with the same name as your model

#TODO where to save the file, etc

pm_family_alpha_form_lambdacov_form_alphacov_form <- function(par,
                                                              fitness,
                                                              neigh_matrix,
                                                              covariates,
                                                              fixed_parameters){
  
  
  # retrieve parameters -----------------------------------------------------
  # parameters to fit are aLL in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,lambda_cov,alpha,alpha_cov,sigma}
  
  # YOU SHOULD NOT NEED TO CHANGE THIS SECTION
  pos <- 1
  
  # if a parameter is passed within the "par" vector,
  # it should be NULL in the "fixed_parameters" list
  if(is.null(fixed_parameters$lambda)){
    lambda <- par[pos]
    pos <- pos + 1
  }else{
    lambda <- fixed_parameters[["lambda"]]
  }
  
  if(is.null(fixed_parameters$lambda_cov)){
    lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
    pos <- pos + ncol(covariates)
  }else{
    lambda_cov <- fixed_parameters[["lambda_cov"]]
  }
  
  if(is.null(fixed_parameters$alpha)){
    alpha <- par[pos:(pos+ncol(neigh_matrix)-1)]
    pos <- pos + ncol(neigh_matrix)
  }else{
    alpha <- fixed_parameters[["alpha"]]
  }
  
  if(is.null(fixed_parameters$alpha_cov)){
    alpha.cov <- par[pos:(pos+(ncol(covariates)*ncol(neigh_matrix))-1)]
    pos <- pos + (ncol(covariates)*ncol(neigh_matrix))
  }else{
    alpha.cov <- fixed_parameters[["alpha.cov"]]
  }
  
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  
  
  # the model should return a "pred" value
  pred <- 0# a function of lambda, alpha, lambda_cov, alpha_cov 
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(log.fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
}