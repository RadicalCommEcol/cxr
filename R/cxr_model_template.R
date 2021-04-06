
# 1 - MODEL FAMILY
# Model names are of the form shown below. "family" is the acronym of general family of your model,
# as e.g. BH for Beverton-Holt, LW for Law-Wilkinson, etc.
# these acronyms may also be useful to name non-linear versions of the models, e.g. "BH-NL"

# 2 - PARAMETER FORMS
# every model has, at least, a lambda parameter. Other potential parameters are 
# intraspecific (alpha_intra) and interspecific interactions (alpha_inter), 
# as well as the effects of covariates over lambda (lambda_cov)
# and over all alphas (alpha_cov).

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
# through the argument "fixed_parameters", which is a list with five components 
# (lambda, alpha_intra, alpha_inter, lambda_cov, alpha_cov),
# where each component should be either NULL if the parameter is to be fit, 
# or the fixed value of the parameter in question.

# 4 - Saving your model
# name the R file with the same name as your model, and for using it within cxr, 
# make it available in the global environment. This is as simple as sourcing the file.

# 5 - adding your model to cxr
# document your model, and file a pull request to the Github repository.

pm_family_alpha_form_lambdacov_form_alphacov_form <- function(par,
                                                              fitness,
                                                              neigh_intra_matrix = NULL,
                                                              neigh_inter_matrix,
                                                              covariates,
                                                              fixed_parameters){
  
  
  # retrieve parameters -----------------------------------------------------
  # parameters to fit are all in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,lambda_cov,alpha,alpha_cov,sigma}
  
  # comment or uncomment sections for the different parameters
  # depending on whether your model includes them
  # note that the section on alpha_inter includes two
  # possibilities, depending on whether a single alpha is 
  # fitted for all interactions (global) or each pairwise alpha is 
  # different (pairwise)
  # both are commented, you need to uncomment the appropriate one
  
  # likewise for the section on alpha_cov
  
  # --------------------------------------------------------------------------
  
  pos <- 1
  
  # if a parameter is passed within the "par" vector,
  # it should be NULL in the "fixed_parameters" list
  
  # lambda
  if(is.null(fixed_parameters$lambda)){
    lambda <- par[pos]
    pos <- pos + 1
  }else{
    lambda <- fixed_parameters[["lambda"]]
  }
  
  # lambda_cov
  if(is.null(fixed_parameters$lambda_cov)){
    lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
    pos <- pos + ncol(covariates)
  }else{
    lambda_cov <- fixed_parameters[["lambda_cov"]]
  }
  
  # alpha_intra
  if(!is.null(neigh_intra_matrix)){
    # intra
    if(is.null(fixed_parameters[["alpha_intra"]])){
      alpha_intra <- par[pos]
      pos <- pos + 1
    }else{
      alpha_intra <- fixed_parameters[["alpha_intra"]]
    }
  }else{
    alpha_intra <- NULL
  }
  
  # alpha_inter
  if(is.null(fixed_parameters[["alpha_inter"]])){
    # uncomment for alpha_global
    # alpha_inter <- par[pos]
    # pos <- pos + 1
    
    # uncomment for alpha_pairwise
    # alpha_inter <- par[pos:(pos+ncol(neigh_inter_matrix)-1)]
    # pos <- pos + ncol(neigh_inter_matrix)
  }else{
    alpha_inter <- fixed_parameters[["alpha_inter"]]
  }
  
  # alpha_cov
  if(is.null(fixed_parameters$alpha_cov)){
    # uncomment for alpha_cov_global
    # alpha_cov <- par[pos:(pos+ncol(covariates)-1)]
    # pos <- pos + ncol(covariates)
    
    # uncomment for alpha_cov_pairwise
    # alpha_cov <- par[pos:(pos+(ncol(covariates)*
    # (ncol(neigh_inter_matrix)+ncol(neigh_intra_matrix)))-1)]
    # pos <- pos + (ncol(covariates)*(ncol(neigh_inter_matrix)+ncol(neigh_intra_matrix)))
  }else{
    alpha_cov <- fixed_parameters[["alpha_cov"]]
  }
  
  # sigma - this is always necessary
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  # the model should return a "pred" value
  # a function of lambda, alpha_intra, alpha_inter, lambda_cov, alpha_cov 
  # and neigh_intra_matrix, neigh_inter_matrix, and covariates
  pred <- 0
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
}