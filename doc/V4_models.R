## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  pm_family_alpha_form_lambdacov_form_alphacov_form <- function(par,
#                                                                fitness,
#                                                                neigh_intra_matrix = NULL,
#                                                                neigh_inter_matrix,
#                                                                covariates,
#                                                                fixed_parameters)

## ----eval=FALSE----------------------------------------------------------
#  pm_UM_alpha_pairwise_lambdacov_none_alphacov_none <- function(par,
#                                                                fitness,
#                                                                neigh_intra_matrix = NULL,
#                                                                neigh_inter_matrix,
#                                                                covariates,
#                                                                fixed_parameters)

## ----eval=FALSE----------------------------------------------------------
#  
#  # retrieve parameters -----------------------------------------------------
#  # parameters to fit are all in the "par" vector,
#  # so we need to retrieve them one by one
#  # order is {lambda,lambda_cov,alpha,alpha_cov,sigma}
#  
#  # comment or uncomment sections for the different parameters
#  # depending on whether your model includes them
#  # note that the section on alpha_inter includes two
#  # possibilities, depending on whether a single alpha is
#  # fitted for all interactions (global) or each pairwise alpha is
#  # different (pairwise)
#  # both are commented, you need to uncomment the appropriate one
#  
#  # likewise for the section on alpha_cov
#  
#  # --------------------------------------------------------------------------
#  
#  pos <- 1
#  
#  # if a parameter is passed within the "par" vector,
#  # it should be NULL in the "fixed_parameters" list
#  
#  # lambda
#  if(is.null(fixed_parameters$lambda)){
#    lambda <- par[pos]
#    pos <- pos + 1
#  }else{
#    lambda <- fixed_parameters[["lambda"]]
#  }
#  
#  # lambda_cov
#  # if(is.null(fixed_parameters$lambda_cov)){
#  #   lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
#  #   pos <- pos + ncol(covariates)
#  # }else{
#  #   lambda_cov <- fixed_parameters[["lambda_cov"]]
#  # }
#  
#  # alpha_intra
#  if(!is.null(neigh_intra_matrix)){
#    # intra
#    if(is.null(fixed_parameters[["alpha_intra"]])){
#      alpha_intra <- par[pos]
#      pos <- pos + 1
#    }else{
#      alpha_intra <- fixed_parameters[["alpha_intra"]]
#    }
#  }else{
#    alpha_intra <- NULL
#  }
#  
#  # alpha_inter
#  if(is.null(fixed_parameters[["alpha_inter"]])){
#    # uncomment for alpha_global
#    # alpha_inter <- par[pos]
#    # pos <- pos + 1
#  
#    # uncomment for alpha_pairwise
#    alpha_inter <- par[pos:(pos+ncol(neigh_inter_matrix)-1)]
#    pos <- pos + ncol(neigh_inter_matrix)
#  }else{
#    alpha_inter <- fixed_parameters[["alpha_inter"]]
#  }
#  
#  # alpha_cov
#  # if(is.null(fixed_parameters$alpha_cov)){
#  #   # uncomment for alpha_cov_global
#  #   # alpha_cov <- par[pos:(pos+ncol(covariates)-1)]
#  #   # pos <- pos + ncol(covariates)
#  #
#  #   # uncomment for alpha_cov_pairwise
#  #   # alpha_cov <- par[pos:(pos+(ncol(covariates)*
#  #   # (ncol(neigh_inter_matrix)+ncol(neigh_intra_matrix)))-1)]
#  #   # pos <- pos + (ncol(covariates)*(ncol(neigh_inter_matrix)+ncol(neigh_intra_matrix)))
#  # }else{
#  #   alpha_cov <- fixed_parameters[["alpha_cov"]]
#  # }
#  
#  # sigma - this is always necessary
#  sigma <- par[length(par)]

## ----eval=FALSE----------------------------------------------------------
#  
#  # MODEL CODE HERE ---------------------------------------------------------
#  
#  # the model should return a "pred" value
#  # a function of lambda, alpha_intra, alpha_inter, lambda_cov, alpha_cov
#  # and neigh_intra_matrix, neigh_inter_matrix, and covariates
#  
#  # we do not differentiate alpha_intra from alpha_inter in this model
#  # so, put together alpha_intra and alpha_inter, and the observations
#  # with intraspecific ones at the beginning
#  if(!is.null(alpha_intra)){
#    alpha <- c(alpha_intra,alpha_inter)
#    all_neigh_matrix <- cbind(neigh_intra_matrix,neigh_inter_matrix)
#  }else{
#    alpha <- alpha_inter
#    all_neigh_matrix <- neigh_inter_matrix
#  }
#  
#  term = 1 #create the denominator term for the model
#  for(z in 1:ncol(all_neigh_matrix)){
#    term <- term + alpha[z]*all_neigh_matrix[,z]
#  }
#  pred <- lambda/ term
#  
#  # MODEL CODE ENDS HERE ----------------------------------------------------

## ----eval=FALSE----------------------------------------------------------
#  # the routine returns the sum of negative log-likelihoods of the data and model:
#  # DO NOT CHANGE THIS
#  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
#  return(sum(-1*llik))

## ----eval=FALSE----------------------------------------------------------
#  # load your model into the global environment
#  source("./pm_UM_alpha_pairwise_lambdacov_none_alphacov_none.R")
#  # fit your data
#  custom_fit <- cxr_pm_fit(data = custom_data, # assuming custom_data is already set...
#                           focal_column = my_focal, # assuming my_focal is already set...
#                           model_family = "UM",
#                           covariates = NULL, # as we have no covariate effects
#                           alpha_form = "pairwise",
#                           lambda_cov_form = "none",
#                           alpha_cov_form = "none")

## ----eval=FALSE----------------------------------------------------------
#  BH_demographic_ratio <- function(pair_lambdas){
#    (pair_lambdas[1]-1)/(pair_lambdas[2]-1)
#  }

## ----eval=FALSE----------------------------------------------------------
#  BH_competitive_ability <- function(lambda, pair_matrix){
#    if(all(pair_matrix >= 0)){
#      (lambda - 1)/sqrt(pair_matrix[1,1] * pair_matrix[1,2])
#    }else{
#      NA_real_
#    }
#  }

## ----eval=FALSE----------------------------------------------------------
#  pred <- lambda.part/ (1+ e.part*r.part )

## ----eval=FALSE----------------------------------------------------------
#  BH_species_fitness <- function(lambda, competitive_response){
#    (lambda-1)/competitive_response
#  }

## ----eval=FALSE----------------------------------------------------------
#  pm_family_alpha_form_lambdacov_form_alphacov_form <- function(par,
#                                                                fitness,
#                                                                neigh_intra_matrix = NULL,
#                                                                neigh_inter_matrix,
#                                                                covariates,
#                                                                fixed_parameters){
#  
#  
#    # retrieve parameters -----------------------------------------------------
#    # parameters to fit are all in the "par" vector,
#    # so we need to retrieve them one by one
#    # order is {lambda,lambda_cov,alpha_intra,alpha_inter,alpha_cov,sigma}
#  
#    # comment or uncomment sections for the different parameters
#    # depending on whether your model includes them
#    # note that the section on alpha_inter includes two
#    # possibilities, depending on whether a single alpha is
#    # fitted for all interactions (global) or each pairwise alpha is
#    # different (pairwise)
#    # both are commented, you need to uncomment the appropriate one
#  
#    # likewise for the section on alpha_cov
#  
#    # --------------------------------------------------------------------------
#  
#    pos <- 1
#  
#    # if a parameter is passed within the "par" vector,
#    # it should be NULL in the "fixed_parameters" list
#  
#    # lambda
#    if(is.null(fixed_parameters$lambda)){
#      lambda <- par[pos]
#      pos <- pos + 1
#    }else{
#      lambda <- fixed_parameters[["lambda"]]
#    }
#  
#    # lambda_cov
#    if(is.null(fixed_parameters$lambda_cov)){
#      lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
#      pos <- pos + ncol(covariates)
#    }else{
#      lambda_cov <- fixed_parameters[["lambda_cov"]]
#    }
#  
#    # alpha_intra
#    if(!is.null(neigh_intra_matrix)){
#      # intra
#      if(is.null(fixed_parameters[["alpha_intra"]])){
#        alpha_intra <- par[pos]
#        pos <- pos + 1
#      }else{
#        alpha_intra <- fixed_parameters[["alpha_intra"]]
#      }
#    }else{
#      alpha_intra <- NULL
#    }
#  
#    # alpha_inter
#    if(is.null(fixed_parameters[["alpha_inter"]])){
#      # uncomment for alpha_global
#      # alpha_inter <- par[pos]
#      # pos <- pos + 1
#  
#      # uncomment for alpha_pairwise
#      # alpha_inter <- par[pos:(pos+ncol(neigh_inter_matrix)-1)]
#      # pos <- pos + ncol(neigh_inter_matrix)
#    }else{
#      alpha_inter <- fixed_parameters[["alpha_inter"]]
#    }
#  
#    # alpha_cov
#    if(is.null(fixed_parameters$alpha_cov)){
#      # uncomment for alpha_cov_global
#      # alpha_cov <- par[pos:(pos+ncol(covariates)-1)]
#      # pos <- pos + ncol(covariates)
#  
#      # uncomment for alpha_cov_pairwise
#      # alpha_cov <- par[pos:(pos+(ncol(covariates)*
#      # (ncol(neigh_inter_matrix)+ncol(neigh_intra_matrix)))-1)]
#      # pos <- pos + (ncol(covariates)*(ncol(neigh_inter_matrix)+ncol(neigh_intra_matrix)))
#    }else{
#      alpha_cov <- fixed_parameters[["alpha_cov"]]
#    }
#  
#    # sigma - this is always necessary
#    sigma <- par[length(par)]
#  
#    # now, parameters have appropriate values (or NULL)
#    # next section is where your model is coded
#  
#    # MODEL CODE HERE ---------------------------------------------------------
#  
#    # the model should return a "pred" value
#    # a function of lambda, alpha_intra, alpha_inter, lambda_cov, alpha_cov
#    # and neigh_intra_matrix, neigh_inter_matrix, and covariates
#    pred <- 0
#  
#    # MODEL CODE ENDS HERE ----------------------------------------------------
#  
#    # the routine returns the sum of log-likelihoods of the data and model:
#    # DO NOT CHANGE THIS
#    llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
#    return(sum(-1*llik))
#  }

## ----eval=FALSE----------------------------------------------------------
#  er_family_lambdacov_form_effectcov_form_responsecov_form <- function(par,
#                                                                       fitness,
#                                                                       target,
#                                                                       density,
#                                                                       covariates,
#                                                                       fixed_parameters){
#  
#    num.sp <- nrow(target)
#  
#    # parameters to fit are all in the "par" vector,
#    # so we need to retrieve them one by one
#    # order is {lambda,lambda_cov,effect,effect_cov,response,response_cov,sigma}
#  
#    # comment or uncomment sections for the different parameters
#    # depending on whether your model includes them
#    # note that effect and response models must always include
#    # lambda, effect, and response, at least.
#  
#    pos <- 1
#  
#    # if a parameter is passed within the "par" vector,
#    # it should be NULL in the "fixed_parameters" list
#  
#    # lambda
#    if(is.null(fixed_parameters[["lambda"]])){
#      lambda <- par[pos:(pos + num.sp - 1)]
#      pos <- pos + num.sp
#    }else{
#      lambda <- fixed_parameters[["lambda"]]
#    }
#  
#    # lambda_cov
#    if(is.null(fixed_parameters$lambda_cov)){
#      # the covariate effects are more efficient in a matrix form
#      # with species in rows (hence byrow = T, because by default
#      # the vector is sorted first by covariates)
#      lambda_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],
#                           nrow = num.sp,
#                           byrow = TRUE)
#      pos <- pos + ncol(covariates)*num.sp
#    }else{
#      lambda_cov <- fixed_parameters[["lambda_cov"]]
#    }
#  
#    # effect
#    if(is.null(fixed_parameters[["effect"]])){
#      effect <- par[pos:(pos + num.sp - 1)]
#      pos <- pos + num.sp
#    }else{
#      effect <- fixed_parameters[["effect"]]
#    }
#  
#    # effect_cov
#    if(is.null(fixed_parameters$effect_cov)){
#      effect_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],
#                           nrow = num.sp,
#                           byrow = TRUE)
#      pos <- pos + ncol(covariates)*num.sp
#    }else{
#      effect_cov <- fixed_parameters[["effect_cov"]]
#    }
#  
#    # response
#    if(is.null(fixed_parameters[["response"]])){
#      response <- par[pos:(pos + num.sp - 1)]
#      pos <- pos + num.sp
#    }else{
#      response <- fixed_parameters[["response"]]
#    }
#  
#    # response_cov
#    if(is.null(fixed_parameters[["response_cov"]])){
#      response_cov <- matrix(par[pos:(pos+(ncol(covariates)*num.sp)-1)],
#                             nrow = num.sp,
#                             byrow = TRUE)
#      pos <- pos + ncol(covariates)*num.sp
#    }else{
#      response_cov <- fixed_parameters[["response_cov"]]
#    }
#  
#    sigma <- par[length(par)]
#  
#    # now, parameters have appropriate values (or NULL)
#    # next section is where your model is coded
#  
#    # MODEL CODE HERE ---------------------------------------------------------
#  
#    # the model should return a "pred" value
#    # a function of lambda, effect, response, lambda_cov, effect_cov, response_cov
#    pred <- 0
#  
#    # MODEL CODE ENDS HERE ----------------------------------------------------
#  
#    # the routine returns the sum of log-likelihoods of the data and model:
#    # DO NOT CHANGE THIS
#    llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
#    return(sum(-1*llik))
#  
#  }

