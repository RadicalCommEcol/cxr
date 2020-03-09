abundance_projection <- function(cxr_fit = NULL, covariates = NULL, timesteps = 0,initial_abundances = 0){
  
# 0 - data ----------------------------------------------------------------

  # beware columns have the same order in alpha_matrix, initial_abundances, etc.
  
numsp <- NULL
spnames <- NULL
lambda <- NULL
alpha_matrix <- NULL
lambda_cov <- NULL
alpha_cov_matrix <- NULL

# 1 - get projection function ---------------------------------------------

  # obtain model family
  mf <- substr(cxr_fit$model_name,4,5)
  # get model
  # character string giving the name of the model
  model_name <- paste(mf,
                      "_project_alpha",alpha_form,
                      "_lambdacov_",lambda_cov_form,
                      "_alphacov_",alpha_cov_form,sep="")
  
  # try to retrieve the function from its name
  # using function "get"
  projection_model <- try(get(model_name),silent = TRUE)
  if(class(projection_model) == "try-error"){
    message(paste("abundance_projection ERROR: model '",model_name,"' could not be retrieved. 
  Make sure it is defined and available in the cxr package or in the global environment.\n"
                  ,sep=""))
    return(NULL)
  }

# 2 - project all timesteps -----------------------------------------------

  # must be multispecies
  
  expected_abund <- matrix(nrow = timesteps,ncol = numsp)
  expected_abund[1,] <- initial_abundances
  colnames(expected_abund) <- spnames
  
  for(it in 2:timesteps){
    for(i.sp in 1:numsp){
      expected_abund[it,i.sp] <- projection_model(lambda = lambda[i.sp],
                                            alpha_intra = alpha_matrix[i.sp,i.sp],
                                            alpha_inter = alpha_matrix[i.sp,-i.sp],
                                            lambda_cov = lambda_cov[i.sp,],
                                            alpha_cov_intra = NULL,
                                            alpha_cov_inter = NULL,
                                            abundance_intra = expected_abund[it-1,i.sp],
                                            abundance_inter = expected_abund[it-1,-i.sp],
                                            covariates = covariates)
    }
  }
  
# return ------------------------------------------------------------------

  tt <- matrix(0, nrow = 3, ncol = 2)
  
}