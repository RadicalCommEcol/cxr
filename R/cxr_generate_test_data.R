#' Generate simulated interaction data
#'
#' Model fitness responses to neighbours and covariates
#' using a Beverton-Holt functional form. This function
#' is fairly restricted and under development, but can be used
#' to generate simple test data to run the main functions of \code{cxr}.
#'
#' @param focal_sp number of focal species, defaults to 1.
#' @param neigh_sp number of neighbour species, defaults to 1.
#' @param covariates number of covariates, defaults to 0.
#' @param observations number of observations, defaults to 10.
#' @param alpha_form what form does the alpha parameter take? one of "none" (no alpha in the model), 
#' "global" (a single alpha for all pairwise interactions), or "pairwise" (one alpha value for every interaction).
#' @param lambda_cov_form form of the covariate effects on lambda. Either "none" (no covariate effects) or "global" (one estimate per covariate).
#' @param alpha_cov_form form of the covariate effects on alpha. One of "none" (no covariate effects), "global" (one estimate per covariate on every alpha),
#' or "pairwise" (one estimate per covariate and pairwise alpha).
#' @param focal_lambda optional 1d vector with lambdas of the focal sp.
#' @param min_lambda if no focal_lambda is provided, lambdas are taken from a uniform distribution
#' with min_lambda and max_lambda as minimum and maximum values.
#' @param max_lambda if no focal_lambda is provided, lambdas are taken from a uniform distribution
#' with min_lambda and max_lambda as minimum and maximum values.
#' @param alpha optional interaction matrix, neigh_sp x neigh_sp
#' @param min_alpha if no focal_alpha is provided, alphas are taken from a uniform distribution
#' with min_alpha and max_alpha as minimum and maximum values.
#' @param max_alpha if no focal_alpha is provided, alphas are taken from a uniform distribution
#' with min_alpha and max_alpha as minimum and maximum values.
#' @param lambda_cov optional matrix of neigh_sp x covariates
#' giving the effect of each covariate over the fecundity (lambda) of each species.
#' @param min_lambda_cov if no focal_lambda_cov is provided, lambda_covs are taken from a uniform distribution
#' with min_lambda_cov and max_lambda_cov as minimum and maximum values.
#' @param max_lambda_cov if no focal_lambda_cov is provided, lambda_covs are taken from a uniform distribution
#' with min_lambda and max_lambda as minimum and maximum values.
#' @param alpha_cov ----------Under development-------------
#' @param min_alpha_cov if no focal_alpha_cov is provided, alpha_covs are taken from a uniform distribution
#' with min_alpha_cov and max_alpha_cov as minimum and maximum values.
#' @param max_alpha_cov if no focal_alpha_cov is provided, alpha_covs are taken from a uniform distribution
#' with min_alpha and max_alpha as minimum and maximum values.
#' @return list with two components: 'observations' is a list with as many components as focal species. 
#' Each component of 'observations' is a dataframe with stochastic number of neighbours and associated fitness.
#' The second component, 'covariates', is again a list with one component per focal species. 
#' Each component of 'covariates' is a dataframe with the values of each covariate for each associated observation.
#' @import stats
#' @export
#' 
#' @example 
#' example_obs <- cxr_generate_test_data(focal_sp = 2,
#'                                       neigh_sp = 2,
#'                                       alpha_form = "pairwise",
#'                                       lambda_cov_form = "global",
#'                                       alpha_cov_form = "global",
#'                                       covariates = 1)
#' 
cxr_generate_test_data <- function(focal_sp = 1,
                             neigh_sp = 1,
                             covariates = 0,
                             observations = 10, 
                             # fitness.model = 1,
                             # model_family = "BH",
                             alpha_form = c("pairwise","none","global"),
                             lambda_cov_form = c("none","global"),
                             alpha_cov_form = c("none","global","pairwise"),
                             
                             focal_lambda = NULL,
                             min_lambda = 0,
                             max_lambda = 10,
                             alpha = NULL,
                             min_alpha = 0,
                             max_alpha = 1,
                             alpha_cov = NULL,
                             min_alpha_cov = -1,
                             max_alpha_cov = 1,
                             lambda_cov = NULL,
                             min_lambda_cov = -1,
                             max_lambda_cov = 1,
                             min_cov = 0,
                             max_cov = 1){
  
  
  
  # check arguments ---------------------------------------------------------
  
  alpha_form <- match.arg(alpha_form)
  lambda_cov_form <- match.arg(lambda_cov_form)
  alpha_cov_form <- match.arg(alpha_cov_form)
  
  # retrieve parameters -----------------------------------------------------
  
  if(is.null(focal_lambda)){
    local_lambda <- runif(focal_sp,min_lambda,max_lambda)
  }else{
    local_lambda <- focal_lambda
  }
  
  if(is.null(alpha)){
    if(alpha_form == "global"){
      local_alpha <- runif(focal_sp,min_alpha,max_alpha)
    }else if(alpha_form == "pairwise"){
      local_alpha <- matrix(runif(focal_sp*neigh_sp,
                                  min_alpha,
                                  max_alpha),
                            nrow = focal_sp,
                            ncol = neigh_sp)
    }
    
  }else{
    local_alpha <- alpha
  }
  
  if(covariates > 0){
    
    if(is.null(lambda_cov)){
      if(lambda_cov_form == "global"){
        local_lambda_cov <- matrix(runif(focal_sp*covariates,min_lambda_cov,max_lambda_cov),
                                   nrow = focal_sp,
                                   ncol = covariates) 
      }
    }else{
      local_lambda_cov <- lambda_cov
    }
    
    if(is.null(alpha_cov)){
      if(alpha_cov_form == "global"){
        local_alpha_cov <- matrix(runif(focal_sp*covariates,min_alpha_cov,max_alpha_cov),
                                  nrow = focal_sp,
                                  ncol = covariates)  
      }else if(alpha_cov_form == "pairwise"){
        # TODO
        # local_alpha_cov <- matrix()
      }
    }else{
      local_alpha_cov <- alpha_cov
    }
  }
  
  # generate neighbour data -----------------------------------------------------------
  
    local_obs <- list()
    for(i.focal in 1:focal_sp){
      local_obs[[i.focal]] <- as.data.frame(matrix(data = round(rpois(neigh_sp*observations,3)),nrow = observations))
      colnames(local_obs[[i.focal]]) <- paste("neigh_",as.character(1:neigh_sp),sep="")
    }
    names(local_obs) <- as.character(paste("focal_",1:focal_sp,sep=""))
    
    if(covariates > 0){
      cov_obs <- list()
      for(i.focal in 1:focal_sp){
        cov_obs[[i.focal]] <- matrix(data = runif(covariates*observations,min_cov,max_cov),nrow = observations)
        colnames(cov_obs[[i.focal]]) <- paste("cov_",1:covariates,sep="")
      }
      names(cov_obs) <- as.character(paste("focal_",1:focal_sp,sep=""))
    }else{
      cov_obs <- NULL
    }
    
  # generate fitness --------------------------------------------------------
  # TODO if important, this can be polished to allow more flexibility
  # by providing separate fitness-generating functions
  # of different model families, etc
  
    for(i.sp in 1:length(local_obs)){
      
      # what type of alpha,lambda_cov,alpha_cov
      if(alpha_form == "none"){
        local_obs[[i.sp]]$fitness <- rep(local_lambda,observations)
      }else if(alpha_form == "global"){
        for(i.obs in 1:nrow(local_obs[[i.sp]])){
          den <- 1+(local_alpha[i.sp,]*sum(local_obs[[i.sp]][i.obs,1:neigh_sp]))
          local_obs[[i.sp]]$fitness[i.obs] <- local_lambda[i.sp]/den
        }
      }else if(alpha_form == "pairwise"){
        if(lambda_cov_form == "none"){
          for(i.obs in 1:nrow(local_obs[[i.sp]])){
            den <- 1+(sum(local_alpha[i.sp,]*local_obs[[i.sp]][i.obs,1:neigh_sp]))
            local_obs[[i.sp]]$fitness[i.obs] <- local_lambda[i.sp]/den
          }
        }else if(lambda_cov_form == "global" & 
                 alpha_cov_form == "global" & !is.null(covariates)){
          for(i.obs in 1:nrow(local_obs[[i.sp]])){
            num <- local_lambda[i.sp] + sum(cov_obs[[i.sp]][i.obs,]*local_lambda_cov[i.sp,])
            cov_term <- sum(cov_obs[[i.sp]][i.obs,]*local_alpha_cov[i.sp,])
            den <- 1+(sum(local_alpha[i.sp,]*local_obs[[i.sp]][i.obs,1:neigh_sp] + cov_term))
            local_obs[[i.sp]]$fitness[i.obs] <- num/den
          }
        }else if(lambda_cov_form == "global" &
                 alpha_cov_form == "pairwise" & !is.null(covariates)){
          for(i.obs in 1:nrow(local_obs[[i.sp]])){
            num <- local_lambda[i.sp] + sum(cov_obs[[i.sp]][i.obs,]*local_lambda_cov[i.sp,])
            # TODO decide on the form of alpha_cov
            # cov_term <- sum(cov_obs[i.obs,]*local_alpha_cov[])
            cov_term <- 0
            den <- 1+(sum(local_alpha[i.sp]*local_obs[[i.sp]][i.obs,1:neigh_sp] + cov_term))
            local_obs[[i.sp]]$fitness[i.obs] <- num/den
          }
        }
      }else{
        stop("cxr_generate_test_data ERROR: the combination of 
    alpha_form, lambda_cov_form, alpha_cov_form provided 
         is not currently supported")
      }
      
    }# for i.sp
 
  list(observations = local_obs,covariates = cov_obs)
   
}
