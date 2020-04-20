#' Title
#' Project abundances from population dynamics models
#' 
#' The function projects a number of steps of a time-discrete model, 
#' with model parameters taken from a `cxr_pm_multifit` object or as
#' function arguments.
#'
#' @param cxr_fit object of type `cxr_pm_multifit`. If this is not specified,
#' all parameters below are needed.
#' @param model_family acronym for model family. Included by default in `cxr` are
#' 'BH' (Beverton-Holt), 'RK' (Ricker), 'LW' (Law-Watkinson), 'LV' (Lotka-Volterra).
#' @param alpha_form character, either "none","global", or "pairwise".
#' @param lambda_cov_form character, either "none" or "global". 
#' @param alpha_cov_form character, either "none","global", or "pairwise".
#' @param lambda named vector with lambda values for all taxa to be projected.
#' @param alpha_matrix square matrix with taxa names in rows and columns.
#' @param lambda_cov optional named matrix with covariates in columns and taxa in rows,
#' representing the effect of each covariate on the lambda parameter of each taxa.
#' @param alpha_cov optional list. Each element of the named list represents the effects of
#' a covariate over alpha values. Thus, each list element contains either 
#' a single element (alpha_cov_form "global").
#' or a square matrix of same dimensions as `alpha_matrix` (alpha_cov_form "pairwise").
#' @param covariates matrix or dataframe with covariates in columns and timesteps in rows.
#' @param timesteps number of timesteps to project.
#' @param initial_abundances named vector of initial abundances for all taxa.
#'
#' @return named matrix with projected abundance values for each taxa at each timestep.
#' @export
#'
abundance_projection <- function(cxr_fit = NULL, 
                                 model_family = NULL,
                                 alpha_form = NULL,
                                 lambda_cov_form = NULL,
                                 alpha_cov_form = NULL,
                                 lambda = NULL,
                                 alpha_matrix = NULL,
                                 lambda_cov = NULL,
                                 alpha_cov = NULL,
                                 covariates = NULL, 
                                 timesteps = 2,
                                 initial_abundances = 0){
  
  # 0 - data ----------------------------------------------------------------
  
  user.par <- c(model_family,
                alpha_form,
                lambda_cov_form,
                alpha_cov_form,
                lambda,
                alpha_matrix,
                lambda_cov,
                alpha_cov)
  
  if(all(!is.null(user.par))){
    numsp <- length(lambda)
    spnames <- names(lambda)
    if(!alpha_form == "none"){
      alpha_rn <- rownames(alpha_matrix)
      alpha_cn <- colnames(alpha_matrix)
    }
    if(!lambda_cov_form == "none"){
      lambda_cov_rn <- rownames(lambda_cov)
      lambda_cov_cn <- colnames(lambda_cov)
    }
    if(!alpha_cov_form == "none"){
      alpha_cov_covn <- names(alpha_cov)
      alpha_cov_rn <- rownames(alpha_cov[[1]])
      alpha_cov_cn <- colnames(alpha_cov[[1]])
      covnames <- colnames(covariates)
    }
    
    abund_names <- names(initial_abundances)
    
    if(!lambda_cov_form == "none"){
      
      spnames.ok <- all(sapply(list(alpha_rn, 
                                    alpha_cn, 
                                    lambda_cov_rn, 
                                    alpha_cov_rn, 
                                    alpha_cov_cn,
                                    abund_names), FUN = identical, spnames))
      covnames.ok <- all(sapply(list(alpha_cov_covn, 
                                     lambda_cov_cn), FUN = identical, covnames))
    }else{
      if(!alpha_form == "none"){
      spnames.ok <- all(sapply(list(alpha_rn, 
                                    alpha_cn,
                                    abund_names), FUN = identical, spnames))
      }else{
        spnames.ok <- all(sapply(list(abund_names), FUN = identical, spnames))
      }
      covnames.ok <- TRUE
    }
    if(!spnames.ok | !covnames.ok){
      message("cxr abundance_projection ERROR: column and/or row names are inconsistent
              across function arguments. Please provide named vectors and matrices
              in all arguments, in order to avoid missmatches.")
      return(NULL)
    }
    
    if(timesteps < 2){
      message("cxr abundance_projection ERROR: number of timesteps cannot be < 2.")
      return(NULL)
    }
    
  }else if(!all(is.null(user.par))){
    message("cxr abundance_projection ERROR: not all parameters were specified.")
    return(NULL)
  }else if(!is.null(cxr_fit)){
    
    if(is.null(covariates) & !is.null(cxr_fit$covariates)){
      message("cxr abundance_projection ERROR: covariates not specified, but cxr_pm_multifit object
              was calculated with covariates.")
      return(NULL)
    }
    
    # obtain model family
    model_family <- substr(cxr_fit$model_name,1,2)
    
    # data
    lambda <- cxr_fit$lambda
    alpha_matrix <- cxr_fit$alpha_matrix
    lambda_cov <- cxr_fit$lambda_cov
    alpha_cov <- cxr_fit$alpha_cov
    
    numsp <- length(lambda)
    spnames <- names(lambda)
    if(is.null(alpha_matrix)){
      alpha_form <- "none"
    }else{
      una <- apply(alpha_matrix,MARGIN = 1,FUN = function(x)length(unique(x)))
      if(any(una > 1)){
        alpha_form <- "pairwise"
      }else{
        alpha_form <- "global"
      } 
    }
    
    lambda_cov_form <- ifelse(is.null(lambda_cov),"none","global")
    
    if(!is.null(alpha_cov)){
      for(i.cov in 1:length(alpha_cov)){
        alpha_cov[[i.cov]] <- alpha_cov[[i.cov]][,spnames]
      }
      un <- sapply(alpha_cov,FUN = function(i){apply(as.matrix(i),
                                                     MARGIN = 1,
                                                     FUN = function(x)length(unique(x)))})
      if(any(un != 1)){
        alpha_cov_form <- "pairwise"
      }else{
        alpha_cov_form <- "global"
      }
    }else{
      alpha_cov_form <- "none"
    }# alpha_cov_form
  }# if cxr_fit
  
  # 1 - get projection function ---------------------------------------------
  
  # get model
  # character string giving the name of the model
  model_name <- paste(model_family,
                      "_project_alpha_",alpha_form,
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

  expected_abund <- matrix(NA_real_,nrow = timesteps,ncol = numsp)
  expected_abund[1,] <- initial_abundances
  colnames(expected_abund) <- spnames
  
  for(it in 2:timesteps){
    for(i.sp in 1:length(spnames)){
      
      other.names <- spnames[which(spnames != spnames[i.sp])]
      
      if(alpha_form == "global"){
        alpha_intra <- NULL
        alpha_inter <- alpha_matrix[spnames[i.sp],]
      }else if(alpha_form == "pairwise"){
        alpha_intra <- alpha_matrix[spnames[i.sp],spnames[i.sp]]
        names(alpha_intra) <- spnames[i.sp]
        alpha_inter <- alpha_matrix[spnames[i.sp],other.names]
      }else{
        alpha_intra <- NULL
        alpha_inter <- NULL
      }
      
      if(lambda_cov_form == "none"){
        ab <- try(projection_model(lambda = lambda[spnames[i.sp]],
                                   alpha_intra = alpha_intra,
                                   alpha_inter = alpha_inter,
                                   lambda_cov = NULL,
                                   alpha_cov = NULL,
                                   abundance = expected_abund[it-1,],
                                   covariates = NULL))
        if(class(ab) != "try-error"){
          expected_abund[it,i.sp] <- ab
        }
      }else{
        alpha_cov_sp <- list()
        
        for(i.cov in 1:length(alpha_cov)){
          alpha_cov_sp[[i.cov]] <- alpha_cov[[i.cov]][spnames[i.sp],]
          if(alpha_cov_form == "global"){
            alpha_cov_sp[[i.cov]] <- unique(alpha_cov_sp[[i.cov]])
          }
        }
        names(alpha_cov_sp) <- names(alpha_cov)
        
        ab <- try(projection_model(lambda = lambda[spnames[i.sp]],
                                   alpha_intra = alpha_intra,
                                   alpha_inter = alpha_inter,
                                   lambda_cov = lambda_cov[spnames[i.sp],],
                                   alpha_cov = alpha_cov_sp,
                                   abundance = expected_abund[it-1,],
                                   covariates = covariates[it]))
        if(class(ab) != "try-error"){
          expected_abund[it,i.sp] <- ab
        }
      }# if-else covariates
    }# for each sp
  }# for each timestep
  
  # return ------------------------------------------------------------------
  
  expected_abund
  
}