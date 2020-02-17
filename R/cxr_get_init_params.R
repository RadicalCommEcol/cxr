
#' sort initial values from cxr_pm/er_fit in an appropriate way
#'
#' @inheritParams cxr_pm_fit
#' @param model_type either 'pm' or 'er'
#' @param neigh_intra name of neigh_intra column
#' @param neigh_inter names of neigh_inter_columns
#'
#' @return list with named initial values and fixed terms
#' @noRd
#'
cxr_get_init_params <- function(initial_values,
                                fixed_terms,
                                alpha_form,
                                lambda_cov_form,
                                alpha_cov_form,
                                model_type = "pm",
                                neigh_intra,
                                neigh_inter,
                                covariates){
  
  # here I store values for fixed parameters
  fixed_parameters <- list()
  
  # here I store initial values for fitted parameters
  init_lambda <- NULL
  init_alpha_intra <- NULL
  init_alpha_inter <- NULL
  init_lambda_cov <- NULL
  init_alpha_cov <- NULL
  
  # initialize params with appropriate length,
  # depending or not on whether they should be fixed
  
  if("lambda" %in% fixed_terms){
    fixed_parameters[["lambda"]] <- initial_values$lambda
  }else{
    init_lambda <- initial_values$lambda
    names(init_lambda) <- "lambda"
  }
  
  if(!is.null(initial_values$sigma)){
    init_sigma <- initial_values$sigma
  }else{
    init_sigma <- 0.1
  }
  names(init_sigma) <- "sigma"
  
  # return_init_length is an auxiliary function
  # in case initial values are not of the same length of the expected parameters
  # e.g. if we want to fit pairwise alphas but only provide a single initial value
  
  if(alpha_form != "none"){
    if("alpha" %in% fixed_terms){
      if(!is.null(neigh_intra)){
        fixed_parameters[["alpha_intra"]] <- initial_values$alpha_intra
      }
      
      fixed_parameters[["alpha_inter"]] <- cxr_return_init_length(alpha_form,
                                                                  initial_values$alpha_inter,
                                                                  neigh_inter,"pm")
    }else{
      
      if(!is.null(neigh_intra)){
        init_alpha_intra <- initial_values$alpha_intra
        names(init_alpha_intra) <- neigh_intra 
      }
      
      init_alpha_inter <- cxr_return_init_length(alpha_form,
                                                 initial_values$alpha_inter,
                                                 neigh_inter,"pm")
    }
  }
  
  
  if(lambda_cov_form != "none" & !is.null(covariates)){
    if("lambda_cov" %in% fixed_terms){
      fixed_parameters[["lambda_cov"]] <- cxr_return_init_length(lambda_cov_form,
                                                                 initial_values$lambda_cov,
                                                                 colnames(covariates),"pm")
    }else{
      init_lambda_cov <- cxr_return_init_length(lambda_cov_form,
                                                initial_values$lambda_cov,
                                                colnames(covariates),"pm")
      names(init_lambda_cov) <- paste("lambda_cov_",colnames(covariates),sep="")
    }
  }
  
  if(alpha_cov_form != "none" & !is.null(covariates)){
    if(alpha_cov_form == "global"){
      name.alpha.cov <- paste("alpha_cov_",colnames(covariates),sep="")
    }else{
      all_neigh <- c(neigh_intra,neigh_inter)
      name.alpha.cov <- paste("alpha_cov",rep(colnames(covariates),
                                              each = length(all_neigh)),
                              rep(all_neigh,ncol(covariates)),sep="_")
    }
    if("alpha_cov" %in% fixed_terms){
      fixed_parameters[["alpha_cov"]] <- cxr_return_init_length(alpha_cov_form,
                                                                initial_values$alpha_cov,
                                                                name.alpha.cov,"pm")
    }else{
      init_alpha_cov <- cxr_return_init_length(alpha_cov_form,
                                               initial_values$alpha_cov,
                                               name.alpha.cov,"pm")
    }
  }
  
  list(init_lambda = init_lambda,
       init_alpha_intra = init_alpha_intra,
       init_alpha_inter = init_alpha_inter,
       init_lambda_cov = init_lambda_cov,
       init_alpha_cov = init_alpha_cov,
       init_sigma = init_sigma,
       fixed_parameters = fixed_parameters)
  
} # function