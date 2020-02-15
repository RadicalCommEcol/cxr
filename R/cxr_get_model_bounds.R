
cxr_get_model_bounds <- function(lower_bounds,upper_bounds,fixed_terms){
  
  lower_lambda <- NULL
  upper_lambda <- NULL
  lower_sigma <- NULL
  upper_sigma <- NULL
  lower_alpha_intra <- NULL
  upper_alpha_intra <- NULL
  lower_alpha_inter <- NULL
  upper_alpha_inter <- NULL
  lower_lambda_cov <- NULL
  upper_lambda_cov <- NULL
  lower_alpha_cov <- NULL
  upper_alpha_cov <- NULL
  
  if(!is.null(lower_bounds$lambda) & 
     !is.null(upper_bounds$lambda) &
     !"lambda" %in% fixed_terms){
    lower_lambda <- lower_bounds$lambda
    upper_lambda <- upper_bounds$lambda
  }
  
  # sigma bounds can be hidden from the user
  # only set if there are bounds for other params
  if(!is.null(lower_bounds) & 
     !is.null(upper_bounds)){
    lower_sigma <- 1e-5
    upper_sigma <- 1e5
  }
  
  if(!is.null(lower_bounds$alpha_intra) & 
     !is.null(upper_bounds$alpha_intra) &
     !"alpha_intra" %in% fixed_terms){
    lower_alpha_intra <- lower_bounds$alpha_intra
    upper_alpha_intra <- upper_bounds$alpha_intra
  }
  
  if(!is.null(lower_bounds$alpha_inter) & 
     !is.null(upper_bounds$alpha_inter) &
     !"alpha_inter" %in% fixed_terms){
    lower_alpha_inter <- lower_bounds$alpha_inter
    upper_alpha_inter <- upper_bounds$alpha_inter
  }

  
  if(!is.null(lower_bounds$lambda_cov) &
     !is.null(upper_bounds$lambda_cov) &
     !"lambda_cov" %in% fixed_terms){
    lower_lambda_cov <- lower_bounds$lambda_cov
    upper_lambda_cov <- upper_bounds$lambda_cov
  }
  
  if(!is.null(lower_bounds$alpha_cov) &
     !is.null(upper_bounds$alpha_cov) &
     !"alpha_cov" %in% fixed_terms){
    lower_alpha_cov <- lower_bounds$alpha_cov
    upper_alpha_cov <- upper_bounds$alpha_cov
  }
  
  list(lower_lambda = lower_lambda,
       lower_alpha_intra = lower_alpha_intra,
       lower_alpha_inter = lower_alpha_inter,
       lower_lambda_cov = lower_lambda_cov,
       lower_alpha_cov = lower_alpha_cov,
       upper_lambda = upper_lambda,
       upper_alpha_intra = upper_alpha_intra,
       upper_alpha_inter = upper_alpha_inter,
       upper_lambda_cov = upper_lambda_cov,
       upper_alpha_cov = upper_alpha_cov,
       lower_sigma = lower_sigma,
       upper_sigma = upper_sigma)
  
}
