
#' Beverton-Holt model for projecting abundances,
#' with specific alpha values and global covariate effects on alpha and lambda
#'
#' @param lambda named numeric lambda value.
#' @param alpha_intra single numeric value.
#' @param alpha_inter numeric vector with interspecific alpha values.
#' @param lambda_cov numeric vector with effects of covariates over lambda.
#' @param alpha_cov named list of named numeric vectors 
#' with effects of each covariate over alpha values.
#' @param abundance named numeric vector of abundances in the previous timestep.
#' @param covariates matrix with observations in rows and covariates in named columns. 
#' Each cell is the value of a covariate in a given observation.
#'
#' @return numeric abundance projected one timestep
#' @export
BH_project_alpha_pairwise_lambdacov_global_alphacov_pairwise <- function(lambda,
                                                               alpha_intra,
                                                               alpha_inter,
                                                               lambda_cov,
                                                               alpha_cov,
                                                               abundance,
                                                               covariates){
  
  # put together intra and inter coefficients,
  # be sure names match
  
  spnames <- names(abundance)
  
  alpha <- c(alpha_intra,alpha_inter)
  alpha <- alpha[spnames]
  alpha_covs <- list()
  for(ia in 1:length(alpha_cov)){
    alpha_covs[[ia]] <- alpha_cov[[ia]][spnames]
  }
  
  numsp <- length(abundance)
  expected_abund <- NA_real_
  
  # model
  num = 1
  focal.cov.matrix <- as.matrix(covariates)
  for(v in 1:ncol(focal.cov.matrix)){
    num <- num + lambda_cov[v]*focal.cov.matrix[,v] 
  }
  cov_term_x <- list()
  for(v in 1:ncol(focal.cov.matrix)){
    cov_temp <- focal.cov.matrix[,v]
    for(z in 1:length(abundance)){
      #create  alpha_cov_i*cov_i vector
      cov_term_x[[z+(length(abundance)*(v-1))]] <- 
        # alpha_cov[z+(ncol(abund)*(v-1))] 
      alpha_cov[[v]][z] * cov_temp  
    }
  }
  cov_term <- list()
  for(z in 0:(length(abundance)-1)){
    cov_term_x_sum <- cov_term_x[[z+1]]
    if(ncol(focal.cov.matrix) > 1){
      for(v in 2:ncol(focal.cov.matrix)){
        cov_term_x_sum <- cov_term_x_sum + 
          cov_term_x[[v + length(abundance)]]
      } 
    }
    cov_term[[z+1]] <- cov_term_x_sum
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:length(abundance)){
    term <- term + (alpha[z] + cov_term[[z]]) * abundance[z]  
  }
  expected_abund <- (lambda * (num) / term) * abundance[names(lambda)]
  expected_abund
}

