
#' Beverton-Holt model for projecting abundances,
#' with specific alpha values and global covariate effects on alpha and lambda
#'
#' @param lambda numeric lambda value.
#' @param alpha_intra single numeric value.
#' @param alpha_inter numeric vector with interspecific alpha values.
#' @param lambda_cov numeric vector with effects of covariates over lambda.
#' @param alpha_cov numeric vector with effects of covariates over alpha. ---------------expand--
#' @param intra_abundance numeric abundance of the focal species in the previous timestep.
#' @param inter_abundances numeric vector of neighbour abundances in the previous timestep.
#' @param covariates matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation.
#'
#' @return numeric abundance projected one timestep
#' @export
#'
#' @examples
BH_project_alpha_pairwise_lambdacov_global_alphacov_pairwise <- function(lambda,
                                                               alpha_intra,
                                                               alpha_inter,
                                                               lambda_cov,
                                                               alpha_cov,
                                                               abundance_intra,
                                                               abundance_inter,
                                                               covariates){
  
  alpha <- c(alpha_intra,alpha_inter)
  abund <- c(abundance_intra,abundance_inter)
  numsp <- length(abund)
  expected_abund <- NA_real_
  
  # model
  num = 1
  focal.cov.matrix <- as.matrix(covariates)
  for(v in 1:ncol(covariates)){
    num <- num + lambda_cov[v]*focal.cov.matrix[,v] 
  }
  cov_term_x <- list()
  for(v in 1:ncol(covariates)){
    cov_temp <- focal.cov.matrix[,v]
    for(z in 1:ncol(abund)){
      #create  alpha_cov_i*cov_i vector
      cov_term_x[[z+(ncol(abund)*(v-1))]] <- 
        alpha_cov[z+(ncol(abund)*(v-1))] * cov_temp  
    }
  }
  cov_term <- list()
  for(z in 0:(ncol(abund)-1)){
    cov_term_x_sum <- cov_term_x[[z+1]]
    if(ncol(covariates) > 1){
      for(v in 2:ncol(covariates)){
        cov_term_x_sum <- cov_term_x_sum + 
          cov_term_x[[v + ncol(abund)]]
      } 
    }
    cov_term[[z+1]] <- cov_term_x_sum
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:ncol(abund)){
    term <- term + (alpha[z] + cov_term[[z]]) * abund[z]  
  }
  expected_abund <- (lambda * (num) / term) * abundance_intra 
  expected_abund
}

