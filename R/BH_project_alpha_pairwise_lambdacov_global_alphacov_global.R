
#' Beverton-Holt model for projecting abundances,
#' with specific alpha values and global covariate effects on alpha and lambda
#'
#' @param lambda numeric lambda value.
#' @param alpha_intra single numeric value.
#' @param alpha_inter numeric vector with interspecific alpha values.
#' @param lambda_cov numeric vector with effects of covariates over lambda.
#' @param alpha_cov_intra included for compatibility, not used in this model.
#' @param alpha_cov_inter list of numeric vectors with effects of each covariate over every alpha.
#' @param intra_abundance numeric abundance of the focal species in the previous timestep.
#' @param inter_abundances 1d vector of neighbour abundances in the previous timestep.
#' @param covariates matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation.
#'
#' @return numeric abundance projected one timestep
#' @export
#'
#' @examples
BH_project_alpha_pairwise_lambdacov_global_alphacov_global <- function(lambda,
                                                               alpha_intra,
                                                               alpha_inter,
                                                               lambda_cov,
                                                               alpha_cov_intra,
                                                               alpha_cov_inter,
                                                               abundance_intra,
                                                               abundance_inter,
                                                               covariates){
  
  alpha <- c(alpha_intra,alpha_inter)
  abund <- c(abundance_intra,abundance_inter)
  alpha_cov <- alpha_cov_inter
  numsp <- length(abund)
  expected_abund <- NA_real_
  
  # numerator
  num = 1
  focal.cov.matrix <- as.matrix(covariates)
  for(z in 1:ncol(focal.cov.matrix)){
    num <- num + lambda_cov[z]*focal.cov.matrix[,z]
  }
  # denominator
  cov_term <- 0 
  for(v in 1:ncol(focal.cov.matrix)){
    cov_term <- cov_term + alpha_cov[[v]] * focal.cov.matrix[,v]
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:length(abund)){
    term <- term + (alpha[z] + cov_term) * abund[z] 
  }
  
  expected_abund <- (lambda * (num / term)) * abundance_intra
  expected_abund
}

