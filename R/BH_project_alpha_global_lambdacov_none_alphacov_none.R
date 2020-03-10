
#' Beverton-Holt model for projecting abundances,
#' with a global alpha and no covariate effects
#'
#' @param lambda numeric lambda value.
#' @param alpha_intra included for compatibility, not used in this model.
#' @param alpha_inter single numeric value.
#' @param lambda_cov included for compatibility, not used in this model.
#' @param alpha_cov included for compatibility, not used in this model.
#' @param abundance named numeric vector of abundances in the previous timestep.
#' @param covariates included for compatibility, not used in this model.
#'
#' @return numeric abundance projected one timestep
#' @export
BH_project_alpha_global_lambdacov_none_alphacov_none <- function(lambda,
                                                                 alpha_intra,
                                                                 alpha_inter,
                                                                 lambda_cov,
                                                                 alpha_cov,
                                                                 abundance,
                                                                 covariates){
  
  spnames <- names(abundance)
  alpha <- alpha_inter
  expected_abund <- NA_real_
  
  num <- lambda
  
  den <- 1
  for(i.sp in 1:length(abundance)){
    den <- den + alpha*abundance[i.sp]
  }# for each sp
  
  expected_abund <- (num/den) * abundance[names(lambda)]
  expected_abund
}

