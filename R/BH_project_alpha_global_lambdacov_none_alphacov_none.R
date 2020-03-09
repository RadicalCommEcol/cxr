
#' Beverton-Holt model for projecting abundances,
#' with a global alpha and no covariate effects
#'
#' @param lambda numeric lambda value.
#' @param alpha_intra included for compatibility, not used in this model.
#' @param alpha_inter single numeric value.
#' @param lambda_cov included for compatibility, not used in this model.
#' @param alpha_cov_intra included for compatibility, not used in this model.
#' @param alpha_cov_inter included for compatibility, not used in this model.
#' @param intra_abundance numeric abundance of the focal species in the previous timestep.
#' @param inter_abundances 1d vector of neighbour abundances in the previous timestep.
#' @param covariates included for compatibility, not used in this model.
#'
#' @return numeric abundance projected one timestep
#' @export
#'
#' @examples
BH_project_alpha_global_lambdacov_none_alphacov_none <- function(lambda,
                                                                   alpha_intra,
                                                                   alpha_inter,
                                                                   lambda_cov,
                                                                   alpha_cov,
                                                                   abundance_intra,
                                                                   abundance_inter,
                                                                   covariates){
  
  expected_abund <- NA_real_
  
  num <- lambda
  
  den <- 1 + alpha*abundance_intra
  for(i.sp in 1:length(abundance_inter)){
    den <- den + alpha*abundance_inter[i.sp]
  }# for each sp
  
  expected_abund <- (num/den) * abundance_intra
  expected_abund
}

