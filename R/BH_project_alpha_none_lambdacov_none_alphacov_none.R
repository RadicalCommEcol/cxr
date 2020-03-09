
#' Beverton-Holt model for projecting abundances,
#' with no alpha and no covariate effects
#'
#' @param lambda numeric lambda value.
#' @param alpha_intra included for compatibility, not used in this model.
#' @param alpha_inter included for compatibility, not used in this model.
#' @param lambda_cov included for compatibility, not used in this model.
#' @param alpha_cov_intra included for compatibility, not used in this model.
#' @param alpha_cov_inter included for compatibility, not used in this model.
#' @param intra_abundance numeric abundance of the focal species in the previous timestep.
#' @param inter_abundances included for compatibility, not used in this model.
#' @param covariates included for compatibility, not used in this model.
#'
#' @return numeric abundance projected one timestep
#' @export
#'
#' @examples
BH_project_alpha_pairwise_lambdacov_none_alphacov_none <- function(lambda,
                                                                   alpha_intra,
                                                                   alpha_inter,
                                                                   lambda_cov,
                                                                   alpha_cov_intra,
                                                                   alpha_cov_inter,
                                                                   abundance_intra,
                                                                   abundance_inter,
                                                                   covariates){
  expected_abund <- lambda * abundance_intra
  expected_abund
}

