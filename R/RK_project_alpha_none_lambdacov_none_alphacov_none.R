
#' Model for projecting abundances,
#' with no alpha and no covariate effects
#'
#' @param lambda numeric lambda value.
#' @param alpha_intra included for compatibility, not used in this model.
#' @param alpha_inter included for compatibility, not used in this model.
#' @param lambda_cov included for compatibility, not used in this model.
#' @param alpha_cov included for compatibility, not used in this model.
#' @param abundance named numeric vector of abundances in the previous timestep.
#' @param covariates included for compatibility, not used in this model.
#'
#' @return numeric abundance projected one timestep
#' @export
RK_project_alpha_none_lambdacov_none_alphacov_none <- function(lambda,
                                                                   alpha_intra,
                                                                   alpha_inter,
                                                                   lambda_cov,
                                                                   alpha_cov,
                                                                   abundance,
                                                                   covariates){
  expected_abund <- lambda * abundance[names(lambda)]
  expected_abund
}

