
#' Ricker model for projecting abundances,
#' with specific alpha values and no covariate effects
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
RK_project_alpha_pairwise_lambdacov_none_alphacov_none <- function(lambda,
                                                                   alpha_intra,
                                                                   alpha_inter,
                                                                   lambda_cov,
                                                                   alpha_cov,
                                                                   abundance,
                                                                   covariates){
  
  spnames <- names(abundance)
  
  alpha <- c(alpha_intra,alpha_inter)
  alpha <- alpha[spnames]
  
  numsp <- length(abundance)
  expected_abund <- NA_real_
  
  num <- lambda
  term <- 0
  for(i.sp in 1:numsp){
    term <- term - abundance[i.sp]*alpha[i.sp]
  }# for each sp
  
  expected_abund <- (lambda * exp(term))* abundance[names(lambda)]
  expected_abund
}

