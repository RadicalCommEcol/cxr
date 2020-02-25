#' Ricker model with no alphas and no covariate effects
#'
#' This model, in all families, is simply given by lambda.
#'
#' @param par 1d vector of initial parameters: lambda and sigma
#' @param fitness 1d vector of fitness observations, in log scale
#' @param neigh_intra_matrix included for compatibility, not used in this model.
#' @param neigh_inter_matrix included for compatibility, not used in this model.
#' @param covariates included for compatibility, not used in this model
#' @param fixed_parameters included for compatibility, not used in this model
#'
#' @return log-likelihood value
#' @export
pm_RK_alpha_none_lambdacov_none_alphacov_none <- function(par,
                                                          fitness,
                                                          neigh_intra_matrix = NULL,
                                                          neigh_inter_matrix,
                                                          covariates,
                                                          fixed_parameters){
  
  
  # retrieve parameters -----------------------------------------------------
  # parameters to fit are all in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,lambda_cov,alpha,alpha_cov,sigma}
  
  # comment or uncomment sections for the different parameters
  # depending on whether your model includes them
  pos <- 1
  
  # if a parameter is passed within the "par" vector,
  # it should be NULL in the "fixed_parameters" list
  if(is.null(fixed_parameters[["lambda"]])){
    lambda <- par[pos]
    pos <- pos + 1
  }else{
    lambda <- fixed_parameters[["lambda"]]
  }
  
  # if(is.null(fixed_parameters$lambda_cov)){
  #   lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
  #   pos <- pos + ncol(covariates)
  # }else{
  #   lambda_cov <- fixed_parameters[["lambda_cov"]]
  # }
  
  # if(is.null(fixed_parameters$alpha)){
  #   alpha <- par[pos:(pos+ncol(neigh_matrix)-1)]
  #   pos <- pos + ncol(neigh_matrix)
  # }else{
  #   alpha <- fixed_parameters[["alpha"]]
  # }
  
  # if(is.null(fixed_parameters$alpha_cov)){
  #   alpha.cov <- par[pos:(pos+(ncol(covariates)*ncol(neigh_matrix))-1)]
  #   pos <- pos + (ncol(covariates)*ncol(neigh_matrix))
  # }else{
  #   alpha.cov <- fixed_parameters[["alpha.cov"]]
  # }
  
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  pred <- rep(lambda, times=length(fitness)) 
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
}