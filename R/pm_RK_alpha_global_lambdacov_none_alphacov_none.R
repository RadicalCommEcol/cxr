#' Ricker model with a global alpha and no covariate effects
#'
#' @param par 1d vector of initial parameters: lambda, alpha, and sigma.
#' @param fitness 1d vector of fitness observations, in log scale.
#' @param neigh_intra_matrix included for compatibility, not used in this model.
#' @param neigh_inter_matrix matrix of arbitrary columns, number of neighbours for each observation.
#' As in this model there is a single alpha argument, do not distinguish neighbour identity
#' @param covariates included for compatibility, not used in this model.
#' @param fixed_parameters optional list specifying values of fixed parameters, 
#' with components "lambda","alpha_inter".
#'
#' @return log-likelihood value
#' @export
pm_RK_alpha_global_lambdacov_none_alphacov_none <- function(par,
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
  
  
  # if(!is.null(neigh_intra_matrix)){
  #   # intra
  #   if(is.null(fixed_parameters[["alpha_intra"]])){
  #     alpha_intra <- par[pos]
  #     pos <- pos + 1
  #   }else{
  #     alpha_intra <- fixed_parameters[["alpha_intra"]]
  #   }
  # }else{
  #   alpha_intra <- NULL
  # }
  
  # inter
  if(is.null(fixed_parameters[["alpha_inter"]])){
    alpha_inter <- par[pos]
    pos <- pos + 1
  }else{
    alpha_inter <- fixed_parameters[["alpha_inter"]]
  }
  
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
  
  background <- rowSums(neigh_inter_matrix)
  pred <- lambda * exp(-alpha_inter*(background))  
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
}