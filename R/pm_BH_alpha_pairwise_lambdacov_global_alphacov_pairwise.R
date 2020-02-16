#' Beverton-Holt model with pairwise alphas, covariate effects on lambda, 
#' and pairwise covariate effects on alpha
#'
#' @param par 1d vector of initial parameters: lambda, lambda_cov, alpha, alpha_cov, and sigma
#' @param fitness 1d vector of fitness observations, in log scale
#' @param neigh_intra_matrix optional matrix of one column, number of intraspecific neighbours for each observation
#' @param neigh_inter_matrix matrix of arbitrary columns, number of interspecific neighbours for each observation
#' @param covariates optional matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation
#' @param fixed_parameters optional list specifying values of fixed parameters, 
#' with components "lambda","alpha_intra","alpha_inter","lambda_cov","alpha_cov".
#'
#' @return log-likelihood value
#' @export
pm_BH_alpha_pairwise_lambdacov_global_alphacov_pairwise <- function(par,
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
  
  if(is.null(fixed_parameters[["lambda_cov"]])){
    lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
    pos <- pos + ncol(covariates)
  }else{
    lambda_cov <- fixed_parameters[["lambda_cov"]]
  }
  
  
  if(!is.null(neigh_intra_matrix)){
    # intra
    if(is.null(fixed_parameters[["alpha_intra"]])){
      alpha_intra <- par[pos]
      pos <- pos + 1
    }else{
      alpha_intra <- fixed_parameters[["alpha_intra"]]
    }
  }else{
    alpha_intra <- NULL
  }
  
  # inter
  if(is.null(fixed_parameters[["alpha_inter"]])){
    alpha_inter <- par[pos:(pos+ncol(neigh_inter_matrix)-1)]
    pos <- pos + ncol(neigh_inter_matrix) -1
  }else{
    alpha_inter <- fixed_parameters[["alpha_inter"]]
  }
  
  if(is.null(fixed_parameters[["alpha_cov"]])){
    alpha_cov <- par[pos:(pos+(ncol(covariates)*ncol(neigh_matrix))-1)]
    pos <- pos + (ncol(covariates)*ncol(neigh_matrix))
  }else{
    alpha_cov <- fixed_parameters[["alpha_cov"]]
  }
  
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  # we do not differentiate alpha_intra from alpha_inter in this model
  # so, put together alpha_intra and alpha_inter, and the observations
  # with intraspecific ones at the beginning
  if(!is.null(alpha_intra)){
    alpha <- c(alpha_intra,alpha_inter)
    all_neigh_matrix <- cbind(neigh_intra_matrix,neigh_inter_matrix)
  }else{
    alpha <- alpha_inter
    all_neigh_matrix <- neigh_inter_matrix
  }
  
  # model
  num = 1
  focal.cov.matrix <- as.matrix(covariates)
  for(v in 1:ncol(covariates)){
    num <- num + lambda_cov[v]*focal.cov.matrix[,v] 
  }
  cov_term_x <- list()
  for(v in 1:ncol(covariates)){
    cov_temp <- focal.cov.matrix[,v]
    for(z in 1:ncol(all_neigh_matrix)){
      #create  alpha_cov_i*cov_i vector
      cov_term_x[[z+(ncol(all_neigh_matrix)*(v-1))]] <- 
        alpha_cov[z+(ncol(all_neigh_matrix)*(v-1))] * cov_temp  
    }
  }
  cov_term <- list()
  for(z in 0:(ncol(all_neigh_matrix)-1)){
    cov_term_x_sum <- cov_term_x[[z+1]]
    if(ncol(covariates) > 1){
      for(v in 2:ncol(covariates)){
        cov_term_x_sum <- cov_term_x_sum + 
          cov_term_x[[v + ncol(all_neigh_matrix)]]
      } 
    }
    cov_term[[z+1]] <- cov_term_x_sum
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:ncol(all_neigh_matrix)){
    term <- term + (alpha[z] + cov_term[[z]]) * all_neigh_matrix[,z]  
  }
  pred <- lambda * (num) / term 
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
}