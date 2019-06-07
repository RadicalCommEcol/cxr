#' Title ML estimation of competitive effects and responses
#' 
#' Calculates the log-likelihood of a Beverton-Holt model parameterized with given values
#' with respect to a fitness metric. The function for calculating fecundity given 
#' effect and response values is taken from Godoy et al. (2014). 
#' Note that, as e is not pair-specific, all species parameters are fit in the same function.
#' In this version, lambda values are fixed.
#'
#' @param init.par 1d vector of initial parameters: r values, e values, and single sigma term.
#' If covariates are given, this vector must include also lambda.cov, response.cov, and effect.cov terms, 
#' after the e values and before the sigma term.
#' @param lambda 1d vector of lambda values
#' @param target_all matrix giving which species is calculated with which values. See ER_optimize
#' @param density_all matrix giving the densities of each species at each observation. See ER_optimize
#' @param log.fitness log of the fitness metric
#' @param covariates if present, it is a dataframe/matrix with as many rows as observationes, and one column for each covariate.
#' Values are covariate values for every observation.
#' @return single numeric value giving the sum of negative log-likelihoods
#' @export
EffectResponse <- function(init.par, lambda, target_all, density_all, log.fitness, covariates = NULL){
  
  num.sp <- nrow(target_all) # same as length(lambda)
  
  r.vector <- init.par[1:num.sp]
  
  if(is.null(covariates)){
    e.vector <- init.par[(num.sp+1):(length(init.par)-1)]
  }else{
    e.vector <- init.par[(num.sp+1):(num.sp + num.sp)]
    lambda.cov <- matrix(init.par[(num.sp + num.sp + 1):(num.sp + num.sp + (num.sp*ncol(covariates)))],
                         nrow = num.sp,ncol = ncol(covariates))
    response.cov <- matrix(init.par[(num.sp + num.sp + (num.sp*ncol(covariates)) + 1):(num.sp + 
                                                                                      num.sp + 
                                                                                      (num.sp*ncol(covariates)) + 
                                                                                      (num.sp*ncol(covariates)))],
                        nrow = num.sp,ncol = ncol(covariates))
    effect.cov <- matrix(init.par[(num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + 1):(num.sp + 
                                                                                         num.sp + 
                                                                                         (num.sp*ncol(covariates)) +
                                                                                         (num.sp*ncol(covariates)) +
                                                                                         (num.sp*ncol(covariates)))],
                           nrow = num.sp,ncol = ncol(covariates))
  }
  
  sigma <- init.par[length(init.par)]
  
  if(is.null(covariates)){
    lambda.part <- colSums(lambda*target_all)
    r.part <- colSums(r.vector*target_all)
    e.part <- colSums(e.vector*density_all)
  }else{
    lambda.cov_all <- matrix(0,nrow = num.sp,ncol = length(log.fitness))
    response.cov_all <- matrix(0,nrow = num.sp,ncol = length(log.fitness))
    effect.cov_all <- matrix(0,nrow = num.sp,ncol = length(log.fitness))
    for(i.sp in 1:num.sp){
      for(i.obs in 1:length(log.fitness)){
        lambda.cov_all[i.sp,i.obs] <- sum(lambda.cov[i.sp,]*covariates[i.obs,])
        response.cov_all[i.sp,i.obs] <- sum(response.cov[i.sp,]*covariates[i.obs,])
        effect.cov_all[i.sp,i.obs] <- sum(effect.cov[i.sp,]*covariates[i.obs,])
        
      }# for i.cov
    }# for i.sp
    
    lambda.part <- colSums(lambda*(1+lambda.cov_all)*target_all)
    r.part <- colSums(r.vector*(1+response.cov_all)*target_all)
    e.part <- colSums(e.vector*(1+effect.cov_all)*density_all)
  }# if-else null covariates
  
  pred <- lambda.part/ (1+ e.part*r.part )
  
  # likelihood as before:
  llik<-dnorm(log.fitness,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

#' Title ML estimation of competitive effects and responses
#' 
#' Calculates the log-likelihood of a Beverton-Holt model parameterized with given values
#' with respect to a fitness metric. The function for calculating fecundity given 
#' effect and response values is taken from Godoy et al. (2014). 
#' Note that, as e is not pair-specific, all species parameters are fit in the same function.
#' In this version, lambda values are also fit.
#'
#' @param init.par 1d vector of initial parameters: lambda values, r values, e values, and single sigma term.
#' If covariates are given, this vector must include also lambda.cov, response.cov, and effect.cov terms, 
#' after the e values and before the sigma term.
#' @param target_all matrix giving which species is calculated with which values. See ER_optimize
#' @param density_all matrix giving the densities of each species at each observation. See ER_optimize
#' @param log.fitness log of the fitness metric
#' @param covariates if present, it is a dataframe/matrix with as many rows as observationes, and one column for each covariate.
#' Values are covariate values for every observation.
#' 
#' @return single numeric value giving the sum of negative log-likelihoods
#' @export
EffectResponse_lambda <- function(init.par, target_all, density_all, log.fitness, covariates = NULL){
  
  num.sp <- nrow(target_all)
  
  lambda.vector <- init.par[1:num.sp]
  r.vector <- init.par[(num.sp+1):(num.sp+num.sp)]
  if(is.null(covariates)){
    e.vector <- init.par[(num.sp+1+num.sp):(length(init.par)-1)]
  }else{
    e.vector <- init.par[(num.sp+1+num.sp):(num.sp+num.sp+num.sp)]
    lambda.cov <- matrix(init.par[(num.sp + num.sp + num.sp + 1):(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)))],
                         nrow = num.sp,ncol = ncol(covariates))
    response.cov <- matrix(init.par[(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + 1):(num.sp + 
                                                                                         num.sp + 
                                                                                         num.sp + 
                                                                                         (num.sp*ncol(covariates)) + 
                                                                                         (num.sp*ncol(covariates)))],
                           nrow = num.sp,ncol = ncol(covariates))
    effect.cov <- matrix(init.par[(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + 1):(num.sp + 
                                                                                                                    num.sp + 
                                                                                                                    num.sp +
                                                                                                                    (num.sp*ncol(covariates)) +
                                                                                                                    (num.sp*ncol(covariates)) +
                                                                                                                    (num.sp*ncol(covariates)))],
                          nrow = num.sp,ncol = ncol(covariates))
  }
  sigma <- init.par[length(init.par)]
  
  if(is.null(covariates)){
    lambda.part <- colSums(lambda.vector*target_all)
    r.part <- colSums(r.vector*target_all)
    e.part <- colSums(e.vector*density_all)
  }else{
    lambda.cov_all <- matrix(0,nrow = num.sp,ncol = length(log.fitness))
    response.cov_all <- matrix(0,nrow = num.sp,ncol = length(log.fitness))
    effect.cov_all <- matrix(0,nrow = num.sp,ncol = length(log.fitness))
    for(i.sp in 1:num.sp){
      for(i.obs in 1:length(log.fitness)){
        lambda.cov_all[i.sp,i.obs] <- sum(lambda.cov[i.sp,]*covariates[i.obs,])
        response.cov_all[i.sp,i.obs] <- sum(response.cov[i.sp,]*covariates[i.obs,])
        effect.cov_all[i.sp,i.obs] <- sum(effect.cov[i.sp,]*covariates[i.obs,])
        
      }# for i.cov
    }# for i.sp
    
    lambda.part <- colSums(lambda.vector*(1+lambda.cov_all)*target_all)
    r.part <- colSums(r.vector*(1+response.cov_all)*target_all)
    e.part <- colSums(e.vector*(1+effect.cov_all)*density_all)
  }# if-else null covariates
  
  pred <- lambda.part/ (1+ e.part*r.part )
  
  # likelihood as before:
  llik<-dnorm(log.fitness,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}