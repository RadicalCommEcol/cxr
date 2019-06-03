#' Title ML estimation of effect-response function for annual plants
#' 
#' Calculates the log-likelihood of a Beverton-Holt model parameterized with given values
#' with respect to a fitness metric. The function for calculating fecundity given 
#' effect and response values is taken from Godoy et al. (2014). 
#' Note that, as e is not pair-specific, all species parameters are fit in the same function.
#' In this version, lambda values are fixed.
#'
#' @param init.par 1d vector of initial parameters: r values, e values, and single sigma term
#' @param lambda 1d vector of lambda values
#' @param target_all matrix giving which species is calculated with which values. See ER_optimize
#' @param density_all matrix giving the densities of each species at each observation. See ER_optimize
#' @param log.fitness log of the fitness metric
#' @param covariates dataframe/matrix with as many rows as observationes, and one column for each covariate.
#' Values are covariate values for every observation.
#'
#' @return single numeric value giving the sum of negative log-likelihoods
#' @export
#'
#' @examples
EffectResponse <- function(init.par, lambda, target_all, density_all, log.fitness){
  
  r.vector <- init.par[1:length(lambda)]
  # r.vector <- par[(num.focal+1):(num.focal+num.focal)]
  e.vector <- init.par[(length(lambda)+1):(length(init.par)-1)]
  
  sigma <- init.par[length(init.par)]
  
  lambda.part <- colSums(lambda*target_all)
  r.part <- colSums(r.vector*target_all)
  e.part <- colSums(e.vector*density_all)
  
  pred <- lambda.part/ (1+ e.part*r.part )
  
  # likelihood as before:
  llik<-dnorm(log.fitness,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

#' Title ML estimation of effect-response function for annual plants
#' 
#' Calculates the log-likelihood of a Beverton-Holt model parameterized with given values
#' with respect to a fitness metric. The function for calculating fecundity given 
#' effect and response values is taken from Godoy et al. (2014). 
#' Note that, as e is not pair-specific, all species parameters are fit in the same function.
#' In this version, lambda values are also fit.
#'
#' @param init.par 1d vector of initial parameters: lambda values, r values, e values, and single sigma term
#' @param target_all matrix giving which species is calculated with which values. See ER_optimize
#' @param density_all matrix giving the densities of each species at each observation. See ER_optimize
#' @param log.fitness log of the fitness metric
#'
#' @return single numeric value giving the sum of negative log-likelihoods
#' @export
#'
#' @examples
EffectResponse_lambda <- function(init.par, target_all, density_all, log.fitness){
  
  lambda.vector <- init.par[1:nrow(target_all)]
  r.vector <- init.par[(nrow(target_all)+1):(nrow(target_all)+nrow(target_all))]
  e.vector <- init.par[(nrow(target_all)+1+nrow(target_all)):(length(init.par)-1)]
  
  sigma<-init.par[length(init.par)]
  
  lambda.part <- colSums(lambda.vector*target_all)
  r.part <- colSums(r.vector*target_all)
  e.part <- colSums(e.vector*density_all)
  
  pred <- lambda.part/ (1+ e.part*r.part )
  
  # likelihood as before:
  llik<-dnorm(log.fitness,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}