#' Title Beverton-Holt fecundity, first model
#' 
#' These functions return the negative log-likelihood of the data
#' given the model and parameters. model_BH1 is \eqn{F_i = \lambda_i}
#'
#' @param par vector containing lambda of focal sp and sigma value
#' @param param.list not used in model_BH1
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix not used in model_BH1
#' @param num.covariates not used in model_BH1
#' @param num.competitors not used in model_BH1
#' @param focal.covariates not used in model_BH1
#' @param fixed.terms not used in model_BH1
#'
#' @return log-likelihood value
#' @import stats 
#' @export
model_BH1 <- function(par, 
                 param.list = NULL, 
                 log.fitness, 
                 focal.comp.matrix = NULL, 
                 num.covariates = NULL, 
                 num.competitors = NULL, 
                 focal.covariates = NULL, 
                 fixed.terms = NULL){
  #lambda and sigma parameters for the normal distribution
  #(assuming lognormal error- seed data are logged) 
  lambda <- par[1]
  sigma <- par[2]
  #this is the predictive model- here is just fitting a horizontal
  #line through the data:
  pred <- rep(lambda, times=length(log.fitness)) 
  #these are the log likelihoods of the data given the model + parameters
  llik <- dnorm(log.fitness, mean = log(pred), sd = (sigma), log = TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}
