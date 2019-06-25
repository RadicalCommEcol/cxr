#' Title Beverton-Holt fecundity, second model
#'
#' These functions return the negative log-likelihood of the data
#' given the model and parameters. BH_2 is \eqn{F_i = \frac{\lambda_i}{1+\alpha \sum_jN_{j}}}
#' 
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second alpha, single interaction coefficient; 
#' last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the "fixed.terms" list
#' @param param.list string listing parameters to optimize. Possible elements are \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}.
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates not used in BH_2
#' @param num.competitors not used in BH_2
#' @param focal.covariates not used in BH_2.
#' @param fixed.terms list with elements \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}. It contains parameters not to be optimized.
#' Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
#' to look for it in \code{par}, and will not consider it here. In this model, \code{lambda.cov} and \code{alpha.cov} are not considered.
#'
#' @return log-likelihood value
#' @import stats 
#' @export
BH_2 <- function(par, 
                 param.list = c("lambda","alpha"), 
                 log.fitness, 
                 focal.comp.matrix, 
                 num.covariates = NULL, 
                 num.competitors = NULL, 
                 focal.covariates = NULL, 
                 fixed.terms = NULL){
  pos <- 1
  if("lambda" %in% param.list){
    lambda <- par[pos] ## same as model 1
    pos <- pos + 1
  }else{
    lambda <- fixed.terms[["lambda"]]
  }
  
  if("alpha" %in% param.list){
    alpha <- par[pos]
    pos <- pos + 1
  }else{
    alpha <- fixed.terms[["alpha"]]
  }
  
  sigma <- par[length(par)] ## same as model 1
  background <- rowSums(focal.comp.matrix)
  # predictive model:
  pred <- lambda/(1+alpha*(background))  
  # log likelihoods of data given the model + parameters:
  llik <- dnorm(log.fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}