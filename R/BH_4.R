#' Title Beverton-Holt fecundity, fourth model
#'
#' These functions return the negative log-likelihood of the data
#' given the model and parameters. TODO: BH_4 is $F_i = \lambda_i/(1+\alpha_{ij}*N_{j})$
#' 
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second lambda.cov, effects of every covariate on lambda; 
#' third alpha, interaction coefficients with every species; 
#' fourth alpha.cov, effects of every covariate on alpha values (single effect for each covariate); 
#' last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the \code{fixed.terms} list
#' @param param.list string listing parameters to optimize. Possible elements are \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}.
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates number of covariates
#' @param num.competitors number of competitor species
#' @param focal.covariates dataframe/matrix with as many rows as observationes, and one column for each covariate.
#' Values of the dataframe are covariate values for every observation.
#' @param fixed.terms list with elements \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}. It contains parameters not to be optimized.
#' Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
#' to look for it in \code{par}, and will not consider it here. In this model, \code{lambda.cov} and \code{alpha.cov} are not considered.
#'
#' @return log-likelihood value
#' @export
BH_4 <- function(par, param.list, log.fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates, fixed.terms){
  
  pos <- 1
  if("lambda" %in% param.list){
    lambda <- par[pos] ## same as model 1
    pos <- pos + 1
  }else{
    lambda <- fixed.terms[["lambda"]]
  }
  
  if("lambda.cov" %in% param.list){
    lambda.cov <- par[pos:(pos+num.covariates-1)]
    pos <- pos + num.covariates
  }else{
    lambda.cov <- fixed.terms[["lambda.cov"]]
  }
  
  if("alpha" %in% param.list){
    alpha <- par[pos:(pos+num.competitors-1)]
    pos <- pos + num.competitors
  }else{
    alpha <- fixed.terms[["alpha"]]
  }
  
  if("alpha.cov" %in% param.list){
    alpha.cov <- par[pos:(pos+num.covariates-1)]
    pos <- pos + num.covariates
  }else{
    alpha.cov <- fixed.terms[["alpha.cov"]]
  }
  
  sigma <- par[length(par)]
  
  num = 1
  focal.cov.matrix <- as.matrix(focal.covariates)
  for(z in 1:num.covariates){
    num <- num + lambda.cov[z]*focal.cov.matrix[,z]
  }
  cov_term <- 0 
  for(v in 1:num.covariates){
    cov_term <- cov_term + alpha.cov[v] * focal.cov.matrix[,v]
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:ncol(focal.comp.matrix)){
    term <- term + (alpha[z] + cov_term) * focal.comp.matrix[,z] 
  }
  pred <- lambda * (num) / term 
  # likelihood as before:
  
  llik<-dnorm(log.fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik))
}