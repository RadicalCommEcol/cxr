#' Title Beverton-Holt fecundity, fifth model
#' 
#' These functions return the negative log-likelihood of the data
#' given the model and parameters. BH_5 is 
#' \eqn{F_i = \frac{\lambda_i + \sum_k b_k \Theta_{i,k}}{1+\sum_j(\alpha_{ij} +\sum_k b_k \phi_{i,j,k})N_{j}}}
#' 
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second lambda.cov, effects of every covariate on lambda; 
#' third alpha, interaction coefficients with every species; 
#' fourth alpha.cov, effects of every covariate on alpha values (varying effect of every covariate over every interaction coefficient); 
#' last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the \code{fixed.terms} list
#' @param param.list string listing parameters to optimize. Possible elements are \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}.
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates number of covariates
#' @param num.competitors number of competitor species
#' @param focal.covariates dataframe with as many rows as observationes, and one column for each covariate.
#' Values of the dataframe are covariate values for every observation.
#' @param fixed.terms list with elements \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}. It contains parameters not to be optimized.
#' Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
#' to look for it in \code{par}, and will not consider it here. In this model, \code{lambda.cov} and \code{alpha.cov} are not considered.
#'
#' @return log-likelihood value
#' @import stats 
#' @export
BH_5 <- function(par, 
                 param.list, 
                 log.fitness, 
                 focal.comp.matrix, 
                 num.covariates, 
                 num.competitors, 
                 focal.covariates, 
                 fixed.terms){
  
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
    alpha.cov <- par[pos:(pos+(num.covariates*num.competitors)-1)]
    pos <- pos + (num.covariates*num.competitors)
  }else{
    alpha.cov <- fixed.terms[["alpha.cov"]]
  }
  
  sigma <- par[length(par)]
  
  num = 1
  focal.cov.matrix <- as.matrix(focal.covariates)
  for(v in 1:num.covariates){
    num <- num + lambda.cov[v]*focal.cov.matrix[,v] 
  }
  cov_term_x <- list()
  for(v in 1:num.covariates){
    cov_temp <- focal.cov.matrix[,v]
    for(z in 1:num.competitors){
      cov_term_x[[z+(num.competitors*(v-1))]] <- alpha.cov[z+(num.competitors*(v-1))] * cov_temp  #create  alpha.cov_i*cov_i vector
    }
  }
  cov_term <- list()
  for(z in 0:(num.competitors-1)){
    cov_term_x_sum <- cov_term_x[[z+1]]
    if(num.covariates > 1){
      for(v in 2:num.covariates){
        cov_term_x_sum <- cov_term_x_sum + cov_term_x[[v + num.competitors]]
      } 
    }
    cov_term[[z+1]] <- cov_term_x_sum
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:num.competitors){
    term <- term + (alpha[z] + cov_term[[z]]) * focal.comp.matrix[,z]  
  }
  pred <- lambda * (num) / term 
  # likelihood as before:
  llik <- dnorm(log.fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik))
}
