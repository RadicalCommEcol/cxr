
################
# Set of Beverton-Holt fecundity models  
# all functions accept the same parameters, even though 
# some of them are not used, for easing replication
# These functions return the negative log-likelihood of the data
# given the model and parameters

#' Title Beverton-Holt fecundity, first model
#'
#' @param par vector containing lambda of focal sp and sigma value
#' @param param.list not used in BH_1
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix not used in BH_1
#' @param num.covariates not used in BH_1
#' @param num.competitors not used in BH_1
#' @param focal.covariates not used in BH_1
#' @param fixed.terms not used in BH_1
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_1 <- function(par, param.list, log.fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates, fixed.terms){
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

################

#' Title Beverton-Holt fecundity, second model
#'
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second alpha, single interaction coefficient; 
#' last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the "fixed.terms" list
#' @param param.list string listing parameters to optimize. Possible elements are "lambda", "lambda.cov", "alpha", "alpha.cov".
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates not used in BH_2
#' @param num.competitors not used in BH_2
#' @param focal.covariates not used in BH_2.
#' @param fixed.terms list with elements "lambda", "lambda.cov", "alpha", "alpha.cov". It contains parameters not to be optimized.
#' Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
#' to look for it in "par", and will not consider it here. In this model, "lambda.cov" and "alpha.cov" are not considered.
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_2 <- function(par, param.list, log.fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates, fixed.terms){
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

################

#' Title Beverton-Holt fecundity, third model
#'
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second alpha, interaction coefficients with every species; 
#' last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the "fixed.terms" list
#' @param param.list string listing parameters to optimize. Possible elements are "lambda", "lambda.cov", "alpha", "alpha.cov".
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates not used in BH_3
#' @param num.competitors number of competitor species
#' @param focal.covariates not used in BH_3.
#' @param fixed.terms list with elements "lambda", "lambda.cov", "alpha", "alpha.cov". It contains parameters not to be optimized.
#' Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
#' to look for it in "par", and will not consider it here. In this model, "lambda.cov" and "alpha.cov" are not considered.
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_3 <- function(par, param.list, log.fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates, fixed.terms){
  
  pos <- 1
  if("lambda" %in% param.list){
    lambda <- par[pos] ## same as model 1
    pos <- pos + 1
  }else{
    lambda <- fixed.terms[["lambda"]]
  }
  
  if("alpha" %in% param.list){
    alpha <- par[pos:(pos+num.competitors-1)]
    pos <- pos + num.competitors
  }else{
    alpha <- fixed.terms[["alpha"]]
  }
  
  # lambda <- par[1] #same as model 2
  # alpha.vector <- par[2:(length(par)-1)] # new parameters- use alpha estimate from model 2 as start 
  # # value for fitting
  sigma <- par[length(par)] ## same as model 2
  # predictive model:
  term = 1 #create the denominator term for the model
  for(z in 1:ncol(focal.comp.matrix)){
    term <- term + alpha[z]*focal.comp.matrix[,z] 
  }
  pred <- lambda/ term
  # likelihood as before:
  llik <- dnorm(log.fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik)) 
}

################

#' Title Beverton-Holt fecundity, fourth model
#'
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second lambda.cov, effects of every covariate on lambda; 
#' third alpha, interaction coefficients with every species; 
#' fourth alpha.cov, effects of every covariate on alpha values (single effect for each covariate); 
#' last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the "fixed.terms" list
#' @param param.list string listing parameters to optimize. Possible elements are "lambda", "lambda.cov", "alpha", "alpha.cov".
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates number of covariates
#' @param num.competitors number of competitor species
#' @param focal.covariates dataframe/matrix with as many rows as observationes, and one column for each covariate.
#' Values of the dataframe are covariate values for every observation.
#' @param fixed.terms list with elements "lambda", "lambda.cov", "alpha", "alpha.cov". It contains parameters not to be optimized.
#' Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
#' to look for it in "par", and will not consider it here.
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
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

################

#' Title Beverton-Holt fecundity, fifth model
#'
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second lambda.cov, effects of every covariate on lambda; 
#' third alpha, interaction coefficients with every species; 
#' fourth alpha.cov, effects of every covariate on alpha values (varying effect of every covariate over every interaction coefficient); 
#' last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the "fixed.terms" list
#' @param param.list string listing parameters to optimize. Possible elements are "lambda", "lambda.cov", "alpha", "alpha.cov".
#' @param log.fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates number of covariates
#' @param num.competitors number of competitor species
#' @param focal.covariates dataframe with as many rows as observationes, and one column for each covariate.
#' Values of the dataframe are covariate values for every observation.
#' @param fixed.terms list with elements "lambda", "lambda.cov", "alpha", "alpha.cov". It contains parameters not to be optimized.
#' Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
#' to look for it in "par", and will not consider it here.
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_5 <- function(par, param.list, log.fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates, fixed.terms){
  
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