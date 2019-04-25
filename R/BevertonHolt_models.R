# functions to optimize over
# all functions accept the same parameters, even though 
# some of them are not used, for easing replication

################
# Set of Beverton-Holt models for fecundity

#' Title
#'
#' @param par vector of length 2 with lambda and sigma values
#' @param log_fitness log of the fitness value 
#' @param focal.comp.matrix not used in BH_1
#' @param num.covariates not used in BH_1
#' @param num.competitors not used in BH_1
#' @param focal.covariates not used in BH_1
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_1 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  #lambda and sigma parameters for the normal distribution
  #(assuming lognormal error- seed data are logged) 
  lambda <- par[1]
  sigma <- par[2]
  #this is the predictive model- here is just fitting a horizontal
  #line through the data:
  pred <- rep(lambda, times=length(log_fitness)) 
  #these are the log likelihoods of the data given the model + parameters
  llik <- dnorm(log_fitness, mean = log(pred), sd = (sigma), log = TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}

################

#' Title
#'
#' @param par vector of length 3, with values for lambda, alpha, and sigma
#' @param log_fitness log of the fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates not used in BH_2
#' @param num.competitors not used in BH_2 
#' @param focal.covariates not used in BH_2
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_2 <-  function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  background <- rowSums(focal.comp.matrix)
  # predictive model:
  pred <- lambda/(1+alpha*(background))  
  # log likelihoods of data given the model + parameters:
  llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}

################

#' Title
#'
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second, interaction coefficients with every species; last, sigma value
#' @param log_fitness log of the fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates not used in BH_3
#' @param num.competitors not used in BH_3
#' @param focal.covariates not used in BH_3
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_3 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] #same as model 2
  alpha.vector <- par[2:(length(par)-1)] # new parameters- use alpha estimate from model 2 as start 
  # value for fitting
  sigma <- par[length(par)] ## same as model 2
  # predictive model:
  term = 1 #create the denominator term for the model
  for(z in 1:ncol(focal.comp.matrix)){
    term <- term + alpha.vector[z]*focal.comp.matrix[,z] 
  }
  pred <- lambda/ term
  # likelihood as before:
  llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik)) 
}

################

#' Title
#'
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second, effects of every covariate on lambda; 
#' third, interaction coefficients with every species; 
#' fourth, effects of every covariate on alpha values (same effect on all interaction coefficients); 
#' last, sigma value
#' @param log_fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates number of covariates
#' @param num.competitors number of competitor species
#' @param focal.covariates dataframe with as many rows as observationes, and one column for each covariate.
#' Values of the dataframe are covariate values for every observation.
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_4 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] 
  lambda.cov <- par[(1+1):(1+num.covariates)] #effect of cov 1, 2, ... on lambda
  alpha.vector <- par[(1+num.covariates+1):(1+num.covariates+num.competitors)] # alpha_ij
  alpha.cov <- par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+num.covariates)] # common effect of cov 1, 2... on alphas
  sigma <- par[length(par)]
  num = 1
  for(z in 1:num.covariates){
    num <- num + lambda.cov[z]*focal.covariates[,z] 
  }
  cov_term <- 0 
  for(v in 1:num.covariates){
    cov_term <- cov_term + alpha.cov[v] * focal.covariates[,v]
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:ncol(focal.comp.matrix)){
    term <- term + (alpha.vector[z] + cov_term) * focal.comp.matrix[,z] 
  }
  pred <- lambda * (num) / term 
  # likelihood as before:
  llik<-dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik))
}

################

#' Title
#'
#' @param par vector of variable length, with the following order: first, lambda of focal sp; 
#' second, effects of every covariate on lambda; 
#' third, interaction coefficients with every species; 
#' fourth, effects of every covariate on alpha values (varying effect of every covariate over every interaction coefficient); 
#' last, sigma value
#' @param log_fitness log of fitness value
#' @param focal.comp.matrix dataframe with as many rows as observations, and one column for each competitor sp. 
#' Values of the dataframe are number of competitors of each sp per observation.
#' @param num.covariates number of covariates
#' @param num.competitors number of competitor species
#' @param focal.covariates dataframe with as many rows as observationes, and one column for each covariate.
#' Values of the dataframe are covariate values for every observation.
#'
#' @return log-likelihood value
#' @export
#'
#' @examples
BH_5 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] 
  lambda.cov <- par[(1+1):(1+num.covariates)] #effect of cov 1, 2... on lambda
  alpha.vector <- par[(1+num.covariates+1):(1+num.covariates+num.competitors)] #alpha_ij
  alpha.cov <- par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.covariates*num.competitors))] #effects of cov 1, 2... on alpha_i
  sigma <- par[length(par)]
  num = 1
  for(v in 1:num.covariates){
    num <- num + lambda.cov[v]*focal.covariates[,v] 
  }
  cov_term_x <- list()
  for(v in 1:num.covariates){
    cov_temp <- focal.covariates[,v]
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
    term <- term + (alpha.vector[z] + cov_term[[z]]) * focal.comp.matrix[,z]  
  }
  pred <- lambda * (num) / term 
  # likelihood as before:
  llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik))
}