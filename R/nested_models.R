# function to optimize over
################

# model1 <- function(par, log_fitness){ 
model1 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  #lambda and sigma parameters for the normal distribution
  #(assuming lognormal error- seed data are logged) 
  lambda <- par[1]
  sigma <- par[2]
  #this is the predictive model- here is just fitting a horizontal
  #line through the data:
  pred <- rep(lambda, times=length(log_fitness)) 
  #these are the log likelihoods of the data given the model + parameters
  llik <- dnorm(log_fitness, mean = log(pred), sd = (sigma), log = TRUE) #CHEK PRED IS LOG!
  #llik <- dnorm(log_fitness, mean = mean(log(pred)), sd = (sigma), log = TRUE) 
  #llik <- dnorm(log_fitness, mean = mean(pred), sd = (sigma), log = TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}

################

# model2 <- function(par, log_fitness, background){ 
model2 <-  function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  background <- rowSums(focal.comp.matrix)
  ## apply is used to get a single vector of densities rather than a matrix
  ## there is only one competitor at a time here- just adding each density to a 
  ## column of other zeros:
  ## predictive model:
  pred <- lambda/(1+alpha*(background))  
  ## log likelihoods of data given the model + parameters:
  llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
  #hist(llik)
  #llik <- dnorm(log_fitness, mean = mean(log(pred)), sd = (sigma), log = TRUE)
  ## return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}

################

# model3 <- function(par, log_fitness, focal.comp.matrix){
model3 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] #same as model 2
  alpha.vector <- par[2:(length(par)-1)] # new parameters- use alpha estimate from model 2 as start 
  #value for fitting
  sigma <- par[length(par)] ## same as model 2
  ## predictive model:
  term = 1 #create the denominator term for the model
  for(z in 1:ncol(focal.comp.matrix)){
    term <- term + alpha.vector[z]*focal.comp.matrix[,z] 
    #THIS IS AS DONE IN LINCX, CHECK WITH OSCAR
  }
  pred <- lambda/ term
  # likelihood as before:
  llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

################

model4 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] 
  lambda.cov <- par[(1+1):(1+num.covariates)] #effect of cov 1, 2, ... on lambda
  alpha.vector <- par[(1+num.covariates+1):(1+num.covariates+num.competitors)] # alfas_ij
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
  return(sum(-1*llik)) #sum of negative log likelihoods
}

################

model5 <- function(par, log_fitness, focal.comp.matrix, num.covariates, num.competitors, focal.covariates){
  lambda <- par[1] 
  lambda.cov <- par[(1+1):(1+num.covariates)] #effect of cov 1, 2... on lambda
  alpha.vector <- par[(1+num.covariates+1):(1+num.covariates+num.competitors)] #alfas_ij
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
  #here I need to reformat cov_term_x to sumatories of the form alpha.cov_i* cov_i + alpha.covj* cov_j + ...
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
  return(sum(-1*llik)) #sum of negative log likelihoods
}