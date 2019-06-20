source("R/BH_1.R")
source("R/BH_2.R")
source("R/BH_3.R")
source("R/BH_4.R")
source("R/BH_5.R")
source("R/AIC.R")
source("R/SEbootstrap_caracoles.R")
source("R/cxr_optimize.R")
source("R/InitParams.R")
source("R/RetrieveParams.R")
# optimization methods
library(nloptr)
library(GenSA)
library(hydroPSO)
library(DEoptimR)
#other packages
library(tidyverse)
###########################
# Caracoles competition data
competition.data <- readr::read_delim(file = "./data/competition.csv",delim = ";")

# spread the data from long to wide format
competition.data <- spread(competition.data,competitor,number,fill = 0)
# how many focal species
focal.sp <- unique(competition.data$focal)

# competition matrix
comp.matrix <- as.matrix(competition.data[,10:ncol(competition.data)])

# covariate: salinity
salinity <- readr::read_delim(file = "./data/salinity.csv",delim = ";")

# one observation per row of competition.data
salinity <- salinity[,c("plot","subplot","year","sum_salinity")]
full.data <- left_join(competition.data,salinity)

# in case we have independent estimates of lambda and/or do not want to optimize it 
# init.lambda <- readr::read_delim(file = "./results/lambda_estimates.csv",delim = ";")
# init.lambda <- arrange(subset(init.lambda, model == max(init.lambda$model) & optim.method == "optim_NM"),focal.sp)
# init.lambda <- init.lambda[,c("lambda")]

# same with other parameters, read them here, modify param.list accordingly, and comment out the appropriate lines below.
# init.alpha <- matrix(rnorm(ncol(comp.matrix)*ncol(comp.matrix),1,0.01),nrow = ncol(comp.matrix),ncol = ncol(comp.matrix))

#############################
# simulation parameters

# models to parameterize
# be aware of including 4 and 5 ONLY if there are covariates
# otherwise it makes no sense (see equations in Lanuza et al. 2018)
models <- 3:4

# which values do we optimize for each model?
param.list <- list(c("lambda","alpha"),c("lambda","alpha","lambda.cov","alpha.cov","lambda.cov_NL","alpha.cov_NL"))

#Choose the non-linearity function
function1 <-function(a,b,x){
  vector <-vector("numeric",length(x))
  for(i in 1:length(x)){
  if (x[i]!=0){vector[i]<- (a*x[i]^2/(b+x[i]^2))}
  }
  return(vector)
  }

function2 <-function(a,b,x){
  for(i in 1:length(x)){
    vector <-vector("numeric",length(x))
    if (x[i]!=0){vector[i]<- (a*x[i]^2/(b+x[i]^2))}
  }
  return(vector)
}
function_NL <- function1
# keep the model definitions in a list, for ease
fitness.models <- list(BH_1 = BH_1,BH_2 = BH_2,BH_3 = BH_3,BH_4 = BH_4,BH_5 = BH_5)

# environmental covariates
covariates <- full.data[,"sum_salinity"]
# if no covariates, comment above and uncomment here
# covariates <- 0

# optimization methods to use
optim.methods <- c(#"optim_NM",
                   "optim_L-BFGS-B"
                  #"nloptr_CRS2_LM"
                   # "nloptr_ISRES"
                   # "nloptr_DIRECT_L_RAND"
                    # "GenSA"
                   # "hydroPSO"
                   # "DEoptimR"
)

# from which method are we taking initial estimates for the next model?
# preliminary observations suggest optim_NM or DEoptimR
init.par.method <- optim.methods[1]
init.method.num <- which(optim.methods == init.par.method)

# values for initial estimates of parameters. 
# Overwrite if necessary
# init.lambda is calculated after log.fitness for each focal species
lower.lambda <- 1
upper.lambda <- 1e4
# sigma
lower.sigma <- 0.000001
upper.sigma <- 1
# alpha
init.alpha <- 1e-4
lower.alpha <- 1e-4
upper.alpha <- 1e5
# lambda.cov
init.lambda.cov <- 1
lower.lambda.cov <- 1e-4
upper.lambda.cov <- 1e4
# alpha.cov
init.alpha.cov <- 1
lower.alpha.cov <- 1e-4
upper.alpha.cov <- 1e4
# lambda.cov_NL
init.lambda.cov_NL <- 1
lower.lambda.cov_NL <- 1e-4
upper.lambda.cov_NL <- 1e2
# alpha.cov
init.alpha.cov_NL <- 1
lower.alpha.cov_NL <- 1e-4
upper.alpha.cov_NL <- 1e4


# if we want quicker calculations, we can disable 
# the bootstrapping for the standard errors
generate.errors <- FALSE
bootstrap.samples <- 3

# store results?
write.results <- TRUE

##############################
# initialize data structures
# elements are matrices/vectors within nested lists, of the form matrix[[focal.sp]][[model]][[method]]

param.matrices <- list()
for(i.sp in 1:length(focal.sp)){
  param.matrices[[i.sp]] <- list()
  for(i.model in 1:length(models)){
    
    param.matrices[[i.sp]][[i.model]] <- list()
    for(i.method in 1:length(optim.methods)){
      param.matrices[[i.sp]][[i.model]][[i.method]] <- list(lambda = 0,
                                                            lambda.lower.error = 0,
                                                            lambda.upper.error = 0,
                                                            sigma = 0,
                                                            alpha = 0,
                                                            alpha.lower.error = 0,
                                                            alpha.upper.error = 0,
                                                            lambda.cov = 0,
                                                            lambda.cov.lower.error = 0,
                                                            lambda.cov.upper.error = 0,
                                                            alpha.cov = 0,
                                                            alpha.cov.lower.error = 0,
                                                            alpha.cov.upper.error = 0,
                                                            lambda.cov_NL = 0,
                                                            lambda.cov_NL.lower.error = 0,
                                                            lambda.cov_NL.upper.error = 0,
                                                            alpha.cov_NL = 0,
                                                            alpha.cov_NL.lower.error = 0,
                                                            alpha.cov_NL.upper.error = 0,
                                                            log.likelihood = 0)
    }
    names(param.matrices[[i.sp]][[i.model]]) <- optim.methods
  }
  names(param.matrices[[i.sp]]) <- names(fitness.models)[models]
}
names(param.matrices) <- focal.sp

###############################
# main loop

for(i.sp in 1:length(focal.sp)){
  
  # subset and prepare the data
  
  focal.sp.data <- subset(competition.data, focal == focal.sp[i.sp])
  # current focal species
  focal <- unique(focal.sp.data$focal)
  # fitness metric...
  # subset >0 records, for calculating logarithms
  focal.sp.data <- subset(focal.sp.data, seed > 0)
  fitness <- focal.sp.data$seed #fitness
  log.fitness <- log(fitness)
  # competition matrix: number of competitors
  focal.comp.matrix <- comp.matrix[which(competition.data$focal == focal.sp[i.sp]),]
  # number of competitors
  num.competitors <- dim(focal.comp.matrix)[2]
  # number of covariates
  num.covariates <- ifelse(is.null(ncol(covariates)),0,ncol(covariates))
  # covariates for the focal species
  if(num.covariates > 0){
    focal.covariates <- covariates[which(competition.data$focal == focal.sp[i.sp]),,drop = FALSE]
  }else{
    focal.covariates <- 0
  }
  
  # generate initial values for the different parameters
  # or gather them from data if they are not to be optimized
  
  # lambda
  if("lambda" %in% param.list[[i.model]]){
    current.init.lambda <- mean(log.fitness)
  }else{
    current.init.lambda <- init.lambda[i.sp]
  }
  # sigma
  current.init.sigma <- sd(log.fitness)
  if(current.init.sigma > upper.sigma){
    current.init.sigma <- upper.sigma
  }
  # alpha
  if("alpha" %in% param.list[[1]]){
    if(models[1]<=2){
      alpha.length <- 1
    }else{
      alpha.length <- num.competitors
    }
    if(length(init.alpha) != alpha.length){
      current.init.alpha <- rep(init.alpha[1],alpha.length) 
    }else{
      current.init.alpha <- init.alpha
    }
  }
  # lambda.cov
  if("lambda.cov" %in% param.list[[2]]){
    if(length(init.lambda.cov) != num.covariates){
      current.init.lambda.cov <- rep(init.lambda.cov[1],num.covariates)
    }else{
      current.init.lambda.cov <- init.lambda.cov  
    }
  }else{
    current.init.lambda.cov <- init.lambda.cov[i.sp]  
  }
  

  #lamda.cov_NL
  if("lambda.cov_NL" %in% param.list[[2]]){
    if(length(init.lambda.cov_NL) != num.covariates){
      current.init.lambda.cov_NL <- rep(init.lambda.cov_NL[1],num.covariates)
    }else{
      current.init.lambda.cov_NL <- init.lambda.cov_NL  
    }
  }
  
  # alpha.cov
  if("alpha.cov" %in% param.list[[2]]){
    if(models[1]<=4){
      length.alpha.cov <- num.covariates
    }else if(models[1]>4){
      length.alpha.cov <- num.covariates*num.competitors
    }
    if(length(init.alpha.cov) != length.alpha.cov){
      current.init.alpha.cov <- rep(init.alpha.cov[1],length.alpha.cov)
    }else{
      current.init.alpha.cov <- init.alpha.cov  
    }  
  }else{
    current.init.alpha.cov <- init.alpha.cov[i.sp]  
  }
  
  # alpha.cov_NL
  if("alpha.cov_NL" %in% param.list[[2]]){
    if(models[1]<=4){
      length.alpha.cov_NL <- num.covariates
    }else if(models[1]>4){
      length.alpha.cov_NL <- num.covariates*num.competitors
    }
    if(length(init.alpha.cov_NL) != length.alpha.cov){
      current.init.alpha.cov_NL <- rep(init.alpha.cov_NL[1],length.alpha.cov)
    }else{
      current.init.alpha.cov_NL <- init.alpha.cov_NL  
    }  
  }

  # model to optimize  
  for(i.model in 1:length(models)){
    
    print("*********************************")
    print(paste(date()," - starting focal sp ",focal.sp[i.sp],", model ",models[i.model],sep=""))
    print("*********************************")

    # obtain initial estimates from either previous model, or from given values
    # also, beware if the initial estimates are single values or vectors.
    
    ######################
    # compute each method
    
    for(i.method in 1:length(optim.methods)){
      
      temp.results <- cxr_optimize(fitness.model = fitness.models[[models[i.model]]],
                                   optim.method = optim.methods[i.method],
                                   param.list = param.list[[i.model]],
                                   log.fitness = log.fitness,
                                   init.lambda = current.init.lambda,
                                   lower.lambda = lower.lambda,
                                   upper.lambda = upper.lambda,
                                   init.sigma = current.init.sigma,
                                   lower.sigma = lower.sigma,
                                   upper.sigma = upper.sigma,
                                   init.alpha = current.init.alpha,
                                   lower.alpha = lower.alpha,
                                   upper.alpha = upper.alpha,
                                   init.lambda.cov = current.init.lambda.cov,
                                   lower.lambda.cov = lower.lambda.cov,
                                   upper.lambda.cov = upper.lambda.cov,
                                   init.alpha.cov = current.init.alpha.cov,
                                   lower.alpha.cov = lower.alpha.cov,
                                   upper.alpha.cov = upper.alpha.cov,
                                   init.lambda.cov_NL = current.init.lambda.cov_NL,
                                   lower.lambda.cov_NL = lower.lambda.cov_NL,
                                   upper.lambda.cov_NL = upper.lambda.cov_NL,
                                   init.alpha.cov_NL = current.init.alpha.cov_NL,
                                   lower.alpha.cov_NL = lower.alpha.cov_NL,
                                   upper.alpha.cov_NL = upper.alpha.cov_NL,
                                   focal.comp.matrix = focal.comp.matrix,
                                   focal.covariates = focal.covariates,
                                   generate.errors = generate.errors,
                                   bootstrap.samples = bootstrap.samples,
                                   function_NL=function_NL,
                                   number.model=models[i.model])
      ###############
      # clean up results
      
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda <- temp.results$lambda
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.lower.error <- temp.results$lambda.lower.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.upper.error <- temp.results$lambda.upper.error
      
      param.matrices[[i.sp]][[i.model]][[i.method]]$sigma <- temp.results$sigma
      
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha <- temp.results$alpha
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.upper.error <- temp.results$alpha.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.lower.error <- temp.results$alpha.lower.error

      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov <- temp.results$lambda.cov
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov.upper.error <- temp.results$lambda.cov.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov.lower.error <- temp.results$lambda.cov.lower.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov_NL <- temp.results$lambda.cov_NL
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov_NL.upper.error <- temp.results$lambda.cov_NL.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov_NL.lower.error <- temp.results$lambda.cov_NL.lower.error
      
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov <- temp.results$alpha.cov
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov.upper.error <- temp.results$alpha.cov.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov.lower.error <- temp.results$alpha.cov.lower.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov_NL <- temp.results$alpha.cov_NL
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov_NL.upper.error <- temp.results$alpha.cov_NL.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov_NL.lower.error <- temp.results$alpha.cov_NL.lower.error
      
      param.matrices[[i.sp]][[i.model]][[i.method]]$log.likelihood <- temp.results$log.likelihood
      param.matrices[[i.sp]][[i.model]][[i.method]]$AIC <- temp.results$AIC
    }# for i.method
    
    #######################
    # update initial values for the different parameters
    if(i.model < length(param.list)){
    # lambda
    if("lambda" %in% param.list[[i.model+1]]){
      if(!is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda)){
        current.init.lambda <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda
      }
    }
    # sigma
    if(!is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$sigma)){
      current.init.sigma <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$sigma
      if(current.init.sigma > upper.sigma){
        current.init.sigma <- upper.sigma
      }
    }
    # alpha
    if("alpha" %in% param.list[[i.model+1]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha)) == 0){
        current.init.alpha <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha
        # is the current estimate of the appropriate length?
        if(models[i.model] > 2){
          if(length(current.init.alpha) == 1){
            current.init.alpha <- rep(current.init.alpha,num.competitors)
          }
        }# if model > 2
      }
    }
    
    # lambda.cov
    if("lambda.cov" %in% param.list[[i.model+1]]){
     
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov)) == 0){
        current.init.lambda.cov <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov
         }
    }
    
    # lambda.cov_NL
    if("lambda.cov_NL" %in% param.list[[i.model+1]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov_NL)) == 0){
        current.init.lambda.cov_NL <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov_NL
      }
    }
    
    # alpha.cov
    if("alpha.cov" %in% param.list[[i.model+1]]){
            if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov)) == 0){
        current.init.alpha.cov <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov
        # is the current estimate of the appropriate length?
        if(models[i.model+1] > 4){ 
          if(length(current.init.alpha.cov) == num.covariates){
            current.init.alpha.cov <- rep(current.init.alpha.cov,num.competitors)
          }
        }# if model > 4
      }
    }
    
    # alpha.cov_NL
    if("alpha.cov_NL" %in% param.list[[i.model+1]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov_NL)) == 0){
        current.init.alpha.cov_NL <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov_NL
        # is the current estimate of the appropriate length?
        if(models[i.model+1] > 4){
          if(length(current.init.alpha.cov_NL) == num.covariates){
            current.init.alpha.cov_NL <- rep(current.init.alpha.cov_NL,num.competitors)
          }
        }# if model > 4
      }
    }
    } #last model
    
  }# for i.model
}# for i.sp

if(write.results){
  save(param.matrices,file = "./temp_param.matrices.Rdata")
  
  # also, create and store dataframes
  
  # just in case it is a factor
  focal.sp <- sort(as.character(focal.sp))
  competitors <- sort(names(param.matrices[[1]][[names(fitness.models)[models[models == max(models)]]]][[1]]$alpha))
  my.models <- names(fitness.models)[models]
  my.covariates <- c("sum_salinity")
  
  # lambda
  lambda.values <- expand.grid(focal.sp,my.models,optim.methods)
  names(lambda.values) <- c("species","model","method")
  lambda.values$lambda <- 0
  lambda.values$lambda.lower <- 0
  lambda.values$lambda.upper <- 0
  
  # alpha
  alpha.values <- expand.grid(focal.sp,competitors,my.models,optim.methods)
  names(alpha.values) <- c("focal","competitor","model","method")
  alpha.values$alpha <- 0
  alpha.values$alpha.lower <- 0
  alpha.values$alpha.upper <- 0
  
  # lambda.cov -- covariate included "by hand"
  lambda.cov.values <- expand.grid(focal.sp,my.models,optim.methods,my.covariates)
  names(lambda.cov.values) <- c("species","model","method","covariate")
  lambda.cov.values$lambda.cov <- 0
  lambda.cov.values$lambda.cov.lower <- 0
  lambda.cov.values$lambda.cov.upper <- 0
  
  # lambda.cov_NL -- covariate included "by hand"
  lambda.cov_NL.values <- expand.grid(focal.sp,my.models,optim.methods,my.covariates)
  names(lambda.cov_NL.values) <- c("species","model","method","covariate")
  lambda.cov_NL.values$lambda.cov_NL <- 0
  lambda.cov_NL.values$lambda.cov__NL.lower <- 0
  lambda.cov_NL.values$lambda.cov_NL.upper <- 0
  
  # alpha.cov -- covariate included "by hand"
  alpha.cov.values <- expand.grid(focal.sp,competitors,my.models,optim.methods,my.covariates)
  names(alpha.cov.values) <- c("focal","competitor","model","method","covariate")
  alpha.cov.values$alpha.cov <- 0
  alpha.cov.values$alpha.cov.lower <- 0
  alpha.cov.values$alpha.cov.upper <- 0
  
  # alpha.cov_NL -- covariate included "by hand"
  alpha.cov_NL.values <- expand.grid(focal.sp,competitors,my.models,optim.methods,my.covariates)
  names(alpha.cov_NL.values) <- c("focal","competitor","model","method","covariate")
  alpha.cov_NL.values$alpha.cov <- 0
  alpha.cov_NL.values$alpha.cov.lower <- 0
  alpha.cov_NL.values$alpha.cov.upper <- 0
  
  #log-likelihood
  llik <-list()
  for (i in models){
    llik[[i]]<-vector("numeric",length(focal.sp))
    names(llik[[i]])<-focal.sp
  }
  for(i.sp in 1:length(focal.sp)){
    for(i.model in 1: length(models)){
      llik[[models[i.model]]][i.sp]<- param.matrices[[i.sp]][[i.model]][[1]]$log.likelihood
    }
  }
  save(llik,file="data/llik.RData")
    
  # fill up the dataframes
  for(i.sp in 1:length(focal.sp)){
    for(i.model in 1:length(my.models)){
      for(i.method in 1:length(optim.methods)){
        
        # lambda
        lambda.pos <- which(lambda.values$species == focal.sp[i.sp] &
                              lambda.values$model == my.models[i.model] &
                              lambda.values$method == optim.methods[i.method])
        lambda.values$lambda[lambda.pos] <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda
        lambda.values$lambda.lower[lambda.pos] <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.lower.error
        lambda.values$lambda.upper[lambda.pos] <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.upper.error
        
        # alpha
        my.alpha.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$alpha
        my.alpha.lower.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$alpha.lower.error
        my.alpha.upper.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$alpha.upper.error
        
        alpha.pos <- which(alpha.values$focal == focal.sp[i.sp] &
                             alpha.values$model == my.models[i.model] &
                             alpha.values$method == optim.methods[i.method] &
                             alpha.values$competitor %in% names(my.alpha.vector))
        
        alpha.values$alpha[alpha.pos] <- my.alpha.vector
        alpha.values$alpha.lower[alpha.pos] <- my.alpha.lower.vector
        alpha.values$alpha.upper[alpha.pos] <- my.alpha.upper.vector
        
        # lambda.cov
        my.lambda.cov.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.cov
        my.lambda.cov.lower.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.cov.lower.error
        my.lambda.cov.upper.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.cov.upper.error
        
        # lambda.cov_NL
        my.lambda.cov_NL.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.cov_NL
        my.lambda.cov_NL.lower.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.cov_NL.lower.error
        my.lambda.cov_NL.upper.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$lambda.cov_NL.upper.error
        
        # all covariates at once, should be ok
        lambda.cov.pos <- which(lambda.cov.values$species == focal.sp[i.sp] &
                                  lambda.cov.values$model == my.models[i.model] &
                                  lambda.cov.values$method == optim.methods[i.method])
        
        lambda.cov.values$lambda.cov[lambda.cov.pos] <- my.lambda.cov.vector
        lambda.cov.values$lambda.cov.lower[lambda.cov.pos] <- my.lambda.cov.lower.vector
        lambda.cov.values$lambda.cov.upper[lambda.cov.pos] <- my.lambda.cov.upper.vector
        
        # alpha.cov
        my.alpha.cov.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$alpha.cov
        my.alpha.cov.lower.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$alpha.cov.lower.error
        my.alpha.cov.upper.vector <- param.matrices[[focal.sp[i.sp]]][[my.models[i.model]]][[optim.methods[i.method]]]$alpha.cov.upper.error
        
        my.alpha.cov.comp <- substr(names(my.alpha.cov.vector),str_length(names(my.alpha.cov.vector))-3,str_length(names(my.alpha.cov.vector)))
        
        for(i.covariate in 1:length(my.covariates)){
          my.cov <- my.alpha.cov.vector[which(grepl(my.covariates[i.covariate],names(my.alpha.cov.vector)))]
          my.lower.cov <- my.alpha.cov.lower.vector[which(grepl(my.covariates[i.covariate],names(my.alpha.cov.lower.vector)))]
          my.upper.cov <- my.alpha.cov.upper.vector[which(grepl(my.covariates[i.covariate],names(my.alpha.cov.upper.vector)))]
          
          alpha.cov.pos <- which(alpha.cov.values$focal == focal.sp[i.sp] &
                                   alpha.cov.values$model == my.models[i.model] &
                                   alpha.cov.values$method == optim.methods[i.method] &
                                   alpha.cov.values$competitor %in% my.alpha.cov.comp & 
                                   alpha.cov.values$covariate == my.covariates[i.covariate])
          
          alpha.cov.values$alpha.cov[alpha.cov.pos] <- my.cov
          alpha.cov.values$alpha.cov.lower[alpha.cov.pos] <- my.lower.cov
          alpha.cov.values$alpha.cov.upper[alpha.cov.pos] <- my.upper.cov
        }
        
      }# for each method
    }# for each model
  }# for each sp
  write.table(lambda.values, file = "./data/lambda_values.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
  write.csv(lambda.values,file = "./data/alpha_values.csv",sep = ";")
  write.csv(lambda.values,file = "./data/lambda_cov_values.csv",sep = ";")
  write.csv(lambda.values,file = "./data/alpha_cov_values.csv",sep = ";")
  
}
