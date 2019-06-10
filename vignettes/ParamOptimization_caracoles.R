source("R/BevertonHolt_models.R")
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
salinity <- readr::read_delim(file = "../Caracoles/data/salinity.csv",delim = ";")

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
models <- 3:5

# which values do we optimize for each model?
param.list <- list(c("lambda","alpha"),
                   c("lambda","alpha","lambda.cov","alpha.cov"),
                   c("lambda","alpha","lambda.cov","alpha.cov"))

#Choose the non-linearity function
#function_NL <- function1
function1 <-function(a,b,x){
  return(a*x^2/(b+x^2))
}
function2 <-function(a,b,x){
  return(a*(1-exp(b*x)))
}

# keep the model definitions in a list, for ease
fitness.models <- list(BH_1 = BH_1,BH_2 = BH_2,BH_3 = BH_3,BH_4 = BH_4,BH_5 = BH_5)

# environmental covariates
covariates <- full.data[,"sum_salinity"]
# if no covariates, comment above and uncomment here
# covariates <- 0

# optimization methods to use
optim.methods <- c(#"optim_NM",
                   "optim_L-BFGS-B",
                  "nloptr_CRS2_LM"
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
lower.alpha <- 0
upper.alpha <- 1e5
# lambda.cov
init.lambda.cov <- 1
lower.lambda.cov <- 0
upper.lambda.cov <- 1e4
# alpha.cov
init.alpha.cov <- 1
lower.alpha.cov <- 0
upper.alpha.cov <- 1e4

# if we want quicker calculations, we can disable 
# the bootstrapping for the standard errors
generate.errors <- TRUE
bootstrap.samples <- 3

# store results?
write.results <- FALSE

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
  if("alpha" %in% param.list[[i.model]]){
    if(models[i.model]<=2){
      alpha.length <- 1
    }else{
      alpha.length <- num.competitors
    }
    if(length(init.alpha) != alpha.length){
     current.init.alpha <- rep(init.alpha[1],alpha.length) 
    }else{
      current.init.alpha <- init.alpha
    }
  }else{
    current.init.alpha <- init.alpha[i.sp,]
  }
  # lambda.cov
  if("lambda.cov" %in% param.list[[i.model]]){
    if(length(init.lambda.cov) != num.covariates){
      current.init.lambda.cov <- rep(init.lambda.cov[1],num.covariates)
    }else{
      current.init.lambda.cov <- init.lambda.cov  
    }
  }else{
    current.init.lambda.cov <- init.lambda.cov[i.sp]  
  }
  # alpha.cov
  if("alpha.cov" %in% param.list[[i.model]]){
    if(models[i.model]<=4){
      length.alpha.cov <- num.covariates
    }else if(models[i.model]>4){
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
                                   focal.comp.matrix = focal.comp.matrix,
                                   focal.covariates = focal.covariates,
                                   generate.errors = generate.errors,
                                   bootstrap.samples = bootstrap.samples)
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
      
    }# for i.method
    
    #######################
    # update initial values for the different parameters
    
    # lambda
    if("lambda" %in% param.list[[i.model]]){
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
    if("alpha" %in% param.list[[i.model]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha)) == 0){
        current.init.alpha <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha
        # is the current estimate of the appropriate length?
        if(i.model > 2){
          if(length(current.init.alpha) == 1){
            current.init.alpha <- rep(current.init.alpha,num.competitors)
          }
        }# if model > 2
      }
    
    # lambda.cov
    if("lambda.cov" %in% param.list[[i.model]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov)) == 0){
        current.init.lambda.cov <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov
      }
    }
    
    # lambda.cov_NL
    if("lambda.cov_NL" %in% param.list[[i.model]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov_NL)) == 0){
        current.init.lambda.cov_NL <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$lambda.cov_NL
      }
    }
    
    # alpha.cov
    if("alpha.cov" %in% param.list[[i.model]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov)) == 0){
        current.init.alpha.cov <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov
        # is the current estimate of the appropriate length?
        if(i.model > 4){
          if(length(current.init.alpha.cov) == num.covariates){
            current.init.alpha.cov <- rep(current.init.alpha.cov,num.competitors)
          }
        }# if model > 4
      }
    }
    
    # alpha.cov_NL
    if("alpha.cov_NL" %in% param.list[[i.model]]){
      if(sum(is.na(param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov_NL)) == 0){
        current.init.alpha.cov_NL <- param.matrices[[i.sp]][[i.model]][[init.par.method]]$alpha.cov_NL
        # is the current estimate of the appropriate length?
        if(i.model > 4){
          if(length(current.init.alpha.cov_NL) == num.covariates){
            current.init.alpha.cov_NL <- rep(current.init.alpha.cov_NL,num.competitors)
          }
        }# if model > 4
      }
    }
    
  }# for i.model
}# for i.sp

if(write.results){
  save(param.matrices,file = "./results/param_estimates.Rdata")
}

