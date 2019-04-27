source("R/BevertonHolt_models.R")
source("R/SEbootstrap_caracoles.R")
source("R/cxr_optimize.R")
source("R/InitParams.R")
# optimization methods
library(nloptr)
library(GenSA)
library(hydroPSO)
library(DEoptimR)
#other packages
library(tidyr)

###########################
# Caracoles competition data

competition.data <- readr::read_delim(file = "./data/competition.csv",delim = ";")

# spread the data from long to wide format
competition.data <- spread(competition.data,competitor,number,fill = 0)

# competition matrix
comp.matrix <- as.matrix(competition.data[,10:ncol(competition.data)])

# covariate: salinity
covariates <- readr::read_delim(file = "../Caracoles/data/salinity.csv",delim = ";")

# one observation per row of competition.data
covariates <- covariates[,c("plot","subplot","year","sum_salinity")]
full.data <- left_join(competition.data,covariates)
covariates <- full.data[,"sum_salinity"]

# if no covariates, comment above and uncomment here
# covariates <- 0

#############################
# simulation parameters

# how many focal species
focal.sp <- unique(competition.data$focal)

# models to parameterize
# be aware of including 4 and 5 ONLY if there are covariates
# otherwise it makes no sense (see equations in Lanuza et al. 2018)
models <- 1:5
# keep the model definitions in a list, for ease
fitness.models <- list(BH_1,BH_2,BH_3,BH_4,BH_5)

# optimization methods to use
optim.methods <- c("optim_NM"#,
                   # "optim_L-BFGS-B",
                   # "nloptr_CRS2_LM", 
                   # "nloptr_ISRES", 
                   # "nloptr_DIRECT_L_RAND", 
                   # "GenSA", 
                   # "hydroPSO", 
                   # "DEoptimR"
                   )
# from which method are we taking initial estimates for the next model?
# early observations suggest optim_NM or DEoptimR
init.par.method <- optim.methods[1]
init.method.num <- which(optim.methods == init.par.method)

# if we want quick calculations, we can disable 
# the bootstrapping for the standard errors
generate.errors <- TRUE
bootstrap.samples <- 2

###
write.results <- FALSE

##############################
# initialize data structures
# lambdas and sigmas will be placed in the dataframe "lambda.results"
# others are matrices/vectors within nested lists, of the form matrix[[focal.sp]][[model]][[method]]

lambda.results <- NULL

# initialize the rest
param.matrices <- list()
for(i.sp in 1:length(focal.sp)){
  param.matrices[[i.sp]] <- list()
  for(i.model in 1:length(models)){
    param.matrices[[i.sp]][[i.model]] <- list()
    for(i.method in 1:length(optim.methods)){
      param.matrices[[i.sp]][[i.model]][[i.method]] <- list(alpha.matrix = 0,
                                                             alpha.lower.error.matrix = 0,
                                                             alpha.upper.error.matrix = 0,
                                                             lambda.cov.matrix = 0,
                                                             lambda.cov.lower.error.matrix = 0,
                                                             lambda.cov.upper.error.matrix = 0,
                                                             alpha.cov.matrix = 0,
                                                             alpha.cov.lower.error.matrix = 0,
                                                             alpha.cov.upper.error.matrix = 0)
    }
  }
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
  focal.covariates <- ifelse(num.covariates > 0,covariates[which(competition.data$focal == focal.sp[i.sp]),,drop = FALSE],0)
  # check
  focal.covariates <- as.data.frame(focal.covariates)
  
  # model to optimize  
  for(i.model in models){
    
    print("*********************************")
    print(paste(date()," - starting focal sp ",focal.sp[i.sp],", model ",i.model,sep=""))
    print("*********************************")
    #message()
    
    # initial parameters are gathered in a separate function
    # this function has default values for the upper and lower bounds,
    # modify if necessary
    init.param.list <- InitParams(model.index = models[i.model],
                                  log.fitness = log.fitness,
                                  lambda.results = lambda.results,
                                  num.competitors = num.competitors,
                                  num.covariates = num.covariates,
                                  init.par.method = init.par.method,
                                  focal.sp = focal.sp[i.sp],
                                  param.matrices = param.matrices)

    ######################
    # compute each method
    
    for(i.method in 1:length(optim.methods)){
      
      temp.results <- cxr_optimize(init.par = init.param.list$init.par,
                                   lower.bounds = init.param.list$lower.bounds,
                                   upper.bounds = init.param.list$upper.bounds,
                                   fitness.model = fitness.models[[i.model]],
                                   optim.method = optim.methods[i.method],
                                   log.fitness = log.fitness,
                                   focal.comp.matrix = focal.comp.matrix,
                                   focal.covariates = focal.covariates,
                                   generate.errors = generate.errors,
                                   bootstrap.samples = bootstrap.samples)
      
      ###############
      # clean up results
      
      temp.lambda <- temp.results$lambda.results
      temp.lambda$focal.sp <- focal.sp[i.sp]
      temp.lambda$model <- models[i.model]
      temp.lambda$optim.method <- optim.methods[i.method]
      
      temp.lambda <- temp.lambda[,c("focal.sp",
                                    "model",
                                    "optim.method",
                                    "lambda",
                                    "lambda.lower.error",
                                    "lambda.upper.error",
                                    "sigma",
                                    "log.likelihood")]

      lambda.results <- rbind(lambda.results,temp.lambda)
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.matrix <- temp.results$alpha.matrix
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.upper.error.matrix <- temp.results$alpha.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.lower.error.matrix <- temp.results$alpha.lower.error
      
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov.matrix <- temp.results$lambda.cov.matrix
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov.upper.error.matrix <- temp.results$lambda.cov.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$lambda.cov.lower.error.matrix <- temp.results$lambda.cov.lower.error
      
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov.matrix <- temp.results$alpha.cov.matrix
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov.upper.error.matrix <- temp.results$alpha.cov.upper.error
      param.matrices[[i.sp]][[i.model]][[i.method]]$alpha.cov.lower.error.matrix <- temp.results$alpha.cov.lower.error
      
    }# for i.method
  }# for i.model
}# for i.sp

if(write.results){
  readr::write_delim(lambda.results,"./results/lambda_estimates.csv",delim = ";")
  save(param.matrices,file = "./results/param_estimates.Rdata")
}

