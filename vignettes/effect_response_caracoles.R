
# find best fits for effect and response parameters, 
# using caracoles data
# and lambda estimates previously calculated

source("R/ER_optimize.R")
source("R/EffectResponse.R")
source("R/ER_SEbootstrap.R")
require(tidyverse)

###########################
# optimization methods to use
optim.methods <- c("optim_NM"
                   # "optim_L-BGFS-B",
                   # "nloptr_CRS2_LM", 
                   # "nloptr_ISRES", 
                   # "nloptr_DIRECT_L_RAND", 
                   # "GenSA"
                   # "hydroPSO", 
                   # "DEoptimR"
)

# if we want quick calculations, we can disable 
# the bootstrapping for the standard errors
generate.errors <- TRUE
bootstrap.samples <- 3

optimize.lambda <- FALSE
# model is different...
if(optimize.lambda){
  effect.response.model <- EffectResponse_lambda
}else{
  effect.response.model <- EffectResponse
}

write.results <- FALSE

###########################
# Caracoles data

competition.data <- readr::read_delim(file = "./data/competition.csv",delim = ";")
### TEST
competition.data <- subset(competition.data,year == 2016)

# covariate: salinity
covariates <- readr::read_delim(file = "../Caracoles/data/salinity.csv",delim = ";")

# one observation per row of competition.data
covariates <- covariates[,c("plot","subplot","year","sum_salinity")]
# in the same order as the competition observations
covariates$site <- paste(covariates$year,covariates$plot,covariates$subplot,sep="_")
covariates <- covariates[,c("site","sum_salinity")]

# assume that seed production does not change in a single year, so group observations
# in year x site records
sp.data <- competition.data %>% group_by(year,plot,subplot,focal,competitor) %>% summarise(seed = sum(seed),number = sum(number))
names(sp.data)[which(names(sp.data) == "seed")] <- "fitness"
sp.data$site <- paste(sp.data$year,sp.data$plot,sp.data$subplot,sep="_")
sp.data <- sp.data[,c("site","focal","fitness","competitor","number")]

# this is the set of species we are fitting
sp.names <- sort(unique(sp.data$focal))

# in case the sp.data dataframe does not include explicit missing competitors, 
# here is a somewhat convoluted way for setting zeros to it
# so that for each focal sp, all competitor sp are included
sites <- unique(sp.data$site)

missing.data <- sp.data
missing.data$site <- "0"
missing.data$focal <- "0"
missing.data$fitness <- 0
missing.data$competitor <- "0"
missing.data$number <- 0
count <- 1

for(i.site in 1:length(sites)){
  for(i.sp in 1:length(sp.names)){
    my.competitors <- unique(sp.data$competitor[sp.data$site == sites[i.site] & sp.data$focal == sp.names[i.sp]])
    if(length(my.competitors) > 0 & length(my.competitors) < length(sp.names)){
      my.fitness <- sp.data$fitness[sp.data$site == sites[i.site] & sp.data$focal == sp.names[i.sp]][1]
      
      missing.competitors <- sp.names[which(!sp.names %in% my.competitors)]
      for(i.com in 1:length(missing.competitors)){
        
        missing.data$site[count] <- sites[i.site]
        missing.data$focal[count] <- sp.names[i.sp]
        missing.data$fitness[count] <- my.fitness
        missing.data$competitor[count] <- missing.competitors[i.com]
        missing.data$number[count] <- 0
        
        count <- count + 1
      }# for each missing
    }# if any missing
    # if(length(my.competitors) > 0 & length(my.competitors) < 19){print(paste(i.site,",",i.sp))}
  }# for i.sp
}# for i.site

missing.data <- droplevels(subset(missing.data,competitor != "0"))

sp.data <- rbind(sp.data,missing.data)
sp.data <- arrange(sp.data, focal, site, competitor)

# discard focal sp with fitness 0
sp.data <- droplevels(subset(sp.data, fitness > 0))

# initial parameter estimates
# read lambda values
load("./results/param_estimates.Rdata")
lambda.values <- data.frame(sp = sp.names,lambda = 0, sigma = 0)

# from which model and optimization method are we taking estimates?
estimates.model <- "BH_5" #Beverton-holt model number 5
estimates.method <- "optim_NM"

# gather lambda from fitted data
for(i.sp in 1:length(sp.names)){
  if(!is.null(param.matrices[[sp.names[i.sp]]])){
    lambda.values$lambda[lambda.values$sp == sp.names[i.sp]] <- param.matrices[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$lambda
    lambda.values$sigma[lambda.values$sp == sp.names[i.sp]] <- param.matrices[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$sigma
  }
}

# sanity check
lambda.values <- arrange(subset(lambda.values, sp %in% sp.data$focal),sp.names)
# sigma is also a parameter
sigma <- mean(lambda.values$sigma)

# only species with proper estimates
sp.data <- subset(sp.data, focal %in% lambda.values$sp & competitor %in% lambda.values$sp)

# sort covariates
positions <- match(sp.data$site,covariates$site)
covariates <- covariates[positions,2:ncol(covariates)]

# initial estimates for r, e
r.values <- rep(1,nrow(lambda.values))
e.values <- rep(1,nrow(lambda.values))

r.lower.bound <- rep(0, times=nrow(lambda.values))
r.upper.bound <- rep(1e2, times=nrow(lambda.values))
e.lower.bound <- r.lower.bound
e.upper.bound <- r.upper.bound
sigma.lower.bound <- 0.0000000001
sigma.upper.bound <- 1

# initial values for lambda.cov, r.cov, e.cov
lambda.cov.values <- matrix(1,nrow = nrow(lambda.values),ncol = ncol(covariates))
e.cov.values <- matrix(1,nrow = nrow(lambda.values),ncol = ncol(covariates))
r.cov.values <- matrix(1,nrow = nrow(lambda.values),ncol = ncol(covariates))

lambda.cov.lower.bound <- matrix(0,nrow = nrow(lambda.values),ncol = ncol(covariates))
lambda.cov.upper.bound <- matrix(1e4,nrow = nrow(lambda.values),ncol = ncol(covariates))
e.cov.lower.bound <- matrix(0,nrow = nrow(lambda.values),ncol = ncol(covariates))
e.cov.upper.bound <- matrix(1e4,nrow = nrow(lambda.values),ncol = ncol(covariates))
r.cov.lower.bound <- matrix(0,nrow = nrow(lambda.values),ncol = ncol(covariates))
r.cov.upper.bound <- matrix(1e4,nrow = nrow(lambda.values),ncol = ncol(covariates))

######################
# build results list and compute each method

param.list <- list()
for(i.method in 1:length(optim.methods)){
  param.list[[i.method]] <- list(lambda = NA,
                                 lambda.lower.error = NA,
                                 lambda.upper.error = NA,
                                 response = NA,
                                 response.lower.error = NA,
                                 response.upper.error = NA,
                                 effect = NA,
                                 effect.lower.error = NA,
                                 effect.upper.error = NA,
                                 sigma = NA,
                                 lambda.cov = NA,
                                 lambda.cov.lower.error = NA,
                                 lambda.cov.upper.error = NA,
                                 response.cov = NA,
                                 response.cov.lower.error = NA,
                                 response.cov.upper.error = NA,
                                 effect.cov = NA,
                                 effect.cov.lower.error = NA,
                                 effect.cov.upper.error = NA,
                                 log.likelihood = NA)

}
names(param.list) <- optim.methods

########## TEST
# i.method <- 1
##########

for(i.method in 1:length(optim.methods)){
  
  param.results <- ER_optimize(lambda.vector = lambda.values$lambda,
                               e.vector = e.values,
                               r.vector = r.values,
                               lambda.cov = lambda.cov.values,
                               e.cov = e.cov.values,
                               r.cov = r.cov.values,
                               sigma = sigma,
                               lambda.lower.bound = lambda.lower.bound,
                               lambda.upper.bound = lambda.upper.bound,
                               e.lower.bound = e.lower.bound,
                               e.upper.bound = e.upper.bound,
                               r.lower.bound = r.lower.bound,
                               r.upper.bound = r.upper.bound,
                               lambda.cov.lower.bound = lambda.cov.lower.bound,
                               lambda.cov.upper.bound = lambda.cov.upper.bound,
                               e.cov.lower.bound = e.cov.lower.bound,
                               e.cov.upper.bound = e.cov.upper.bound,
                               r.cov.lower.bound = r.cov.lower.bound,
                               r.cov.upper.bound = r.cov.upper.bound,
                               sigma.lower.bound = sigma.lower.bound,
                               sigma.upper.bound = sigma.upper.bound,
                               effect.response.model = effect.response.model,
                               optim.method = optim.methods[i.method],
                               sp.data = sp.data,
                               covariates = covariates,
                               optimize.lambda = optimize.lambda,
                               generate.errors = generate.errors,
                               bootstrap.samples = bootstrap.samples)
  
  # store results
  param.list[[optim.methods[i.method]]]$lambda <- param.results$lambda
  param.list[[optim.methods[i.method]]]$lambda.lower.error <- param.results$lambda.lower.error
  param.list[[optim.methods[i.method]]]$lambda.upper.error <- param.results$lambda.upper.error
  param.list[[optim.methods[i.method]]]$response <- param.results$response
  param.list[[optim.methods[i.method]]]$response.lower.error <- param.results$response.lower.error
  param.list[[optim.methods[i.method]]]$response.upper.error <- param.results$response.upper.error
  param.list[[optim.methods[i.method]]]$effect <- param.results$effect
  param.list[[optim.methods[i.method]]]$effect.lower.error <- param.results$effect.lower.error
  param.list[[optim.methods[i.method]]]$effect.upper.error <- param.results$effect.upper.error
  param.list[[optim.methods[i.method]]]$sigma <- param.results$sigma
  param.list[[optim.methods[i.method]]]$lambda.cov <- param.results$lambda.cov
  param.list[[optim.methods[i.method]]]$lambda.cov.lower.error <- param.results$lambda.cov.lower.error
  param.list[[optim.methods[i.method]]]$lambda.cov.upper.error <- param.results$lambda.cov.upper.error
  param.list[[optim.methods[i.method]]]$response.cov <- param.results$response.cov
  param.list[[optim.methods[i.method]]]$response.cov.lower.error <- param.results$response.cov.lower.error
  param.list[[optim.methods[i.method]]]$response.cov.upper.error <- param.results$response.cov.upper.error
  param.list[[optim.methods[i.method]]]$effect.cov <- param.results$effect.cov
  param.list[[optim.methods[i.method]]]$effect.cov.lower.error <- param.results$effect.cov.lower.error
  param.list[[optim.methods[i.method]]]$effect.cov.upper.error <- param.results$effect.cov.upper.error
  param.list[[optim.methods[i.method]]]$log.likelihood <- param.results$log.likelihood
  
}

if(write.results){
  save(param.matrices,file = "./results/effect_response_estimates.Rdata")
}
