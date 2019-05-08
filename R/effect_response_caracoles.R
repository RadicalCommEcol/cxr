
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

# assume that seed production does not change in a single year, so group observations
# in year x site records
sp.data <- competition.data %>% group_by(year,plot,subplot,focal,competitor) %>% summarise(seed = sum(seed),number = sum(number))
names(sp.data)[which(names(sp.data) == "seed")] <- "fitness"
sp.data$site <- paste(sp.data$year,sp.data$plot,sp.data$subplot,sep="_")
sp.data <- sp.data[,c("site","focal","fitness","competitor","number")]

# Initial parameter estimates
lambda.values <- readr::read_delim(file = "./results/lambda_estimates.csv",delim = ";")

# take estimates from the most complex model parameterized
max.model <- max(lambda.values$model)

# stick with those values fitted with max.model
lambda.values <- subset(lambda.values, model == max.model)

# select method that minimises the overall loglikelihood across all lambdas
loglik <- lambda.values %>% group_by(optim.method) %>% summarise(sum.loglik = sum(log.likelihood))
lambda.values <- subset(lambda.values, optim.method == loglik$optim.method[loglik$sum.loglik == min(loglik$sum.loglik)])
sigma <- mean(lambda.values$sigma)

# sanity check
lambda.values <- arrange(subset(lambda.values, focal.sp %in% sp.data$focal),focal.sp)

# initial estimates for r, e
r.values <- rep(1,nrow(lambda.values))
e.values <- rep(1,nrow(lambda.values))

r.lower.bound <- rep(0, times=nrow(lambda.values))
r.upper.bound <- rep(1e2, times=nrow(lambda.values))
e.lower.bound <- r.lower.bound
e.upper.bound <- r.upper.bound
sigma.lower.bound <- 0.0000000001
sigma.upper.bound <- 1

# only species with proper estimates
sp.data <- subset(sp.data, focal %in% lambda.values$focal.sp & competitor %in% lambda.values$focal.sp)

######################
# compute each method

full.results <- NULL

for(i.method in 1:length(optim.methods)){
  
  param.results <- ER_optimize(lambda.vector = lambda.values$lambda,
                               e.vector = e.values,
                               r.vector = r.values,
                               sigma = sigma,
                               e.lower.bound = e.lower.bound,
                               e.upper.bound = e.upper.bound,
                               r.lower.bound = r.lower.bound,
                               r.upper.bound = r.upper.bound,
                               sigma.lower.bound = sigma.lower.bound,
                               sigma.upper.bound = sigma.upper.bound,
                               effect.response.model = effect.response.model,
                               optim.method = optim.methods[i.method],
                               sp.data = sp.data,
                               optimize.lambda = optimize.lambda,
                               generate.errors = generate.errors,
                               bootstrap.samples = bootstrap.samples)
  
  param.results$optim.method <- optim.methods[i.method]
  
  full.results <- rbind(full.results, param.results)
}

if(write.results){
  readr::write_delim(full.results,"./results/effect_response_estimates.csv",delim = ";")
}
