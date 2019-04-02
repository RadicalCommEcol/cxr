
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
                   # "GenSA", 
                   # "hydroPSO", 
                   # "DEoptimR"
)

# if we want quick calculations, we can disable 
# the bootstrapping for the standard errors
generate.errors <- TRUE
bootstrap.samples <- 10

###
write.results <- FALSE

###########################
# Caracoles data

competition.data <- readr::read_delim(file = "./data/competition.csv",delim = ";")

# subset for quick tests only!
competition.data <- subset(competition.data, year == 2018)

# assume that seed production does not change in a single year, so group observations
# in year x site records
sp.data <- competition.data %>% group_by(year,plot,subplot,focal,competitor) %>% summarise(seed = sum(seed),number = sum(number))
names(sp.data)[which(names(sp.data) == "seed")] <- "fitness"
sp.data$site <- paste(sp.data$year,sp.data$plot,sp.data$subplot,sep="_")
sp.data <- sp.data[,c("site","focal","fitness","competitor","number")]

###########################
# Initial parameter estimates
lambda.values <- readr::read_delim(file = "./results/lambda_estimates.csv",delim = ";")

# we will take estimates from the most complex model parameterized
max.model <- max(lambda.values$model)

# delete extremely unlikely lambda values, e.g. > 1000
# and stick with those values fitted with max.model
lambda.values <- subset(lambda.values, lambda < 1e3 & model == max.model)

############## TEST
lambda.values <- subset(lambda.values, focal.sp %in% sp.data$focal)
##############

sigma <- as.numeric(lambda.values %>% summarise(mean(sigma)))

# beware that these values should be sorted alphabetically
# so that BEMA = 1, CETE = 2,...
lambda.values <- lambda.values %>% group_by(focal.sp) %>% summarise(lambda = mean(lambda))

r.values <- rep(1,nrow(lambda.values))
e.values <- rep(1,nrow(lambda.values))

init.par <- c(lambda.values$lambda,r.values,e.values, sigma)
lower.bounds <- c(rep(1, times=nrow(lambda.values)), rep(0, times=nrow(lambda.values)),rep(0, times=nrow(lambda.values)),0.0000000001)
upper.bounds <- c(rep(1e4, times=nrow(lambda.values)), rep(1e2, times=nrow(lambda.values)),rep(1e2, times=nrow(lambda.values)),1)

############
# only species with proper estimates
# TODO: group species without lambda/e in a general category?
sp.data <- subset(sp.data, focal %in% lambda.values$focal.sp & competitor %in% lambda.values$focal.sp)

######################
# compute each method

full.results <- NULL

for(i.method in 1:length(optim.methods)){
  
  param.results <- ER_optimize(init.par = init.par,
                               lower.bounds = lower.bounds,
                               upper.bounds = upper.bounds,
                               effect.response.model = EffectResponse,
                               optim.method = optim.methods[i.method],
                               sp.data = sp.data,
                               generate.errors = generate.errors,
                               bootstrap.samples = bootstrap.samples)
  
  param.results$optim.method <- optim.methods[i.method]
  
  full.results <- rbind(full.results, param.results)
}

if(write.results){
  readr::write_delim(full.results,"./results/effect_response_estimates.csv",delim = ";")
}
