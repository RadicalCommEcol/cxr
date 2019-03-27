
# find best fits for effect and response parameters, 
# from caracoles data

source("R/effect_response_model.R")
require(tidyverse)

###########################
# Caracoles data

competition.data <- readr::read_delim(file = "./data/competition.csv",delim = ";")

# assume that seed production does not change in a single year, so group observations
# in year x site records
sp.data <- competition.data %>% group_by(year,plot,subplot,focal,competitor) %>% summarise(seed = sum(seed),number = sum(number))
names(sp.data)[which(names(sp.data) == "seed")] <- "fitness"
sp.data$site <- paste(sp.data$year,sp.data$plot,sp.data$subplot,sep="_")
sp.data <- sp.data[,c("site","focal","fitness","competitor","number")]

# lambda initial estimates
lambda.values <- readr::read_delim(file = "./results/lambda_estimates.csv",delim = ";")

# we will take estimates from the most complex model parameterized
max.model <- max(lambda.values$model)

# delete extremely unlikely lambda values, e.g. > 1000
# and stick with those values fitted with max.model
lambda.values <- subset(lambda.values, lambda < 1e3 & model == max.model)

sigma <- as.numeric(lambda.values %>% summarise(mean(sigma)))

# beware that these values should be sorted alphabetically
# so that BEMA = 1, CETE = 2,...
lambda.values <- lambda.values %>% group_by(focal.sp) %>% summarise(lambda = mean(lambda))

r.values <- rep(1,nrow(lambda.values))
e.values <- rep(1,nrow(lambda.values))

############
# only species with proper estimates
# TODO: group species without lambda/e in a general category?
sp.data <- subset(sp.data, focal %in% lambda.values$focal.sp & competitor %in% lambda.values$focal.sp)

# set zeros to missing data
# so that for each focal sp, all competitors are included
sites <- unique(sp.data$site)
focal.sp <- sort(unique(sp.data$focal))
missing.data <- sp.data
missing.data$site <- "0"
missing.data$focal <- "0"
missing.data$fitness <- 0
missing.data$competitor <- "0"
missing.data$number <- 0
count <- 1

for(i.site in 1:length(sites)){
  for(i.focal in 1:length(focal.sp)){
    my.competitors <- unique(sp.data$competitor[sp.data$site == sites[i.site] & sp.data$focal == focal.sp[i.focal]])
    if(length(my.competitors) > 0 & length(my.competitors) < length(focal.sp)){
      my.fitness <- sp.data$fitness[sp.data$site == sites[i.site] & sp.data$focal == focal.sp[i.focal]][1]
      
      missing.competitors <- focal.sp[which(!focal.sp %in% my.competitors)]
      for(i.com in 1:length(missing.competitors)){
      
      missing.data$site[count] <- sites[i.site]
      missing.data$focal[count] <- focal.sp[i.focal]
      missing.data$fitness[count] <- my.fitness
      missing.data$competitor[count] <- missing.competitors[i.com]
      missing.data$number[count] <- 0
      
      count <- count + 1
      }# for each missing
    }# if any missing
    # if(length(my.competitors) > 0 & length(my.competitors) < 19){print(paste(i.site,",",i.focal))}
  }# for i.focal
}# for i.site

missing.data <- droplevels(subset(missing.data,competitor != "0"))

sp.data <- rbind(sp.data,missing.data)
sp.data <- arrange(sp.data, focal, site, competitor)

# discard focal sp with fitness 0
sp.data <- droplevels(subset(sp.data, fitness > 0))

############
par <- c(lambda.values$lambda,r.values,e.values, sigma)

############
# optim

optim.par <- optim(par, 
                   effect_response_model, 
                   gr = NULL, 
                   method = "Nelder-Mead", 
                   # lower = lower.bounds,
                   # upper = upper.bounds,
                   control = list(), 
                   hessian = F,
                   sp.data = sp.data)
















