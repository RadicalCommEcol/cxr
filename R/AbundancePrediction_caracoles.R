######
library(tidyverse)

######
# timesteps
timesteps <- 50

######
source("./R/nested_models.R")
source("./R/PredictAbundances.R")

####################
competition.data <- readr::read_delim(file = "./data/competition.csv",delim = ";")

# spread the data from long to wide format
competition.data <- spread(competition.data,competitor,number,fill = 0)

# competition matrix
comp.matrix <- as.matrix(competition.data[,10:ncol(competition.data)])
colnames(comp.matrix) <- names(competition.data)[10:ncol(competition.data)]

# covariates
covariates <- 0

abundances <- readr::read_delim("./data/abundances.csv",delim = ";")
# complete missing sp-year combinations
abundances <- complete(abundances, year,plot,subplot,species, fill = list(individuals = 0, month = 0, day = 0, order = 0))

###########
# read lambda,s,g, and alpha values
lambda.values <- readr::read_delim(file = "./results/lambda_estimates.csv",delim = ";")
load("./results/param_estimates")

############ TODO:CHECK
# we will take estimates from the most complex model parameterized
max.model <- max(lambda.values$model)

# delete extremely unlikely lambda values, e.g. > 1000
# and stick with those values fitted with max.model
lambda.values <- subset(lambda.values, lambda < 1e3 & model == max.model)

# ideally we will have a complete parameterization of all species present
# otherwise, the dynamics will only predict those parameterized,
# and the effect of other competitors will dissappear from time 2 onwards
# UNLESS we decide to take as a constant the abundance of other, non-focal, competitors
# check with nacho and oscar

# in any case, for now, the projected species need to have:
# initial abundance data
focal.sp.abund <- unique(abundances$species)
# lambda values
focal.sp.lambda <- unique(lambda.values$focal.sp)
# and be present in the competition matrices
competition.sp <- unique(colnames(comp.matrix))

# the set of projected species is the intersection of these three sets
focal.species <- Reduce(intersect, list(focal.sp.abund,focal.sp.lambda,competition.sp))

# order projected focal species alphabetically
focal.species <- sort(focal.species)

# first take on lambda values: average of the different optim methods
sp.par <- lambda.values %>% filter(focal.sp %in% focal.species) %>% group_by(focal.sp,model) %>% summarise(lambda = mean(lambda))

# TODO:update with proper estimates
sp.par$germ.rate <- runif(nrow(sp.par),0,0.5)
sp.par$survival.rate <- runif(nrow(sp.par),0,0.5)

# for now, select focal species and only use them. thus, I need to subset the alpha matrices appropriately
# indexes of the focal species in the alpha matrices
focal.indices <- which(colnames(comp.matrix) %in% focal.species)

###########
num.sp <- length(focal.species)
num.cov <- 0

# gather the complete competition matrix from the fitted data
# TODO: check whether averaging across optim methods is appropriate

alpha.matrix <- matrix(0,nrow=num.sp,ncol=num.sp)
rownames(alpha.matrix) <- focal.species
colnames(alpha.matrix) <- focal.species

for(i.sp in 1:num.sp){
  my.sp.alpha <- rep(0,num.sp)
  for(i.method in 1:length(unique(lambda.values$optim.method))){
    temp.alpha <- param.matrices[[focal.species[i.sp]]][[max.model]][[i.method]]$alpha.matrix
    # subset to include only interactions with other focal species
    my.sp.alpha <- my.sp.alpha + temp.alpha[focal.indices]
  }
  # average over all optim methods
  my.sp.alpha <- my.sp.alpha/length(unique(lambda.values$optim.method))
  
  alpha.matrix[focal.species[i.sp],] <- my.sp.alpha
}

#####################

# initial abundances
init.abund <- abundances %>% filter(species %in% focal.species) %>% group_by(year,plot,subplot,species) %>% summarise(abundance = sum(individuals))
init.abund$site <- paste(init.abund$plot,"_",init.abund$subplot,sep="")
init.abund <- init.abund[,c("year","site","species","abundance")]

# init.abund <- expand.grid(1:num.obs,1:num.sp)
# names(init.abund) <- c("site","sp")
# init.abund$abundance <- rnorm(nrow(init.abund),100,50)

# environmental heterogeneity
if(num.cov > 0){
cov.time <- expand.grid(1:num.obs,1:timesteps,1:num.cov)
names(cov.time) <- c("site","timestep","covariate")
cov.time$value <- runif(nrow(cov.time),0,10)

# effect of covariates on alpha and lambda
lambda.cov.matrix <- matrix(unlist(lambda.cov.orig),nrow = num.sp)

alpha.cov.matrix <- list()
for(i.cov in 1:num.cov){
  alpha.cov.matrix[[i.cov]] <- matrix(nrow = num.sp,ncol = num.sp)
  for(i.sp in 1:num.sp){
    alpha.cov.matrix[[i.cov]][i.sp,] <- alpha.cov.orig[[i.sp]][,i.cov]
  }
}

}else{
  cov.time <- 0
  lambda.cov.matrix <- 0
  alpha.cov.matrix <- 0
}

###################
# repeat for each year
year.abund <- subset(init.abund, year == 2017)

par <- list(sp.par = sp.par, initial.values = year.abund, 
            covariates = cov.time, other.par = list(alpha.matrix = alpha.matrix, 
                                                    lambda.cov.matrix = lambda.cov.matrix, 
                                                    alpha.cov.matrix = alpha.cov.matrix))


abundance.model <- abund.fun.3
predicted.abundances <- PredictAbundances(par = par,timesteps = timesteps,abundance.model = abundance.model)
predicted.abundances$timestep <- as.factor(predicted.abundances$timestep)
predicted.abundances$site <- as.factor(predicted.abundances$site)
predicted.abundances$sp <- as.factor(predicted.abundances$sp)
# TODO:why NAs?
predicted.abundances$abundance[is.na(predicted.abundances$abundance)] <- 0

# some quick summarising for plotting
plot.data <- predicted.abundances %>% group_by(timestep,sp) %>% summarise(mean.abund = mean(abundance), sd.abund = sd(abundance))

abund.plot <- ggplot(plot.data,aes(x = timestep,y = mean.abund, group = sp)) + 
  geom_point(aes(color = sp)) + 
  geom_errorbar(aes(ymin = mean.abund - sd.abund, ymax = mean.abund + sd.abund))+
  # facet_wrap(site~.,ncol = 4)+
  # ylim(-10,10)+
  NULL
abund.plot

