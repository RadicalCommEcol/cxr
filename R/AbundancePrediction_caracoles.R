######
library(tidyverse)

######
# timesteps
timesteps <- 4

######
source("./R/BevertonHolt_abundance_models.R")
source("./R/PredictAbundances.R")

####################
# covariates
num.cov <- 1
covariates <- readr::read_delim("../Caracoles/data/salinity.csv",delim = ";")
covariates$site <- paste(covariates$plot,"_",covariates$subplot,sep="")
covariates <- covariates[,c("year","site","sum_salinity")]
covariates <- gather(covariates,key = "covariate",value = "value",-year,-site)

# years are the different timesteps; rescale them
names(covariates)[1] <- "timestep"
covariates$timestep <- covariates$timestep - 2014

# initial abundances
abundances <- readr::read_delim("./data/abundances.csv",delim = ";")
# complete missing sp-year combinations
abundances <- complete(abundances, year,plot,subplot,species, fill = list(individuals = 0, month = 0, day = 0, order = 0))

###########
# read lambda,s,g, and alpha values
load("./results/param_estimates.Rdata")
species.rates <- readr::read_delim("../Caracoles/raw_data/seed_germination_survival.txt",delim = "\t")

# only species with germ/survival rates
sp.names <- sort(unique(species.rates$code))

estimates.model <- "BH_5" #Beverton-holt model number 5
estimates.method <- "optim_NM"

# gather lambda from fitted data
species.rates$lambda <- 0
for(i.sp in 1:length(sp.names)){
  if(!is.null(param.matrices[[sp.names[i.sp]]])){
    species.rates$lambda[species.rates$code == sp.names[i.sp]] <- param.matrices[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$lambda
  }
}

# subset species set to these with valid s,g, and lambda estimates
species.rates <- subset(species.rates,lambda != 0)
# clean up the species rates dataframe
sp.par <- species.rates[,c("lambda","germination","seed survival")]
names(sp.par) <- c("lambda","germ.rate","survival.rate")

# update sp.names
sp.names <- sort(unique(species.rates$code))
num.sp <- length(sp.names)

# gather the complete competition matrix from the fitted data
alpha.matrix <- matrix(0,nrow=num.sp,ncol=num.sp)
rownames(alpha.matrix) <- sp.names
colnames(alpha.matrix) <- sp.names

estimated.names <- names(param.matrices)
valid.positions <- match(estimated.names,sp.names)

for(i.sp in 1:num.sp){
  my.sp.alpha <- rep(0,num.sp)
  
  # including a temporary hack for valid sorting of species positions
  my.sp.alpha <- param.matrices[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$alpha[valid.positions]
  my.sp.alpha <- my.sp.alpha[!is.na(my.sp.alpha)]
  
  alpha.matrix[sp.names[i.sp],] <- my.sp.alpha
}

#####################

# initial abundances
init.abund <- abundances %>% filter(species %in% sp.names) %>% group_by(year,plot,subplot,species) %>% summarise(abundance = sum(individuals))
init.abund$site <- paste(init.abund$plot,"_",init.abund$subplot,sep="")
init.abund <- init.abund[,c("year","site","species","abundance")]

# environmental heterogeneity
if(num.cov > 0){
  
  lambda.cov <- matrix(0,nrow = num.sp,ncol = num.cov)
  # lambda.cov
  for(i.sp in 1:num.sp){
    for(i.cov in 1:num.cov){
      lambda.cov[i.sp,i.cov] <- param.matrices[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$lambda.cov[i.cov]
    }
  }
  
  alpha.cov <- list()
  for(i.cov in 1:num.cov){
    alpha.cov[[i.cov]] <- matrix(nrow = num.sp,ncol = num.sp)
    for(i.sp in 1:num.sp){
      alpha.cov[[i.cov]][i.sp,] <- param.matrices[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$alpha.cov[i.sp+(num.sp*(i.cov-1))]
    }
  }
  
}else{
  covariates <- 0
  lambda.cov <- 0
  alpha.cov <- 0
}

###################
# which year to take as a starting point
year.abund <- subset(init.abund, year == 2015)

par <- list(sp.par = sp.par, initial.values = year.abund, 
            covariates = covariates, other.par = list(alpha.matrix = alpha.matrix, 
                                                    lambda.cov.matrix = lambda.cov, 
                                                    alpha.cov.matrix = alpha.cov))


abundance.model <- BH_abundance_5
predicted.abundances <- PredictAbundances(par = par,timesteps = timesteps,abundance.model = abundance.model)
predicted.abundances$timestep <- as.factor(predicted.abundances$timestep)
predicted.abundances$site <- as.factor(predicted.abundances$site)
predicted.abundances$sp <- as.factor(predicted.abundances$sp)

# some quick summarising for plotting
plot.data <- predicted.abundances %>% group_by(timestep,sp) %>% summarise(mean.abund = mean(abundance), sd.abund = sd(abundance))

abund.mean.plot <- ggplot(plot.data,aes(x = timestep,y = mean.abund, group = sp)) + 
  geom_point(aes(color = sp)) + 
  geom_line(aes(color = sp)) +
  # geom_errorbar(aes(ymin = mean.abund - sd.abund, ymax = mean.abund + sd.abund))+
  # facet_wrap(site~.,ncol = 4)+
  # ylim(-10,10)+
  NULL
# abund.mean.plot

# prepare data for comparing observed vs predicted
obs.pred <- abundances[,c("year","plot","subplot","species","individuals")]
# timestep to real years
predicted.abundances$timestep <- as.numeric(as.character(predicted.abundances$timestep)) + 2015
# get back plots and subplots
predicted.abundances$plot <- as.numeric(substr(predicted.abundances$site,1,1))
predicted.abundances$subplot <- substr(predicted.abundances$site,3,4)
predicted.abundances <- predicted.abundances[,c("timestep","plot","subplot","sp","abundance")]
names(predicted.abundances) <- c("year","plot","subplot","species","predicted")
obs.pred <- left_join(obs.pred,predicted.abundances)
obs.pred <- subset(obs.pred,year > 2015)
obs.pred <- obs.pred[which(!is.na(obs.pred$predicted)),]
obs.pred$year <- as.factor(obs.pred$year)

obs.pred.plot <- ggplot(obs.pred,aes(x = individuals,y = predicted,group = species)) + 
  geom_point(aes(color = species, shape = year)) + 
  geom_abline(slope = 1) +
  facet_wrap(plot~., scales = "free_y") +
  xlim(0,100) + ylim(0,100) +
  NULL
# obs.pred.plot


