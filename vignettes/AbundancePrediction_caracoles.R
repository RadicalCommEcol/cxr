######
library(tidyverse)

######
# timesteps
timesteps <- 4

######
source("./R/BH_abundance_5.R")
source("./R/PredictAbundances.R")

####################
# covariates
num.cov <- 1
##############
# TODO: change when salinity is part of the package
covariates <- readr::read_delim("../Caracoles/data/salinity.csv",delim = ";")
##############
covariates$site <- paste(covariates$plot,"_",covariates$subplot,sep="")
covariates <- covariates[,c("year","site","sum_salinity")]
covariates <- gather(covariates,key = "covariate",value = "value",-year,-site)

# years are the different timesteps; rescale them
names(covariates)[1] <- "timestep"
covariates$timestep <- covariates$timestep - 2014

# initial abundances
##############
# TODO: change when abundances is part of the package
##############
# abundances <- readr::read_delim("./data/abundances.csv",delim = ";")
# complete missing sp-year combinations
# abundances <- tidyr::complete(abundance, year,plot,subplot,species, fill = list(individuals = 0, month = 0, day = 0, order = 0))
abundances <- read.table(file = "../Caracoles/data/abundances.csv",header = T,sep = ";",stringsAsFactors = F)
###########
# read lambda,s,g, and alpha values
# load("./results/param_estimates.Rdata")
##############
# TODO: change when species_rates is part of the package
##############
# species_rates <- readr::read_delim("../Caracoles/raw_data/seed_germination_survival.txt",delim = "\t")
data(species_rates)
data(param_estimates)

# only species with germination/survival rates
sp.names <- sort(as.character(unique(species_rates$code)))

estimates.model <- "BH_5" #Beverton-holt model number 5
estimates.method <- "optim_L-BFGS-B"

# gather lambda from fitted data
species_rates$lambda <- 0
for(i.sp in 1:length(sp.names)){
  if(!is.null(param_estimates[[sp.names[i.sp]]])){
    species_rates$lambda[species_rates$code == sp.names[i.sp]] <- param_estimates[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$lambda
  }
}

# subset species set to these with valid s,g, and lambda estimates
species_rates <- subset(species_rates,lambda != 0)
# clean up the species rates dataframe
sp.par <- species_rates[,c("lambda","germination","seed.survival")]
names(sp.par) <- c("lambda","germ.rate","survival.rate")

# update sp.names
sp.names <- sort(as.character(unique(species_rates$code)))
num.sp <- length(sp.names)

# gather the complete competition matrix from the fitted data
alpha.matrix <- matrix(0,nrow=num.sp,ncol=num.sp)
rownames(alpha.matrix) <- sp.names
colnames(alpha.matrix) <- sp.names

# estimated.names <- names(param_estimates)
# valid.positions <- match(estimated.names,sp.names)

for(i.sp in 1:num.sp){
  # my.sp.alpha <- rep(0,num.sp)
  
  # including a temporary hack for valid sorting of species positions
  my.sp.alpha <- param_estimates[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$alpha
  
  my.alpha.pos <- my.sp.alpha[which(!is.na(match(names(my.sp.alpha),sp.names)))]
  
  alpha.matrix[sp.names[i.sp],] <- my.alpha.pos
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
      lambda.cov[i.sp,i.cov] <- param_estimates[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$lambda.cov[i.cov]
    }
  }
  
  alpha.cov <- list()
  for(i.cov in 1:num.cov){
    alpha.cov[[i.cov]] <- matrix(nrow = num.sp,ncol = num.sp)
    for(i.sp in 1:num.sp){
      alpha.cov[[i.cov]][i.sp,] <- param_estimates[[sp.names[i.sp]]][[estimates.model]][[estimates.method]]$alpha.cov[i.sp+(num.sp*(i.cov-1))]
    }
  }
  
}else{
  covariates <- 0
  lambda.cov <- 0
  alpha.cov <- 0
}

###################
# which year to take as a starting point
# year.abund <- init.abund
year.abund <- subset(init.abund, year == 2015)

par <- list(sp.par = sp.par, initial.values = year.abund, 
            covariates = covariates, other.par = list(alpha.matrix = alpha.matrix, 
                                                    lambda.cov.matrix = lambda.cov, 
                                                    alpha.cov.matrix = alpha.cov))


abundance.model <- BH_abundance_5
predicted.abundances <- PredictAbundances(par = par,timesteps = timesteps,abundance.model = abundance.model,return.seeds = F)
predicted.abundances$timestep <- as.factor(predicted.abundances$timestep)
predicted.abundances$site <- as.factor(predicted.abundances$site)
predicted.abundances$sp <- as.factor(predicted.abundances$sp)

# some quick summarising for plotting
# plot.data <- predicted.abundances %>% group_by(timestep,sp) %>% summarise(mean.abund = mean(abundance), sd.abund = sd(abundance))
abund.mean.plot <- ggplot(plot.data,aes(x = timestep,y = mean.abund, group = sp)) +
  geom_point(aes(color = sp)) +
  geom_line(aes(color = sp)) +
  # geom_errorbar(aes(ymin = mean.abund - sd.abund, ymax = mean.abund + sd.abund))+
  # facet_wrap(site~.,ncol = 4)+
   # ylim(-10,10)+
  NULL
abund.mean.plot

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

# observed and predicted abundances
obs.pred.plot <- ggplot(obs.pred,aes(x = individuals,y = predicted,group = species)) + 
  geom_point(aes(color = species))+ #shape = year)) + 
  geom_abline(slope = 1, color = "lightgrey") +
  facet_wrap(year~plot, scales = "free", ncol = 9) +
  # xlim(0,10) + ylim(0,10) +
  NULL
obs.pred.plot

# only predicted
pred.mean <- obs.pred %>% group_by(year,plot,species) %>% summarise(all.ind = sum(individuals),all.pred = sum(predicted))

pred.plot <- ggplot(pred.mean,aes(x = year,y = all.pred,group = species)) + 
  geom_point(aes(color = species))+ #shape = year)) + 
  geom_line(aes(color = species)) +
  
  # geom_abline(slope = 1, color = "lightgrey") +
  facet_grid(plot~., scales = "free_y") +
  # xlim(0,10) + ylim(0,10) +
  NULL
pred.plot

pred.mean %>% group_by(species) %>% summarise(mean.total = mean(all.pred))
