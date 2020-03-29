## ----setup,echo=FALSE----------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## ------------------------------------------------------------------------
library(cxr)
data("neigh_list") 

## ------------------------------------------------------------------------
my.sp <- "HOMA" 
# get data from the list
obs_homa <- neigh_list[[my.sp]]
# no need for ID column
obs_homa <- subset(obs_homa,select = -c(obs_ID))
# For each observation, we need the individual plant fitness and the number of neighbours per species (in columns).
head(obs_homa)

## ------------------------------------------------------------------------
#?cxr_pm_fit #check the help file for a description of the arguments
fit_homa <- cxr_pm_fit(data = obs_homa,
                       focal_column = my.sp,
                       model_family = "RK",
                       covariates = NULL,
                       optimization_method = "Nelder-Mead",
                       alpha_form = "pairwise",
                       lambda_cov_form = "none",
                       alpha_cov_form = "none",
                       initial_values = list(lambda = 1,
                                             alpha_intra = .1,
                                             alpha_inter = .1),
                       #not aplicable to this optimazation method
                       # lower_bounds = list(lambda = 0, 
                       #                     alpha_intra = 0,
                       #                     alpha_inter = 0),
                       # upper_bounds = list(lambda = 10,
                       #                     alpha_intra = 1,
                       #                     alpha_inter = 1),
                       fixed_terms = NULL,
                       bootstrap_samples = 3) # a low number for demonstration purposes, increase it for robust results.

## ------------------------------------------------------------------------
summary(fit_homa)

## ------------------------------------------------------------------------
names(fit_homa) #list of all available elements.

#reproduction success in the absence of neighbors
fit_homa$lambda
# intraspecific interaction
fit_homa$alpha_intra
# interspecific interactions
fit_homa$alpha_inter

## ------------------------------------------------------------------------
my.sp <- c("BEMA","CETE","LEMA")
obs_3sp <- neigh_list[my.sp]
# discard ID column
for(i in 1:length(obs_3sp)){
  obs_3sp[[i]] <- obs_3sp[[i]][,2:length(obs_3sp[[i]])]
}
# load covariates: salinity
data("salinity_list")
salinity <- salinity_list[my.sp]
# keep only salinity column
for(i in 1:length(salinity)){
  salinity[[i]] <- as.matrix(salinity[[i]][,2:length(salinity[[i]])])
  colnames(salinity[[i]]) <- "salinity"
}

## ------------------------------------------------------------------------
names(obs_3sp)
head(obs_3sp[[1]])
nrow(obs_3sp[[1]])

head(salinity[[1]])
nrow(salinity[[1]])

## ------------------------------------------------------------------------
fit_3sp <- cxr_pm_multifit(data = obs_3sp,
                           focal_column = my.sp,
                           model_family = "RK",
                           # here we use a bounded method for demonstration purposes
                           optimization_method = "L-BFGS-B", 
                           covariates = salinity,
                           alpha_form = "pairwise",
                           lambda_cov_form = "global", # effect of covariates over lambda
                           alpha_cov_form = "pairwise", # effect of covariates over alpha
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,
                                                 alpha_inter = 0.1,
                                                 lambda_cov = 0.1,
                                                 alpha_cov = 0.1),
                           lower_bounds = list(lambda = 0,
                                               alpha_intra = 0,
                                               alpha_inter = -1,
                                               lambda_cov = 0,
                                               alpha_cov = 0),
                           upper_bounds = list(lambda = 100,
                                               alpha_intra = 1,
                                               alpha_inter = 1,
                                               lambda_cov = 1,
                                               alpha_cov = 1),
                           bootstrap_samples = 3)

## ------------------------------------------------------------------------
summary(fit_3sp)

## ------------------------------------------------------------------------
fit_3sp$log_likelihood

## ------------------------------------------------------------------------
fit_3sp$lambda_cov
fit_3sp$alpha_cov

