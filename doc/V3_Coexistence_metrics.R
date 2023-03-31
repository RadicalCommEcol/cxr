## ----setup,echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
#note - we disable warnings and other output for the vignette, otherwise the functions are quite "chatty". 

## -----------------------------------------------------------------------------
library(cxr)
data("neigh_list")
data <- neigh_list
# keep only fitness and neighbours columns
for(i in 1:length(data)){
  data[[i]] <- data[[i]][,2:length(data[[i]])]
}

focal_column <- names(data)
model_family <- "RK" 
optimization_method <- "bobyqa" # we use a bounded method
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"
initial_values = list(lambda = 1,
                      alpha_intra = 0.1,
                      alpha_inter = 0.1)
lower_bounds = list(lambda = 0,
                    alpha_intra = 0.01,
                    alpha_inter = 0.01)
upper_bounds = list(lambda = 100,
                    alpha_intra = 1,
                    alpha_inter = 1)
fixed_terms <- NULL
bootstrap_samples <- 0

## -----------------------------------------------------------------------------
all.sp.fit <- cxr_pm_multifit(data = data,
                              model_family = model_family,
                              focal_column = focal_column,
                              optimization_method = optimization_method,
                              alpha_form = alpha_form,
                              lambda_cov_form = lambda_cov_form,
                              alpha_cov_form = alpha_cov_form,
                              initial_values = initial_values,
                              lower_bounds = lower_bounds,
                              upper_bounds = upper_bounds,
                              fixed_terms = fixed_terms,
                              bootstrap_samples = bootstrap_samples)


## -----------------------------------------------------------------------------
niche_overlap_all_pairs <- niche_overlap(cxr_multifit = all.sp.fit)
# as not all species combinations occur, 
# several interaction coefficients are NA
# which, in turn, return NA values for niche overlap
# so, remove them and check the first computed values
niche_overlap_all_pairs <- niche_overlap_all_pairs[complete.cases(niche_overlap_all_pairs),]
head(niche_overlap_all_pairs)

## -----------------------------------------------------------------------------
niche_overlap_all_pairs$MCT_niche_diff <- 1 - niche_overlap_all_pairs$niche_overlap_MCT
niche_overlap_all_pairs$SA_niche_diff <- 1 - niche_overlap_all_pairs$niche_overlap_SA

## -----------------------------------------------------------------------------
avg_fitness_diff_all_pairs <- avg_fitness_diff(cxr_multifit = all.sp.fit)

# as with niche overlap, remove NA values
avg_fitness_diff_all_pairs <- avg_fitness_diff_all_pairs[complete.cases(
  avg_fitness_diff_all_pairs),]

# average fitness ratio of sp1 over sp2
# if < 1, sp2 is the superior competitor, 
# and the average fitness difference is the inverse ratio,
# i.e. sp2 over sp1.
head(avg_fitness_diff_all_pairs)

## -----------------------------------------------------------------------------
competitive_ability_all_pairs <- competitive_ability(cxr_multifit = all.sp.fit)
competitive_ability_all_pairs <- competitive_ability_all_pairs[complete.cases(
  competitive_ability_all_pairs),]

head(competitive_ability_all_pairs)

## -----------------------------------------------------------------------------
data("neigh_list")

# For obtaining effect and responses, all species 
# need to have the same number of observations. 
# We selct 3 species that have >250 observations
names(neigh_list)
sapply(neigh_list,nrow)
# BEMA, HOMA, LEMA, SASO, have > 250 observations.
example_sp <- c(1,5,6) #corresponds to c("BEMA","HOMA","LEMA") 
n.obs <- 250
data <- neigh_list[example_sp]

# use a bounded optimization method
optimization_method <- "L-BFGS-B"

# no fixed terms, i.e. we fit all parameters
fixed_terms <- NULL

# according to a Ricker model (for consistency with previous examples)
model_family <- "RK"

# no standard error calculation in this example
bootstrap_samples <- 0

# keep only fitness and neighbours columns
# and subset to 'n.obs' rows
for(i in 1:length(data)){
  data[[i]] <- data[[i]][1:n.obs,c(2,example_sp+2)]#2:length(data[[i]])]
}

# set initial values and bounds
initial_values_er = list(lambda = 10, 
                         effect = 1, 
                         response = 1)
lower_bounds_er = list(lambda = 1, 
                       effect = 0.1, 
                       response = 0.1)
upper_bounds_er = list(lambda = 100, 
                       effect = 10, 
                       response = 10)

## -----------------------------------------------------------------------------
er.fit <- cxr_er_fit(data = data,
                          model_family = model_family,
                          optimization_method = optimization_method,
                          initial_values = initial_values_er,
                          lower_bounds = lower_bounds_er,
                          upper_bounds = upper_bounds_er,
                          fixed_terms = fixed_terms,
                          bootstrap_samples = bootstrap_samples)

## -----------------------------------------------------------------------------
spfitness <- species_fitness(er.fit)
spfitness

