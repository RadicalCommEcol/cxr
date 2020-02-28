library(cxr)
source("R/cxr_check_initial_values.R")
source("R/cxr_check_input_data.R")
source("R/cxr_check_method_boundaries.R")
source("R/cxr_check_pm_input.R")
source("R/cxr_get_init_params.R")
source("R/cxr_get_model_bounds.R")
source("R/cxr_pm_bootstrap.R")
source("R/cxr_retrieve_params.R")
source("R/cxr_return_init_length.R")
source("R/cxr_sort_params.R")
source("R/pm_BH_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("R/pm_BH_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("R/cxr_pm_fit.R")
source("R/cxr_pm_multifit.R")


 data("neigh_list")
 data <- neigh_list[1:3]
 # keep only fitness and neighbours columns
 for(i in 1:length(data)){
   data[[i]] <- data[[i]][,2:length(data[[i]])]
 }
 focal_column <- names(data)
 data("salinity_list")
salinity <- salinity_list[1:3]
 # keep only salinity column
for(i in 1:length(salinity)){
  salinity[[i]] <- as.matrix(salinity[[i]][,2:length(salinity[[i]])])
  colnames(salinity[[i]]) <- "salinity"
}
model_family <- "BH"
covariates <- salinity
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "pairwise"
initial_values = list(lambda = 1,
                      alpha_intra = 0.1,
                      alpha_inter = 0.1,
                      lambda_cov = 0.1,
                      alpha_cov = 0.1)
lower_bounds = list(lambda = 0,
                      alpha_intra = 0,
                      alpha_inter = -1,
                      lambda_cov = -1,
                      alpha_cov = -1)
upper_bounds = list(lambda = 100,
                      alpha_intra = 1,
                      alpha_inter = 1,
                      lambda_cov = 1,
                      alpha_cov = 1)
fixed_terms <- NULL
bootstrap_samples <- 3

ttt <- cxr_pm_multifit(data = data,
                      focal_column = focal_column,
                      model_family = model_family,
                      covariates = covariates,
                      optimization_method = optimization_method,
                      alpha_form = alpha_form,
                      lambda_cov_form = lambda_cov_form,
                      alpha_cov_form = alpha_cov_form,
                      initial_values = initial_values,
                      lower_bounds = lower_bounds,
                      upper_bounds = upper_bounds,
                      fixed_terms = fixed_terms,
                      bootstrap_samples = bootstrap_samples)
