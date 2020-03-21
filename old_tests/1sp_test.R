library(cxr)
source("R/cxr_check_pm_input.R")
source("R/cxr_get_init_params.R")
source("R/cxr_get_model_bounds.R")
source("R/cxr_sort_params.R")
source("R/cxr_return_init_length.R")

source("R/pm_LV_alpha_pairwise_lambdacov_global_alphacov_pairwise.R")
source("R/pm_LV_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("R/pm_LV_alpha_pairwise_lambdacov_global_alphacov_global.R")

source("R/cxr_retrieve_params.R")
source("R/cxr_pm_bootstrap.R")

source("R/cxr_check_initial_values.R")
source("R/cxr_check_input_data.R")
source("R/cxr_check_method_boundaries.R")

source("R/cxr_pm_fit.R")
source("R/summary.R")

data("neigh_list")

model_family = "LV"
# model_family = "RK"
# 
# model_family = "LW"

my.sp <- "BEMA"
data <- neigh_list[[my.sp]][2:ncol(neigh_list[[1]])]
focal_column = my.sp

data("salinity_list")
salinity <- salinity_list[[my.sp]][2]
covariates <- salinity
optimization_method <- "bobyqa"
alpha_form <- "pairwise" # "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "pairwise"
initial_values = list(lambda = 1,
                      alpha_intra = 0,
                      alpha_inter = 0,
                      lambda_cov = -0.1,
                      alpha_cov = -0.1)
lower_bounds = list(lambda = 0,
                      alpha_intra = -1,
                      alpha_inter = -1,
                      lambda_cov = -1,
                      alpha_cov = -1)
upper_bounds = list(lambda = 100,
                      alpha_intra = 0,
                      alpha_inter = 0,
                      lambda_cov = 0,
                      alpha_cov = 0)
fixed_terms <- NULL
bootstrap_samples <- 3
tt <- cxr_pm_fit(data = data,
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
summary(tt)
tt$alpha_inter
