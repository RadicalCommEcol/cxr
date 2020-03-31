context("niche overlap and species fitness functions")
# test niche overlap and species fitness

# data --------------------------------------------------------------------

data("neigh_list")
data <- neigh_list
# keep only fitness and neighbours columns
for(i in 1:length(data)){
  data[[i]] <- data[[i]][,2:length(data[[i]])]
}

focal_column <- names(data)
data("salinity_list")
salinity <- salinity_list
# keep only salinity column
for(i in 1:length(salinity)){
  salinity[[i]] <- as.matrix(salinity[[i]][,2:length(salinity[[i]])])
  colnames(salinity[[i]]) <- "salinity"
}
model_families <- c("BH","LV","RK","LW")
covariates <- salinity
optimization_method <- "nlminb"
alpha_form <- "pairwise"
lambda_cov_form <- "none"
alpha_cov_form <- "none"
initial_values = list(lambda = 1,
                      alpha_intra = 0.1,
                      alpha_inter = 0.1,
                      lambda_cov = 0.1,
                      alpha_cov = 0.1)
lower_bounds = list(lambda = 0,
                    alpha_intra = 0.01,
                    alpha_inter = 0,
                    lambda_cov = -1,
                    alpha_cov = -1)
upper_bounds = list(lambda = 100,
                    alpha_intra = 1,
                    alpha_inter = 1,
                    lambda_cov = 1,
                    alpha_cov = 1)

initial_values_er = list(lambda = 1,effect = 0.01,response = 0.01)
lower_bounds_er = list(lambda = 0,effect = 0,response = 0)
upper_bounds_er = list(lambda = 100,effect = 1,response = 1)

fixed_terms <- NULL
bootstrap_samples <- 3


# 1 - two cxr objects ------------------------------------------------

sp1.fit <- cxr::cxr_pm_fit(data = data[[1]],
                           focal_column = focal_column[1],
                           model_family = model_families[1], # irrelevant
                           covariates = covariates[[1]],
                           optimization_method = optimization_method,
                           alpha_form = alpha_form,
                           lambda_cov_form = lambda_cov_form,
                           alpha_cov_form = alpha_cov_form,
                           initial_values = initial_values,
                           lower_bounds = lower_bounds,
                           upper_bounds = upper_bounds,
                           fixed_terms = fixed_terms,
                           bootstrap_samples = bootstrap_samples)

sp2.fit <- cxr::cxr_pm_fit(data = data[[2]],
                           focal_column = focal_column[2],
                           model_family = model_families[1],
                           covariates = covariates[[2]],
                           optimization_method = optimization_method,
                           alpha_form = alpha_form,
                           lambda_cov_form = lambda_cov_form,
                           alpha_cov_form = alpha_cov_form,
                           initial_values = initial_values,
                           lower_bounds = lower_bounds,
                           upper_bounds = upper_bounds,
                           fixed_terms = fixed_terms,
                           bootstrap_samples = bootstrap_samples)


# 2 - generalization to n species, across all pairs ---------

cxr_multifit <- cxr_pm_multifit(data = data,
                                model_family = model_families[1], # irrelevant
                                focal_column = focal_column,
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


# 3 - pairwise matrix --------------------------------

# make sure intra>inter
posmatrix <- matrix(c(runif(1,0.5,1),runif(2,0,0.5),runif(1,0.5,1)),nrow = 2)
poslambdas <- runif(2,1,10)
model_family <- "BH"

# test --------------------------------------------------------------------
test_that("niche overlaps are correctly calculated", {
  # multispecies
  expect_s3_class(niche_overlap(cxr_multifit = cxr_multifit),"data.frame")
  # 2 cxr objects
  # expect_vector(object = niche_overlap(cxr_sp1 = sp1.fit,cxr_sp2 = sp2.fit),ptype = numeric,size = 2)
  expect_length(niche_overlap(cxr_sp1 = sp1.fit,cxr_sp2 = sp2.fit),2)
  # pairwise matrix
  # expect_vector(niche_overlap(pair_matrix = posmatrix))
  expect_length(niche_overlap(pair_matrix = posmatrix),2)
})

# test --------------------------------------------------------------------
test_that("average fitness differenes are correctly calculated", {
  # multispecies
  expect_s3_class(avg_fitness_diff(cxr_multifit = cxr_multifit),"data.frame")
  # 2 cxr objects
  expect_s3_class(avg_fitness_diff(cxr_sp1 = sp1.fit,
                                   cxr_sp2 = sp2.fit),"data.frame")
  # pairwise matrix
  expect_s3_class(avg_fitness_diff(pair_lambdas = poslambdas,
                                   pair_matrix = posmatrix,
                                   model_family = model_family),"data.frame")
  
})

# test --------------------------------------------------------------------
test_that("competitive ability is correctly calculated", {
  # multispecies
  expect_s3_class(competitive_ability(cxr_multifit = cxr_multifit),"data.frame")
  # 2 cxr objects
  expect_s3_class(competitive_ability(cxr_sp1 = sp1.fit,
                                   cxr_sp2 = sp2.fit),"data.frame")
  # pairwise matrix
  expect_s3_class(competitive_ability(lambda = poslambdas[1],
                                   pair_matrix = posmatrix,
                                   model_family = model_family),"data.frame")
  
})

# 4 - n-sp fitness, and difference, for all model families --------------

# fit three species at once
data("neigh_list")
# these species all have >250 observations
example_sp <- c(1,5,6)
n.obs <- 250
data <- neigh_list[example_sp]
optimization_method <- "nlminb"
fixed_terms <- NULL
model_families <- c("BH","LV","RK","LW")
bootstrap_samples <- 0
# keep only fitness and neighbours columns
for(i in 1:length(data)){
  data[[i]] <- data[[i]][1:n.obs,c(2,example_sp+2)]#2:length(data[[i]])]
}
initial_values_er = list(lambda = 1, 
                         effect = 1, 
                         response = 1)
lower_bounds_er = list(lambda = 0, 
                       effect = 0, 
                       response = 0)
upper_bounds_er = list(lambda = 100, 
                       effect = 10, 
                       response = 10)
iv_LV <- list(lambda = 1, 
              effect = 0, 
              response = 0)
lb_LV = list(lambda = 0, 
                       effect = -1, 
                       response = -1)
ub_LV = list(lambda = 100, 
                       effect = 1, 
                       response = 1)
# test --------------------------------------------------------------------
test_that("species fitness are correctly calculated", {
  
  # skip this on CRAN, 
  # as it may take long
  skip_on_cran()
  
  for(i.model in 1:length(model_families)){
    
    # lotka-volterra are trickier to fit,
    # need different sets of initial values and bounds
    if(model_families[i.model] == "LV"){
      iv <- iv_LV
      lb <- lb_LV
      ub <- ub_LV
    }else{
      iv <- initial_values_er
      lb <- lower_bounds_er
      ub <- upper_bounds_er
    }
    
    er.fit <- cxr::cxr_er_fit(data = data,
                              model_family = model_families[i.model],
                              # covariates = covariates,
                              optimization_method = optimization_method,
                              initial_values = iv,
                              lower_bounds = lb,
                              upper_bounds = ub,
                              fixed_terms = fixed_terms,
                              bootstrap_samples = bootstrap_samples)
    
    spfitness <- cxr::species_fitness(er.fit)
    # expect_vector(spfitness)
    expect_length(spfitness,length(data))
  }
})