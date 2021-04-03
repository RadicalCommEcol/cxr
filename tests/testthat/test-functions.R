context("parameter fitting functions")

# set data ----------------------------------------------------------------

data("neigh_list")
data("salinity_list")
my.sp <- "BEMA"
model_families <- c("BH","RK","LV","LW")
sp_data <- neigh_list[[my.sp]][2:ncol(neigh_list[[1]])]
covariates <- salinity_list[[my.sp]][2]
all.sp <- names(sp_data)[2:length(names(sp_data))]
neigh.sp <- all.sp[which(all.sp != my.sp)]

# test every configuration of alpha, lambda_cov, alpha_cov
# for every model family included by default

param.conf <- list()
param.conf[[1]] <- c("none","none","none")
param.conf[[2]] <- c("global","none","none")
param.conf[[3]] <- c("pairwise","none","none")
param.conf[[4]] <- c("pairwise","global","global")
param.conf[[5]] <- c("pairwise","global","pairwise")
initial_values <- list(BH = list(lambda = 1,
                                 alpha_intra = 0,
                                 alpha_inter = 0,
                                 lambda_cov = 0,
                                 alpha_cov = 0),
                       LV = list(lambda = 1,
                                 alpha_intra = 0,
                                 alpha_inter = 0,
                                 lambda_cov = 0,
                                 alpha_cov = 0),
                       RK = list(lambda = 1,
                                 alpha_intra = 0,
                                 alpha_inter = 0,
                                 lambda_cov = 0,
                                 alpha_cov = 0),
                       LW = list(lambda = 1,
                                 alpha_intra = 0,
                                 alpha_inter = 0,
                                 lambda_cov = 0,
                                 alpha_cov = 0))

lower_bounds = list(BH = list(lambda = 0,
                              alpha_intra = 0,
                              alpha_inter = -1,
                              lambda_cov = -1,
                              alpha_cov = -1),
                    LV = list(lambda = 1,
                              alpha_intra = -1,
                              alpha_inter = -1,
                              lambda_cov = -1,
                              alpha_cov = -1),
                    RK = list(lambda = 0,
                              alpha_intra = 0,
                              alpha_inter = -1,
                              lambda_cov = -1,
                              alpha_cov = -1),
                    LW = list(lambda = 0,
                              alpha_intra = 0,
                              alpha_inter = -1,
                              lambda_cov = -1,
                              alpha_cov = -1))

upper_bounds = list(BH = list(lambda = 1000,
                              alpha_intra = 1,
                              alpha_inter = 1,
                              lambda_cov = 1,
                              alpha_cov = 1),
                    LV = list(lambda = 100,
                              alpha_intra = 0,
                              alpha_inter = 0,
                              lambda_cov = 0,
                              alpha_cov = 0),
                    RK = list(lambda = 100,
                              alpha_intra = 1,
                              alpha_inter = 1,
                              lambda_cov = 1,
                              alpha_cov = 1),
                    LW = list(lambda = 100,
                              alpha_intra = 1,
                              alpha_inter = 1,
                              lambda_cov = 1,
                              alpha_cov = 1))
# when run with nonsensical parameters/initial values
# it returns NULL
test_that("invalid arguments return NULL",{
  wrong.data <- data.frame(c1 = c(1,2,3),c2 = c("ab","a","a"))
  wrong_fit1 <- cxr_pm_fit(data = wrong.data,
                           model_family = "BH",
                           alpha_form = "none",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none")
  
  wrong_fit2 <- cxr_pm_fit(data = sp_data,
                           model_family = "BH",
                           alpha_form = "none",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none")
  
  wrong_fit3 <- cxr_pm_fit(data = sp_data,
                           model_family = "BH",
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,# here is the error
                                                 alpha_inter = 0.1, 
                                                 lambda_cov = 0.1, 
                                                 alpha_cov = 0.1),
                           alpha_form = "none", 
                           lambda_cov_form = "none",
                           alpha_cov_form = "none")
  expect_null(wrong_fit1)
  expect_null(wrong_fit2)
  expect_null(wrong_fit3)
  
})

test_that("parameter form 1 returns valid object", {
  
  for(im in 1:length(model_families)){
    sp_fit <- cxr_pm_fit(data = sp_data,
                         model_family = model_families[im],
                         focal_column = my.sp,
                         optimization_method = "bobyqa",
                         alpha_form = param.conf[[1]][1],
                         lambda_cov_form = param.conf[[1]][2],
                         alpha_cov_form = param.conf[[1]][3],
                         initial_values = initial_values[[model_families[im]]],
                         lower_bounds = lower_bounds[[model_families[im]]],
                         upper_bounds = upper_bounds[[model_families[im]]],
                         bootstrap_samples = 0)
    
    lambda <- sp_fit$lambda
    alpha_intra <- sp_fit$alpha_intra
    alpha_inter <- sp_fit$alpha_inter
    lambda_cov <- sp_fit$lambda_cov
    alpha_cov <- sp_fit$alpha_cov
    
    lambda_standard_error <- sp_fit$lambda_standard_error
    alpha_intra_standard_error <- sp_fit$alpha_intra_standard_error
    alpha_inter_standard_error <- sp_fit$alpha_inter_standard_error
    lambda_cov_standard_error <- sp_fit$lambda_cov_standard_error
    alpha_cov_standard_error <- sp_fit$alpha_cov_standard_error
    
    log_likelihood <- sp_fit$log_likelihood
    
    expect_s3_class(sp_fit, "cxr_pm_fit")
    expect_type(lambda, "double")
    
    expect_null(alpha_intra)
    expect_null(alpha_inter)
    expect_null(lambda_cov)
    expect_null(alpha_cov)
    expect_null(lambda_standard_error)
    expect_null(alpha_intra_standard_error)
    expect_null(alpha_inter_standard_error)
    expect_null(lambda_cov_standard_error)
    expect_null(alpha_cov_standard_error)
    expect_type(log_likelihood, "double")
  }# for each model family
})

test_that("parameter form 2 returns valid object", {
  
  for(im in 1:length(model_families)){
    sp_fit <- cxr_pm_fit(data = sp_data,
                         model_family = model_families[im],
                         focal_column = my.sp,
                         optimization_method = "bobyqa",
                         alpha_form = param.conf[[2]][1],
                         lambda_cov_form = param.conf[[2]][2],
                         alpha_cov_form = param.conf[[2]][3],
                         initial_values = initial_values[[model_families[im]]],
                         lower_bounds = lower_bounds[[model_families[im]]],
                         upper_bounds = upper_bounds[[model_families[im]]],
                         bootstrap_samples = 0)
  
    lambda <- sp_fit$lambda
    alpha_intra <- sp_fit$alpha_intra
    alpha_inter <- sp_fit$alpha_inter
    lambda_cov <- sp_fit$lambda_cov
    alpha_cov <- sp_fit$alpha_cov
    
    lambda_standard_error <- sp_fit$lambda_standard_error
    alpha_intra_standard_error <- sp_fit$alpha_intra_standard_error
    alpha_inter_standard_error <- sp_fit$alpha_inter_standard_error
    lambda_cov_standard_error <- sp_fit$lambda_cov_standard_error
    alpha_cov_standard_error <- sp_fit$alpha_cov_standard_error
    
    log_likelihood <- sp_fit$log_likelihood
    
    expect_s3_class(sp_fit, "cxr_pm_fit")
    expect_type(lambda, "double")
    expect_null(alpha_intra)
    expect_type(alpha_inter, "double")
    
    expect_null(lambda_cov)
    expect_null(alpha_cov)
    expect_null(lambda_standard_error)
    expect_null(alpha_intra_standard_error)
    expect_null(alpha_inter_standard_error)
    expect_null(lambda_cov_standard_error)
    expect_null(alpha_cov_standard_error)
    expect_type(log_likelihood, "double")
  }
})

test_that("parameter form 3 returns valid object", {
  
  for(im in 1:length(model_families)){
    sp_fit <- cxr_pm_fit(data = sp_data,
                         model_family = model_families[im],
                         focal_column = my.sp,
                         optimization_method = "bobyqa",
                         alpha_form = param.conf[[3]][1],
                         lambda_cov_form = param.conf[[3]][2],
                         alpha_cov_form = param.conf[[3]][3],
                         initial_values = initial_values[[model_families[im]]],
                         lower_bounds = lower_bounds[[model_families[im]]],
                         upper_bounds = upper_bounds[[model_families[im]]],
                         bootstrap_samples = 0)
    
    lambda <- sp_fit$lambda
    alpha_intra <- sp_fit$alpha_intra
    alpha_inter <- sp_fit$alpha_inter
    lambda_cov <- sp_fit$lambda_cov
    alpha_cov <- sp_fit$alpha_cov
    
    lambda_standard_error <- sp_fit$lambda_standard_error
    alpha_intra_standard_error <- sp_fit$alpha_intra_standard_error
    alpha_inter_standard_error <- sp_fit$alpha_inter_standard_error
    lambda_cov_standard_error <- sp_fit$lambda_cov_standard_error
    alpha_cov_standard_error <- sp_fit$alpha_cov_standard_error
    
    log_likelihood <- sp_fit$log_likelihood
    
    expect_s3_class(sp_fit, "cxr_pm_fit")
    expect_type(lambda, "double")
    expect_type(alpha_intra, "double")
    expect_type(alpha_inter, "double")
    expect_length(alpha_inter,length(neigh.sp))
    
    expect_null(lambda_cov)
    expect_null(alpha_cov)
    expect_null(lambda_standard_error)
    expect_null(alpha_intra_standard_error)
    expect_null(alpha_inter_standard_error)
    expect_null(lambda_cov_standard_error)
    expect_null(alpha_cov_standard_error)
    expect_type(log_likelihood, "double")
  }
})

test_that("parameter form 4 returns valid object", {
  
  for(im in 1:length(model_families)){
    sp_fit <- cxr_pm_fit(data = sp_data,
                         model_family = model_families[im],
                         focal_column = my.sp,
                         covariates = covariates,
                         optimization_method = "bobyqa",
                         alpha_form = param.conf[[4]][1],
                         lambda_cov_form = param.conf[[4]][2],
                         alpha_cov_form = param.conf[[4]][3],
                         initial_values = initial_values[[model_families[im]]],
                         lower_bounds = lower_bounds[[model_families[im]]],
                         upper_bounds = upper_bounds[[model_families[im]]],
                         bootstrap_samples = 0)
  
    lambda <- sp_fit$lambda
    alpha_intra <- sp_fit$alpha_intra
    alpha_inter <- sp_fit$alpha_inter
    lambda_cov <- sp_fit$lambda_cov
    alpha_cov <- sp_fit$alpha_cov
    
    lambda_standard_error <- sp_fit$lambda_standard_error
    alpha_intra_standard_error <- sp_fit$alpha_intra_standard_error
    alpha_inter_standard_error <- sp_fit$alpha_inter_standard_error
    lambda_cov_standard_error <- sp_fit$lambda_cov_standard_error
    alpha_cov_standard_error <- sp_fit$alpha_cov_standard_error
    
    log_likelihood <- sp_fit$log_likelihood
    
    expect_s3_class(sp_fit, "cxr_pm_fit")
    expect_type(lambda, "double")
    expect_type(alpha_intra, "double")
    expect_type(alpha_inter, "double")
    expect_length(alpha_inter,length(neigh.sp))
    expect_type(lambda_cov, "double")
    expect_type(alpha_cov, "list")
    for(i.cov in 1:length(covariates)){
      expect_type(alpha_cov[[i.cov]], "double")
    }
    
    expect_null(lambda_standard_error)
    expect_null(alpha_intra_standard_error)
    expect_null(alpha_inter_standard_error)
    expect_null(lambda_cov_standard_error)
    expect_null(alpha_cov_standard_error)
    expect_type(log_likelihood, "double")
  }
})

test_that("parameter form 5 returns valid object", {
  
  for(im in 1:length(model_families)){
    sp_fit <- cxr_pm_fit(data = sp_data,
                         model_family = model_families[im],
                         focal_column = my.sp,
                         covariates = covariates,
                         optimization_method = "bobyqa",
                         alpha_form = param.conf[[5]][1],
                         lambda_cov_form = param.conf[[5]][2],
                         alpha_cov_form = param.conf[[5]][3],
                         initial_values = initial_values[[model_families[im]]],
                         lower_bounds = lower_bounds[[model_families[im]]],
                         upper_bounds = upper_bounds[[model_families[im]]],
                         bootstrap_samples = 0)
  
    lambda <- sp_fit$lambda
    alpha_intra <- sp_fit$alpha_intra
    alpha_inter <- sp_fit$alpha_inter
    lambda_cov <- sp_fit$lambda_cov
    alpha_cov <- sp_fit$alpha_cov
    
    lambda_standard_error <- sp_fit$lambda_standard_error
    alpha_intra_standard_error <- sp_fit$alpha_intra_standard_error
    alpha_inter_standard_error <- sp_fit$alpha_inter_standard_error
    lambda_cov_standard_error <- sp_fit$lambda_cov_standard_error
    alpha_cov_standard_error <- sp_fit$alpha_cov_standard_error
    
    log_likelihood <- sp_fit$log_likelihood
    
    expect_s3_class(sp_fit, "cxr_pm_fit")
    expect_type(lambda, "double")
    expect_type(alpha_intra, "double")
    expect_type(alpha_inter, "double")
    expect_length(alpha_inter,length(neigh.sp))
    expect_type(lambda_cov, "double")
    expect_type(alpha_cov, "list")
    for(i.cov in 1:length(covariates)){
      expect_type(alpha_cov[[i.cov]], "double")
      expect_length(alpha_cov[[i.cov]],length(all.sp))
    }
    
    expect_null(lambda_standard_error)
    expect_null(alpha_intra_standard_error)
    expect_null(alpha_inter_standard_error)
    expect_null(lambda_cov_standard_error)
    expect_null(alpha_cov_standard_error)
    expect_type(log_likelihood, "double")
  }
})

test_that("errors are correctly calculated", {
  
  for(im in 1:length(model_families)){
    sp_fit <- cxr_pm_fit(data = sp_data,
                         model_family = model_families[im],
                         focal_column = my.sp,
                         covariates = covariates,
                         optimization_method = "bobyqa",
                         alpha_form = param.conf[[5]][1],
                         lambda_cov_form = param.conf[[5]][2],
                         alpha_cov_form = param.conf[[5]][3],
                         initial_values = initial_values[[model_families[im]]],
                         lower_bounds = lower_bounds[[model_families[im]]],
                         upper_bounds = upper_bounds[[model_families[im]]],
                         bootstrap_samples = 3)
    
    lambda <- sp_fit$lambda
    alpha_intra <- sp_fit$alpha_intra
    alpha_inter <- sp_fit$alpha_inter
    lambda_cov <- sp_fit$lambda_cov
    alpha_cov <- sp_fit$alpha_cov
    
    lambda_standard_error <- sp_fit$lambda_standard_error
    alpha_intra_standard_error <- sp_fit$alpha_intra_standard_error
    alpha_inter_standard_error <- sp_fit$alpha_inter_standard_error
    lambda_cov_standard_error <- sp_fit$lambda_cov_standard_error
    alpha_cov_standard_error <- sp_fit$alpha_cov_standard_error
    
    log_likelihood <- sp_fit$log_likelihood
    
    expect_s3_class(sp_fit, "cxr_pm_fit")
    expect_type(lambda, "double")
    expect_type(alpha_intra, "double")
    expect_type(alpha_inter, "double")
    expect_length(alpha_inter,length(neigh.sp))
    expect_type(lambda_cov, "double")
    expect_type(alpha_cov, "list")
    for(i.cov in 1:length(covariates)){
      expect_type(alpha_cov[[i.cov]], "double")
      expect_length(alpha_cov[[i.cov]],length(all.sp))
    }
    expect_type(lambda_standard_error, "double")
    expect_type(alpha_intra_standard_error, "double")
    expect_type(alpha_inter_standard_error, "double")
    expect_length(alpha_inter_standard_error,length(neigh.sp))
    expect_type(lambda_cov_standard_error, "double")
    expect_type(alpha_cov_standard_error, "list")
    for(i.cov in 1:length(covariates)){
      expect_type(alpha_cov_standard_error[[i.cov]], "double")
      expect_length(alpha_cov_standard_error[[i.cov]],length(all.sp))
    }
    expect_type(log_likelihood, "double")
  }
})

# test multifit -----------------------------------------------------------

three_sp <- c("BEMA","LEMA","HOMA")
sp.pos <- which(names(neigh_list) %in% three_sp)

data <- neigh_list[sp.pos]
# keep only fitness and neighbours columns
for(i in 1:length(data)){
  data[[i]] <- data[[i]][,2:length(data[[i]])]
}

focal_column <- names(data)

salinity <- salinity_list[sp.pos]

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
fixed_terms <- NULL
bootstrap_samples <- 3

test_that("multiple species are correctly fitted", {
  # skip this on CRAN, 
  # as it may take long
  skip_on_cran()
  multifit <- cxr_pm_multifit(data = data,
                              focal_column = focal_column,
                              model_family = model_family,
                              covariates = covariates,
                              optimization_method = optimization_method,
                              alpha_form = alpha_form,
                              lambda_cov_form = lambda_cov_form,
                              alpha_cov_form = alpha_cov_form,
                              initial_values = initial_values[[model_family]],
                              lower_bounds = lower_bounds[[model_family]],
                              upper_bounds = upper_bounds[[model_family]],
                              fixed_terms = fixed_terms,
                              bootstrap_samples = bootstrap_samples)
  
  lambda <- multifit$lambda
  alpha_matrix <- multifit$alpha_matrix
  lc <- multifit$lambda_cov
  ac <- multifit$alpha_cov
  lse <- multifit$lambda_standard_error
  ase <- multifit$alpha_matrix_standard_error
  lcse <- multifit$lambda_cov_standard_error
  acse <- multifit$alpha_cov_standard_error
  llik <- multifit$log_likelihood
  
  expect_s3_class(multifit, "cxr_pm_multifit")
  expect_type(lambda, "double")
  expect_type(alpha_matrix, "double")
  expect_type(lc, "double")
  expect_type(ac, "list")
  for(i.cov in 1:ncol(covariates[[1]])){
    expect_type(ac[[i.cov]], "double")
  }
  expect_type(lse, "double")
  expect_type(ase, "double")
  expect_type(lcse, "double")
  expect_type(acse, "list")
  for(i.cov in 1:ncol(covariates[[1]])){
    expect_type(acse[[i.cov]], "double")
  }
  expect_type(llik, "double")
  
})

