context("abundance projections")
# test abundance projections

# data --------------------------------------------------------------------

model_family = c("BH")
alpha_form = "pairwise"
lambda_cov_form = "global"
alpha_cov_form = "pairwise"
lambda = c(sp1 = 1,sp2 = 2,sp3 = 3)
alpha_matrix = matrix(runif(9,0,1),nrow = 3, dimnames = list(c("sp1","sp2","sp3"),
                                                             c("sp1","sp2","sp3")))
lambda_cov = matrix(runif(3,0,1),nrow = 3, dimnames = list(c("sp1","sp2","sp3"),
                                                           c("c1")))
alpha_cov = list(c1 = matrix(runif(9,0,0.5),
                                   nrow = 3,
                                   dimnames = list(c("sp1","sp2","sp3"),
                                                   c("sp1","sp2","sp3"))))

covariates_proj <- matrix(runif(10,0,1),nrow = 10,dimnames = list(c(paste("t",1:10,sep="")),c("c1")))
timesteps <- 10
initial_abundances <- c(sp1 = 10,sp2 = 12,sp3 = 13)


# projection function -----------------------------------------------------

for(i.m in 1:length(model_family)){
  ab <- abundance_projection(cxr_fit = NULL,
                             model_family = model_family[i.m],
                             alpha_form = alpha_form,
                             lambda_cov_form = lambda_cov_form,
                             alpha_cov_form = alpha_cov_form,
                             lambda = lambda,
                             alpha_matrix = alpha_matrix,
                             lambda_cov = lambda_cov,
                             alpha_cov = alpha_cov,
                             covariates = covariates_proj,
                             timesteps = timesteps,
                             initial_abundances = initial_abundances)
  
  test_that("abundances are correctly projected", {

    expect_is(ab, "matrix")
    expect_type(ab, "double")
    
    })
}

# cxr object --------------------------------------------------------------

data("neigh_list")
data("salinity_list")
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
param.conf <- list()
param.conf[[1]] <- c("none","none","none")
param.conf[[2]] <- c("global","none","none")
param.conf[[3]] <- c("pairwise","none","none")
param.conf[[4]] <- c("pairwise","global","global")
param.conf[[5]] <- c("pairwise","global","pairwise")
fixed_terms <- NULL
bootstrap_samples <- 3
timesteps <- 10
initial_abundances <- c("BEMA" = 10,"HOMA" = 12,"LEMA" = 13)
covariates_proj <- matrix(runif(10,0,1),nrow = 10,dimnames = list(c(paste("t",1:10,sep="")),c("c1")))

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

# function ----------------------------------------------------------------
for(i.conf in 1:length(param.conf)){
for(i.m in 1:length(model_family)){
  
  cxr_fit <- cxr_pm_multifit(data = data,
                             focal_column = focal_column,
                             model_family = model_family[i.m],
                             covariates = salinity,
                             optimization_method = optimization_method,
                             alpha_form = param.conf[[i.conf]][1],
                             lambda_cov_form = param.conf[[i.conf]][2],
                             alpha_cov_form = param.conf[[i.conf]][3],
                             initial_values = initial_values[[model_family[i.m]]],
                             lower_bounds = lower_bounds[[model_family[i.m]]],
                             upper_bounds = upper_bounds[[model_family[i.m]]],
                             fixed_terms = fixed_terms,
                             bootstrap_samples = bootstrap_samples)
  
  ab <- abundance_projection(cxr_fit = cxr_fit,
                             covariates = covariates_proj,
                             timesteps = timesteps,
                             initial_abundances = initial_abundances)
  
  test_that("abundances are correctly projected", {
    
    expect_is(ab, "matrix")
    expect_type(ab, "double")
    
  })
}# for each family
}# for each model
