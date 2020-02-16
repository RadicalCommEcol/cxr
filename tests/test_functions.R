context("functions")

# set data ----------------------------------------------------------------

data("neigh_list")
my.sp <- "BEMA"
sp_data <- neigh_list[[my.sp]][2:ncol(neigh_list[[1]])]

# test every configuration of alpha, lambda_cov, alpha_cov
param.conf <- list()
param.conf[[1]] <- c("none","none","none")
param.conf[[2]] <- c("global","none","none")
param.conf[[3]] <- c("pairwise","none","none")
param.conf[[4]] <- c("pairwise","global","global")
param.conf[[5]] <- c("pairwise","global","pairwise")

for(i.conf in 1:length(param.conf)){

  test_that("Expected classes", {
  sp_fit <- cxr_pm_fit(data = sp_data,
                       focal_column = my.sp,
                       optimization_method = "bobyqa",
                       alpha_form = param.conf[[i.conf]][1],
                       lambda_cov_form = param.conf[[i.conf]][2],
                       alpha_cov_form = param.conf[[i.conf]][3],
                       initial_values = list(lambda = 1,alpha_intra = 0.1,alpha_inter = 0.1, lambda_cov = 0.1, alpha_cov = 0.1),
                       lower_bounds = list(lambda = 0,alpha_intra = 0,alpha_inter = 0,lambda_cov = 0, alpha_cov = 0),
                       upper_bounds = list(lambda = 100,alpha_intra = 1,alpha_inter = 1, lambda_cov = 1, alpha_cov = 1),
                       bootstrap_samples = 3)
  

  expect_equal(class(sp_fit), "cxr_pm_fit")
  expect_equal(class(sp_fit$lambda), "numeric")
  expect_equal(class(sp_fit$lambda.lower.error), "numeric")
  expect_equal(class(sp_fit$lambda.upper.error), "numeric")
  expect_equal(class(sp_fit$sigma), "numeric")
  expect_equal(class(sp_fit$alpha), "numeric")
  expect_equal(class(sp_fit$alpha.lower.error), "numeric")
  expect_equal(class(sp_fit$alpha.upper.error), "numeric")
  expect_equal(class(sp_fit$lambda.cov), "numeric")
  expect_equal(class(sp_fit$lambda.cov.lower.error), "numeric")
  expect_equal(class(sp_fit$lambda.cov.upper.error), "numeric")
  expect_equal(class(sp_fit$alpha.cov), "numeric")
  expect_equal(class(sp_fit$alpha.cov.lower.error), "numeric")
  expect_equal(class(sp_fit$alpha.cov.upper.error), "numeric")
  expect_equal(class(sp_fit$log.likelihood), "numeric")
})
}