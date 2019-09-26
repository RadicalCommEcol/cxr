context("models")

skip_on_cran()

# 1 - BH fecundity models----
# generate test data
focal.sp <- c(1,2)
num.sp <- 5
num.cov <- 2
num.obs <- 10 # sites, per focal species

focal.lambda <- c(100,200)
alpha.matrix.orig <- matrix(data = runif(num.sp*num.sp,0.001,0.1),nrow = num.sp, ncol = num.sp)
lambda.cov.orig <- matrix(runif(num.sp*num.cov,0.001,0.1),nrow = num.sp, ncol = num.cov)
alpha.cov.orig <- list()
for(i.cov in 1:num.cov){
  alpha.cov.orig[[i.cov]] <- matrix(data = runif(num.sp*num.cov,0.001,0.1),nrow = num.sp, ncol = num.sp) 
}

test.data <- GenerateTestData(focal.sp = focal.sp,
                              num.sp = num.sp,
                              num.cov = num.cov,
                              num.obs = num.obs,
                              fitness.model = 5,
                              focal.lambda = focal.lambda,
                              alpha = alpha.matrix.orig,
                              alpha.cov = alpha.cov.orig,
                              lambda.cov = lambda.cov.orig)
# select a focal sp
test.focal <- test.data[test.data$focal == 1,]
focal.comp.matrix <- test.focal[,2:6]
focal.covariates <- test.focal[,7:8]
my.log.fitness <- log(test.focal$fitness)

# param sets for the five models
my.lambda <- focal.lambda[1]
my.sigma <- sd(log(test.focal$fitness))
my.mean.alpha <- mean(alpha.matrix.orig[1,])
my.full.alpha <- alpha.matrix.orig[1,]
my.lambda.cov <- lambda.cov.orig[1,]
my.mean.alpha.cov <- c(mean(alpha.cov.orig[[1]][1,]),mean(alpha.cov.orig[[2]][1,]))
my.full.alpha.cov <- c(alpha.cov.orig[[1]][1,],alpha.cov.orig[[2]][1,])

param.list <- list(c("lambda","alpha"),
                c("lambda","alpha"),
                c("lambda","alpha","lambda.cov","alpha.cov"),
                c("lambda","alpha","lambda.cov","alpha.cov"))

par.1 <- c(my.lambda,my.sigma)
par.2 <- c(my.lambda,my.mean.alpha,my.sigma)
par.3 <- c(my.lambda,my.full.alpha,my.sigma)
par.4 <- c(my.lambda,my.lambda.cov,my.full.alpha,my.mean.alpha.cov,my.sigma)
par.5 <- c(my.lambda,my.lambda.cov,my.full.alpha,my.full.alpha.cov,my.sigma)

results_BH_1 <- model_BH1(par = par.1,log.fitness = my.log.fitness)
results_BH_2 <- model_BH2(par = par.2,
                     log.fitness = my.log.fitness,
                     param.list = param.list[[1]],
                     focal.comp.matrix = focal.comp.matrix,
                     num.competitors = ncol(focal.comp.matrix))
results_BH_3 <- model_BH3(par = par.3,
                     log.fitness = my.log.fitness,
                     param.list = param.list[[2]],
                     focal.comp.matrix = focal.comp.matrix,
                     num.competitors = ncol(focal.comp.matrix))
results_BH_4 <- model_BH4(par = par.4,
                     log.fitness = my.log.fitness,
                     param.list = param.list[[3]],
                     focal.comp.matrix = focal.comp.matrix,
                     num.competitors = ncol(focal.comp.matrix),
                     num.covariates = ncol(focal.covariates),
                     focal.covariates = focal.covariates)
results_BH_5 <- model_BH5(par = par.5,
                     log.fitness = my.log.fitness,
                     param.list = param.list[[4]],
                     focal.comp.matrix = focal.comp.matrix,
                     num.competitors = ncol(focal.comp.matrix),
                     num.covariates = ncol(focal.covariates),
                     focal.covariates = focal.covariates)
# test----
test_that("Expected classes", {
  expect_equal(class(results_BH_1), "numeric")
  expect_equal(class(results_BH_2), "numeric")
  expect_equal(class(results_BH_3), "numeric")
  expect_equal(class(results_BH_4), "numeric")
  expect_equal(class(results_BH_5), "numeric")
})

# 2 - BH abundance models----
# use the two focal species above
sp.par <- data.frame(lambda = focal.lambda,
                     germ.rate = runif(2,0.2,0.9),
                     survival.rate = runif(2,0.2,0.9))
init.abund <- runif(2,1,100)
alpha.matrix.sub <- alpha.matrix.orig[1:2,1:2]
mean.alpha.sub <- mean(alpha.matrix.sub)
cov.values <- c(test.focal$cov1[1],test.focal$cov2[1])
lambda.cov.sub <- lambda.cov.orig[1:2,]
alpha.cov.sub <- list(alpha.cov.orig[[1]][1:2,1:2],alpha.cov.orig[[2]][1:2,1:2])
mean.alpha.cov.sub <- list(mean(alpha.cov.sub[[1]]),mean(alpha.cov.sub[[2]]))

results_abund_1 <- model_abundBH1(sp.par = sp.par,
                                  init.abund = init.abund,
                                  cov.values = cov.values,
                                  alpha.matrix = alpha.matrix.sub,
                                  lambda.cov.matrix = lambda.cov.sub,
                                  alpha.cov.matrix = alpha.cov.sub)
results_abund_2 <- model_abundBH2(sp.par = sp.par,
                                  init.abund = init.abund,
                                  cov.values = cov.values,
                                  alpha.matrix = mean.alpha.sub,
                                  lambda.cov.matrix = lambda.cov.sub,
                                  alpha.cov.matrix = alpha.cov.sub)
results_abund_3 <- model_abundBH3(sp.par = sp.par,
                                  init.abund = init.abund,
                                  cov.values = cov.values,
                                  alpha.matrix = alpha.matrix.sub,
                                  lambda.cov.matrix = lambda.cov.sub,
                                  alpha.cov.matrix = alpha.cov.sub)
results_abund_4 <- model_abundBH4(sp.par = sp.par,
                                  init.abund = init.abund,
                                  cov.values = cov.values,
                                  alpha.matrix = alpha.matrix.sub,
                                  lambda.cov.matrix = lambda.cov.sub,
                                  alpha.cov.matrix = mean.alpha.cov.sub)
results_abund_5 <- model_abundBH5(sp.par = sp.par,
                                  init.abund = init.abund,
                                  cov.values = cov.values,
                                  alpha.matrix = alpha.matrix.sub,
                                  lambda.cov.matrix = lambda.cov.sub,
                                  alpha.cov.matrix = alpha.cov.sub)
# test----
test_that("Expected classes", {
  expect_equal(class(results_abund_1), "numeric")
  expect_equal(class(results_abund_2), "numeric")
  expect_equal(class(results_abund_3), "numeric")
  expect_equal(class(results_abund_4), "numeric")
  expect_equal(class(results_abund_5), "numeric")
})

