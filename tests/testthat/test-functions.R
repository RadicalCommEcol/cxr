context("functions")

skip_on_cran()

# test optimizer functions----

# set params, models, methods...

# generate test data
focal.sp <- c(1,2)
num.sp <- 5
num.cov <- 1
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

# function params
fitness.model <- BH_5
optim.method <- "optim_L-BFGS-B"
param.list <- c("lambda","alpha","lambda.cov","alpha.cov")
log.fitness <- log(test.focal$fitness)
init.lambda <- focal.lambda[1]
lower.lambda <- 1e-5
upper.lambda <- 1e4
init.sigma <- sd(log(test.focal$fitness))
lower.sigma <- 1e-5
upper.sigma <- 1e2
init.alpha <- alpha.matrix.orig[1,]
lower.alpha <- 1e-5
upper.alpha <- 1e4
init.lambda.cov <- lambda.cov.orig[1,]
lower.lambda.cov <- 1e-5
upper.lambda.cov <- 1e4
init.alpha.cov <- alpha.cov.orig[[1]][1,]
lower.alpha.cov <- 1e-5
upper.alpha.cov <- 1e4
focal.comp.matrix <- test.focal[,2:6]
focal.covariates <- as.matrix(test.focal[,7])
generate.errors <- TRUE
bootstrap.samples <- 3

results_optimize <- cxr_optimize(fitness.model = fitness.model,
                                 optim.method = optim.method,
                                 param.list = param.list,
                                 log.fitness = log.fitness,
                                 init.lambda = init.lambda,
                                 lower.lambda = lower.lambda,
                                 upper.lambda = upper.lambda,
                                 init.sigma = init.sigma,
                                 lower.sigma = lower.sigma,
                                 upper.sigma = upper.sigma,
                                 init.alpha = init.alpha,
                                 lower.alpha = lower.alpha,
                                 upper.alpha = upper.alpha,
                                 init.lambda.cov = init.lambda.cov,
                                 lower.lambda.cov = lower.lambda.cov,
                                 upper.lambda.cov = upper.lambda.cov,
                                 init.alpha.cov = init.alpha.cov,
                                 lower.alpha.cov = lower.alpha.cov,
                                 upper.alpha.cov = upper.alpha.cov,
                                 focal.comp.matrix = focal.comp.matrix,
                                 focal.covariates = focal.covariates,
                                 generate.errors = generate.errors,
                                 bootstrap.samples = bootstrap.samples)

# test----
test_that("Expected classes", {
  expect_equal(class(results_optimize), "list")
  expect_equal(class(results_optimize$lambda), "numeric")
  expect_equal(class(results_optimize$lambda.lower.error), "numeric")
  expect_equal(class(results_optimize$lambda.upper.error), "numeric")
  expect_equal(class(results_optimize$sigma), "numeric")
  expect_equal(class(results_optimize$alpha), "numeric")
  expect_equal(class(results_optimize$alpha.lower.error), "numeric")
  expect_equal(class(results_optimize$alpha.upper.error), "numeric")
  expect_equal(class(results_optimize$lambda.cov), "numeric")
  expect_equal(class(results_optimize$lambda.cov.lower.error), "numeric")
  expect_equal(class(results_optimize$lambda.cov.upper.error), "numeric")
  expect_equal(class(results_optimize$alpha.cov), "numeric")
  expect_equal(class(results_optimize$alpha.cov.lower.error), "numeric")
  expect_equal(class(results_optimize$alpha.cov.upper.error), "numeric")
  expect_equal(class(results_optimize$log.likelihood), "numeric")
})

# effect-response function----

lambda.vector <- focal.lambda
e.vector <- runif(2,1,5)
r.vector <- runif(2,1,5)
lambda.cov <- lambda.cov.orig[1:2,]
e.cov <- matrix(runif(2,0.01,0.1),nrow = 2)
r.cov <- matrix(runif(2,0.01,0.1),nrow = 2)
sigma <- sd(log(test.focal$fitness))
lower.e <- 1e-5
upper.e <- 1e2
lower.r <- 1e-5
upper.r <- 1e2
lower.e.cov <- 1e-5
upper.e.cov <- 1e2
lower.r.cov <- 1e-5
upper.r.cov <- 1e2
optim.method <- "optim_L-BFGS-B"
sp.data <- test.data
sp.data$site <- rep((1:num.obs),2)
# sp.data <- sp.data[,c("site","focal","1","2","fitness")]
sp.data.long <- tidyr::gather(sp.data,key ="competitor",value = "number","1","2")
ER.covariates <- as.matrix(sp.data.long[,c("cov1")])
sp.data.long <- sp.data.long[,c("site","focal","fitness","competitor","number")]
generate.errors <- TRUE
bootstrap.samples <- 3

# first function
optimize.lambda <- TRUE
effect.response.model <- EffectResponse_lambda

results_ER <- ER_optimize(lambda.vector = lambda.vector,
                          e.vector = e.vector,
                          r.vector = r.vector,
                          lambda.cov = lambda.cov,
                          e.cov = e.cov,
                          r.cov = r.cov,
                          sigma = init.sigma,
                          lower.lambda = lower.lambda,
                          upper.lambda = upper.lambda,
                          lower.e = lower.e,
                          upper.e = upper.e,
                          lower.r = lower.r,
                          upper.r = upper.r,
                          lower.lambda.cov = lower.lambda.cov,
                          upper.lambda.cov = upper.lambda.cov,
                          lower.e.cov = lower.e.cov,
                          upper.e.cov = upper.e.cov,
                          lower.r.cov = lower.r.cov,
                          upper.r.cov = upper.r.cov,
                          lower.sigma = lower.sigma,
                          upper.sigma = upper.sigma,
                          effect.response.model = effect.response.model,
                          optim.method = optim.method,
                          sp.data = sp.data.long,
                          covariates = ER.covariates,
                          optimize.lambda = optimize.lambda,
                          generate.errors = generate.errors,
                          bootstrap.samples = bootstrap.samples)

# test----
test_that("Expected classes", {
  expect_equal(class(results_ER), "list")
  expect_equal(class(results_ER$lambda), "numeric")
  expect_equal(class(results_ER$lambda.lower.error), "numeric")
  expect_equal(class(results_ER$lambda.upper.error), "numeric")
  expect_equal(class(results_ER$response), "numeric")
  expect_equal(class(results_ER$response.lower.error), "numeric")
  expect_equal(class(results_ER$response.upper.error), "numeric")
  expect_equal(class(results_ER$effect), "numeric")
  expect_equal(class(results_ER$effect.lower.error), "numeric")
  expect_equal(class(results_ER$effect.upper.error), "numeric")
  expect_equal(class(results_ER$sigma), "numeric")
  expect_equal(class(results_ER$lambda.cov), "numeric")
  expect_equal(class(results_ER$lambda.cov.lower.error), "numeric")
  expect_equal(class(results_ER$lambda.cov.upper.error), "numeric")
  expect_equal(class(results_ER$response.cov), "numeric")
  expect_equal(class(results_ER$response.cov.lower.error), "numeric")
  expect_equal(class(results_ER$response.cov.upper.error), "numeric")
  expect_equal(class(results_ER$effect.cov), "numeric")
  expect_equal(class(results_ER$effect.cov.lower.error), "numeric")
  expect_equal(class(results_ER$effect.cov.upper.error), "numeric")
  expect_equal(class(results_ER$log.likelihood), "numeric")
})
