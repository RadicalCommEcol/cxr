
# timesteps
timesteps <- 100

# test data
focal.sp <- c(1,2)
num.sp <- 5
num.cov <- 2
num.obs <- 25 # per focal species

focal.lambda <- c(10,20)
alpha.matrix.orig <- matrix(data = runif(num.sp*num.sp,-0.001,0),nrow = num.sp, ncol = num.sp)
alpha.cov.orig <- list()
lambda.cov.orig <- list()

for(i.sp in 1:num.sp){
  alpha.cov.orig[[i.sp]] <- matrix(data = rnorm(num.sp*num.cov,0,0.005),nrow = num.sp, ncol = num.cov) # rows: competitors, columns: covariates
  lambda.cov.orig[[i.sp]] <- rnorm(num.cov,0,0.005)#rep(0,num.cov)
}

test.data <- GenerateTestData(focal.sp = focal.sp,
                              num.sp = num.sp,
                              num.cov = num.cov,
                              num.obs = num.obs,
                              fitness.model = 5,
                              focal.lambda = focal.lambda,
                              alpha.matrix = alpha.matrix.orig,
                              alpha.cov = alpha.cov.orig,
                              lambda.cov = lambda.cov.orig)

# initial abundances
init.abund <- expand.grid(1:num.obs,focal.sp)
names(init.abund) <- c("site","sp")
init.abund$abundance <- rnorm(nrow(init.abund),100,50)

# lambda, s, g
sp.par <- data.frame(sp = focal.sp,lambda = focal.lambda,
                     germ.rate = runif(length(focal.sp),0,0.5),
                     survival.rate = runif(length(focal.sp),0,0.5))

# environmental heterogeneity
cov.time <- expand.grid(1:num.obs,1:timesteps,1:num.cov)
names(cov.time) <- c("site","timestep","covariate")
cov.time$value <- runif(nrow(cov.time),0,10)

PredictAbundances <- function(par,timesteps,abundance.model){
  
}

