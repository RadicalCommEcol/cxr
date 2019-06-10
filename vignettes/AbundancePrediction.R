source("./R/GenerateTestData.R")
source("./R/BevertonHolt_abundance_models.R")
source("./R/PredictAbundances.R")

# timesteps
timesteps <- 50

# test data
focal.sp <- c(1,2)
num.sp <- 5
num.cov <- 2
num.obs <- 3 # per focal species

focal.lambda <- c(100,200)
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
init.abund <- expand.grid(1:num.obs,1:num.sp)
names(init.abund) <- c("site","species")
init.abund$abundance <- rnorm(nrow(init.abund),100,50)

# lambda, s, g
# values for not focal species are given zero
# ideally we will have a complete parameterization of all species present
# otherwise, the dynamics will only predict those parameterized,
# and the effect of other competitors will dissappear from time 2 onwards
sp.par <- data.frame(species = 1:num.sp,lambda = 0,germ.rate = 0, survival.rate = 0)
sp.par$lambda[focal.sp] <- focal.lambda
sp.par$germ.rate[focal.sp] <- runif(length(focal.sp),0,0.5)
sp.par$survival.rate[focal.sp] <- runif(length(focal.sp),0,0.5)

# environmental heterogeneity
cov.time <- expand.grid(1:num.obs,1:timesteps,1:num.cov)
names(cov.time) <- c("site","timestep","covariate")
cov.time$value <- runif(nrow(cov.time),0,10)

# alpha matrix
alpha.matrix <- alpha.matrix.orig

# effect of covariates on alpha and lambda
lambda.cov.matrix <- matrix(unlist(lambda.cov.orig),nrow = num.sp)

alpha.cov.matrix <- list()
for(i.cov in 1:num.cov){
  alpha.cov.matrix[[i.cov]] <- matrix(nrow = num.sp,ncol = num.sp)
  for(i.sp in 1:num.sp){
    alpha.cov.matrix[[i.cov]][i.sp,] <- alpha.cov.orig[[i.sp]][,i.cov]
  }
}

par <- list(sp.par = sp.par, initial.values = init.abund, 
            covariates = cov.time, other.par = list(alpha.matrix = alpha.matrix, 
                                                    lambda.cov.matrix = lambda.cov.matrix, 
                                                    alpha.cov.matrix = alpha.cov.matrix))

abundance.model <- BH_abundance_5
predicted.abundances <- PredictAbundances(par = par,timesteps = timesteps,abundance.model = abundance.model)
predicted.abundances$timestep <- as.factor(predicted.abundances$timestep)
predicted.abundances$site <- as.factor(predicted.abundances$site)
predicted.abundances$sp <- as.factor(predicted.abundances$sp)

abund.plot <- ggplot2::ggplot(predicted.abundances,aes(x = timestep,y = abundance, group = sp)) + 
  ggplot2::geom_line(aes(color = sp)) + 
  ggplot2::facet_grid(site~.)+
  # ylim(-10,10)+
  NULL
abund.plot

