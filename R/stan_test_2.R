library(rstan)

########
# generate some data
focal.sp <- c(1)
num.sp <- 5
num.cov <- 2
num.obs <- 25 # per focal species

focal.lambda <- c(10)
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

###################
# stan model

my.model <- "
data{
  
  int<lower = 0> n; // number of observations

  vector[n] n_comp; //number of competitors per observation
  
  vector[n] fitness; // vector of fitness (log??)
  
}

parameters{
  
  real<lower = 0> lambda; 
  
  real<lower = 0.0001, upper = 1> alpha; //single alpha
  
}
transformed parameters{
  
  vector[n] Fhat; //computed fitness
  
  Fhat = (lambda ./ (1+(alpha * n_comp))); // model 2
  
}
model{
  
  lambda ~ normal(mean(fitness),sd(fitness)); // prior expectation and distribution for lambda
  alpha ~ uniform(0.0001,1);
}
"
############
# compile stan model
stan.model <- stan_model(model_code = my.model) 

############
# initial data
fitness <- test.data$fitness#log(test.data$fitness)
N <- nrow(test.data)
n_comp <- rowSums(test.data[,2:(num.sp+1)])
dat <- list(fitness = fitness, n = N, n_comp = n_comp); 

############
# fit the model
fit <- sampling(stan.model, data = dat, iter = 1000, warmup=200) 
fit
