

#####################

# test1 - no effect of covariates on lambdas and alphas (models 1-3)

test1 <- compete(focal = test.data$focal, 
                fitness = test.data$fitness, 
                comp_matrix = comp_matrix
                , lower1 = c(1,0.0001)
                , lower2 = c(1,0,0.0001)
                , lower3 = c(1, rep(0, times=ncol(comp_matrix)),0.0000000001)
                , hessian3 = TRUE) #works

# test2 - full set of models (1-5)

test2 <- compete(test.data$focal, test.data$fitness, comp_matrix, covariates#,method = "Nelder-Mead"
                     , lower1 = c(1,0.0001)
                     , lower2 = c(1,0,0.0001)
                     , lower3 = c(1, rep(0, times=ncol(comp_matrix)),0.0000000001)
                 # , hessian3 = F, hessian4 = F
                     , lower4 = c(1, rep(0.001, times=ncol(covariates)), #n_cov
                                  rep(0, times=ncol(comp_matrix)), #alfas
                                  rep(0.001, times=ncol(covariates)), #n_cov
                                  0.0000000001)
                     , lower5 = c(1, rep(0.001, times=ncol(covariates)), #n_cov
                                  rep(0, times=ncol(comp_matrix)), #alfas
                                  rep(0.001, times=(ncol(covariates)*ncol(comp_matrix))), #n_cov*n_bg
                                  0.0000000001)
                     , hessian5 = F)

########## alternative hessian calculation, etc

focal <- test.data$focal
fitness <- test.data$fitness
log_fitness <- log(fitness)
comp_matrix_i <- comp_matrix
n_cov <- ncol(covariates) 
covariates_i <- covariates
alphas <- test1$alpha_matrix3
sigma <- test1$`lambdas$co`$sigma_est3
method <- "L-BFGS-B"
n_bg <- dim(comp_matrix_i)[2]

# model 4!
par <- c(487,  #lambda 1
          rep(0.0001, times = n_cov), #l_cov n_cov
          alphas,#testcomp3$par[2:(n_bg+1)], #alfas n_bg
          rep(0.0001, times = n_cov), #a_cov n_cov
          sigma)#testcomp3$par[length(par3)]) # sigma 1

lower4 <- c(1, rep(0.001, times=ncol(covariates)), #n_cov
            rep(0, times=ncol(comp_matrix)), #alfas
            rep(0.001, times=ncol(covariates)), #n_cov
            0.0000000001)
upper4 <- rep(1e5,length(lower4))

compmodel4 <- function(par, log_fitness, comp_matrix_i, n_cov, n_bg, covariates_i){
  lambda <- par[1] 
  l_cov <- par[(1+1):(1+n_cov)] #effect of cov 1, 2, ... on lambda
  a_comp <- par[(1+n_cov+1):(1+n_cov+n_bg)] # alfas_ij
  a_cov <- par[(1+n_cov+n_bg+1):(1+n_cov+n_bg+n_cov)] # common effect of cov 1, 2... on alphas
  sigma <- par[length(par)]
  num = 1
  for(z in 1:n_cov){
    num <- num + l_cov[z]*covariates_i[,z] 
  }
  cov_term <- 0 
  for(v in 1:n_cov){
    cov_term <- cov_term + a_cov[v] * covariates_i[,v]
  }
  term <- 1 #create the denominator term for the model
  for(z in 1:ncol(comp_matrix_i)){
    term <- term + (a_comp[z] + cov_term) * comp_matrix_i[,z] 
  }
  pred <- lambda * (num) / term 
  # likelihood as before:
  llik<-dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  testcomp4 <- optim(par, 
                     compmodel4, 
                     gr = NULL, 
                     method = method, 
                     lower = lower4, 
                     # upper = Inf,
                     control = list(), 
                     hessian = F,
                     log_fitness = log_fitness, 
                     comp_matrix_i = comp_matrix_i,
                     n_cov = n_cov, 
                     n_bg = n_bg, 
                     covariates_i = covariates_i) #check hessian can be done with all methods.
  par4 <- testcomp4$par
  
  if(testcomp4$convergence == 0){
    message(paste(splist[i], "model 4 converged on rep", k, sep = " "))
    break
  }
}

#################################
# methods comparison



#########################
# other optim functions
# genSA
library(GenSA)

par4.gensa <- GenSA(par = par,fn = compmodel4,lower = lower4,upper = upper4, control = list(maxit = 25), 
                    log_fitness = log_fitness, 
                    comp_matrix_i = comp_matrix_i,
                    n_cov = n_cov, 
                    n_bg = n_bg, 
                    covariates_i = covariates_i)

# nloptr: algorithm must be restricted to derivative-free ones, the other options need a gradient function
#"NLOPT_GN_CRS2_LM"
# NLOPT_GN_ISRES
# NLOPT_GN_DIRECT_L
# NLOPT_GN_DIRECT
library(nloptr)
par4.nloptr <- nloptr(x0 = par,eval_f = compmodel4,opts = list("algorithm"="NLOPT_GN_ISRES"),
                      lb = lower4,
                      ub = upper4,
                      log_fitness = log_fitness, 
                      comp_matrix_i = comp_matrix_i,
                      n_cov = n_cov, 
                      n_bg = n_bg, 
                      covariates_i = covariates_i)

# psoptim no
# DEoptim
library(DEoptimR)
par4.deoptimr <- DEoptimR::JDEoptim(lower = lower4,upper = upper4,fn = compmodel4,
                                    log_fitness = log_fitness, 
                                    comp_matrix_i = comp_matrix_i,
                                    n_cov = n_cov, 
                                    n_bg = n_bg, 
                                    covariates_i = covariates_i)

# hydroPSO
library(hydroPSO)
par4.hydroPSO <- hydroPSO::hydroPSO(par = par,fn = compmodel4,lower = lower4,upper = upper4,
                                    log_fitness = log_fitness, 
                                    comp_matrix_i = comp_matrix_i,
                                    n_cov = n_cov, 
                                    n_bg = n_bg, 
                                    covariates_i = covariates_i)

##########################
# hessian
library(numDeriv)
hessian.par <- par4.deoptimr$par
my.hessian <- hessian(compmodel4,x = hessian.par,log_fitness = log_fitness, 
                      comp_matrix_i = comp_matrix_i,
                      n_cov = n_cov, 
                      n_bg = n_bg, 
                      covariates_i = covariates_i)

#####################
# diagnostic plots
library(tidyverse)

# predicted.values <- test1$`lambdas$co`
# 
# # this value cannot be compared to a real one, as it merges all effects into one alpha
# 
# alpha2.plot <- ggplot(predicted.values, aes(x = splist, y = alpha_est2)) + 
#   geom_errorbar(aes(ymin=alpha.lower.error.2, ymax = alpha.upper.error.2), width = 0.05, size = 0.3) +
#   geom_point(size = 2) + 
#   ggtitle("alpha from model 2") + 
#   NULL
# alpha2.plot

##########
# models 3,4,5
# convert to tidy format
alpha.values <- NULL
for(i.model in c(3,4,5)){
  
  if(i.model == 3){
    my.predicted.values <- test2$alpha_matrix3
    my.upper.errors <- test2$alpha_upper_error3
    my.lower.errors <- test2$alpha_lower_error3
  }else if(i.model == 4){
    my.predicted.values <- test2$alpha_matrix4
    my.upper.errors <- test2$alpha_upper_error4
    my.lower.errors <- test2$alpha_lower_error4
  }else if(i.model == 5){
    my.predicted.values <- test2$alpha_matrix5
    my.upper.errors <- test2$alpha_upper_error5
    my.lower.errors <- test2$alpha_lower_error5
  }
  
  predicted.values <- as.data.frame(my.predicted.values)
  predicted.values <- gather(predicted.values,key = "sp",value = "alpha")
  
  upper.errors <- as.data.frame(my.upper.errors)
  upper.errors <- gather(upper.errors,key = "sp",value = "alpha.upper.error")
  
  lower.errors <- as.data.frame(my.lower.errors)
  lower.errors <- gather(lower.errors,key = "sp",value = "alpha.lower.error")
  
  temp.alpha.values <- left_join(predicted.values,upper.errors)
  temp.alpha.values <- left_join(temp.alpha.values,lower.errors)
  
  temp.alpha.values$model <- i.model
  
  alpha.values <- rbind(alpha.values, temp.alpha.values)
}

alpha.values$focal.sp <- 1
alpha.values$model <- as.factor(alpha.values$model)

# append true values
true.values <- data.frame(sp = c("sp1","sp2","sp3"), alpha = c(a_sp1,a_sp2,a_sp3), alpha.upper.error = 0, alpha.lower.error = 0,
                          model = "truth", focal.sp = 1)

alpha.values <- rbind(alpha.values, true.values)

dodge <- position_dodge(.3)
alpha.plot <- ggplot(alpha.values, aes(x = sp, y = alpha, group = model)) + 
  geom_errorbar(aes(ymin=alpha.lower.error, ymax = alpha.upper.error), width = 0.05, size = 0.3, position = dodge) +
  geom_point(aes(color = model), size = 2, position = dodge) + 
  facet_grid(focal.sp~.) + 
  NULL
alpha.plot


