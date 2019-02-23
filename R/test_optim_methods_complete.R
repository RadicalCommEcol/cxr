source("R/generate_test_data.R")
# source("R/compete.R")
source("R/nested_models.R")
# hessian calc
library(pracma)

# optimization methods
library(nloptr)
library(GenSA)
library(hydroPSO)
library(DEoptimR)

# factors to loop over
focal.sp <- unique(test.data$focal)
optim.methods <- c("optim_NM", "optim_L-BGFS-B","nloptr_CRS2_LM", "nloptr_ISRES", "nloptr_DIRECT_L_RAND", "GenSA", "hydroPSO", "DEoptimR")
models <- 1:5

# initialize data structures
# lambdas and sigmas will be placed in the dataframe "lambda.results"
# others are matrices/vectors within nested lists, of the form matrix[[focal.sp]][[model]][[method]]

lambda.results <- NULL

# general purpose list
temp.list <- list()
for(i.sp in 1:length(focal.sp)){
  temp.list[[i.sp]] <- list()
  for(i.model in 1:length(models)){
    temp.list[[i.sp]][[i.model]] <- list()
    for(i.method in 1:length(optim.methods)){
      temp.list[[i.sp]][[i.model]][[i.method]] <- 0
    }
  }
}

# initialize the rest
alpha.matrix <- temp.list
alpha.upper.error.matrix <- temp.list
alpha.lower.error.matrix <- temp.list
lambda.cov.matrix <- temp.list
lambda.cov.upper.error.matrix <- temp.list
lambda.cov.lower.error.matrix <- temp.list
alpha.cov.matrix <- temp.list
alpha.cov.upper.error.matrix <- temp.list
alpha.cov.lower.error.matrix <- temp.list

# i.sp <- 1
# i.model <- 1
# i.method <- 1

# main loop
for(i.sp in 1:length(focal.sp)){
  
  # subset and prepare the data
  
  focal.sp.data <- subset(test.data, focal == focal.sp[i.sp])
  # current focal species
  focal <- unique(focal.sp.data$focal)
  # fitness metric
  fitness = focal.sp.data$fitness
  log_fitness <- log(fitness)
  # competition matrix: number of competitors
  focal.comp.matrix <- comp_matrix[which(test.data$focal == focal.sp[i.sp]),]
  # competitors at each observation
  background <- rowSums(focal.comp.matrix)
  # number of competitors
  num.competitors <- dim(focal.comp.matrix)[2]
  # number of covariates
  num.covariates <- ncol(covariates)
  # covariates for the focal species
  focal.covariates <- covariates[which(test.data$focal == focal.sp[i.sp]),,drop = FALSE]

  # model to optimize  
  for(i.model in models){
    
    # initial data should be taken from the previous model
    # except for model1

    if(models[i.model] == 1){
      my.model <- model1
      
      init.lambda <- mean(log_fitness)
      init.sigma <- sd(log_fitness)

      init.par <- c(init.lambda, init.sigma)
      lower.bounds <- c(1, # lambda 
                        0.0000000001) # sigma
      upper.bounds <- rep(1e5,length(lower.bounds))
      
    }else if(models[i.model] == 2){
      my.model <- model2
      
      init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp[i.sp] & 
                                             lambda.results$model == 1 & 
                                             lambda.results$optim.method == "DEoptimR"]
      init.alphas <- 0.0001
      init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp[i.sp] & 
                                           lambda.results$model == 1 & 
                                           lambda.results$optim.method == "DEoptimR"]
      
      init.par <- c(init.lambda,init.alphas,init.sigma)
      lower.bounds <- c(1, # lambda 
                        0,
                        0.0000000001) # sigma
      upper.bounds <- rep(1e5,length(lower.bounds))
      
    }else if(models[i.model] == 3){
      my.model <- model3
      
      init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp[i.sp] & 
                                             lambda.results$model == 2 & 
                                             lambda.results$optim.method == "DEoptimR"]
      init.alphas <- rep(alpha.matrix[[focal.sp[i.sp]]][[2]][[8]], times=ncol(comp_matrix))
      init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp[i.sp] & 
                                           lambda.results$model == 2 & 
                                           lambda.results$optim.method == "DEoptimR"]
      
      init.par <- c(init.lambda,init.alphas,init.sigma)
      lower.bounds <- c(1, # lambda 
                        rep(0, times=ncol(comp_matrix)),
                        0.0000000001) # sigma
      upper.bounds <- rep(1e5,length(lower.bounds))
      
    }else if(models[i.model] == 4){
      my.model <- model4
      
      init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp[i.sp] & 
                                             lambda.results$model == 3 & 
                                             lambda.results$optim.method == "DEoptimR"]
      init.alphas <- alpha.matrix[[focal.sp[i.sp]]][[3]][[8]]
      init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp[i.sp] & 
                                           lambda.results$model == 3 & 
                                           lambda.results$optim.method == "DEoptimR"]
      init.l_cov <- rep(0.001, times = num.covariates)
      init.a_cov <- rep(0.001, times = num.covariates)
      
      init.par <- c(init.lambda,init.l_cov,init.alphas,init.a_cov,init.sigma)
      lower.bounds <- c(1, rep(0.001, times=ncol(covariates)), #n_cov
                        rep(0, times=ncol(comp_matrix)), #alfas
                        rep(0.0001, times=ncol(covariates)), #n_cov
                        0.0000000001)
      upper.bounds <- rep(1e5,length(lower.bounds))
      
    }else if(models[i.model] == 5){
      my.model <- model5
      
      init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp[i.sp] & 
                                             lambda.results$model == 4 & 
                                             lambda.results$optim.method == "DEoptimR"]
      init.alphas <- alpha.matrix[[focal.sp[i.sp]]][[4]][[8]]
      init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp[i.sp] & 
                                           lambda.results$model == 4 & 
                                           lambda.results$optim.method == "DEoptimR"]
      
      init.l_cov <- lambda.cov.matrix[[focal.sp[i.sp]]][[4]][[8]]
      a_cov4 <- alpha.cov.matrix[[focal.sp[i.sp]]][[4]][[8]]
      
      init.a_cov <- c()
      for(w in 1:length(a_cov4)){
        init.a_cov <- c(init.a_cov, rep(a_cov4[w], times = num.competitors))
      }
      
      init.par <- c(init.lambda,init.l_cov,init.alphas,init.a_cov,init.sigma)
      lower.bounds <- c(1, rep(0.001, times=ncol(covariates)), #n_cov
                        rep(0, times=ncol(comp_matrix)), #alfas
                        rep(0.0001, times=(ncol(covariates)*ncol(comp_matrix))), #n_cov*n_bg
                        0.0000000001)
      upper.bounds <- rep(1e5,length(lower.bounds))
      
    }
    
    ######################
    # compute each method

    for(i.method in 1:length(optim.methods)){
      
      # optimization
      optim.par <- list(par = rep(NA,length(init.par)), value = NA)
      temp.results <- data.frame(focal.sp = focal.sp[i.sp],
                                 model = models[i.model],
                                 optim.method = optim.methods[i.method],
                                 lambda = NA,
                                 lambda.lower.error = NA,
                                 lambda.upper.error = NA,
                                 sigma = NA,
                                 log.likelihood = NA)
      
      if(optim.methods[i.method] == "optim_NM"){
        
        optim.par <- optim(init.par, 
                           my.model, 
                           gr = NULL, 
                           method = "Nelder-Mead", 
                           # lower = lower.bounds,
                           # upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           log_fitness = log_fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates)
        
      }else if(optim.methods[i.method] == "optim_L-BGFS-B"){
        
        optim.par <- optim(init.par, 
                           my.model, 
                           gr = NULL, 
                           method = "L-BFGS-B", 
                           lower = lower.bounds, 
                           upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           log_fitness = log_fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates)
        
      }else if(optim.methods[i.method] == "nloptr_CRS2_LM"){
        
        optim.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log_fitness = log_fitness, 
                            focal.comp.matrix = focal.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
        
      }else if(optim.methods[i.method] == "nloptr_ISRES"){
        
        optim.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log_fitness = log_fitness, 
                            focal.comp.matrix = focal.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
        
      }else if(optim.methods[i.method] == "nloptr_DIRECT_L_RAND"){
        
        optim.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log_fitness = log_fitness, 
                            focal.comp.matrix = focal.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
        
      }else if(optim.methods[i.method] == "GenSA"){
        
        optim.par <- GenSA(par = init.par,fn = my.model,lower = lower.bounds,upper = upper.bounds, control = list(maxit = 1e3), 
                           log_fitness = log_fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates)
        
      }else if(optim.methods[i.method] == "hydroPSO"){
        
        optim.par <- hydroPSO::hydroPSO(par = init.par,fn = my.model,lower = lower.bounds,upper = upper.bounds, control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                        log_fitness = log_fitness, 
                                        focal.comp.matrix = focal.comp.matrix,
                                        num.covariates = num.covariates, 
                                        num.competitors = num.competitors, 
                                        focal.covariates = focal.covariates)
        
      }else if(optim.methods[i.method] == "DEoptimR"){
        
        optim.par <- DEoptimR::JDEoptim(lower = lower.bounds,upper = upper.bounds,fn = my.model,
                                        log_fitness = log_fitness, 
                                        focal.comp.matrix = focal.comp.matrix,
                                        num.covariates = num.covariates, 
                                        num.competitors = num.competitors, 
                                        focal.covariates = focal.covariates)
      }
      
      ##################################
      # tidy the output from the method
      # if-else the method outputs optim-like values
      if(optim.methods[i.method] %in% c("optim_NM","optim_L-BGFS-B","DEoptimR","hydroPSO","GenSA")){
        
        temp.results$lambda <- optim.par$par[1]
        temp.results$sigma <- optim.par$par[length(optim.par$par)]
        temp.results$log.likelihood <- optim.par$value
        temp.alpha.matrix <- NA
        temp.lambda.cov.matrix <- NA
        temp.alpha.cov.matrix <- NA
        temp.alpha.lower.error <- NA 
        temp.alpha.upper.error <- NA
        
        # obtain the different matrices
        if(models[i.model] == 2){
          temp.alpha.matrix <- optim.par$par[2]
        }else if(models[i.model] == 3){
          temp.alpha.matrix <- optim.par$par[(1+num.covariates+1):(1+num.covariates+num.competitors)]
        }else if(models[i.model] == 4){
          temp.alpha.matrix <- optim.par$par[(1+num.covariates+1):(1+num.covariates+num.competitors)]
          temp.lambda.cov.matrix <- optim.par$par[(1+1):(1+num.covariates)]
          temp.alpha.cov.matrix <- optim.par$par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+num.covariates)]
        }else if(models[i.model] == 5){
          temp.alpha.matrix <- optim.par$par[(1+num.covariates+1):(1+num.covariates+num.competitors)]
          temp.lambda.cov.matrix <- optim.par$par[(1+1):(1+num.covariates)]
          temp.alpha.cov.matrix <- optim.par$par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.covariates*num.competitors))] #effects of cov 1, 2... on alpha_i
        }
        my.par <- optim.par$par
        
      }else{ # methods with different nomenclature
        temp.results$lambda <- optim.par$solution[1]
        temp.results$sigma <- optim.par$solution[length(optim.par$solution)]
        temp.results$log.likelihood <- optim.par$objective
        
        # obtain the different matrices
        if(models[i.model] > 2){
          temp.alpha.matrix <- optim.par$solution[(1+num.covariates+1):(1+num.covariates+num.competitors)]
          if(models[i.model] == 4){
            temp.lambda.cov.matrix <- optim.par$solution[(1+1):(1+num.covariates)]
            temp.alpha.cov.matrix <- optim.par$solution[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+num.covariates)]
          }else if(models[i.model] == 5){
            temp.lambda.cov.matrix <- optim.par$solution[(1+1):(1+num.covariates)]
            temp.alpha.cov.matrix <- optim.par$solution[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.covariates*num.competitors))] #effects of cov 1, 2... on alpha_i
          }
        }
        my.par <- optim.par$solution
      }  
      
      #####################
      # hessian calculation
      hessian.par <- numDeriv::hessian(my.model,my.par,
                                       log_fitness = log_fitness,
                                       focal.comp.matrix = focal.comp.matrix,
                                       num.covariates = num.covariates,
                                       num.competitors = num.competitors,
                                       focal.covariates = focal.covariates)
      # errors
      inverse <- solve(hessian.par)
      errors <- sqrt(diag(inverse))
      
      # save estimates errors, if hessian does not fail
      if(sum(is.na(errors)) == 0){
        
        temp.results$lambda.lower.error <- my.par[1]-1.96*errors[1]
        temp.results$lambda.upper.error <- my.par[1]+1.96*errors[1]
        
        # l_cov_error5 <- c(my.par[2:(1+num.covariates)]-1.96*errors[2:(1+num.covariates)],
        #                       my.par[2:(1+num.covariates)]+1.96*errors[2:(1+num.covariates)])
        # temp.alpha.upper.error <- my.par[(1+num.covariates+1):(1+num.covariates+num.competitors)]-1.96*errors[(1+num.covariates+1):(1+num.covariates+num.competitors)]
        # temp.alpha.lower.error <- my.par[(1+num.covariates+1):(1+num.covariates+num.competitors)]+1.96*errors[(1+num.covariates+1):(1+num.covariates+num.competitors)]
        # temp.alpha.cov <- c(my.par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.competitors*num.covariates))]-1.96*errors[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.competitors*num.covariates))],
        #                       my.par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.competitors*num.covariates))]+1.96*errors[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.competitors*num.covariates))])
      }
      ###############
      # merge results
      
      lambda.results <- rbind(lambda.results,temp.results)
      alpha.matrix[[i.sp]][[i.model]][[i.method]] <- temp.alpha.matrix
      alpha.upper.error.matrix[[i.sp]][[i.model]][[i.method]] <- temp.alpha.upper.error
      alpha.lower.error.matrix[[i.sp]][[i.model]][[i.method]] <- temp.alpha.lower.error
      
      lambda.cov.matrix[[i.sp]][[i.model]][[i.method]] <- temp.lambda.cov.matrix
      # lambda.cov.upper.error.matrix <- temp.list
      # lambda.cov.lower.error.matrix <- temp.list
      alpha.cov.matrix[[i.sp]][[i.model]][[i.method]] <- temp.alpha.cov.matrix
      # alpha.cov.upper.error.matrix <- temp.list
      # alpha.cov.lower.error.matrix <- temp.list
      
    }# for i.method
  }# for i.model
}# for i.sp

