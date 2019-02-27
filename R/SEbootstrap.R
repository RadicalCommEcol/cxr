####
# standard error estimates from bootstrap samples

SEbootstrap <- function(optim.method,
                              fitness.model,
                              lower.bounds,
                              upper.bounds,
                              init.data,
                              init.par,
                              num.sp,
                              num.cov,
                              nsamples){
  
  boot.results <- matrix(nrow = nsamples, ncol =  length(init.par))
  
  for(i.sample in 1:nsamples){
    
    my.sample <- sample(nrow(init.data),nrow(init.data),replace = T)
    
    boot.data <- init.data[my.sample,]
    
    boot.fitness <- boot.data$fitness
    boot.log.fitness <- log(boot.fitness)
    boot.comp.matrix <- boot.data[,2:(num.sp+1)]
    boot.covariates <- boot.data[,(num.sp+2):(num.sp+2+num.cov-1), drop = FALSE]
    
    if(optim.method == "optim_NM"){
      
      my.boot.par <- optim(init.par, 
                           my.model, 
                           gr = NULL, 
                           method = "Nelder-Mead", 
                           # lower = lower.bounds,
                           # upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           log_fitness = boot.log.fitness, 
                           focal.comp.matrix = boot.comp.matrix,
                           num.covariates = num.cov, 
                           num.competitors = num.sp, 
                           focal.covariates = boot.covariates)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "optim_L-BGFS-B"){
      
      my.boot.par <- optim(init.par, 
                         my.model, 
                         gr = NULL, 
                         method = "L-BFGS-B", 
                         lower = lower.bounds, 
                         upper = upper.bounds,
                         control = list(), 
                         hessian = F,
                         log_fitness = boot.log.fitness, 
                         focal.comp.matrix = boot.comp.matrix,
                         num.covariates = num.covariates, 
                         num.competitors = num.competitors, 
                         focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "nloptr_CRS2_LM"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          log_fitness = boot.log.fitness, 
                          focal.comp.matrix = boot.comp.matrix,
                          num.covariates = num.covariates, 
                          num.competitors = num.competitors, 
                          focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_ISRES"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          log_fitness = boot.log.fitness, 
                          focal.comp.matrix = boot.comp.matrix,
                          num.covariates = num.covariates, 
                          num.competitors = num.competitors, 
                          focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_DIRECT_L_RAND"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          log_fitness = boot.log.fitness, 
                          focal.comp.matrix = boot.comp.matrix,
                          num.covariates = num.covariates, 
                          num.competitors = num.competitors, 
                          focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "GenSA"){
      
      my.boot.par <- GenSA(par = init.par,fn = my.model,
                         lower = lower.bounds,
                         upper = upper.bounds, 
                         control = list(maxit = 1e2), 
                         log_fitness = boot.log.fitness, 
                         focal.comp.matrix = boot.comp.matrix,
                         num.covariates = num.covariates, 
                         num.competitors = num.competitors, 
                         focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "hydroPSO"){
      
      # suppress annoying output
      sink("/dev/null")
      my.boot.par <- hydroPSO::hydroPSO(par = init.par,fn = my.model,
                                      lower = lower.bounds,
                                      upper = upper.bounds, 
                                      control=list(write2disk=FALSE, maxit = 1e2, MinMax = "min", verbose = F),
                                      log_fitness = boot.log.fitness, 
                                      focal.comp.matrix = boot.comp.matrix,
                                      num.covariates = num.covariates, 
                                      num.competitors = num.competitors, 
                                      focal.covariates = focal.covariates)
      sink()
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "DEoptimR"){
      
      my.boot.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                      upper = upper.bounds,
                                      maxiter = 100,
                                      fn = my.model,
                                      log_fitness = boot.log.fitness, 
                                      focal.comp.matrix = boot.comp.matrix,
                                      num.covariates = num.covariates, 
                                      num.competitors = num.competitors, 
                                      focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$par
      
    }
    
    boot.results[i.sample,] <- my.boot.par
    
  }  
  
  boot.se <- apply(boot.results,2,sd)
  boot.se
}


