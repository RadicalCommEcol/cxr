####
# standard error estimates from bootstrap samples

SEbootstrap <- function(optim.method,
                        fitness.model,
                        lower.bounds,
                        upper.bounds,
                        init.par,
                        log.fitness,
                        focal.comp.matrix,
                        focal.covariates,
                        num.competitors,
                        num.covariates,
                        nsamples){
  
  boot.results <- matrix(nrow = nsamples, ncol = length(init.par))
  
  for(i.sample in 1:nsamples){
    
    my.sample <- sample(length(log.fitness),length(log.fitness),replace = T)
    
    # sample fitness, competition matrix, and covariates matrix
    boot.fitness <- log.fitness[my.sample]
    boot.comp.matrix <- focal.comp.matrix[my.sample,]
    boot.covariates <- ifelse(is.data.frame(focal.covariates),focal.covariates[my.sample,],0)
    
    num.competitors <- dim(focal.comp.matrix)[2]
    num.covariates <- ifelse(is.null(ncol(focal.covariates)),0,ncol(focal.covariates))
    
    # boot.data <- init.data[my.sample,]
    
    ############
    # TODO: include fitness as a param, independent, so that I will not
    # need to get it from the dataframe
    # boot.fitness <- boot.data$seed
    # ############
    # 
    # boot.fitness <- log(boot.fitness)
    
    ############
    # TODO: same with the competitors and covariates matrix
    # boot.comp.matrix <- as.matrix(boot.data[,10:(num.competitors+9)])
    # 
    # if(num.covariates > 0){
    #   boot.covariates <- boot.data[,(num.competitors+10):(num.competitors+10+num.covariates-1), drop = FALSE]
    # }else{
    #   boot.covariates <- 0
    # }
    
    ############
    if(optim.method == "optim_NM"){
      
      my.boot.par <- optim(init.par, 
                           my.model, 
                           gr = NULL, 
                           method = "Nelder-Mead", 
                           # lower = lower.bounds,
                           # upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           log.fitness = boot.fitness, 
                           focal.comp.matrix = boot.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
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
                           log.fitness = boot.fitness, 
                           focal.comp.matrix = boot.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "nloptr_CRS2_LM"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_ISRES"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_DIRECT_L_RAND"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = my.model,opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log.fitness = boot.fitness, 
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
                           log.fitness = boot.fitness, 
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
                                        log.fitness = boot.fitness, 
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
                                        log.fitness = boot.fitness, 
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


