####
# standard error estimates from bootstrap samples

SEbootstrap <- function(fitness.model,
                        optim.method,
                        param.list,
                        fixed.terms,
                        log.fitness,
                        init.par,
                        lower.bounds,
                        upper.bounds,
                        focal.comp.matrix,
                        focal.covariates,
                        nsamples){
  
  if(nsamples<2){
    print("SEbootstrap: number of bootstrap samples cannot be < 2. Setting bootstrap samples to 2.")
    nsamples <- 2
  }
  
  num.competitors <- dim(focal.comp.matrix)[2]
  num.covariates <- ifelse(is.null(ncol(focal.covariates)),0,ncol(focal.covariates))
  
  boot.results <- matrix(nrow = nsamples, ncol = length(init.par))
  
  for(i.sample in 1:nsamples){
    
    my.sample <- sample(length(log.fitness),length(log.fitness),replace = T)
    
    # sample fitness, competition matrix, and covariates matrix
    boot.fitness <- log.fitness[my.sample]
    boot.comp.matrix <- focal.comp.matrix[my.sample,]
    # boot.covariates <- ifelse(is.data.frame(focal.covariates),focal.covariates[my.sample,],0)
    if(is.data.frame(focal.covariates)){
      boot.covariates <- as.data.frame(focal.covariates[my.sample,])
    }else{
      boot.covariates <- 0
    }
    
    ############
    if(optim.method == "optim_NM"){
      
      my.boot.par <- optim(init.par, 
                           fitness.model, 
                           gr = NULL, 
                           method = "Nelder-Mead", 
                           # lower = lower.bounds,
                           # upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           param.list = param.list,
                           log.fitness = boot.fitness, 
                           focal.comp.matrix = boot.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = boot.covariates,
                           fixed.terms = fixed.terms)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "optim_L-BGFS-B"){
      
      my.boot.par <- optim(init.par, 
                           fitness.model, 
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
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_ISRES"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_DIRECT_L_RAND"){
      
      my.boot.par <- nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = focal.covariates)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "GenSA"){
      
      my.boot.par <- GenSA(par = init.par,fn = fitness.model,
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
      my.boot.par <- hydroPSO::hydroPSO(par = init.par,fn = fitness.model,
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
                                        fn = fitness.model,
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


