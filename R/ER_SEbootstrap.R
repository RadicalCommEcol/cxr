####
# standard error estimates from bootstrap samples

ER_SEbootstrap <- function(optim.method,
                        effect.response.model,
                        lower.bounds,
                        upper.bounds,
                        init.par,
                        sp.data,
                        nsamples){
  
  boot.results <- matrix(nrow = nsamples, ncol =  length(init.par))
  
  for(i.sample in 1:nsamples){
    
    boot.data <- sp.data[0,]
    focal.sp <- unique(sp.data$focal)
    # in this case we have to make sure there are samples of all species,
    # because the parameters are estimated in parallel for all of them
    # furthermore, we cannot sample a site only partially, as all competitors
    # are needed for building the matrix. So, we have two conditions:
    # 1 - all focal sp are sampled
    # 2 - when a site is selected, all its data are in the sample
    # also, remember that a bootstrap sample has the same size
    # as the original dataset. In this case, that might not be exact,
    # as we are sampling complete sites. Thus, the bootstrap data might become
    # slightly bigger, but not much.
    
    # dummy variables
    all.sp.sampled <- FALSE
    # this is for creating a new site index
    # in case a site is repeatedly sampled
    count <- 1
    
    while(nrow(boot.data) < nrow(sp.data) | !all.sp.sampled){
      # my.sp <- sample(focal.sp,1)
      # my.sites <- unique(sp.data$site[sp.data$focal == my.sp])
      # boot.data <- rbind(boot.data,sp.data[sp.data$site == sample(my.sites,1) & sp.data$focal == my.sp,])
      sample.data <- sp.data[sp.data$site == sample(unique(sp.data$site),1),]
      sample.data$site <- count
      boot.data <- rbind(boot.data,sample.data)
      if(length(unique(boot.data$focal)) == length(focal.sp)){
        all.sp.sampled <- TRUE
      }# if
      count <- count + 1
    }# while
    boot.data$site <- as.character(boot.data$site)
    
    ############
    boot.fitness <- boot.data$fitness
    boot.log.fitness <- log(boot.fitness)
    
    if(optim.method == "optim_NM"){
      
      my.boot.par <- optim(init.par, 
                           effect.response.model, 
                           gr = NULL, 
                           method = "Nelder-Mead", 
                           # lower = lower.bounds,
                           # upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           sp.data = boot.data)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "optim_L-BGFS-B"){
      
      my.boot.par <- optim(init.par, 
                           effect.response.model, 
                           gr = NULL, 
                           method = "L-BFGS-B", 
                           lower = lower.bounds,
                           upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           sp.data = boot.data)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "nloptr_CRS2_LM"){
      
      my.boot.par <- nloptr(x0 = init.par,
                            eval_f = effect.response.model,
                            opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            sp.data = boot.data)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_ISRES"){
      
      my.boot.par <- nloptr(x0 = init.par,
                            eval_f = effect.response.model,
                            opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            sp.data = boot.data)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "nloptr_DIRECT_L_RAND"){
      
      my.boot.par <- nloptr(x0 = init.par,
                            eval_f = effect.response.model,
                            opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            sp.data = boot.data)
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.methods[i.method] == "GenSA"){
      
      my.boot.par <- GenSA(par = init.par,
                           fn = effect.response.model,
                           lower = lower.bounds,
                           upper = upper.bounds, 
                           control = list(maxit = 1e3), 
                           sp.data = boot.data)
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "hydroPSO"){
      
      # suppress annoying output
      # sink("/dev/null")
      my.boot.par <- hydroPSO::hydroPSO(par = init.par,
                                        fn = effect.response.model,
                                        lower = lower.bounds,
                                        upper = upper.bounds, 
                                        control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                        sp.data = boot.data)
      sink()
      my.boot.par <- my.boot.par$par
      
    }else if(optim.methods[i.method] == "DEoptimR"){
      
      my.boot.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                        upper = upper.bounds,
                                        fn = effect.response.model,
                                        sp.data = boot.data)
      my.boot.par <- my.boot.par$par
      
    }
    
    boot.results[i.sample,] <- my.boot.par
    
  }  
  
  boot.se <- apply(boot.results,2,sd)
  boot.se
}


