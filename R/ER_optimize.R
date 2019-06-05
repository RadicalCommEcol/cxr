
#' Estimate competition effect and response for a set of species
#' 
#' This function is similar in spirit to cxr_optimize, in that it optimizes a set of parameters
#' via maximum likelihood.
#'
#' @param lambda.vector 1d vector of lambda estimates/initial values (depending on whether lambda values are optimized or not)
#' @param e.vector 1d vector of competitive effect initial values
#' @param r.vector 1d vector of competitive response initial values
#' @param lambda.cov numeric matrix of num.sp x num.covariates, effect of every covariate on species' lambda. 
#' Discarded if covariates are not passed. 
#' @param e.cov numeric matrix of num.sp x num.covariates, effect of every covariate on species' competitive effect.
#' Discarded if covariates are not passed.
#' @param r.cov numeric matrix of num.sp x num.covariates, effect of every covariate on species' competitive response. 
#' Discarded if covariates are not passed.
#' @param sigma initial value for variation estimate.
#' @param lambda.lower.bound lower bounds for lambda, in case it is optimized. Same dimensions as lambda.vector
#' @param lambda.upper.bound upper bounds for lambda, in case it is optimized. Same dimensions as lambda.vector
#' @param e.lower.bound lower bounds for competitive effects. Same dimensions as e.vector
#' @param e.upper.bound upper bounds for competitive effects. Same dimensions as e.vector
#' @param r.lower.bound lower bounds for competitive responses. Same dimensions as r.vector
#' @param r.upper.bound upper bounds for competitive responses. Same dimensions as r.vector
#' @param lambda.cov.lower.bound lower bounds for covariate effects on lambda. Same dimensions as lambda.cov. 
#' Discarded if covariates are not passed.
#' @param lambda.cov.upper.bound upper bounds for covariate effects on lambda. Same dimensions as lambda.cov. 
#' Discarded if covariates are not passed. 
#' @param e.cov.lower.bound lower bounds for covariate effects on e. Same dimensions as e.cov. 
#' Discarded if covariates are not passed.
#' @param e.cov.upper.bound upper bounds for covariate effects on e. Same dimensions as e.cov. 
#' Discarded if covariates are not passed. 
#' @param r.cov.lower.bound lower bounds for covariate effects on r. Same dimensions as r.cov. 
#' Discarded if covariates are not passed.
#' @param r.cov.upper.bound upper bounds for covariate effects on r. Same dimensions as r.cov. 
#' Discarded if covariates are not passed.
#' @param sigma.lower.bound lower bound for sigma.
#' @param sigma.upper.bound upper bound for sigma.
#' @param effect.response.model function returning a value to optimize over, e.g. maximum likelihood
#' @param optim.method one of the following: "optim_NM","optim_L-BFGS-B","nloptr_CRS2_LM", 
#' "nloptr_ISRES","nloptr_DIRECT_L_RAND","GenSA","hydroPSO","DEoptimR".
#' @param sp.data dataframe with all the necessary information in long format. It should have the following columns:
#' - site: character ID
#' - focal: character ID of the focal species. Any number of focal species is allowed, but the number of focal species
#' must match the number of initial parameters (one lambda, e, and r per species).
#' - fitness: numeric, a fitness metric
#' - competitor: character, ID of a competitor for that observation. The set of competitors must be, for now, the same
#' as the set of focal species.
#' - number: number of competitor individuals from the associated species. Observations without competitors of a given species
#' can be explicit, i.e. setting number to zero, or implicit, not included in the dataframe.
#' @param covariates optional matrix/dataframe with as many rows as observations and covariates in columns.
#' @param optimize.lambda boolean, whether we want to optimize lambda values or not
#' @param generate.errors whether we want to compute bootstrap standard errors for the parameters
#' @param bootstrap.samples number of bootstrap samples
#' @return list with estimated species values for e, r, lambda (optional), and if covariates are given, the effects of covariates on lambda, r, and e.
#' @export
#'
ER_optimize <- function(lambda.vector,
                        e.vector,
                        r.vector,
                        lambda.cov = NULL,
                        e.cov = NULL,
                        r.cov = NULL,
                        sigma,
                        lambda.lower.bound = 0,
                        lambda.upper.bound = 1e3,
                        e.lower.bound,
                        e.upper.bound,
                        r.lower.bound,
                        r.upper.bound,
                        lambda.cov.lower.bound = NULL,
                        lambda.cov.upper.bound = NULL,
                        e.cov.lower.bound = NULL,
                        e.cov.upper.bound = NULL,
                        r.cov.lower.bound = NULL,
                        r.cov.upper.bound = NULL,
                        sigma.lower.bound,
                        sigma.upper.bound,
                        effect.response.model,
                        optim.method,
                        sp.data,
                        covariates = NULL,
                        optimize.lambda = FALSE,
                        generate.errors = FALSE,
                        bootstrap.samples = 0){
  
  if(!is.null(covariates)){
    if(is.null(lambda.cov) | is.null(e.cov) | is.null(r.cov) | is.null(lambda.cov.lower.bound) | is.null(lambda.cov.upper.bound) | is.null(e.cov.lower.bound) | is.null(e.cov.upper.bound)
       | is.null(r.cov.lower.bound) | is.null(r.cov.upper.bound)){
      stop("ER_optimize ERROR: Covariates are given but initial values/bounds for lambda.cov, r.cov, or e.cov are missing")
    }
  }
  
  # fill up matrices
  # num.sp x num.observations. 1 if species is focal in a given observation, 0 otherwise
  target_all <- NULL
  # num.sp x num.observations. density of each species in each observation
  density_all <- NULL

  # fitness metric of the focal sp at each observation
  log.fitness <- log(sp.data$fitness)
  
  # species names and number
  sp.list <- unique(sp.data$focal)
  num.sp <- length(sp.list)
  
  for(i.sp in 1:num.sp){
    
    target.my.sp <- integer(nrow(sp.data))
    target.my.sp <- ifelse(sp.data$focal == sp.list[i.sp],1,0)
    
    density.my.sp <- integer(nrow(sp.data))
    for(i.obs in 1:nrow(sp.data)){
      if(sp.data$competitor[i.obs] == sp.list[i.sp]){
        density.my.sp[i.obs] <- sp.data$number[i.obs]
      }
    }
    
    target_all <- rbind(target_all,target.my.sp)
    density_all <- rbind(density_all,density.my.sp)
  }
  
  # which model to use depending on whether we optimize lambda or not
  # also depending on whether there are covariates
  if(optimize.lambda){
    if(is.null(covariates)){
      init.par <- c(lambda.vector,r.vector,e.vector,sigma)
      lower.bounds <- c(lambda.lower.bound,r.lower.bound,e.lower.bound,sigma.lower.bound)
      upper.bounds <- c(lambda.upper.bound,r.upper.bound,e.upper.bound,sigma.upper.bound)
    }else{
      init.par <- c(lambda.vector,r.vector,e.vector,lambda.cov,r.cov,e.cov,sigma)
      lower.bounds <- c(lambda.lower.bound,r.lower.bound,e.lower.bound,lambda.cov.lower.bound,r.cov.lower.bound,e.cov.lower.bound,sigma.lower.bound)
      upper.bounds <- c(lambda.upper.bound,r.upper.bound,e.upper.bound,lambda.cov.upper.bound,r.cov.upper.bound,e.cov.upper.bound,sigma.upper.bound)
    }
  }else{
    if(is.null(covariates)){
      init.par <- c(r.vector,e.vector,sigma)
      lower.bounds <- c(r.lower.bound,e.lower.bound,sigma.lower.bound)
      upper.bounds <- c(r.upper.bound,e.upper.bound,sigma.upper.bound)
    }else{
      init.par <- c(r.vector,e.vector,lambda.cov,r.cov,e.cov,sigma)
      lower.bounds <- c(r.lower.bound,e.lower.bound,lambda.cov.lower.bound,r.cov.lower.bound,e.cov.lower.bound,sigma.lower.bound)
      upper.bounds <- c(r.upper.bound,e.upper.bound,lambda.cov.upper.bound,r.cov.upper.bound,e.cov.upper.bound,sigma.upper.bound)
    }# if-else covariates
  }# if-else optimize lambda
  
  optim.par <- list(par = rep(NA,length(init.par)), value = NA)
  
  #############
  # initialize result vectors with proper names
  fit.lambda <- rep(NA, num.sp)
  names(fit.lambda) <- sp.list
  
  fit.response <- rep(NA, num.sp)
  names(fit.response) <- sp.list
  
  fit.effect <- rep(NA, num.sp)
  names(fit.effect) <- sp.list
  
  fit.sigma <- NA
  fit.log.likelihood <- NA
  
  fit.lambda.lower.error <- rep(NA, num.sp)
  names(fit.lambda.lower.error) <- sp.list
  
  fit.lambda.upper.error <- rep(NA, num.sp)
  names(fit.lambda.upper.error) <- sp.list
  
  fit.response.lower.error <- rep(NA, num.sp)
  names(fit.response.lower.error) <- sp.list
  
  fit.response.upper.error <- rep(NA, num.sp)
  names(fit.response.upper.error) <- sp.list
  
  fit.effect.lower.error <- rep(NA, num.sp)
  names(fit.effect.lower.error) <- sp.list
  
  fit.effect.upper.error <- rep(NA, num.sp)
  names(fit.effect.upper.error) <- sp.list
  
  if(!is.null(covariates)){
    fit.lambda.cov <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.r.cov <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.e.cov <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.lambda.cov.lower.error <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.lambda.cov.upper.error <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.r.cov.lower.error <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.r.cov.upper.error <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.e.cov.lower.error <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
    fit.e.cov.upper.error <- matrix(NA,nrow = num.sp,ncol = ncol(covariates))
  }else{
    fit.lambda.cov <- NA
    fit.r.cov <- NA
    fit.e.cov <- NA
    fit.lambda.cov.lower.error <- NA
    fit.lambda.cov.upper.error <- NA
    fit.r.cov.lower.error <- NA
    fit.r.cov.upper.error <- NA
    fit.e.cov.lower.error <- NA
    fit.e.cov.upper.error <- NA
  }
  
  # optimization methods
  if(optim.method == "optim_NM"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- optim(init.par, 
                         effect.response.model, 
                         gr = NULL, 
                         method = "Nelder-Mead", 
                         # lower = lower.bounds,
                         # upper = upper.bounds,
                         control = list(), 
                         hessian = F,
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness,
                         covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- optim(init.par, 
                         effect.response.model, 
                         gr = NULL, 
                         method = "Nelder-Mead", 
                         # lower = lower.bounds,
                         # upper = upper.bounds,
                         control = list(), 
                         hessian = F,
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness,
                         lambda = lambda.vector,
                         covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "optim_L-BFGS-B"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- optim(init.par, 
                         effect.response.model, 
                         gr = NULL, 
                         method = "L-BFGS-B", 
                         lower = lower.bounds,
                         upper = upper.bounds,
                         control = list(), 
                         hessian = F,
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness) 
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- optim(init.par, 
                         effect.response.model, 
                         gr = NULL, 
                         method = "L-BFGS-B", 
                         lower = lower.bounds,
                         upper = upper.bounds,
                         control = list(), 
                         hessian = F,
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness,
                         lambda = lambda.vector)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "nloptr_CRS2_LM"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness,
                          lambda = lambda.vector)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "nloptr_ISRES"){
    
    if(optimize.lambda){
      tryCatch({
    optim.par <- nloptr(x0 = init.par,
                        eval_f = effect.response.model,
                        opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                        lb = lower.bounds,
                        ub = upper.bounds,
                        target_all = target_all,
                        density_all = density_all,
                        log.fitness = log.fitness)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness,
                          lambda = lambda.vector)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "nloptr_DIRECT_L_RAND"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness,
                          lambda = lambda.vector)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "GenSA"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- GenSA(par = init.par,
                         fn = effect.response.model,
                         lower = lower.bounds,
                         upper = upper.bounds, 
                         control = list(maxit = 1e3), 
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- GenSA(par = init.par,
                         fn = effect.response.model,
                         lower = lower.bounds,
                         upper = upper.bounds, 
                         control = list(maxit = 1e3), 
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness,
                         lambda = lambda.vector)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "hydroPSO"){
    
    if(optimize.lambda){
      # suppress annoying output??
      # sink("/dev/null")
      tryCatch({
      optim.par <- hydroPSO::hydroPSO(par = init.par,
                                      fn = effect.response.model,
                                      lower = lower.bounds,
                                      upper = upper.bounds, 
                                      control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                      target_all = target_all,
                                      density_all = density_all,
                                      log.fitness = log.fitness)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      # suppress annoying output??
      # sink("/dev/null")
      tryCatch({
      optim.par <- hydroPSO::hydroPSO(par = init.par,
                                      fn = effect.response.model,
                                      lower = lower.bounds,
                                      upper = upper.bounds, 
                                      control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                      target_all = target_all,
                                      density_all = density_all,
                                      log.fitness = log.fitness,
                                      lambda = lambda.vector)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "DEoptimR"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                      upper = upper.bounds,
                                      fn = effect.response.model,
                                      target_all = target_all,
                                      density_all = density_all,
                                      log.fitness = log.fitness)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                      upper = upper.bounds,
                                      fn = effect.response.model,
                                      target_all = target_all,
                                      density_all = density_all,
                                      log.fitness = log.fitness,
                                      lambda = lambda.vector)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
  }# if method
  
  ##################################
  # tidy the output from the method
  # if-else the method outputs optim-like values
  if(optim.method %in% c("optim_NM","optim_L-BFGS-B","DEoptimR","hydroPSO","GenSA")){
    
    print(paste("Effect-Response:",optim.method," finished with convergence status ",optim.par$convergence,sep=""))
    
    if(optimize.lambda){
      fit.lambda <- optim.par$par[1:num.sp]
      fit.response <- optim.par$par[(num.sp+1):(num.sp+num.sp)]
      
      if(!is.null(covariates)){
        fit.effect <- optim.par$par[(num.sp+1+num.sp):num.sp+num.sp+num.sp]
        fit.lambda.cov <- optim.par$par[(num.sp + num.sp + num.sp + 1):(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)))]
        fit.r.cov <- optim.par$par[(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
        fit.e.cov <- optim.par$par[(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
      }else{
        fit.effect <- optim.par$par[(num.sp+1+num.sp):(length(optim.par$par)-1)]
      }
    }else{
      fit.lambda <- lambda.vector
      fit.response <- optim.par$par[1:num.sp]
      
      if(!is.null(covariates)){
        fit.effect <- optim.par$par[(num.sp+1):num.sp+num.sp]
        fit.lambda.cov <- optim.par$par[(num.sp + num.sp + 1):(num.sp + num.sp + (num.sp*ncol(covariates)))]
        fit.r.cov <- optim.par$par[(num.sp + num.sp + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
        fit.e.cov <- optim.par$par[(num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
      }else{
        fit.effect <- optim.par$par[(num.sp+1):(length(optim.par$par)-1)]
      }
    }
    
    fit.sigma <- optim.par$par[length(optim.par$par)]
    fit.log.likelihood <- optim.par$value
    
  }else{ # methods with different nomenclature
    
    print(paste("Effect-Response:",optim.method," finished with convergence status ",optim.par$status,sep=""))
    
    if(optimize.lambda){
      fit.lambda <- optim.par$solution[1:num.sp]
      fit.response <- optim.par$solution[(num.sp+1):(num.sp+num.sp)]
      
      if(!is.null(covariates)){
        fit.effect <- optim.par$solution[(num.sp+1+num.sp):num.sp+num.sp+num.sp]
        fit.lambda.cov <- optim.par$solution[(num.sp + num.sp + num.sp + 1):(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)))]
        fit.r.cov <- optim.par$solution[(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
        fit.e.cov <- optim.par$solution[(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
      }else{
        fit.effect <- optim.par$solution[(num.sp+1+num.sp):(length(optim.par$solution)-1)]
      }
    }else{
      fit.lambda <- lambda.vector
      fit.response <- optim.par$solution[1:num.sp]
      
      if(!is.null(covariates)){
        fit.effect <- optim.par$solution[(num.sp+1):num.sp+num.sp]
        fit.lambda.cov <- optim.par$solution[(num.sp + num.sp + 1):(num.sp + num.sp + (num.sp*ncol(covariates)))]
        fit.r.cov <- optim.par$solution[(num.sp + num.sp + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
        fit.e.cov <- optim.par$solution[(num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + 1):(num.sp + num.sp + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)) + (num.sp*ncol(covariates)))]
      }else{
        fit.effect <- optim.par$solution[(num.sp+1):(length(optim.par$solution)-1)]
      }
    }
    
    fit.sigma <- optim.par$solution[length(optim.par$solution)]
    fit.log.likelihood <- optim.par$objective
  }  
  
  #####################
  # standard errors via bootstrapping
  if(generate.errors){
    errors <- ER_SEbootstrap(lambda.vector = lambda.vector,
                             e.vector = e.vector,
                             r.vector = r.vector,
                             sigma = sigma,
                             e.lower.bound = e.lower.bound,
                             e.upper.bound = e.upper.bound,
                             r.lower.bound = r.lower.bound,
                             r.upper.bound = r.upper.bound,
                             sigma.lower.bound = sigma.lower.bound,
                             sigma.upper.bound = sigma.upper.bound,
                             sp.data = sp.data,
                             effect.response.model = effect.response.model,
                             optim.method = optim.method,
                             optimize.lambda = optimize.lambda,
                             nsamples = bootstrap.samples)
  }else{
    errors <- NA
  }
  
  if(sum(is.na(errors)) == 0){
    
    if(optimize.lambda){
      fit.lambda.lower.error <- fit.lambda - 1.96*errors[1:num.sp]
      fit.lambda.upper.error <- fit.lambda + 1.96*errors[1:num.sp]
      fit.response.lower.error <- fit.response - 1.96*errors[(num.sp+1):(num.sp+num.sp)]
      fit.response.upper.error <- fit.response + 1.96*errors[(num.sp+1):(num.sp+num.sp)]
      
      if(is.null(covariates)){
        fit.effect.lower.error <- fit.effect - 1.96*errors[(num.sp+1+num.sp):(length(init.par)-1)]
        fit.effect.upper.error <- fit.effect + 1.96*errors[(num.sp+1+num.sp):(length(init.par)-1)]
      }else{
        fit.effect.lower.error <- fit.effect - 1.96*errors[(num.sp+1+num.sp):(num.sp+num.sp+num.sp)]
        fit.effect.upper.error <- fit.effect + 1.96*errors[(num.sp+1+num.sp):(num.sp+num.sp+num.sp)]
        fit.lambda.cov.lower.error <- fit.lambda.cov - 1.96*errors[(num.sp+num.sp+num.sp+1):(num.sp+num.sp+num.sp+(num.sp*ncol(covariates)))]
        fit.lambda.cov.upper.error <- fit.lambda.cov + 1.96*errors[(num.sp+num.sp+num.sp+1):(num.sp+num.sp+num.sp+(num.sp*ncol(covariates)))]
        fit.response.cov.lower.error <- fit.response.cov - 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.response.cov.upper.error <- fit.response.cov + 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.effect.cov.lower.error <- fit.effect.cov - 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
        fit.effect.cov.upper.error <- fit.effect.cov + 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
      }      
    }else{
      fit.response.lower.error <- fit.response - 1.96*errors[1:num.sp]
      fit.response.upper.error <- fit.response + 1.96*errors[1:num.sp]
      
      if(is.null(covariates)){
        fit.effect.lower.error <- fit.effect - 1.96*errors[(num.sp+1):(length(init.par)-1)]
        fit.effect.upper.error <- fit.effect + 1.96*errors[(num.sp+1):(length(init.par)-1)]
      }else{
        fit.effect.lower.error <- fit.effect - 1.96*errors[(num.sp+1):(num.sp+num.sp)]
        fit.effect.upper.error <- fit.effect + 1.96*errors[(num.sp+1):(num.sp+num.sp)]
        fit.lambda.cov.lower.error <- fit.lambda.cov - 1.96*errors[(num.sp+num.sp+1):(num.sp+num.sp+(num.sp*ncol(covariates)))]
        fit.lambda.cov.upper.error <- fit.lambda.cov + 1.96*errors[(num.sp+num.sp+1):(num.sp+num.sp+(num.sp*ncol(covariates)))]
        fit.response.cov.lower.error <- fit.response.cov - 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.response.cov.upper.error <- fit.response.cov + 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.effect.cov.lower.error <- fit.effect.cov - 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
        fit.effect.cov.upper.error <- fit.effect.cov + 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
      }

    }# if-else
  }# if errors
  
  return.list <- list(lambda = fit.lambda,
                      lambda.lower.error = fit.lambda.lower.error,
                      lambda.upper.error = fit.lambda.upper.error,
                      response = fit.response,
                      response.lower.error = fit.response.lower.error,
                      response.upper.error = fit.response.upper.error,
                      effect = fit.effect,
                      effect.lower.error = fit.effect.lower.error,
                      effect.upper.error = fit.effect.upper.error,
                      sigma = sigma,
                      lambda.cov = fit.lambda.cov,
                      lambda.cov.lower.error = fit.lambda.cov.lower.error,
                      lambda.cov.upper.error = fit.lambda.cov.upper.error,
                      response.cov = fit.r.cov,
                      response.cov.lower.error = fit.r.cov.lower.error,
                      response.cov.upper.error = fit.r.cov.upper.error,
                      effect.cov = fit.e.cov,
                      effect.cov.lower.error = fit.e.cov.lower.error,
                      effect.cov.upper.error = fit.e.cov.upper.error,
                      log.likelihood = fit.log.likelihood)
  
  return.list[lengths(return.list) == 0] <- NA_character_
  
  return.list
  
}





