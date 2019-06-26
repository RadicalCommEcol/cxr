
#' Estimate competition effects and responses for a set of species
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
#' @param lower.lambda lower bound for lambda, in case it is optimized. Either length 1 or same length as lambda.vector
#' @param upper.lambda upper bound for lambda, in case it is optimized. Either length 1 or same length as lambda.vector
#' @param lower.e lower bound for competitive effects. Either length 1 or same length as e.vector
#' @param upper.e upper bound for competitive effects. Either length 1 or same length as e.vector
#' @param lower.r lower bound for competitive responses. Either length 1 or same length as r.vector
#' @param upper.r upper bound for competitive responses. Either length 1 or same length as r.vector
#' @param lower.lambda.cov lower bound for covariate effects on lambda. Either length 1 or same length as lambda.cov. 
#' Discarded if covariates are not passed.
#' @param upper.lambda.cov upper bound for covariate effects on lambda. Either length 1 or same length as lambda.cov. 
#' Discarded if covariates are not passed. 
#' @param lower.e.cov lower bound for covariate effects on e. Either length 1 or same length as e.cov. 
#' Discarded if covariates are not passed.
#' @param upper.e.cov upper bound for covariate effects on e. Either length 1 or same length as e.cov. 
#' Discarded if covariates are not passed. 
#' @param lower.r.cov lower bound for covariate effects on r. Either length 1 or same length as r.cov. 
#' Discarded if covariates are not passed.
#' @param upper.r.cov upper bound for covariate effects on r. Either length 1 or same length as r.cov. 
#' Discarded if covariates are not passed.
#' @param lower.sigma lower bound for sigma. Length 1.
#' @param upper.sigma upper bound for sigma. Length 1.
#' @param effect.response.model function returning a value to optimize over, e.g. maximum likelihood
#' @param optim.method optimization method to use. One of the following: "optim_NM","optim_L-BFGS-B","nloptr_CRS2_LM", 
#' "nloptr_ISRES","nloptr_DIRECT_L_RAND","GenSA","hydroPSO","DEoptimR".
#' @param sp.data dataframe with all the necessary information in long format. It should have the following columns:
#' - site: character ID
#' - focal: character ID of the focal species. Any number of focal species is allowed, but the number of focal species
#' must match the number of initial parameters (one lambda, e, and r per species).
#' - fitness: numeric, a fitness metric
#' - competitor: character, ID of a competitor for that observation. The set of competitors must be, for now, the same
#' as the set of focal species.
#' - number: number of neighbouring/competitor individuals from the associated species. Observations without competitors of a given species
#' must be explicit, i.e. setting number to zero.
#' @param covariates optional matrix/dataframe with as many rows as the sp.data dataframe, and covariates in columns.
#' @param optimize.lambda boolean, whether we want to optimize lambda values or not.
#' @param generate.errors boolean, whether to compute bootstrap errors for the fitted parameters. Note that, depending on 
#' the data and optimization method, this may be computationally expensive.
#' @param bootstrap.samples how many bootstrap samples to compute.
#' @return list with estimated species values for e, r, lambda (optional), and if covariates are given, the effects of covariates on lambda, r, and e.
#' @import stats
#' @export
#'
ER_optimize <- function(lambda.vector,
                        e.vector,
                        r.vector,
                        lambda.cov = NULL,
                        e.cov = NULL,
                        r.cov = NULL,
                        sigma,
                        lower.lambda = 0,
                        upper.lambda = 1e3,
                        lower.e,
                        upper.e,
                        lower.r,
                        upper.r,
                        lower.lambda.cov = NULL,
                        upper.lambda.cov = NULL,
                        lower.e.cov = NULL,
                        upper.e.cov = NULL,
                        lower.r.cov = NULL,
                        upper.r.cov = NULL,
                        lower.sigma,
                        upper.sigma,
                        effect.response.model,
                        optim.method,
                        sp.data,
                        covariates = NULL,
                        optimize.lambda = FALSE,
                        generate.errors = FALSE,
                        bootstrap.samples = 0){
  
  # some sanity checks
  if(!is.null(covariates)){
    if(is.null(lambda.cov) | is.null(e.cov) | is.null(r.cov) | is.null(lower.lambda.cov) | is.null(upper.lambda.cov) | is.null(lower.e.cov) | is.null(upper.e.cov)
       | is.null(lower.r.cov) | is.null(upper.r.cov)){
      stop("ER_optimize ERROR: Covariates are given but initial values/bounds for lambda.cov, r.cov, or e.cov are missing")
    }
  }
  
  if (optim.method %in% c("nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND") & !requireNamespace("nloptr", quietly = TRUE)) {
    stop("ER_optimize ERROR: Package \"nloptr\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method == "GenSA" & !requireNamespace("GenSA", quietly = TRUE)) {
    stop("ER_optimize ERROR: Package \"GenSA\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method == "hydroPSO" & !requireNamespace("hydroPSO", quietly = TRUE)) {
    stop("ER_optimize ERROR: Package \"hydroPSO\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method == "DEoptimR" & !requireNamespace("DEoptimR", quietly = TRUE)) {
    stop("ER_optimize ERROR: Package \"DEoptimR\" needed for the method selected to work.",
         call. = FALSE)
  }
  # fill up matrices
  # num.sp x num.observations. 1 if species is focal in a given observation, 0 otherwise
  target_all <- NULL
  # num.sp x num.observations. density of each species in each observation
  density_all <- NULL

  # fitness metric of the focal sp at each observation
  log.fitness <- log(sp.data$fitness)
  
  # species names and number
  sp.list <- as.character(unique(sp.data$focal))
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
  rownames(target_all) <- sp.list
  rownames(density_all) <- sp.list
  
  # check the dimensions of lower and upper bounds for the optimization vectors
  if(length(lower.lambda) == 1){
    lower.lambda.vector <- rep(lower.lambda,length(lambda.vector))
  }
  if(length(upper.lambda) == 1){
    upper.lambda.vector <- rep(upper.lambda,length(lambda.vector))
  }
  if(length(lower.e) == 1){
    lower.e.vector <- rep(lower.e,length(e.vector))
  }
  if(length(upper.e) == 1){
    upper.e.vector <- rep(upper.e,length(e.vector))
  }
  if(length(lower.r) == 1){
    lower.r.vector <- rep(lower.r,length(r.vector))
  }
  if(length(upper.r) == 1){
    upper.r.vector <- rep(upper.r,length(r.vector))
  }
  if(!is.null(covariates)){
    if(length(lower.lambda.cov) == 1){
      lower.lambda.cov.vector <- rep(lower.lambda.cov,length(lambda.cov))
    }
    if(length(upper.lambda.cov) == 1){
      upper.lambda.cov.vector <- rep(upper.lambda.cov,length(lambda.cov))
    }
    if(length(lower.e.cov) == 1){
      lower.e.cov.vector <- rep(lower.e.cov,length(e.cov))
    }
    if(length(upper.e.cov) == 1){
      upper.e.cov.vector <- rep(upper.e.cov,length(e.cov))
    }
    if(length(lower.r.cov) == 1){
      lower.r.cov.vector <- rep(lower.r.cov,length(r.cov))
    }
    if(length(upper.r.cov) == 1){
      upper.r.cov.vector <- rep(upper.r.cov,length(r.cov))
    }
  }
  
  # which model to use depending on whether we optimize lambda or not
  # also depending on whether there are covariates
  if(optimize.lambda){
    if(is.null(covariates)){
      init.par <- c(lambda.vector,r.vector,e.vector,sigma)
      lower.bounds <- c(lower.lambda.vector,lower.r.vector,lower.e.vector,lower.sigma)
      upper.bounds <- c(upper.lambda.vector,upper.r.vector,upper.e.vector,upper.sigma)
    }else{
      init.par <- c(lambda.vector,r.vector,e.vector,lambda.cov,r.cov,e.cov,sigma)
      lower.bounds <- c(lower.lambda.vector,lower.r.vector,lower.e.vector,lower.lambda.cov.vector,lower.r.cov.vector,lower.e.cov.vector,lower.sigma)
      upper.bounds <- c(upper.lambda.vector,upper.r.vector,upper.e.vector,upper.lambda.cov.vector,upper.r.cov.vector,upper.e.cov.vector,upper.sigma)
    }
  }else{
    if(is.null(covariates)){
      init.par <- c(r.vector,e.vector,sigma)
      lower.bounds <- c(lower.r.vector,lower.e.vector,lower.sigma)
      upper.bounds <- c(upper.r.vector,upper.e.vector,upper.sigma)
    }else{
      init.par <- c(r.vector,e.vector,lambda.cov,r.cov,e.cov,sigma)
      lower.bounds <- c(lower.r.vector,lower.e.vector,lower.lambda.cov.vector,lower.r.cov.vector,lower.e.cov.vector,lower.sigma)
      upper.bounds <- c(upper.r.vector,upper.e.vector,upper.lambda.cov.vector,upper.r.cov.vector,upper.e.cov.vector,upper.sigma)
    }# if-else covariates
  }# if-else optimize lambda
  
  optim.par <- list(par = rep(NA,length(init.par)), value = NA)
  
  #############
  # initialize result vectors with proper names
  fit.lambda <- rep(NA, num.sp)
  fit.response <- rep(NA, num.sp)
  fit.effect <- rep(NA, num.sp)
  
  fit.sigma <- NA
  fit.log.likelihood <- NA
  
  fit.lambda.lower.error <- rep(NA, num.sp)
  fit.lambda.upper.error <- rep(NA, num.sp)
  fit.response.lower.error <- rep(NA, num.sp)
  fit.response.upper.error <- rep(NA, num.sp)
  fit.effect.lower.error <- rep(NA, num.sp)
  fit.effect.upper.error <- rep(NA, num.sp)
  
  if(!is.null(covariates)){
    fit.lambda.cov <- rep(NA, num.sp*ncol(covariates))
    fit.r.cov <- rep(NA, num.sp*ncol(covariates))
    fit.e.cov <- rep(NA, num.sp*ncol(covariates))
    fit.lambda.cov.lower.error <- rep(NA, num.sp*ncol(covariates))
    fit.lambda.cov.upper.error <- rep(NA, num.sp*ncol(covariates))
    fit.r.cov.lower.error <- rep(NA, num.sp*ncol(covariates))
    fit.r.cov.upper.error <- rep(NA, num.sp*ncol(covariates))
    fit.e.cov.lower.error <- rep(NA, num.sp*ncol(covariates))
    fit.e.cov.upper.error <- rep(NA, num.sp*ncol(covariates))
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
                         log.fitness = log.fitness,
                         covariates = covariates) 
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
                         lambda = lambda.vector,
                         covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "nloptr_CRS2_LM"){
    
    if(optimize.lambda){
      tryCatch({
        optim.par <- nloptr::nloptr(x0 = init.par,
                                    eval_f = effect.response.model,
                                    opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                                    lb = lower.bounds,
                                    ub = upper.bounds,
                                    target_all = target_all,
                                    density_all = density_all,
                                    log.fitness = log.fitness,
                                    covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- nloptr::nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness,
                          lambda = lambda.vector,
                          covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "nloptr_ISRES"){
    
    if(optimize.lambda){
      tryCatch({
    optim.par <- nloptr::nloptr(x0 = init.par,
                        eval_f = effect.response.model,
                        opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                        lb = lower.bounds,
                        ub = upper.bounds,
                        target_all = target_all,
                        density_all = density_all,
                        log.fitness = log.fitness,
                        covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- nloptr::nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness,
                          lambda = lambda.vector,
                          covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "nloptr_DIRECT_L_RAND"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- nloptr::nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness,
                          covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- nloptr::nloptr(x0 = init.par,
                          eval_f = effect.response.model,
                          opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                          lb = lower.bounds,
                          ub = upper.bounds,
                          target_all = target_all,
                          density_all = density_all,
                          log.fitness = log.fitness,
                          lambda = lambda.vector,
                          covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
    
  }else if(optim.method == "GenSA"){
    
    if(optimize.lambda){
      tryCatch({
      optim.par <- GenSA::GenSA(par = init.par,
                         fn = effect.response.model,
                         lower = lower.bounds,
                         upper = upper.bounds, 
                         control = list(maxit = 1e3), 
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness,
                         covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- GenSA::GenSA(par = init.par,
                         fn = effect.response.model,
                         lower = lower.bounds,
                         upper = upper.bounds, 
                         control = list(maxit = 1e3), 
                         target_all = target_all,
                         density_all = density_all,
                         log.fitness = log.fitness,
                         lambda = lambda.vector,
                         covariates = covariates)
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
                                      log.fitness = log.fitness,
                                      covariates = covariates)
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
                                      lambda = lambda.vector,
                                      covariates = covariates)
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
                                      log.fitness = log.fitness,
                                      covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }else{
      tryCatch({
      optim.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                      upper = upper.bounds,
                                      fn = effect.response.model,
                                      target_all = target_all,
                                      density_all = density_all,
                                      log.fitness = log.fitness,
                                      lambda = lambda.vector,
                                      covariates = covariates)
      }, error=function(e){cat("ER_optimize ERROR :",conditionMessage(e), "\n")})
    }
  }# if method
  
  ##################################
  # tidy the output from the method
  # if-else the method outputs optim-like values
  if(optim.method %in% c("optim_NM","optim_L-BFGS-B","DEoptimR","hydroPSO","GenSA")){
    
    #print(paste("Effect-Response:",optim.method," finished with convergence status ",optim.par$convergence,sep=""))
    
    if(optimize.lambda){
      fit.lambda <- optim.par$par[1:num.sp]
      fit.response <- optim.par$par[(num.sp+1):(num.sp+num.sp)]
      
      if(!is.null(covariates)){
        fit.effect <- optim.par$par[(num.sp+1+num.sp):(num.sp+num.sp+num.sp)]
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
    
    #print(paste("Effect-Response:",optim.method," finished with convergence status ",optim.par$status,sep=""))
    
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
    errors <- ER_SEbootstrap(effect.response.model = effect.response.model,
                             optim.method = optim.method,
                             sp.data = sp.data,
                             init.par = init.par,
                             lower.bounds = lower.bounds,
                             upper.bounds = upper.bounds,
                             covariates = covariates,
                             optimize.lambda = optimize.lambda,
                             lambda.vector = lambda.vector,
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
        fit.r.cov.lower.error <- fit.r.cov - 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.r.cov.upper.error <- fit.r.cov + 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.e.cov.lower.error <- fit.e.cov - 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
        fit.e.cov.upper.error <- fit.e.cov + 1.96*errors[(num.sp+num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
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
        fit.r.cov.lower.error <- fit.r.cov - 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.r.cov.upper.error <- fit.r.cov + 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+1):(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates)))]
        fit.e.cov.lower.error <- fit.e.cov - 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
        fit.e.cov.upper.error <- fit.e.cov + 1.96*errors[(num.sp+num.sp+(num.sp*ncol(covariates))+(num.sp*ncol(covariates))+1):(length(init.par)-1)]
      }

    }# if-else
  }# if errors
  
  # vectors with appropriate names
  names(fit.lambda) <- sp.list
  names(fit.response) <- sp.list
  names(fit.effect) <- sp.list
  
  names(fit.lambda.lower.error) <- sp.list
  names(fit.lambda.upper.error) <- sp.list
  names(fit.response.lower.error) <- sp.list
  names(fit.response.upper.error) <- sp.list
  names(fit.effect.lower.error) <- sp.list
  names(fit.effect.upper.error) <- sp.list
  
  if(!is.null(covariates)){
    if(length(fit.lambda.cov) == length(sp.list)){
      names(fit.lambda.cov) <- sp.list
    }
    if(length(fit.lambda.cov.lower.error) == length(sp.list)){
      names(fit.lambda.cov.lower.error) <- sp.list
    }
    if(length(fit.lambda.cov.upper.error) == length(sp.list)){
      names(fit.lambda.cov.upper.error) <- sp.list
    }
    if(length(fit.r.cov) == length(sp.list)){
      names(fit.r.cov) <- sp.list
    }
    if(length(fit.r.cov.lower.error) == length(sp.list)){
      names(fit.r.cov.lower.error) <- sp.list
    }
    if(length(fit.r.cov.upper.error) == length(sp.list)){
      names(fit.r.cov.upper.error) <- sp.list
    }
    if(length(fit.e.cov) == length(sp.list)){
      names(fit.e.cov) <- sp.list
    }
    if(length(fit.e.cov.lower.error) == length(sp.list)){
      names(fit.e.cov.lower.error) <- sp.list
    }
    if(length(fit.e.cov.upper.error) == length(sp.list)){
      names(fit.e.cov.upper.error) <- sp.list
    }
  }
  
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





