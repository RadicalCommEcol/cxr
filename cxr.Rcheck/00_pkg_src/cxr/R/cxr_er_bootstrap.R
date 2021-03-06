#' standard error estimates for effect and response parameters
#' 
#' Computes bootstrap standard errors for a given effect/response function
#'
#' @param effect.response.model effect/response function
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
#' @param init.par 1d vector of initial parameters
#' @param lower.bounds 1d vector of lower bounds
#' @param upper.bounds 1d vector of upper bounds
#' @param covariates dataframe/matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' from an observation.
#' @param optimize.lambda boolean, whether to optimize the values of lambda or not.
#' @param lambda.vector in case lambda is not to be optimized, fixed values for it.
#' @param nsamples how many bootstrap samples to compute.
#'
#' @return 1d vector, the standard error of each parameter in init.par
#' @export
cxr_er_bootstrap <- function(effect.response.model,
                           optim.method,
                           sp.data,
                           init.par,
                           lower.bounds,
                           upper.bounds,
                           covariates,
                           optimize.lambda,
                           lambda.vector,
                           nsamples){
  
  if(nsamples<2){
    print("cxr_er_bootstrap: number of bootstrap samples cannot be < 2. Setting bootstrap samples to 2.")
    nsamples <- 2
  }
  
  boot.results <- matrix(nrow = nsamples, ncol = length(init.par))
  
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
      my.sample.site <- sample(unique(sp.data$site),1)
      # positions <- c(positions,my.sample.site)
      sample.data <- sp.data[sp.data$site == my.sample.site,]
      sample.data$orig.site <- sample.data$site
      sample.data$site <- count
      boot.data <- rbind(boot.data,sample.data)
      if(length(unique(boot.data$focal)) == length(focal.sp)){
        all.sp.sampled <- TRUE
      }# if
      count <- count + 1
    }# while
    
    # are there covariates?
    if(is.null(covariates)){
      boot.covariates <- NULL
    }else{
      positions <- match(boot.data$orig.site,sp.data$site)
      boot.covariates <- as.matrix(covariates[positions,])
    }    
    boot.data$site <- as.character(boot.data$site)
    
    ############
    boot.fitness <- boot.data$fitness
    boot.log.fitness <- log(boot.fitness)
    
    # target and density matrices
    # num.sp x num.observations. 1 if species is focal in a given observation, 0 otherwise
    target_all <- NULL
    # num.sp x num.observations. density of each species in each observation
    density_all <- NULL
    
    sp.list <- unique(boot.data$focal)
    
    for(i.sp in 1:length(sp.list)){
      
      target.my.sp <- integer(nrow(boot.data))
      target.my.sp <- ifelse(boot.data$focal == sp.list[i.sp],1,0)
      
      density.my.sp <- integer(nrow(boot.data))
      for(i.obs in 1:nrow(boot.data)){
        if(boot.data$competitor[i.obs] == sp.list[i.sp]){
          density.my.sp[i.obs] <- boot.data$number[i.obs]
        }
      }
      
      target_all <- rbind(target_all,target.my.sp)
      density_all <- rbind(density_all,density.my.sp)
    }
    
    if(optim.method == "optim_NM"){
      
      if(optimize.lambda){
        tryCatch({
        my.boot.par <- optim(init.par, 
                           effect.response.model, 
                           gr = NULL, 
                           method = "Nelder-Mead", 
                           # lower = lower.bounds,
                           # upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           target_all = target_all,
                           density_all = density_all,
                           log.fitness = boot.log.fitness,
                           covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        tryCatch({
        my.boot.par <- optim(init.par, 
                           effect.response.model, 
                           gr = NULL, 
                           method = "Nelder-Mead", 
                           # lower = lower.bounds,
                           # upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           target_all = target_all,
                           density_all = density_all,
                           log.fitness = boot.log.fitness,
                           lambda = lambda.vector,
                           covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }
      my.boot.par <- my.boot.par$par
      
    }else if(optim.method == "optim_L-BFGS-B"){
      
      if(optimize.lambda){
        tryCatch({
        my.boot.par <- optim(init.par, 
                             effect.response.model, 
                             gr = NULL, 
                             method = "L-BFGS-B", 
                             lower = lower.bounds,
                             upper = upper.bounds,
                             control = list(), 
                             hessian = F,
                             target_all = target_all,
                             density_all = density_all,
                             log.fitness = boot.log.fitness,
                             covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        tryCatch({
        my.boot.par <- optim(init.par, 
                             effect.response.model, 
                             gr = NULL, 
                             method = "L-BFGS-B", 
                             lower = lower.bounds,
                             upper = upper.bounds,
                             control = list(), 
                             hessian = F,
                             target_all = target_all,
                             density_all = density_all,
                             log.fitness = boot.log.fitness,
                             lambda = lambda.vector,
                             covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }
      
      my.boot.par <- my.boot.par$par
      
    }else if(optim.method == "nloptr_CRS2_LM"){
      
      if(optimize.lambda){
        tryCatch({
        my.boot.par <- nloptr::nloptr(x0 = init.par,
                              eval_f = effect.response.model,
                              opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                              lb = lower.bounds,
                              ub = upper.bounds,
                              target_all = target_all,
                              density_all = density_all,
                              log.fitness = boot.log.fitness,
                              covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        tryCatch({
        my.boot.par <- nloptr::nloptr(x0 = init.par,
                              eval_f = effect.response.model,
                              opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                              lb = lower.bounds,
                              ub = upper.bounds,
                              target_all = target_all,
                              density_all = density_all,
                              log.fitness = boot.log.fitness,
                              lambda = lambda.vector,
                              covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }

      my.boot.par <- my.boot.par$solution
      
    }else if(optim.method == "nloptr_ISRES"){
      
      if(optimize.lambda){
        tryCatch({
        my.boot.par <- nloptr::nloptr(x0 = init.par,
                              eval_f = effect.response.model,
                              opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                              lb = lower.bounds,
                              ub = upper.bounds,
                              target_all = target_all,
                              density_all = density_all,
                              log.fitness = boot.log.fitness,
                              covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        tryCatch({
        my.boot.par <- nloptr::nloptr(x0 = init.par,
                              eval_f = effect.response.model,
                              opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                              lb = lower.bounds,
                              ub = upper.bounds,
                              target_all = target_all,
                              density_all = density_all,
                              log.fitness = boot.log.fitness,
                              lambda = lambda.vector,
                              covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }
    
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.method == "nloptr_DIRECT_L_RAND"){
      
      if(optimize.lambda){
        tryCatch({
        my.boot.par <- nloptr::nloptr(x0 = init.par,
                              eval_f = effect.response.model,
                              opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                              lb = lower.bounds,
                              ub = upper.bounds,
                              target_all = target_all,
                              density_all = density_all,
                              log.fitness = boot.log.fitness,
                              covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        tryCatch({
        my.boot.par <- nloptr::nloptr(x0 = init.par,
                              eval_f = effect.response.model,
                              opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                              lb = lower.bounds,
                              ub = upper.bounds,
                              target_all = target_all,
                              density_all = density_all,
                              log.fitness = boot.log.fitness,
                              lambda = lambda.vector,
                              covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }
      
      my.boot.par <- my.boot.par$solution
      
    }else if(optim.method == "GenSA"){
      
      if(optimize.lambda){
        tryCatch({
        my.boot.par <- GenSA::GenSA(par = init.par,
                             fn = effect.response.model,
                             lower = lower.bounds,
                             upper = upper.bounds, 
                             control = list(maxit = 1e3), 
                             target_all = target_all,
                             density_all = density_all,
                             log.fitness = boot.log.fitness,
                             covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        tryCatch({
        my.boot.par <- GenSA::GenSA(par = init.par,
                             fn = effect.response.model,
                             lower = lower.bounds,
                             upper = upper.bounds, 
                             control = list(maxit = 1e3), 
                             target_all = target_all,
                             density_all = density_all,
                             log.fitness = boot.log.fitness,
                             lambda = lambda.vector,
                             covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }
      
      my.boot.par <- my.boot.par$par
      
    }else if(optim.method == "hydroPSO"){
      
      if(optimize.lambda){
        # suppress annoying output
        # sink("/dev/null")
        tryCatch({
        my.boot.par <- hydroPSO::hydroPSO(par = init.par,
                                          fn = effect.response.model,
                                          lower = lower.bounds,
                                          upper = upper.bounds, 
                                          control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                          target_all = target_all,
                                          density_all = density_all,
                                          log.fitness = boot.log.fitness,
                                          covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        # suppress annoying output
        # sink("/dev/null")
        tryCatch({
        my.boot.par <- hydroPSO::hydroPSO(par = init.par,
                                          fn = effect.response.model,
                                          lower = lower.bounds,
                                          upper = upper.bounds, 
                                          control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                          target_all = target_all,
                                          density_all = density_all,
                                          log.fitness = boot.log.fitness,
                                          lambda = lambda.vector,
                                          covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }
      
      my.boot.par <- my.boot.par$par
      
    }else if(optim.method == "DEoptimR"){
      
      if(optimize.lambda){
        tryCatch({
        my.boot.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                          upper = upper.bounds,
                                          fn = effect.response.model,
                                          target_all = target_all,
                                          density_all = density_all,
                                          log.fitness = boot.log.fitness,
                                          covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }else{
        tryCatch({
        my.boot.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                          upper = upper.bounds,
                                          fn = effect.response.model,
                                          target_all = target_all,
                                          density_all = density_all,
                                          log.fitness = boot.log.fitness,
                                          lambda = lambda.vector,
                                          covariates = boot.covariates)
        }, error=function(e){cat("er_optim ERROR :",conditionMessage(e), "\n")})
      }
      
      my.boot.par <- my.boot.par$par
      
    }
    
    boot.results[i.sample,] <- my.boot.par
    
  }  
  
  boot.se <- apply(boot.results,2,sd)
  boot.se
}


