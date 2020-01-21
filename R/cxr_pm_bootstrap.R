
#' Standard error estimates for model parameters
#' 
#' Computes bootstrap standard errors for a given population dynamics model.
#'
#' @param fitness.model function returning a single value to minimize, given a set of parameters and a fitness metric
#' @param optim.method optimization method to use. One of the following: "optim_NM","optim_L-BFGS-B","nloptr_CRS2_LM", 
#' "nloptr_ISRES","nloptr_DIRECT_L_RAND","GenSA","hydroPSO","DEoptimR".
#' @param param.list string vector giving the parameters that are to be optimized for the fitness model.
#' @param fixed.terms string vector giving the parameters that are NOT optimized for the fitness model.
#' @param log.fitness 1d vector, log of the fitness metric for every observation
#' @param init.par 1d vector of initial parameters
#' @param lower.bounds 1d vector of lower bounds
#' @param upper.bounds 1d vector of upper bounds
#' @param focal.comp.matrix matrix with observations in rows and neighbours in columns. Each cell is the number of neighbours
#' of a given species in a given observation.
#' @param focal.covariates optional matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation.
#' @param nsamples how many bootstrap samples to compute.
#'
#' @return 1d vector, the standard error of each parameter in init.par
#' @import stats 
#' @export
cxr_pm_bootstrap <- function(fitness.model,
                             optim.method,
                             param.list, # drop
                             fixed.terms,
                             log.fitness, # change this and focal.comp.matrix to data
                             init.par,
                             lower.bounds,
                             upper.bounds,
                             focal.comp.matrix,
                             focal.covariates, # change to covariates
                             nsamples){ # change to bootstrap_samples
  
  if(nsamples<2){
    print("cxr_pm_bootstrap: number of bootstrap samples cannot be < 2. Setting bootstrap samples to 2.")
    nsamples <- 2
  }
  
  bootres <- matrix(nrow = nsamples, ncol = length(init.par))
  
  for(i.sample in 1:nsamples){
    
    bsample <- sample(nrow(data),nrow(data),replace = TRUE)
    # bsample <- sample(length(log.fitness),length(log.fitness),replace = T)
    
    # sample data
    bdata <- data[bsample,]
    
    bneigh <- subset(bdata, select = -c(fitness))
    bneigh <- as.matrix(bneigh)
    # how many neighbour species?
    nneigh <- ncol(bneigh)
    ncov
    if(is.data.frame(covariates)){
      bcov <- as.data.frame(covariates[bsample,])
    }else if(is.matrix(covariates)){
      bcov <- as.data.frame(covariates[bsample,])
    }else{
      bcov <- 0
    }
    
    bpar <- NULL
    
    ############
    if(optim.method %in% c("BFGS", "CG", "Nelder-Mead", "lbfgsb3", "Rtnmin", "snewton",
                           "snewtonm", "ucminf", "newuoa", "hjn", "lbfgs", "subplex")){
      tryCatch({
        
        bpar <- optimx::optimx(init.par, 
                               fitness.model, 
                               gr = NULL, 
                               method = optim.method, 
                               # lower = lower.bounds,
                               # upper = upper.bounds,
                               control = list(), 
                               hessian = F,
                               param.list = param.list,
                               log.fitness = boot.fitness, 
                               focal.comp.matrix = bneigh,
                               # num.covariates = num.covariates, 
                               # num.competitors = num.competitors, 
                               focal.covariates = bcov, 
                               fixed.terms = fixed.terms)
        
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        bpar <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method %in% c("L-BFGS-B", "nlm", "nlminb", 
                                 "Rcgmin", "Rvmmin", "spg", 
                                 "bobyqa", "nmkb", "hjkb")){
      tryCatch({
        
        bpar <-  optimx::optimx(init.par, 
                                fitness.model, 
                                gr = NULL, 
                                method = optim.method, 
                                lower = lower.bounds, 
                                upper = upper.bounds,
                                control = list(), 
                                hessian = F,
                                param.list = param.list,
                                log.fitness = boot.fitness, 
                                focal.comp.matrix = bneigh,
                                num.covariates = num.covariates, 
                                num.competitors = num.competitors, 
                                focal.covariates = bcov,
                                fixed.terms = fixed.terms)
        
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        my.boot.par <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "bobyqa"){
      tryCatch({
        bpar <- optimx::optimr(init.par, 
                               fitness.model, 
                               gr = NULL, 
                               method = "bobyqa", 
                               lower = lower.bounds, 
                               upper = upper.bounds,
                               control = list(parscale = abs(init.par)), 
                               hessian = F,
                               param.list = param.list,
                               log.fitness = boot.fitness, 
                               focal.comp.matrix = bneigh,
                               num.covariates = num.covariates, 
                               num.competitors = num.competitors, 
                               focal.covariates = bcov,
                               fixed.terms = fixed.terms)
        bpar <- bpar[,1:length(init.par)]
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "nloptr_CRS2_LM"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                               lb = lower.bounds,
                               ub = upper.bounds,
                               param.list = param.list,
                               log.fitness = boot.fitness, 
                               focal.comp.matrix = bneigh,
                               num.covariates = num.covariates, 
                               num.competitors = num.competitors, 
                               focal.covariates = bcov,
                               fixed.terms = fixed.terms)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "nloptr_ISRES"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                               lb = lower.bounds,
                               ub = upper.bounds,
                               param.list = param.list,
                               log.fitness = boot.fitness, 
                               focal.comp.matrix = bneigh,
                               num.covariates = num.covariates, 
                               num.competitors = num.competitors, 
                               focal.covariates = bcov,
                               fixed.terms = fixed.terms)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "nloptr_DIRECT_L_RAND"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                               lb = lower.bounds,
                               ub = upper.bounds,
                               param.list = param.list,
                               log.fitness = boot.fitness, 
                               focal.comp.matrix = bneigh,
                               num.covariates = num.covariates, 
                               num.competitors = num.competitors, 
                               focal.covariates = bcov,
                               fixed.terms = fixed.terms)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "GenSA"){
      tryCatch({
        bpar <- GenSA::GenSA(par = init.par,fn = fitness.model,
                             lower = lower.bounds,
                             upper = upper.bounds, 
                             control = list(maxit = 1e2), 
                             param.list = param.list,
                             log.fitness = boot.fitness, 
                             focal.comp.matrix = bneigh,
                             num.covariates = num.covariates, 
                             num.competitors = num.competitors, 
                             focal.covariates = bcov,
                             fixed.terms = fixed.terms)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "hydroPSO"){
      tryCatch({
        # suppress annoying output
        bpar <- hydroPSO::hydroPSO(par = init.par,fn = fitness.model,
                                   lower = lower.bounds,
                                   upper = upper.bounds, 
                                   control=list(write2disk=FALSE, maxit = 1e2, MinMax = "min", verbose = F),
                                   param.list = param.list,
                                   log.fitness = boot.fitness, 
                                   focal.comp.matrix = bneigh,
                                   num.covariates = num.covariates, 
                                   num.competitors = num.competitors, 
                                   focal.covariates = bcov,
                                   fixed.terms = fixed.terms)
        
        bpar <- bpar$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "DEoptimR"){
      tryCatch({
        bpar <- DEoptimR::JDEoptim(lower = lower.bounds,
                                   upper = upper.bounds,
                                   maxiter = 100,
                                   fn = fitness.model,
                                   param.list = param.list,
                                   log.fitness = boot.fitness, 
                                   focal.comp.matrix = bneigh,
                                   num.covariates = num.covariates, 
                                   num.competitors = num.competitors, 
                                   focal.covariates = bcov,
                                   fixed.terms = fixed.terms)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }
    
    if(!is.null(bpar)){
      if(sum(is.na(bpar)) == 0){
        bootres[i.sample,] <- bpar
      }
    }
    
  }  
  
  bootres <- bootres[which(!is.na(rowSums(bootres))),]
  if(nrow(bootres)>2){
    boot.se <- apply(boot.results,2,sd)
  }else{
    boot.se <- boot.results[1,]
  }
  boot.se
}


