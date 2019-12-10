
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
    print("cxr_pm_bootstrap: number of bootstrap samples cannot be < 2. Setting bootstrap samples to 2.")
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
    }else if(is.matrix(focal.covariates)){
      boot.covariates <- as.data.frame(focal.covariates[my.sample,])
    }else{
      boot.covariates <- 0
    }
    my.boot.par <- NULL
    ############
    if(optim.method %in% c("BFGS", "CG", "Nelder-Mead", "lbfgsb3", "Rtnmin", "snewton",
                           "snewtonm", "ucminf", "newuoa", "hjn", "lbfgs", "subplex")){
      tryCatch({
      my.boot.par <- optimx::optimx(init.par, 
                           fitness.model, 
                           gr = NULL, 
                           method = optim.method, 
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
      par.pos <- which(!names(my.boot.par) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
      my.boot.par <- as.numeric(my.boot.par[,par.pos])
      row.names(my.boot.par) <- NULL
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method %in% c("L-BFGS-B", "nlm", "nlminb", 
                                 "Rcgmin", "Rvmmin", "spg", 
                                 "bobyqa", "nmkb", "hjkb")){
      tryCatch({
      my.boot.par <-  optimx::optimx(init.par, 
                           fitness.model, 
                           gr = NULL, 
                           method = optim.method, 
                           lower = lower.bounds, 
                           upper = upper.bounds,
                           control = list(), 
                           hessian = F,
                           param.list = param.list,
                           log.fitness = boot.fitness, 
                           focal.comp.matrix = boot.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = boot.covariates,
                           fixed.terms = fixed.terms)
      par.pos <- which(!names(my.boot.par) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
      my.boot.par <- as.numeric(my.boot.par[,par.pos])
      row.names(my.boot.par) <- NULL
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "nloptr_CRS2_LM"){
      tryCatch({
      my.boot.par <- nloptr::nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            param.list = param.list,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = boot.covariates,
                            fixed.terms = fixed.terms)
      my.boot.par <- my.boot.par$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "nloptr_ISRES"){
      tryCatch({
      my.boot.par <- nloptr::nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            param.list = param.list,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = boot.covariates,
                            fixed.terms = fixed.terms)
      my.boot.par <- my.boot.par$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "nloptr_DIRECT_L_RAND"){
      tryCatch({
      my.boot.par <- nloptr::nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                            lb = lower.bounds,
                            ub = upper.bounds,
                            param.list = param.list,
                            log.fitness = boot.fitness, 
                            focal.comp.matrix = boot.comp.matrix,
                            num.covariates = num.covariates, 
                            num.competitors = num.competitors, 
                            focal.covariates = boot.covariates,
                            fixed.terms = fixed.terms)
      my.boot.par <- my.boot.par$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "GenSA"){
      tryCatch({
      my.boot.par <- GenSA::GenSA(par = init.par,fn = fitness.model,
                           lower = lower.bounds,
                           upper = upper.bounds, 
                           control = list(maxit = 1e2), 
                           param.list = param.list,
                           log.fitness = boot.fitness, 
                           focal.comp.matrix = boot.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = boot.covariates,
                           fixed.terms = fixed.terms)
      my.boot.par <- my.boot.par$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "hydroPSO"){
      tryCatch({
      # suppress annoying output
      my.boot.par <- hydroPSO::hydroPSO(par = init.par,fn = fitness.model,
                                        lower = lower.bounds,
                                        upper = upper.bounds, 
                                        control=list(write2disk=FALSE, maxit = 1e2, MinMax = "min", verbose = F),
                                        param.list = param.list,
                                        log.fitness = boot.fitness, 
                                        focal.comp.matrix = boot.comp.matrix,
                                        num.covariates = num.covariates, 
                                        num.competitors = num.competitors, 
                                        focal.covariates = boot.covariates,
                                        fixed.terms = fixed.terms)

      my.boot.par <- my.boot.par$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optim.method == "DEoptimR"){
      tryCatch({
      my.boot.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                        upper = upper.bounds,
                                        maxiter = 100,
                                        fn = fitness.model,
                                        param.list = param.list,
                                        log.fitness = boot.fitness, 
                                        focal.comp.matrix = boot.comp.matrix,
                                        num.covariates = num.covariates, 
                                        num.competitors = num.competitors, 
                                        focal.covariates = boot.covariates,
                                        fixed.terms = fixed.terms)
      my.boot.par <- my.boot.par$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }
    
    if(!is.null(my.boot.par)){
      if(sum(is.na(my.boot.par)) == 0){
        boot.results[i.sample,] <- my.boot.par
      }
    }

  }  
  
  boot.results <- boot.results[which(!is.na(rowSums(boot.results))),]
  if(nrow(boot.results)>2){
    boot.se <- apply(boot.results,2,sd)
  }else{
    boot.se <- boot.results[1,]
  }
  boot.se
}


