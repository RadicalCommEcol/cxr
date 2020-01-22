
#' Standard error estimates for model parameters
#' 
#' Computes bootstrap standard errors for a given population dynamics model.
#'
#' @param fitness_model function returning a single value to minimize, given a set of parameters and a fitness metric
#' @param optimization_method numerical optimization method
#' @param data dataframe with observations in rows and two sets of columns:
#' * fitness: fitness metric for the focal individual
#' * neighbours: columns with user-defined names with number of neighbours for each group
#' @param covariates optional matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation.
#' @param init_par 1d vector of initial parameters
#' @param lower_bounds 1d vector of lower bounds
#' @param upper_bounds 1d vector of upper bounds
#' @param fixed_parameters optional list specifying values of fixed parameters, with components "lambda","alpha","lambda_cov", and "alpha_cov".
#' @param bootstrap_samples how many bootstrap samples to compute.
#'
#' @return 1d vector, the standard error of each parameter in init.par
#' @import stats 
#' @md
#' @export
cxr_pm_bootstrap <- function(fitness_model,
                             optimization_method,
                             data,
                             covariates,
                             init_par,
                             lower_bounds,
                             upper_bounds,
                             fixed_parameters,
                             bootstrap_samples){ 
  if(bootstrap_samples<2){
    print("cxr_pm_bootstrap: number of bootstrap samples cannot be < 2. Setting bootstrap samples to 2.")
    bootstrap_samples <- 2
  }
  
  bootres <- matrix(nrow = bootstrap_samples, ncol = length(init_par))
  
  for(i.sample in 1:bootstrap_samples){
    
    bsample <- sample(nrow(data),nrow(data),replace = TRUE)
    # bsample <- sample(length(log.fitness),length(log.fitness),replace = T)
    
    # sample data
    bdata <- data[bsample,]
    
    bneigh <- subset(bdata, select = -c(fitness))
    bneigh <- as.matrix(bneigh)
    # how many neighbour species?
    nneigh <- ncol(bneigh)
    
    if(is.data.frame(covariates)){
      bcov <- as.data.frame(covariates[bsample,])
    }else if(is.matrix(covariates)){
      bcov <- as.data.frame(covariates[bsample,])
    }else{
      bcov <- 0
    }
    
    bpar <- NULL
    
    ############
    if(optimization_method %in% c("BFGS", "CG", "Nelder-Mead", "ucminf")){
      tryCatch({
        bpar <- optimx::optimx(par = init_par, 
                               fn = fitness_model, 
                               gr = NULL, 
                               method = optimization_method,
                               # lower = lower_bounds,
                               # upper = upper_bounds,
                               control = list(), 
                               hessian = F,
                               fitness = log(bdata$fitness), 
                               neigh_matrix = bneigh,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        bpar <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
        
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method %in% c("L-BFGS-B", "nlm", "nlminb", 
                                        "Rcgmin", "Rvmmin", "spg", 
                                        "bobyqa", "nmkb", "hjkb")){
      tryCatch({
        bpar <- optimx::optimx(par = init_par, 
                               fn = fitness_model, 
                               gr = NULL, 
                               method = optimization_method,
                               lower = lower_bounds,
                               upper = upper_bounds,
                               control = list(), 
                               hessian = F,
                               fitness = log(bdata$fitness), 
                               neigh_matrix = bneigh,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        bpar <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_CRS2_LM"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                               eval_f = fitness_model,
                               opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e4),
                               lb = lower_bounds,
                               ub = upper_bounds,
                               fitness = log(bdata$fitness), 
                               neigh_matrix = bneigh,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_ISRES"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                               eval_f = fitness_model,
                               opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e4),
                               lb = lower_bounds,
                               ub = upper_bounds,
                               fitness = log(bdata$fitness), 
                               neigh_matrix = bneigh,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_DIRECT_L_RAND"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                               eval_f = fitness_model,
                               opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e4),
                               lb = lower_bounds,
                               ub = upper_bounds,
                               fitness = log(bdata$fitness), 
                               neigh_matrix = bneigh,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "GenSA"){
      tryCatch({
        bpar <- GenSA::GenSA(par = init_par,
                             fn = fitness_model,
                             lower = lower_bounds,
                             upper = upper_bounds, 
                             control = list(maxit = 1e3), 
                             fitness = log(bdata$fitness), 
                             neigh_matrix = bneigh,
                             covariates = covariates, 
                             fixed_parameters = fixed_parameters)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "hydroPSO"){
      tryCatch({
        # suppress annoying output
        bpar <- hydroPSO::hydroPSO(par = init_par,
                                   fn = fitness_model,
                                   lower = lower_bounds,
                                   upper = upper_bounds, 
                                   control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                   fitness = log(bdata$fitness), 
                                   neigh_matrix = bneigh,
                                   covariates = covariates, 
                                   fixed_parameters = fixed_parameters)
        
        bpar <- bpar$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "DEoptimR"){
      tryCatch({
        bpar <- DEoptimR::JDEoptim(lower = lower_bounds,
                                   upper = upper_bounds,
                                   fn = fitness_model,
                                   fitness = log(bdata$fitness), 
                                   neigh_matrix = bneigh,
                                   covariates = covariates, 
                                   fixed_parameters = fixed_parameters)
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
    boot.se <- apply(bootres,2,sd)
  }else{
    boot.se <- bootres[1,]
  }
  boot.se
}


