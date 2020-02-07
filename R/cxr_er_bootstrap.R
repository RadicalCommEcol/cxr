#' standard error estimates for effect and response parameters
#' 
#' Computes bootstrap standard errors for a given effect/response function
#'
#' @param fitness_model effect/response function, see \code{cxr_er_fit}
#' @param init_par initial values for parameters
#' @param fixed_parameters list with values for fixed parameters, or NULL.
#' @inheritParams cxr_er_fit
#' @export
cxr_er_bootstrap <- function(fitness_model,
                             optimization_method,
                             data,
                             covariates,
                             init_par,
                             lower_bounds,
                             upper_bounds,
                             fixed_parameters,
                             bootstrap_samples){ 
  
  if(bootstrap_samples<2){
    print("cxr_er_bootstrap: number of bootstrap samples cannot be < 2. Setting bootstrap samples to 2.")
    bootstrap_samples <- 2
  }
  
  bootres <- matrix(nrow = bootstrap_samples, ncol = length(init_par))
  sp.list <- unique(data$focal)
  
  # for sampling each focal species equally in base R,
  # a combined split-sample approach is effective.
  # However, this does not preserve the indexes of the sample
  # so if there are covariates, an option is to append them to bdata
  # and, after the sample, recover them
  if(!is.null(covariates)){
    datacov <- cbind(data,covariates) 
  }
  
  for(i.sample in 1:bootstrap_samples){
    
    # sample grouping by focal species
    if(!is.null(covariates)){
      bdata <- do.call(rbind, 
                       lapply(split(datacov, datacov$focal), 
                              function(x) x[sample(nrow(x), nrow(x),replace = TRUE), ]))
      bcov <- bdata[,(ncol(data)+1):(ncol(datacov))]
      bdata <- bdata[,1:ncol(data)]
    }else{
      bdata <- do.call(rbind, 
                       lapply(split(data, data$focal), 
                              function(x) x[sample(nrow(x), nrow(x),replace = TRUE), ]))
      bcov <- NULL
    }
    
    # target and density matrices
    # num.sp x num.observations. 1 if species is focal in a given observation, 0 otherwise
    btarget_all <- NULL
    # num.sp x num.observations. density of each species in each observation
    bdensity_all <- NULL
    
    for(i.sp in 1:length(sp.list)){
      
      target.my.sp <- integer(nrow(bdata))
      target.my.sp <- ifelse(bdata$focal == sp.list[i.sp],1,0)
      density.my.sp <- bdata[[sp.list[i.sp]]]
      
      btarget_all <- rbind(btarget_all,target.my.sp)
      bdensity_all <- rbind(bdensity_all,density.my.sp)
    }
    rownames(btarget_all) <- sp.list
    rownames(bdensity_all) <- sp.list

# parameter optimization --------------------------------------------------

    bpar <- NULL
    
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
                                       target = btarget_all,
                                       density = bdensity_all,
                                       covariates = covariates,  
                                       fixed_parameters = fixed_parameters)
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        bpar <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
        
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
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
                                       target = btarget_all,
                                       density = bdensity_all,
                                       covariates = covariates, 
                                       fixed_parameters = fixed_parameters)
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        bpar <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_CRS2_LM"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                                       eval_f = fitness_model,
                                       opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e4),
                                       lb = lower_bounds,
                                       ub = upper_bounds,
                                       fitness = log(bdata$fitness), 
                                       target = btarget_all,
                                       density = bdensity_all,
                                       covariates = covariates,  
                                       fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_ISRES"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                                       eval_f = fitness_model,
                                       opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e4),
                                       lb = lower_bounds,
                                       ub = upper_bounds,
                                       fitness = log(bdata$fitness), 
                                       target = btarget_all,
                                       density = bdensity_all,
                                       covariates = covariates, 
                                       fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_DIRECT_L_RAND"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                                       eval_f = fitness_model,
                                       opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e4),
                                       lb = lower_bounds,
                                       ub = upper_bounds,
                                       fitness = log(bdata$fitness), 
                                       target = btarget_all,
                                       density = bdensity_all,
                                       covariates = covariates, 
                                       fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "GenSA"){
      tryCatch({
        bpar <- GenSA::GenSA(par = init_par,
                                     fn = fitness_model,
                                     lower = lower_bounds,
                                     upper = upper_bounds, 
                                     control = list(maxit = 1e3), 
                                     fitness = log(bdata$fitness), 
                                     target = btarget_all,
                                     density = bdensity_all,
                                     covariates = covariates, 
                                     fixed_parameters = fixed_parameters)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "hydroPSO"){
      tryCatch({
        # suppress annoying output??
        # sink("/dev/null")
        bpar <- hydroPSO::hydroPSO(par = init_par,
                                           fn = fitness_model,
                                           lower = lower_bounds,
                                           upper = upper_bounds, 
                                           control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                           fitness = log(bdata$fitness), 
                                           target = btarget_all,
                                           density = bdensity_all,
                                           covariates = covariates, 
                                           fixed_parameters = fixed_parameters)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
      
    }else if(optimization_method == "DEoptimR"){
      tryCatch({
        bpar <- DEoptimR::JDEoptim(lower = lower_bounds,
                                           upper = upper_bounds,
                                           fn = fitness_model,
                                           fitness = log(bdata$fitness), 
                                           target = btarget_all,
                                           density = bdensity_all,
                                           covariates = covariates,  
                                           fixed_parameters = fixed_parameters)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_er_bootstrap optimization ERROR :",conditionMessage(e), "\n")})
    }
    
    if(!is.null(bpar)){
      if(sum(is.na(bpar)) == 0){
        bootres[i.sample,] <- bpar
      }
    }
    
  }# for i.sample  
  
  bootres <- bootres[which(!is.na(rowSums(bootres))),]
  if(nrow(bootres)>2){
    boot.se <- apply(bootres,2,sd)
  }else{
    boot.se <- bootres[1,]
  }
  
  if(!is.null(names(init_par))){
    names(boot.se) <- paste(names(init_par),"_se",sep="")
  }
  boot.se
  
}


