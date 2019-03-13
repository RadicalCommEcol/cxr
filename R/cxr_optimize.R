
#' Title optimize parameters according to a given fitness model
#'
#' @param init.par 1d vector of initial parameters estimates (see function InitParams)
#' @param lower.bounds 1d vector of lower bounds for parameter values
#' @param upper.bounds 1d vector of upper bounds for parameter values
#' @param fitness.model function returning a value to optimize over
#' @param optim.method one of the following: "optim_NM","optim_L-BGFS-B","nloptr_CRS2_LM", 
#' "nloptr_ISRES","nloptr_DIRECT_L_RAND","GenSA","hydroPSO","DEoptimR".
#' @param log.fitness fitness values
#' @param focal.comp.matrix competition matrix (one row per observation, columns are abundance of competitors of the different sp)
#' @param focal.covariates covariates matrix (one row per observation, columns are covariates)
#' @param generate.errors whether we want to compute bootstrap standard errors for the parameters
#' @param bootstrap.samples number of bootstrap samples
#'
#' @return list with several named dataframes and parameter matrices. I think names of the elements are self-explanatory.
#' @export
#'
#' @examples
cxr_optimize <- function(init.par,
                        lower.bounds,
                        upper.bounds,
                        fitness.model,
                        optim.method,
                        log.fitness,
                        focal.comp.matrix,
                        focal.covariates,
                        generate.errors = FALSE,
                        bootstrap.samples = 0){

  num.competitors <- dim(focal.comp.matrix)[2]
  num.covariates <- ifelse(is.null(ncol(focal.covariates)),0,ncol(focal.covariates))
  
optim.par <- list(par = rep(NA,length(init.par)), value = NA)
temp.results <- data.frame(lambda = NA,
                           lambda.lower.error = NA,
                           lambda.upper.error = NA,
                           sigma = NA,
                           log.likelihood = NA)

if(optim.method == "optim_NM"){
  
  optim.par <- optim(init.par, 
                     fitness.model, 
                     gr = NULL, 
                     method = "Nelder-Mead", 
                     # lower = lower.bounds,
                     # upper = upper.bounds,
                     control = list(), 
                     hessian = F,
                     log.fitness = log.fitness, 
                     focal.comp.matrix = focal.comp.matrix,
                     num.covariates = num.covariates, 
                     num.competitors = num.competitors, 
                     focal.covariates = focal.covariates)
  
}else if(optim.method == "optim_L-BGFS-B"){
  
  optim.par <- optim(init.par, 
                     fitness.model, 
                     gr = NULL, 
                     method = "L-BFGS-B", 
                     lower = lower.bounds, 
                     upper = upper.bounds,
                     control = list(), 
                     hessian = F,
                     log.fitness = log.fitness, 
                     focal.comp.matrix = focal.comp.matrix,
                     num.covariates = num.covariates, 
                     num.competitors = num.competitors, 
                     focal.covariates = focal.covariates)
  
}else if(optim.method == "nloptr_CRS2_LM"){
  
  optim.par <- nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                      lb = lower.bounds,
                      ub = upper.bounds,
                      log.fitness = log.fitness, 
                      focal.comp.matrix = focal.comp.matrix,
                      num.covariates = num.covariates, 
                      num.competitors = num.competitors, 
                      focal.covariates = focal.covariates)
  
}else if(optim.method == "nloptr_ISRES"){
  
  optim.par <- nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                      lb = lower.bounds,
                      ub = upper.bounds,
                      log.fitness = log.fitness, 
                      focal.comp.matrix = focal.comp.matrix,
                      num.covariates = num.covariates, 
                      num.competitors = num.competitors, 
                      focal.covariates = focal.covariates)
  
}else if(optim.method == "nloptr_DIRECT_L_RAND"){
  
  optim.par <- nloptr(x0 = init.par,eval_f = fitness.model,opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                      lb = lower.bounds,
                      ub = upper.bounds,
                      log.fitness = log.fitness, 
                      focal.comp.matrix = focal.comp.matrix,
                      num.covariates = num.covariates, 
                      num.competitors = num.competitors, 
                      focal.covariates = focal.covariates)
  
}else if(optim.method == "GenSA"){
  
  optim.par <- GenSA(par = init.par,fn = fitness.model,lower = lower.bounds,upper = upper.bounds, control = list(maxit = 1e3), 
                     log.fitness = log.fitness, 
                     focal.comp.matrix = focal.comp.matrix,
                     num.covariates = num.covariates, 
                     num.competitors = num.competitors, 
                     focal.covariates = focal.covariates)
  
}else if(optim.method == "hydroPSO"){
  
  # suppress annoying output??
  # sink("/dev/null")
  optim.par <- hydroPSO::hydroPSO(par = init.par,fn = fitness.model,lower = lower.bounds,upper = upper.bounds, control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                  log.fitness = log.fitness, 
                                  focal.comp.matrix = focal.comp.matrix,
                                  num.covariates = num.covariates, 
                                  num.competitors = num.competitors, 
                                  focal.covariates = focal.covariates)
  
  
}else if(optim.method == "DEoptimR"){
  
  optim.par <- DEoptimR::JDEoptim(lower = lower.bounds,upper = upper.bounds,fn = fitness.model,
                                  log.fitness = log.fitness, 
                                  focal.comp.matrix = focal.comp.matrix,
                                  num.covariates = num.covariates, 
                                  num.competitors = num.competitors, 
                                  focal.covariates = focal.covariates)
}

##################################
# tidy the output from the method
# if-else the method outputs optim-like values
if(optim.method %in% c("optim_NM","optim_L-BGFS-B","DEoptimR","hydroPSO","GenSA")){
  
  print(paste("...",optim.method," finished with convergence status ",optim.par$convergence,sep=""))
  
  temp.results$lambda <- optim.par$par[1]
  temp.results$sigma <- optim.par$par[length(optim.par$par)]
  temp.results$log.likelihood <- optim.par$value
  temp.alpha.matrix <- NA
  temp.lambda.cov.matrix <- NA
  temp.alpha.cov.matrix <- NA
  temp.alpha.lower.error <- NA 
  temp.alpha.upper.error <- NA
  temp.lambda.cov.lower.error <- NA
  temp.lambda.cov.upper.error <- NA
  temp.alpha.cov.lower.error <- NA
  temp.alpha.cov.upper.error <- NA
  
  # obtain the different matrices
  if(models[i.model] == 2){
    temp.alpha.matrix <- optim.par$par[2]
  }else if(models[i.model] == 3){
    temp.alpha.matrix <- optim.par$par[2:(1+num.competitors)]
  }else if(models[i.model] == 4){
    temp.alpha.matrix <- optim.par$par[(1+num.covariates+1):(1+num.covariates+num.competitors)]
    temp.lambda.cov.matrix <- optim.par$par[(1+1):(1+num.covariates)]
    temp.alpha.cov.matrix <- optim.par$par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+num.covariates)]
  }else if(models[i.model] == 5){
    temp.alpha.matrix <- optim.par$par[(1+num.covariates+1):(1+num.covariates+num.competitors)]
    temp.lambda.cov.matrix <- optim.par$par[(1+1):(1+num.covariates)]
    temp.alpha.cov.matrix <- optim.par$par[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.covariates*num.competitors))] #effects of cov 1, 2... on alpha_i
  }
  my.par <- optim.par$par
  
}else{ # methods with different nomenclature
  
  print(paste("...",optim.method," finished with convergence status ",optim.par$status,sep=""))
  
  temp.results$lambda <- optim.par$solution[1]
  temp.results$sigma <- optim.par$solution[length(optim.par$solution)]
  temp.results$log.likelihood <- optim.par$objective
  
  # obtain the different matrices
  if(models[i.model] > 2){
    temp.alpha.matrix <- optim.par$solution[(1+num.covariates+1):(1+num.covariates+num.competitors)]
    if(models[i.model] == 4){
      temp.lambda.cov.matrix <- optim.par$solution[(1+1):(1+num.covariates)]
      temp.alpha.cov.matrix <- optim.par$solution[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+num.covariates)]
    }else if(models[i.model] == 5){
      temp.lambda.cov.matrix <- optim.par$solution[(1+1):(1+num.covariates)]
      temp.alpha.cov.matrix <- optim.par$solution[(1+num.covariates+num.competitors+1):(1+num.covariates+num.competitors+(num.covariates*num.competitors))] #effects of cov 1, 2... on alpha_i
    }
  }
  my.par <- optim.par$solution
}  

#####################
# standard errors via bootstrapping
if(generate.errors){
  errors <- SEbootstrap(optim.method,
                        fitness.model,
                        lower.bounds,
                        upper.bounds,
                        init.par,
                        log.fitness,
                        focal.comp.matrix,
                        focal.covariates,
                        num.competitors,
                        num.covariates,
                        nsamples = bootstrap.samples)
}else{
  errors <- NA
}

if(sum(is.na(errors)) == 0){
  
  temp.results$lambda.lower.error <- my.par[1]-1.96*errors[1]
  temp.results$lambda.upper.error <- my.par[1]+1.96*errors[1]
  
  # different models have different parameters
  # written sequentially for readability
  # as optim models need a vector for multidim optimization,
  # parameters need to be recovered from this vector and from the error matrix
  # hence these long sequences
  if(models[i.model] == 2){
    
    temp.alpha.lower.error <- my.par[2]-1.96*errors[2]
    temp.alpha.upper.error <- my.par[2]+1.96*errors[2]
    
  }else if(models[i.model] == 3){
    
    temp.alpha.lower.error <- my.par[2:(num.competitors+1)]-1.96*errors[2:(num.competitors+1)]
    temp.alpha.upper.error <- my.par[2:(num.competitors+1)]+1.96*errors[2:(num.competitors+1)]
    
  }else if(models[i.model] == 4){
    
    temp.lambda.cov.lower.error <- my.par[2:(1+num.covariates)]-1.96*errors[2:(1+num.covariates)]
    temp.lambda.cov.upper.error <- my.par[2:(1+num.covariates)]+1.96*errors[2:(1+num.covariates)]
    
    temp.alpha.lower.error <- my.par[(1+1+num.covariates):(num.covariates+num.competitors+1)]-1.96*
      errors[(1+1+num.covariates):(num.covariates+num.competitors+1)]
    temp.alpha.upper.error <- my.par[(1+1+num.covariates):(num.covariates+num.competitors+1)]+1.96*
      errors[(1+1+num.covariates):(num.covariates+num.competitors+1)]
    
    temp.alpha.cov.lower.error <- my.par[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+num.covariates)]-1.96*
      errors[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+num.covariates)]
    temp.alpha.cov.upper.error <- my.par[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+num.covariates)]+1.96*
      errors[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+num.covariates)]
    
  }else if(models[i.model] == 5){
    
    temp.lambda.cov.lower.error <- my.par[2:(1+num.covariates)]-1.96*errors[2:(1+num.covariates)]
    temp.lambda.cov.upper.error <- my.par[2:(1+num.covariates)]+1.96*errors[2:(1+num.covariates)]
    
    temp.alpha.lower.error <- my.par[(1+1+num.covariates):(num.covariates+num.competitors+1)]-1.96*
      errors[(1+1+num.covariates):(num.covariates+num.competitors+1)]
    temp.alpha.upper.error <- my.par[(1+1+num.covariates):(num.covariates+num.competitors+1)]+1.96*
      errors[(1+1+num.covariates):(num.covariates+num.competitors+1)]
    
    temp.alpha.cov.lower.error <- my.par[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+(num.covariates*num.competitors))]-1.96*
      errors[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+num.covariates)]
    temp.alpha.cov.upper.error <- my.par[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+(num.covariates*num.competitors))]+1.96*
      errors[(2+num.covariates+num.competitors):(num.covariates+num.competitors+1+num.covariates)]
    
  }
}
    
  return(list(lambda.results = temp.results,
            alpha.matrix = temp.alpha.matrix,
            alpha.upper.error = temp.alpha.upper.error,
            alpha.lower.error = temp.alpha.lower.error,
            alpha.cov.matrix = temp.alpha.cov.matrix,
            lambda.cov.matrix = temp.lambda.cov.matrix,
            lambda.cov.upper.error = temp.lambda.cov.upper.error,
            lambda.cov.lower.error = temp.lambda.cov.lower.error,
            alpha.cov.upper.error = temp.alpha.cov.upper.error,
            alpha.cov.lower.error = temp.alpha.cov.lower.error))

}





