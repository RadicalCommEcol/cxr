
#' Title estimate competition effect and response for a set of species
#' 
#' This function is similar in spirit to cxr_optimize, in that it optimizes a set of parameters
#' via maximum likelihood.
#'
#' @param init.par 1d vector of initial parameters estimates in the order: lambda (for all sp), response (for all sp), effect (for all sp).
#' @param lower.bounds 1d vector of lower bounds for parameter values
#' @param upper.bounds 1d vector of upper bounds for parameter values
#' @param effect.response.model function returning a value to optimize over, e.g. maximum likelihood
#' @param optim.method one of the following: "optim_NM","optim_L-BGFS-B","nloptr_CRS2_LM", 
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
#' @param generate.errors whether we want to compute bootstrap standard errors for the parameters
#' @param bootstrap.samples number of bootstrap samples
#' @return dataframe with estimated values for lambda, e, and r.
#' @export
#'
#' @examples
ER_optimize <- function(init.par,
                        lower.bounds,
                        upper.bounds,
                        effect.response.model,
                        optim.method,
                        sp.data,
                        generate.errors = FALSE,
                        bootstrap.samples = 0){
  
  # in case the sp.data dataframe does not include explicit missing competitors, 
  # here is a somewhat convoluted way for setting zeros to it
  # so that for each focal sp, all competitor sp are included
  sites <- unique(sp.data$site)
  focal.sp <- sort(unique(sp.data$focal))
  missing.data <- sp.data
  missing.data$site <- "0"
  missing.data$focal <- "0"
  missing.data$fitness <- 0
  missing.data$competitor <- "0"
  missing.data$number <- 0
  count <- 1
  
  for(i.site in 1:length(sites)){
    for(i.focal in 1:length(focal.sp)){
      my.competitors <- unique(sp.data$competitor[sp.data$site == sites[i.site] & sp.data$focal == focal.sp[i.focal]])
      if(length(my.competitors) > 0 & length(my.competitors) < length(focal.sp)){
        my.fitness <- sp.data$fitness[sp.data$site == sites[i.site] & sp.data$focal == focal.sp[i.focal]][1]
        
        missing.competitors <- focal.sp[which(!focal.sp %in% my.competitors)]
        for(i.com in 1:length(missing.competitors)){
          
          missing.data$site[count] <- sites[i.site]
          missing.data$focal[count] <- focal.sp[i.focal]
          missing.data$fitness[count] <- my.fitness
          missing.data$competitor[count] <- missing.competitors[i.com]
          missing.data$number[count] <- 0
          
          count <- count + 1
        }# for each missing
      }# if any missing
      # if(length(my.competitors) > 0 & length(my.competitors) < 19){print(paste(i.site,",",i.focal))}
    }# for i.focal
  }# for i.site
  
  missing.data <- droplevels(subset(missing.data,competitor != "0"))
  
  sp.data <- rbind(sp.data,missing.data)
  sp.data <- arrange(sp.data, focal, site, competitor)
  
  # discard focal sp with fitness 0
  sp.data <- droplevels(subset(sp.data, fitness > 0))
  
  
  optim.par <- list(par = rep(NA,length(init.par)), value = NA)
  temp.results <- data.frame(species = focal.sp,
                             lambda = NA,
                             lambda.lower.error = NA,
                             lambda.upper.error = NA,
                             effect.par = NA,
                             effect.lower.error = NA,
                             effect.upper.error = NA,
                             response.par = NA,
                             response.lower.error = NA,
                             response.upper.error = NA,
                             sigma = NA,
                             log.likelihood = NA)
  
  if(optim.method == "optim_NM"){
    
    optim.par <- optim(init.par, 
                       effect.response.model, 
                       gr = NULL, 
                       method = "Nelder-Mead", 
                       # lower = lower.bounds,
                       # upper = upper.bounds,
                       control = list(), 
                       hessian = F,
                       sp.data = sp.data)
    
  }else if(optim.method == "optim_L-BFGS-B"){
    
    optim.par <- optim(init.par, 
                       effect.response.model, 
                       gr = NULL, 
                       method = "L-BFGS-B", 
                       lower = lower.bounds,
                       upper = upper.bounds,
                       control = list(), 
                       hessian = F,
                       sp.data = sp.data)
    
  }else if(optim.method == "nloptr_CRS2_LM"){
    
    optim.par <- nloptr(x0 = init.par,
                        eval_f = effect.response.model,
                        opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                        lb = lower.bounds,
                        ub = upper.bounds,
                        sp.data = sp.data)
    
  }else if(optim.method == "nloptr_ISRES"){
    
    optim.par <- nloptr(x0 = init.par,
                        eval_f = effect.response.model,
                        opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                        lb = lower.bounds,
                        ub = upper.bounds,
                        sp.data = sp.data)
    
  }else if(optim.method == "nloptr_DIRECT_L_RAND"){
    
    optim.par <- nloptr(x0 = init.par,
                        eval_f = effect.response.model,
                        opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                        lb = lower.bounds,
                        ub = upper.bounds,
                        sp.data = sp.data)
    
  }else if(optim.method == "GenSA"){
    
    optim.par <- GenSA(par = init.par,
                       fn = effect.response.model,
                       lower = lower.bounds,
                       upper = upper.bounds, 
                       control = list(maxit = 1e3), 
                       sp.data = sp.data)
    
  }else if(optim.method == "hydroPSO"){
    
    # suppress annoying output??
    # sink("/dev/null")
    optim.par <- hydroPSO::hydroPSO(par = init.par,
                                    fn = effect.response.model,
                                    lower = lower.bounds,
                                    upper = upper.bounds, 
                                    control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                    sp.data = sp.data)
    
    
  }else if(optim.method == "DEoptimR"){
    
    optim.par <- DEoptimR::JDEoptim(lower = lower.bounds,
                                    upper = upper.bounds,
                                    fn = effect.response.model,
                                    sp.data = sp.data)
  }
  
  ##################################
  # tidy the output from the method
  # if-else the method outputs optim-like values
  if(optim.method %in% c("optim_NM","optim_L-BGFS-B","DEoptimR","hydroPSO","GenSA")){
    
    print(paste("Effect-Response:",optim.method," finished with convergence status ",optim.par$convergence,sep=""))
    
    temp.results$lambda <- optim.par$par[1:length(focal.sp)]
    temp.results$response.par <- optim.par$par[(length(focal.sp)+1):(length(focal.sp)+length(focal.sp))]
    temp.results$effect.par <- optim.par$par[(length(focal.sp)+1+length(focal.sp)):(length(optim.par$par)-1)]
    temp.results$sigma <- optim.par$par[length(optim.par$par)]
    temp.results$log.likelihood <- optim.par$value
    
    my.par <- optim.par$par
    
  }else{ # methods with different nomenclature
    
    print(paste("Effect-Response:",optim.method," finished with convergence status ",optim.par$status,sep=""))
    
    temp.results$lambda <- optim.par$solution[1:length(focal.sp)]
    temp.results$response.par <- optim.par$solution[(length(focal.sp)+1):(length(focal.sp)+length(focal.sp))]
    temp.results$effect.par <- optim.par$solution[(length(focal.sp)+1+length(focal.sp)):(length(par)-1)]
    temp.results$sigma <- optim.par$solution[length(optim.par$solution)]
    temp.results$log.likelihood <- optim.par$objective
    
    my.par <- optim.par$solution
  }  
  
  #####################
  # standard errors via bootstrapping
  if(generate.errors){
    errors <- ER_SEbootstrap(optim.method,
                             effect.response.model,
                             lower.bounds,
                             upper.bounds,
                             init.par,
                             sp.data,
                             nsamples = bootstrap.samples)
  }else{
    errors <- NA
  }
  
  if(sum(is.na(errors)) == 0){
    
    temp.results$lambda.lower.error <- my.par[1:length(focal.sp)]-1.96*errors[1:length(focal.sp)]
    temp.results$lambda.upper.error <- my.par[1:length(focal.sp)]+1.96*errors[1:length(focal.sp)]
    
    temp.results$response.lower.error <- my.par[(length(focal.sp)+1):(length(focal.sp)+length(focal.sp))]-1.96*errors[(length(focal.sp)+1):(length(focal.sp)+length(focal.sp))]
    temp.results$response.upper.error <- my.par[(length(focal.sp)+1):(length(focal.sp)+length(focal.sp))]+1.96*errors[(length(focal.sp)+1):(length(focal.sp)+length(focal.sp))]
    
    temp.results$effect.lower.error <- my.par[(length(focal.sp)+1+length(focal.sp)):(length(par)-1)]-1.96*errors[(length(focal.sp)+1+length(focal.sp)):(length(par)-1)]
    temp.results$effect.upper.error <- my.par[(length(focal.sp)+1+length(focal.sp)):(length(par)-1)]+1.96*errors[(length(focal.sp)+1+length(focal.sp)):(length(par)-1)]
    
  }
  
  return(temp.results)
  
}





