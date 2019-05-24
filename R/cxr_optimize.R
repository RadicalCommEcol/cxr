
cxr_optimize <- function(fitness.model,
                         optim.method,
                         param.list, #lambda,sigma,alpha,lambda.cov,alpha.cov
                         log.fitness,
                         init.lambda = NULL,
                         lower.lambda = 1,
                         upper.lambda = 1e5,
                         init.sigma = NULL,
                         lower.sigma = 1e-10,
                         upper.sigma = 1e5,
                         init.alpha = 1e-4,
                         lower.alpha = 0,
                         upper.alpha = 1e5,
                         init.lambda.cov = 1e-3,
                         lower.lambda.cov = 1e-4,
                         upper.lambda.cov = 1e5,
                         init.alpha.cov = 1e-3,
                         lower.alpha.cov = 1e-4,
                         upper.alpha.cov = 1e5,
                         focal.comp.matrix,
                         focal.covariates,
                         generate.errors = FALSE,
                         bootstrap.samples = 0){
  
  num.competitors <- dim(focal.comp.matrix)[2]
  name.competitors <- colnames(focal.comp.matrix)
  num.covariates <- ifelse(is.null(ncol(focal.covariates)),0,ncol(focal.covariates))
  name.covariates <- colnames(focal.covariates)
  
  # generate vector of initial parameters
  # depending on which ones we want to optimize
  
  # initialize parameters to optimize
  my.init.lambda <- NULL
  my.init.alpha <- NULL
  my.init.lambda.cov <- NULL
  my.init.alpha.cov <- NULL
  
  # initialize fixed terms list
  fixed.terms <- list(lambda = NULL, lambda.cov = NULL, alpha = NULL, alpha.cov = NULL)
  
  # only parameters on "param.list" will be optimized
  # otherwise, get their values and put them in the "fixed.term" list
  # lambda
  if("lambda" %in% param.list){
    my.init.lambda <- init.lambda
  }else{
    fixed.terms[["lambda"]] <- init.lambda
  }
  
  # sigma. This one is always present
  if(is.null(init.sigma)){
    my.init.sigma <- sd(log.fitness)
  }else{
    my.init.sigma <- init.sigma
  }
  
  # alpha, single value or matrix
  if("alpha" %in% param.list){
    my.init.alpha <- init.alpha
  }else{
    fixed.terms[["alpha"]] <- init.alpha
  }
  
  # lambda covariates
  if("lambda.cov" %in% param.list){
    my.init.lambda.cov <- init.lambda.cov
  }else{
    fixed.terms[["lambda.cov"]] <- init.lambda.cov
  }
  
  # alpha covariates, single value or matrix
  if("alpha.cov" %in% param.list){
    my.init.alpha.cov <- init.alpha.cov
  }else{
    fixed.terms[["alpha.cov"]] <- init.alpha.cov
  }
  
  # put them all together in a single vector, also lower and upper bounds
  init.par <- InitParams(init.lambda = my.init.lambda,
                         init.sigma = my.init.sigma,
                         init.alpha = my.init.alpha,
                         init.lambda.cov = my.init.lambda.cov,
                         init.alpha.cov = my.init.alpha.cov,
                         lower.lambda = lower.lambda,
                         upper.lambda = upper.lambda,
                         lower.sigma = lower.sigma,
                         upper.sigma = upper.sigma,
                         lower.alpha = lower.alpha,
                         upper.alpha = upper.alpha,
                         lower.lambda.cov = lower.lambda.cov,
                         upper.lambda.cov = upper.lambda.cov,
                         lower.alpha.cov = lower.alpha.cov,
                         upper.alpha.cov = upper.alpha.cov,
                         num.competitors = num.competitors,
                         num.covariates = num.covariates)
  
  # optim functions
  if(optim.method == "optim_NM"){
    
    optim.result <- optim(init.par$init.par, 
                          fitness.model, 
                          gr = NULL, 
                          method = "Nelder-Mead", 
                          # lower = lower.bounds,
                          # upper = upper.bounds,
                          control = list(), 
                          hessian = F,
                          param.list = param.list,
                          log.fitness = log.fitness, 
                          focal.comp.matrix = focal.comp.matrix,
                          num.covariates = num.covariates, 
                          num.competitors = num.competitors, 
                          focal.covariates = focal.covariates,
                          fixed.terms = fixed.terms)
    
  }else if(optim.method == "optim_L-BFGS-B"){
    
    optim.result <- optim(init.par$init.par, 
                          fitness.model, 
                          gr = NULL, 
                          method = "L-BFGS-B", 
                          lower = init.par$lower.bounds, 
                          upper = init.par$upper.bounds,
                          control = list(), 
                          hessian = F,
                          param.list = param.list,
                          log.fitness = log.fitness, 
                          focal.comp.matrix = focal.comp.matrix,
                          num.covariates = num.covariates, 
                          num.competitors = num.competitors, 
                          focal.covariates = focal.covariates,
                          fixed.terms = fixed.terms)
    
  }else if(optim.method == "nloptr_CRS2_LM"){
    
    optim.result <- nloptr(x0 = init.par$init.par,
                           eval_f = fitness.model,
                           opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e3),
                           lb = init.par$lower.bounds,
                           ub = init.par$upper.bounds,
                           param.list = param.list,
                           log.fitness = log.fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates,
                           fixed.terms = fixed.terms)
    
  }else if(optim.method == "nloptr_ISRES"){
    
    optim.result <- nloptr(x0 = init.par$init.par,
                           eval_f = fitness.model,
                           opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e3),
                           lb = init.par$lower.bounds,
                           ub = init.par$upper.bounds,
                           param.list = param.list,
                           log.fitness = log.fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates,
                           fixed.terms = fixed.terms)
    
  }else if(optim.method == "nloptr_DIRECT_L_RAND"){
    
    optim.result <- nloptr(x0 = init.par$init.par,
                           eval_f = fitness.model,
                           opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e3),
                           lb = init.par$lower.bounds,
                           ub = init.par$upper.bounds,
                           param.list = param.list,
                           log.fitness = log.fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates,
                           fixed.terms = fixed.terms)
    
  }else if(optim.method == "GenSA"){
    
    optim.result <- GenSA(par = init.par$init.par,
                          fn = fitness.model,
                          lower = init.par$lower.bounds,
                          upper = init.par$upper.bounds, 
                          control = list(maxit = 1e3), 
                          param.list = param.list,
                          log.fitness = log.fitness, 
                          focal.comp.matrix = focal.comp.matrix,
                          num.covariates = num.covariates, 
                          num.competitors = num.competitors, 
                          focal.covariates = focal.covariates,
                          fixed.terms = fixed.terms)
    
  }else if(optim.method == "hydroPSO"){
    
    # suppress annoying output??
    # sink("/dev/null")
    optim.result <- hydroPSO::hydroPSO(par = init.par$init.par,
                                       fn = fitness.model,
                                       lower = lower.bounds,
                                       upper = upper.bounds, 
                                       control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                       param.list = param.list,
                                       log.fitness = log.fitness, 
                                       focal.comp.matrix = focal.comp.matrix,
                                       num.covariates = num.covariates, 
                                       num.competitors = num.competitors, 
                                       focal.covariates = focal.covariates,
                                       fixed.terms = fixed.terms)
    
    
  }else if(optim.method == "DEoptimR"){
    
    optim.result <- DEoptimR::JDEoptim(lower = init.par$lower.bounds,upper = init.par$upper.bounds,fn = fitness.model,
                                       param.list = param.list,
                                       log.fitness = log.fitness, 
                                       focal.comp.matrix = focal.comp.matrix,
                                       num.covariates = num.covariates, 
                                       num.competitors = num.competitors, 
                                       focal.covariates = focal.covariates,
                                       fixed.terms = fixed.terms)
  }
  
  ##################################
  # gather the output from the method
  # if-else the method outputs optim-like values
  if(optim.method %in% c("optim_NM","optim_L-BGFS-B","DEoptimR","hydroPSO","GenSA")){
    
    print(paste("...",optim.method," finished with convergence status ",optim.result$convergence,sep=""))
    
    optim.params <- RetrieveParams(optim.params = optim.result$par,
                                   param.list = param.list,
                                   alpha.length = length(init.alpha),
                                   alpha.cov.length = length(init.alpha.cov),
                                   num.competitors = num.competitors,
                                   num.covariates = num.covariates)
    
    log.likelihood <- optim.result$value
    
  }else{ # methods with different nomenclature
    
    print(paste("...",optim.method," finished with convergence status ",optim.result$status,sep=""))
    
    optim.params <- RetrieveParams(optim.params = optim.result$solution,
                                   param.list = param.list,
                                   alpha.length = length(init.alpha),
                                   alpha.cov.length = length(init.alpha.cov),
                                   num.competitors = num.competitors,
                                   num.covariates = num.covariates)
    
    log.likelihood <- optim.result$objective
  }  
  
  #####################
  # standard errors via bootstrapping
  if(generate.errors){
    
    errors <- SEbootstrap(fitness.model = fitness.model,
                          optim.method = optim.method,
                          param.list = param.list,
                          fixed.terms = fixed.terms,
                          log.fitness = log.fitness,
                          init.par = init.par$init.par,
                          lower.bounds = init.par$lower.bounds,
                          upper.bounds = init.par$upper.bounds,
                          focal.comp.matrix = focal.comp.matrix,
                          focal.covariates = focal.covariates,
                          nsamples = bootstrap.samples)
  }else{
    errors <- rep(NA,length(init.par))
  }
  
  error.params <- RetrieveParams(optim.params = errors,
                                 param.list = param.list,
                                 alpha.length = length(init.alpha),
                                 alpha.cov.length = length(init.alpha.cov),
                                 num.competitors = num.competitors,
                                 num.covariates = num.covariates)
  
  # append names to the results
  lambda <- optim.params[["lambda"]]
  lambda.lower.error <- optim.params[["lambda"]]-1.96*error.params[["lambda"]]
  lambda.upper.error <- optim.params[["lambda"]]+1.96*error.params[["lambda"]]
  sigma <- ifelse(optim.params[["sigma"]] > 0, optim.params[["sigma"]], 1e-10)
  alpha <- optim.params[["alpha"]]
  alpha.lower.error <- optim.params[["alpha"]]-1.96*error.params[["alpha"]]
  alpha.upper.error <- optim.params[["alpha"]]+1.96*error.params[["alpha"]]
  if(!is.null(name.competitors)){
    names(alpha) <- name.competitors
    if(!is.null(alpha.lower.error)){
      names(alpha.lower.error) <- name.competitors
    }
    if(!is.null(alpha.upper.error)){
      names(alpha.upper.error) <- name.competitors
    }
  }
  lambda.cov <- optim.params[["lambda.cov"]]
  lambda.cov.lower.error <- optim.params[["lambda.cov"]]-1.96*error.params[["lambda.cov"]]
  lambda.cov.upper.error <- optim.params[["lambda.cov"]]+1.96*error.params[["lambda.cov"]]
  if(!is.null(name.covariates)){
    if(length(lambda.cov) == length(name.covariates)){
      names(lambda.cov) <- name.covariates
    }
    if(length(lambda.cov.lower.error) == length(name.covariates)){
      names(lambda.cov.lower.error) <- name.covariates
    }
    if(length(lambda.cov.upper.error) == length(name.covariates)){
      names(lambda.cov.upper.error) <- name.covariates
    }
  }
  alpha.cov <- optim.params[["alpha.cov"]]
  alpha.cov.lower.error <- optim.params[["alpha.cov"]]-1.96*error.params[["alpha.cov"]]
  alpha.cov.upper.error <- optim.params[["alpha.cov"]]+1.96*error.params[["alpha.cov"]]
  if(!is.null(name.competitors) & !is.null(name.covariates)){
    name.alpha.cov <- paste(rep(name.covariates,each = num.competitors),rep(name.competitors,num.covariates),sep="_")
    if(length(alpha.cov) == length(name.alpha.cov)){
      names(alpha.cov) <- name.alpha.cov
    }
    if(length(alpha.cov.lower.error) == length(name.alpha.cov)){
      names(alpha.cov.lower.error) <- name.alpha.cov
    }
    if(length(alpha.cov.upper.error) == length(name.alpha.cov)){
      names(alpha.cov.upper.error) <- name.alpha.cov
    }
  }

  return(list(lambda = lambda,
              lambda.lower.error = lambda.lower.error,
              lambda.upper.error = lambda.upper.error,
              sigma = sigma,
              alpha = alpha,
              alpha.lower.error = alpha.lower.error,
              alpha.upper.error = alpha.upper.error,
              lambda.cov = lambda.cov,
              lambda.cov.lower.error = lambda.cov.lower.error,
              lambda.cov.upper.error = lambda.cov.upper.error,
              alpha.cov = alpha.cov,
              alpha.cov.lower.error = alpha.cov.lower.error,
              alpha.cov.upper.error = alpha.cov.upper.error,
              log.likelihood = log.likelihood
  ))
  
}

