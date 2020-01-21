#' General optimization for population models
#' 
#' Wrapper for optimization procedures. It accepts a population dynamics model, defined as a function, and a series of parameters. 
#' It returns the optimal value for the parameters given a fitness metric and an optimization method. Optionally, bootstrap confidence intervals
#' can also be computed.
#'
#' @param fitness.model function giving the population dynamics model. Any functional form is allowed, but the model must be constrained
#' to free parameters \code{lambda} (fecundity of each sp in absence of competition), \code{alpha} (interaction coefficients),
#' \code{lambda.cov} (effect of covariates on lambda), \code{alpha.cov} (effect of covariates on alpha)  
#' @param optim.method optimization method to use. See vignette "Data and Model formats" for a list of available methods. 
#' @param param.list string vector giving the parameters that are to be optimized for the fitness model 
#' (to choose among "lambda", "alpha", "lambda.cov", and "alpha.cov").
#' @param log.fitness 1d vector, log of the fitness metric for every observation
#' @param init.lambda 1d vector, initial value of lambda
#' @param lower.lambda lower bound for lambda
#' @param upper.lambda upper bound for lambda
#' @param init.sigma initial value for sigma (standard deviation)
#' @param lower.sigma lower bound for sigma
#' @param upper.sigma upper bound for sigma
#' @param init.alpha initial value for the alpha vector/matrix
#' @param lower.alpha lower bound for alpha
#' @param upper.alpha upper bound for alpha
#' @param init.lambda.cov initial value for the lambda.cov matrix. Discarded if no covariates are given.
#' @param lower.lambda.cov lower bound for lambda.cov
#' @param upper.lambda.cov upper bound for lambda.cov
#' @param init.alpha.cov initial value for the alpha.cov matrix. Discarded if no covariates are given.
#' @param lower.alpha.cov lower bound for alpha.cov
#' @param upper.alpha.cov upper bound for alpha.cov
#' @param focal.comp.matrix matrix with observations in rows and neighbours in columns. Each cell is the number of neighbours
#' of a given species in a given observation.
#' @param focal.covariates optional matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation.
#' @param generate.errors boolean, whether to compute bootstrap errors for the fitted parameters. Note that, depending on 
#' the data, model, and optimization method, this may be computationally expensive.
#' @param bootstrap.samples how many bootstrap samples to compute.
#' @param verbose work in progress
#'
#' @return list with the fitted parameters, and the loglikelihood of the fit. If a parameter is taken as a constant, the list will return 
#' the original value given.
#' @import stats
#' @export
pm_optim <- function(fitness.model,
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
                         focal.covariates = NULL,
                         generate.errors = FALSE,
                         bootstrap.samples = 0,
                         verbose = FALSE){
  # some sanity checks
  # first, discard if method selected is not robust in optimx,
  # or if needs a gradient function
  if(optim.method %in% c("lbfgsb3","Rtnmin","snewton","snewtonm","hjn",
                         "lbfgs","newuoa","subplex")){
    stop(paste("pm_optim ERROR: Method ",optim.method," is not currently supported.",sep=""),
         call. = FALSE)
  }
  # also, throw an error if any extra package is needed
  if (optim.method == "spg" & !requireNamespace("BB", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"BB\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method == "ucminf" & !requireNamespace("ucminf", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"ucminf\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method %in% c("nmkb","hjkb") & !requireNamespace("dfoptim", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"dfoptim\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method %in% c("nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND") & !requireNamespace("nloptr", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"nloptr\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method == "GenSA" & !requireNamespace("GenSA", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"GenSA\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method == "hydroPSO" & !requireNamespace("hydroPSO", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"hydroPSO\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optim.method == "DEoptimR" & !requireNamespace("DEoptimR", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"DEoptimR\" needed for the method selected to work.",
         call. = FALSE)
  }
  
  if(verbose){
    #suppressWarnings(message(date()," -- pm_optim: fitting ",length(log.fitness)," observations with method ",optim.method," started"))
  }
  
  num.competitors <- dim(as.matrix(focal.comp.matrix))[2]
  name.competitors <- colnames(focal.comp.matrix)
  
  if(!is.null(focal.covariates)){
    num.covariates <- ncol(as.matrix(focal.covariates))
    name.covariates <- colnames(focal.covariates)
  }else{
    num.covariates <- 0
    name.covariates <- NULL
  }
  
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
    my.init.lambda <- c(lambda = init.lambda)
  }else{
    fixed.terms[["lambda"]] <- init.lambda
  }
  
  # sigma. This one is always present
  if(is.null(init.sigma)){
    my.init.sigma <- c(sigma = sd(log.fitness))
    # double check
    if(my.init.sigma > upper.sigma){
      my.init.sigma <- upper.sigma
    }else if(my.init.sigma < lower.sigma){
      my.init.sigma <- lower.sigma
    }
  }else{
    my.init.sigma <- c(sigma = init.sigma)
    # double check
    if(my.init.sigma > upper.sigma){
      my.init.sigma <- upper.sigma
    }else if(my.init.sigma < lower.sigma){
      my.init.sigma <- lower.sigma
    }
  }
  
  # alpha, single value or matrix
  if("alpha" %in% param.list){
    my.init.alpha <- init.alpha
    if(length(my.init.alpha) == 1){
      names(my.init.alpha) <- "alpha"
    }else{
      names(my.init.alpha) <- name.competitors
    }
  }else{
    fixed.terms[["alpha"]] <- init.alpha
  }
  
  # lambda covariates
  if("lambda.cov" %in% param.list){
    my.init.lambda.cov <- setNames(init.lambda.cov,name.covariates)
  }else{
    fixed.terms[["lambda.cov"]] <- init.lambda.cov
  }
  
  # alpha covariates, single value or matrix
  if("alpha.cov" %in% param.list){
    my.init.alpha.cov <- init.alpha.cov
    if(length(my.init.alpha.cov) == num.covariates){
      names(my.init.alpha.cov) <- name.covariates
    }else{
      names(my.init.alpha.cov) <- paste(rep(name.covariates,each = num.competitors),
                                        rep(name.competitors,num.covariates),sep="_")
    }
  }else{
    fixed.terms[["alpha.cov"]] <- init.alpha.cov
  }
  
  # put them all together in a single vector, 
  # also returning lower and upper bounds
  init.par <- cxr_init_params(init.lambda = my.init.lambda,
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
  
  # in case of error
  optim.result <- NULL
  # optim functions
  if(optim.method %in% c("BFGS", "CG", "Nelder-Mead", "ucminf")){
    tryCatch({
      optim.result <- optimx::optimx(par = init.par$init.par, 
                     fn = fitness.model, 
                     gr = NULL, 
                     method = optim.method,
                     # lower = init.par$lower.bounds,
                     # upper = init.par$upper.bounds,
                     control = list(), 
                     hessian = F,
                     param.list = param.list,
                     log.fitness = log.fitness, 
                     focal.comp.matrix = focal.comp.matrix,
                     num.covariates = num.covariates, 
                     num.competitors = num.competitors, 
                     focal.covariates = focal.covariates,
                     fixed.terms = fixed.terms)
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
  }else if(optim.method %in% c("L-BFGS-B", "nlm", "nlminb", 
                               "Rcgmin", "Rvmmin", "spg", 
                               "bobyqa", "nmkb", "hjkb")){
    tryCatch({
      optim.result <- optimx::optimx(init.par$init.par, 
                                     fitness.model, 
                                     gr = NULL, 
                                     method = optim.method, 
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
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
  }else if(optim.method == "nloptr_CRS2_LM"){
    tryCatch({
    optim.result <- nloptr::nloptr(x0 = init.par$init.par,
                           eval_f = fitness.model,
                           opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e4),
                           lb = init.par$lower.bounds,
                           ub = init.par$upper.bounds,
                           param.list = param.list,
                           log.fitness = log.fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates,
                           fixed.terms = fixed.terms)
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
  }else if(optim.method == "nloptr_ISRES"){
    tryCatch({
    optim.result <- nloptr::nloptr(x0 = init.par$init.par,
                           eval_f = fitness.model,
                           opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e4),
                           lb = init.par$lower.bounds,
                           ub = init.par$upper.bounds,
                           param.list = param.list,
                           log.fitness = log.fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates,
                           fixed.terms = fixed.terms)
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
  }else if(optim.method == "nloptr_DIRECT_L_RAND"){
    tryCatch({
    optim.result <- nloptr::nloptr(x0 = init.par$init.par,
                           eval_f = fitness.model,
                           opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e4),
                           lb = init.par$lower.bounds,
                           ub = init.par$upper.bounds,
                           param.list = param.list,
                           log.fitness = log.fitness, 
                           focal.comp.matrix = focal.comp.matrix,
                           num.covariates = num.covariates, 
                           num.competitors = num.competitors, 
                           focal.covariates = focal.covariates,
                           fixed.terms = fixed.terms)
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
  }else if(optim.method == "GenSA"){
    tryCatch({
    optim.result <- GenSA::GenSA(par = init.par$init.par,
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
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
  }else if(optim.method == "hydroPSO"){
    tryCatch({
    # suppress annoying output??
    # sink("/dev/null")
    optim.result <- hydroPSO::hydroPSO(par = init.par$init.par,
                                       fn = fitness.model,
                                       lower = init.par$lower.bounds,
                                       upper = init.par$upper.bounds, 
                                       control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                       param.list = param.list,
                                       log.fitness = log.fitness, 
                                       focal.comp.matrix = focal.comp.matrix,
                                       num.covariates = num.covariates, 
                                       num.competitors = num.competitors, 
                                       focal.covariates = focal.covariates,
                                       fixed.terms = fixed.terms)
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
    
  }else if(optim.method == "DEoptimR"){
    tryCatch({
    optim.result <- DEoptimR::JDEoptim(lower = init.par$lower.bounds,upper = init.par$upper.bounds,fn = fitness.model,
                                       param.list = param.list,
                                       log.fitness = log.fitness, 
                                       focal.comp.matrix = focal.comp.matrix,
                                       num.covariates = num.covariates, 
                                       num.competitors = num.competitors, 
                                       focal.covariates = focal.covariates,
                                       fixed.terms = fixed.terms)
    }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})
  }

  ##################################
  # gather the output from the method
  # if-else the method outputs optim-like values
  if(optim.method %in% c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", 
                         "nlminb", "Rcgmin", "Rvmmin", "spg", "ucminf", 
                         "bobyqa", "nmkb", "hjkb")){
    if(!is.null(optim.result)){
      row.names(optim.result) <- NULL
      par.pos <- which(!names(optim.result) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
      # rownames(optim.result) <- NULL
      optim.params <- cxr_retrieve_params(optim.params = as.numeric(optim.result[,par.pos]),
                                          param.list = param.list,
                                          alpha.length = length(init.alpha),
                                          alpha.cov.length = length(init.alpha.cov),
                                          num.competitors = num.competitors,
                                          num.covariates = num.covariates)
      
      log.likelihood <- optim.result$value
      
    }else{
      optim.params <- cxr_retrieve_params(optim.params = rep(NA_real_,length(init.par$init.par)),
                                          param.list = param.list,
                                          alpha.length = length(init.alpha),
                                          alpha.cov.length = length(init.alpha.cov),
                                          num.competitors = num.competitors,
                                          num.covariates = num.covariates)
      log.likelihood <- NA_real_
    }
  }else
  if(optim.method %in% c("DEoptimR","hydroPSO","GenSA")){
    
    if(verbose){
      #suppressWarnings(message(date()," -- pm_optim: method ",optim.method," completed with convergence status ",optim.result$convergence))
    }

    if(!is.null(optim.result)){
    
    optim.params <- cxr_retrieve_params(optim.params = optim.result$par,
                                   param.list = param.list,
                                   alpha.length = length(init.alpha),
                                   alpha.cov.length = length(init.alpha.cov),
                                   num.competitors = num.competitors,
                                   num.covariates = num.covariates)

    log.likelihood <- optim.result$value
    
    }else{
      optim.params <- cxr_retrieve_params(optim.params = rep(NA_real_,length(init.par$init.par)),
                                     param.list = param.list,
                                     alpha.length = length(init.alpha),
                                     alpha.cov.length = length(init.alpha.cov),
                                     num.competitors = num.competitors,
                                     num.covariates = num.covariates)
      log.likelihood <- NA_real_
    }
    
  }else{ # methods with different nomenclature
    
    if(verbose){
      #suppressWarnings(message(date()," -- pm_optim: ",length(log.fitness)," observations, method ",optim.method," completed with convergence status ",optim.result$status))
    }    
    if(!is.null(optim.result)){
    
    optim.params <- cxr_retrieve_params(optim.params = optim.result$solution,
                                   param.list = param.list,
                                   alpha.length = length(init.alpha),
                                   alpha.cov.length = length(init.alpha.cov),
                                   num.competitors = num.competitors,
                                   num.covariates = num.covariates)
    
    log.likelihood <- optim.result$objective
    
    }else{
      optim.params <- cxr_retrieve_params(optim.params = rep(NA_real_,length(init.par$init.par)),
                                     param.list = param.list,
                                     alpha.length = length(init.alpha),
                                     alpha.cov.length = length(init.alpha.cov),
                                     num.competitors = num.competitors,
                                     num.covariates = num.covariates)
      log.likelihood <- NA_real_
    }
  }  
  
  #####################
  # standard errors via bootstrapping
  if(generate.errors){
    
    if(verbose){
      #suppressWarnings(message(date()," -- pm_optim: generating standard errors for ",length(log.fitness)," observations and ",bootstrap.samples," bootstrap samples"))
    }
    
    errors <- cxr_pm_bootstrap(fitness.model = fitness.model,
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
    errors <- rep(NA_real_,length(init.par))
  }
  
  if(verbose){
    #suppressWarnings(message(date()," -- pm_optim: standard errors computed"))
  }
  
  error.params <- cxr_retrieve_params(optim.params = errors,
                                 param.list = param.list,
                                 alpha.length = length(init.alpha),
                                 alpha.cov.length = length(init.alpha.cov),
                                 num.competitors = num.competitors,
                                 num.covariates = num.covariates)
  
  # return list of results
  lambda <- optim.params[["lambda"]]
  lambda.lower.error <- optim.params[["lambda"]]-1.96*error.params[["lambda"]]
  lambda.upper.error <- optim.params[["lambda"]]+1.96*error.params[["lambda"]]
  
  if(!is.na(optim.params[["sigma"]])){
    sigma <- ifelse(optim.params[["sigma"]] > 0, optim.params[["sigma"]], 1e-10)
  }else{
    sigma <- NA_real_
  }
  
  alpha <- optim.params[["alpha"]]
  alpha.lower.error <- optim.params[["alpha"]]-1.96*error.params[["alpha"]]
  alpha.upper.error <- optim.params[["alpha"]]+1.96*error.params[["alpha"]]
  
  # named vectors
  if(!is.null(name.competitors) & !is.null(alpha)){
    if(length(alpha) == length(name.competitors)){
      names(alpha) <- name.competitors
      if(!is.null(alpha.lower.error)){
        names(alpha.lower.error) <- name.competitors
      }
      if(!is.null(alpha.upper.error)){
        names(alpha.upper.error) <- name.competitors
      }
    }
  }
  
  lambda.cov <- optim.params[["lambda.cov"]]
  lambda.cov.lower.error <- optim.params[["lambda.cov"]]-1.96*error.params[["lambda.cov"]]
  lambda.cov.upper.error <- optim.params[["lambda.cov"]]+1.96*error.params[["lambda.cov"]]
  
  # named vectors
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
  
  # named vectors
  if(!is.null(name.competitors) & !is.null(name.covariates)){
    
    if(length(alpha.cov) == num.covariates){
      name.alpha.cov <- name.covariates
    }else{
      name.alpha.cov <- paste(rep(name.covariates,each = num.competitors),
                              rep(name.competitors,num.covariates),sep="_")
    }
    # double-check
    if(length(alpha.cov) == length(name.alpha.cov)){
      names(alpha.cov) <- name.alpha.cov
    }
    if(length(alpha.cov.lower.error) == length(name.alpha.cov)){
      names(alpha.cov.lower.error) <- name.alpha.cov
    }
    if(length(alpha.cov.upper.error) == length(name.alpha.cov)){
      names(alpha.cov.upper.error) <- name.alpha.cov
    }
  }# !is.null

  return.list <- list(lambda = lambda,
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
                      log.likelihood = log.likelihood)
  
  return.list[lengths(return.list) == 0] <- NA_real_
  
  return.list

}

# character vector with all methods available, just for reference
all.methods <- c("BFGS","CG","Nelder-Mead","L-BFGS-B","nlm","nlminb","Rcgmin",
                 "Rvmmin","spg","ucminf","bobyqa","nmkb","hjkb",
                 "nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND",
                 "GenSA", "hydroPSO","DEoptimR")


