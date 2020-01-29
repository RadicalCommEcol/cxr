
source("R/cxr_check_input_er.R")
source("R/cxr_init_er_params.R")
source("R/cxr_retrieve_er_params.R")
source("R/cxr_return_init_length.R")
source('R/er_BH_lambdacov_none_effectcov_none_responsecov_none.R')

spdata <- data.frame(fitness = runif(10,0,1),f1 = round(runif(10,1,10)), f2 = round(runif(10,1,5)))
spdata2 <- data.frame(fitness = runif(10,0,1),f1 = round(runif(10,1,10)), f2 = round(runif(10,1,5)))

spdata$focal <- "f1"
spdata2$focal <- "f2"

splist <- list(f1 = spdata,f2 = spdata2)
spdf <- rbind(spdata,spdata2)

c1 <- data.frame(c1 = rnorm(10,1,.1))
# c2 <- data.frame(c2 = rnorm(10,1,.1))
clist <- list(f1 = c1, f2 = c1*2)
cdf <- rbind(c1,c1*2)

optimization_method <- "bobyqa"
model_family <- "BH"
data <- spdf
covariates <- NULL
lambda_cov_form <- effect_cov_form <- response_cov_form <- "none"
initial_values = list(lambda = 1, 
                      effect = 0.1, 
                      response = 0.1, 
                      lambda_cov = 0.1, 
                      effect_cov = 0.1, 
                      response_cov = 0.1,
                      sigma = 0.1)
lower_bounds <- list(lambda = 0, 
                     effect = 0, 
                     response = 0, 
                     lambda_cov = 0, 
                     effect_cov = 0, 
                     response_cov = 0,
                     sigma = 0)
upper_bounds <- list(lambda = 10, 
                     effect = 1, 
                     response = 1, 
                     lambda_cov = 1, 
                     effect_cov = 1, 
                     response_cov = 1,
                     sigma = 1)
bootstrap_samples <- 0
fixed_terms <- NULL


cxr_er_fit <- function(data, 
                       model_family = c("BH"),
                       covariates = NULL, 
                       optimization_method = c("BFGS", "CG", "Nelder-Mead", 
                                               "ucminf","L-BFGS-B", "nlm", "nlminb", 
                                               "Rcgmin", "Rvmmin", "spg", 
                                               "bobyqa", "nmkb", "hjkb",
                                               "nloptr_CRS2_LM","nloptr_ISRES",
                                               "nloptr_DIRECT_L_RAND","DEoptimR",
                                               "hydroPSO","GenSA"), 
                       lambda_cov_form = c("none","global"),
                       effect_cov_form = c("none","global"),
                       response_cov_form = c("none","global"),
                       initial_values = list(lambda = 0, 
                                             effect = 0, 
                                             response = 0, 
                                             lambda_cov = 0, 
                                             effect_cov = 0, 
                                             response_cov = 0,
                                             sigma = 0),
                       lower_bounds = NULL,
                       upper_bounds = NULL,
                       fixed_terms = NULL,
                       # errors
                       bootstrap_samples = 0
){
  
  # TODO add cxr:: to the internal functions once they are added
  
  # argument checks ---------------------------------------------------------
  
  optimization_method <- match.arg(optimization_method)
  lambda_cov_form <- match.arg(lambda_cov_form)
  effect_cov_form <- match.arg(effect_cov_form)
  response_cov_form <- match.arg(response_cov_form)
  
  # TODO fixed terms?
  
  # check installed packages for optimization method
  if (optimization_method %in% c("nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND") & !requireNamespace("nloptr", quietly = TRUE)) {
    stop("cxr_er_fit ERROR: Package \"nloptr\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "GenSA" & !requireNamespace("GenSA", quietly = TRUE)) {
    stop("cxr_er_fit ERROR: Package \"GenSA\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "hydroPSO" & !requireNamespace("hydroPSO", quietly = TRUE)) {
    stop("cxr_er_fit ERROR: Package \"hydroPSO\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "DEoptimR" & !requireNamespace("DEoptimR", quietly = TRUE)) {
    stop("cxr_er_fit ERROR: Package \"DEoptimR\" needed for the method selected to work.",
         call. = FALSE)
  }
  
  # check input data is a single dataframe with all sp
  # or a list with single-sp data
  # check input data
  # TODO focal = neighbours?
  data.ok <- cxr_check_input_er(data,covariates)
  if(!data.ok){
    stop("cxr_er_fit ERROR: check the consistency of your input data: 
    1) 'data' is either a list of dataframes or a single dataframe
    2) if 'data' is a list, each component is a dataframe; all
    dataframes have the same number of observations and same columns
    3) if 'data' is a dataframe, it has a 'focal' column, and all 
    focal species have the same number of observations
    4) 'data' and 'covariates' (if present) are of the same class and
    have the same number of observations
    5) no NAs are allowed either in 'data' or 'covariates'")
  }
  
  # check covariates if alpha_cov or lambda_cov are to be fit
  if(is.null(covariates) & (effect_cov_form != "none" | 
                            response_cov_form != "none" | 
                            lambda_cov_form != "none")){
    stop("cxr_er_fit ERROR: need to specify covariates if any of lambda_cov,
         effect_cov, or response_cov are to be fit")
  }
  
  # retrieve model ----------------------------------------------------------
  # character string giving the name of the model
  model_name <- paste("er_",model_family,"_lambdacov_",lambda_cov_form,"_effectcov_",effect_cov_form,"_responsecov_",response_cov_form,sep="")
  
  # try to retrieve the function from its name
  # using function "get"
  tryCatch({
    fitness_model <- get(model_name)
  }, error=function(e){cat("cxr_er_fit ERROR : model '",model_name,"' 
  could not be retrieved. Make sure it is defined and available in the cxr package 
                           or in the global environment")})
  
  # prepare data ------------------------------------------------------------
  
  # if list, create a single dataframe
  covdf <- NULL
  if(class(data) == "list"){
    spdf <- do.call(rbind,data)
    if(!is.null(covariates)){
      covdf <- do.call(rbind,covariates)
    }
  }else{
    spdf <- data
    covdf <- covariates
  }
  
  # fill up matrices
  # num.sp x num.observations. 1 if species is focal in a given observation, 0 otherwise
  target_all <- NULL
  # num.sp x num.observations. density of each species in each observation
  density_all <- NULL
  
  # fitness metric of the focal sp at each observation
  # log.fitness <- log(spdf$fitness)
  
  # species names and number
  sp.list <- as.character(unique(spdf$focal))
  num.sp <- length(sp.list)
  
  for(i.sp in 1:num.sp){
    
    target.my.sp <- integer(nrow(spdf))
    target.my.sp <- ifelse(spdf$focal == sp.list[i.sp],1,0)
    density.my.sp <- spdf[,sp.list[i.sp]]
    
    target_all <- rbind(target_all,target.my.sp)
    density_all <- rbind(density_all,density.my.sp)
  }
  rownames(target_all) <- sp.list
  rownames(density_all) <- sp.list
  
  # are covariates named?
  if(!is.null(covdf)){
    if(is.null(names(covdf))){
      names(covdf) <- paste("cov",1:ncol(covdf),sep="")
    }
  }

  # initial parameters
  # check wheter each model parameter is to be fitted or is fixed
  # note how initial_values also function as values 
  # for those parameters that are fixed
  
  # here I store values for fixed parameters
  fixed_parameters <- list()
  
  # here I store initial values for fitted parameters
  init_lambda <- NULL
  init_effect <- NULL
  init_response <- NULL
  init_lambda_cov <- NULL
  init_effect_cov <- NULL
  init_response_cov <- NULL
  
  if("lambda" %in% fixed_terms){
    fixed_parameters[["lambda"]] <- initial_values$lambda
  }else{
    if(length(initial_values$lambda) == 1){
      init_lambda <- rep(initial_values$lambda,num.sp)
    }else{
      init_lambda <- initial_values$lambda
    }
    names(init_lambda) <- paste("lambda_",sp.list,sep="")
  }
  
  if(!is.null(initial_values$sigma)){
    init_sigma <- initial_values$sigma
  }else{
    init_sigma <- 0.1
  }
  names(init_sigma) <- "sigma"
  
  if("effect" %in% fixed_terms){
    fixed_parameters[["effect"]] <- initial_values$effect
  }else{
    if(length(initial_values$effect) == 1){
      init_effect <- rep(initial_values$effect,num.sp)
    }else{
      init_effect <- initial_values$effect
    }
    names(init_effect) <- paste("effect_",sp.list,sep="")
  }
  
  if("response" %in% fixed_terms){
    fixed_parameters[["response"]] <- initial_values$response
  }else{
    if(length(initial_values$response) == 1){
      init_response <- rep(initial_values$response,num.sp)
    }else{
      init_response <- initial_values$response
    }
    names(init_response) <- paste("response_",sp.list,sep="")
  }
  
  # return_init_length is an auxiliary function
  # in case initial values are not of the same length of the expected parameters

  # TODO appropriate length is num.sp * num.cov
  if(lambda_cov_form != "none" & !is.null(covdf)){
    if("lambda_cov" %in% fixed_terms){
      fixed_parameters[["lambda_cov"]] <- cxr_return_init_length(lambda_cov_form,initial_values$lambda_cov,names(covdf))
    }else{
      init_lambda_cov <- cxr_return_init_length(lambda_cov_form,initial_values$lambda_cov,num.sp*names(covdf))
      names(init_lambda_cov) <- paste("lambda_cov_",names(covariates),sep="")
    }
  }

  if(effect_cov_form != "none" & !is.null(covdf)){
    if("effect_cov" %in% fixed_terms){
      fixed_parameters[["effect_cov"]] <- cxr_return_init_length(effect_cov_form,initial_values$effect_cov,names(covdf))
    }else{
      init_effect_cov <- cxr_return_init_length(effect_cov_form,initial_values$effect_cov,names(covdf))
      names(init_effect_cov) <- paste("effect_cov_",names(covariates),sep="")
    }
  }
  
  if(response_cov_form != "none" & !is.null(covdf)){
    if("response_cov" %in% fixed_terms){
      fixed_parameters[["response_cov"]] <- cxr_return_init_length(response_cov_form,initial_values$response_cov,names(covdf))
    }else{
      init_response_cov <- cxr_return_init_length(response_cov_form,initial_values$response_cov,names(covariates))
      names(init_response_cov) <- paste("response_cov_",names(covariates),sep="")
    }
  }
  
  # retrieve lower and upper bounds, if present
  
  lower_lambda <- NULL
  upper_lambda <- NULL
  lower_sigma <- NULL
  upper_sigma <- NULL
  lower_effect <- NULL
  upper_effect <- NULL
  lower_response<- NULL
  upper_response <- NULL
  lower_lambda_cov <- NULL
  upper_lambda_cov <- NULL
  lower_effect_cov <- NULL
  upper_effect_cov <- NULL
  lower_response_cov <- NULL
  upper_response_cov <- NULL
  
  if(!is.null(lower_bounds$lambda) & 
     !is.null(upper_bounds$lambda) &
     !"lambda" %in% fixed_terms){
    lower_lambda <- lower_bounds$lambda
    upper_lambda <- upper_bounds$lambda
  }
  
  # sigma bounds can be hidden from the user
  # only set if there are bounds for other params
  if(!is.null(lower_bounds) & 
     !is.null(upper_bounds)){
    lower_sigma <- 1e-10
    upper_sigma <- 1
  }
  
  if(!is.null(lower_bounds$effect) & 
     !is.null(upper_bounds$effect) &
     !"effect" %in% fixed_terms){
    lower_effect <- lower_bounds$effect
    upper_effect <- upper_bounds$effect
  }
  
  if(!is.null(lower_bounds$response) & 
     !is.null(upper_bounds$response) &
     !"response" %in% fixed_terms){
    lower_response <- lower_bounds$response
    upper_response <- upper_bounds$response
  }
  
  if(!is.null(lower_bounds$lambda_cov) &
     !is.null(upper_bounds$lambda_cov) &
     !"lambda_cov" %in% fixed_terms){
    lower_lambda_cov <- lower_bounds$lambda_cov
    upper_lambda_cov <- upper_bounds$lambda_cov
  }
  
  if(!is.null(lower_bounds$effect_cov) &
     !is.null(upper_bounds$effect_cov) &
     !"effect_cov" %in% fixed_terms){
    lower_effect_cov <- lower_bounds$effect_cov
    upper_effect_cov <- upper_bounds$effect_cov
  }
  
  if(!is.null(lower_bounds$response_cov) &
     !is.null(upper_bounds$response_cov) &
     !"response_cov" %in% fixed_terms){
    lower_response_cov <- lower_bounds$response_cov
    upper_response_cov <- upper_bounds$response_cov
  }
  
  # sort parameters for optim routine ---------------------------------------
  
  # TODO cxr_init_er_params
  
  init_par <- cxr_init_er_params(init_lambda = init_lambda,
                                 init_sigma = init_sigma,
                                 init_effect = init_effect,
                                 init_response = init_response,
                                 init_lambda_cov = init_lambda_cov,
                                 init_effect_cov = init_effect_cov,
                                 init_response_cov = init_response_cov,
                                 lower_lambda = lower_lambda,
                                 upper_lambda = upper_lambda,
                                 lower_sigma = lower_sigma,
                                 upper_sigma = upper_sigma,
                                 lower_effect = lower_effect,
                                 upper_effect = upper_effect,
                                 lower_response = lower_response,
                                 upper_response = upper_response,
                                 lower_lambda_cov = lower_lambda_cov,
                                 upper_lambda_cov = upper_lambda_cov,
                                 lower_effect_cov = lower_effect_cov,
                                 upper_effect_cov = upper_effect_cov,
                                 lower_response_cov = lower_response_cov,
                                 upper_response_cov = upper_response_cov)
  
  # fit parameters ----------------------------------------------------------

  optim_result <- NULL

  if(optimization_method %in% c("BFGS", "CG", "Nelder-Mead", "ucminf")){
    tryCatch({
      optim_result <- optimx::optimx(par = init_par$init_par, 
                                     fn = fitness_model, 
                                     gr = NULL, 
                                     method = optimization_method,
                                     # lower = init_par$lower_bounds,
                                     # upper = init_par$upper_bounds,
                                     control = list(), 
                                     hessian = F,
                                     fitness = log(spdf$fitness), 
                                     target = target_all,
                                     density = density_all,
                                     covariates = covdf,  
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method %in% c("L-BFGS-B", "nlm", "nlminb", 
                               "Rcgmin", "Rvmmin", "spg", 
                               "bobyqa", "nmkb", "hjkb")){
    tryCatch({
      optim_result <- optimx::optimx(par = init_par$init_par, 
                                     fn = fitness_model, 
                                     gr = NULL, 
                                     method = optimization_method,
                                     lower = init_par$lower_bounds,
                                     upper = init_par$upper_bounds,
                                     control = list(), 
                                     hessian = F,
                                     fitness = log(spdf$fitness), 
                                     target = target_all,
                                     density = density_all,
                                     covariates = covdf, 
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_CRS2_LM"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = init_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e4),
                                     lb = init_par$lower_bounds,
                                     ub = init_par$upper_bounds,
                                     fitness = log(spdf$fitness), 
                                     target = target_all,
                                     density = density_all,
                                     covariates = covdf,  
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_ISRES"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = init_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e4),
                                     lb = init_par$lower_bounds,
                                     ub = init_par$upper_bounds,
                                     fitness = log(spdf$fitness), 
                                     target = target_all,
                                     density = density_all,
                                     covariates = covdf, 
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_DIRECT_L_RAND"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = init_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e4),
                                     lb = init_par$lower_bounds,
                                     ub = init_par$upper_bounds,
                                     fitness = log(spdf$fitness), 
                                     target = target_all,
                                     density = density_all,
                                     covariates = covdf, 
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "GenSA"){
    tryCatch({
      optim_result <- GenSA::GenSA(par = init_par$init_par,
                                   fn = fitness_model,
                                   lower = init_par$lower_bounds,
                                   upper = init_par$upper_bounds, 
                                   control = list(maxit = 1e3), 
                                   fitness = log(spdf$fitness), 
                                   target = target_all,
                                   density = density_all,
                                   covariates = covdf, 
                                   fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "hydroPSO"){
    tryCatch({
      # suppress annoying output??
      # sink("/dev/null")
      optim_result <- hydroPSO::hydroPSO(par = init_par$init_par,
                                         fn = fitness_model,
                                         lower = init_par$lower_bounds,
                                         upper = init_par$upper_bounds, 
                                         control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                         fitness = log(spdf$fitness), 
                                         target = target_all,
                                         density = density_all,
                                         covariates = covdf, 
                                         fixed_parameters = fixed_parameters)

    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
    
  }else if(optimization_method == "DEoptimR"){
    tryCatch({
      optim_result <- DEoptimR::JDEoptim(lower = init_par$lower_bounds,
                                         upper = init_par$upper_bounds,
                                         fn = fitness_model,
                                         fitness = log(spdf$fitness), 
                                         target = target_all,
                                         density = density_all,
                                         covariates = covdf,  
                                         fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }
  
  # gather output -----------------------------------------------------------
  
  # this is cumbersome and boring because different methods 
  # have different output types
  
  if(is.null(optim_result)){
    optim_params <- cxr_retrieve_er_params(optim_params = rep(NA_real_,length(init_par$init_par)),
                                           lambda_length = length(init_lambda),
                                           effect_length = length(init_effect),
                                           response_length = length(init_response),
                                           lambda_cov_length = length(init_lambda_cov),
                                           effect_cov_length = length(init_effect_cov),
                                           response_cov_length = length(init_response_cov))
    llik <- NA_real_
  }else{
    if(optimization_method %in% c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", 
                           "nlminb", "Rcgmin", "Rvmmin", "spg", "ucminf", 
                           "bobyqa", "nmkb", "hjkb")){
      par.pos <- which(!names(optim_result) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
      outnames <- names(optim_result[,par.pos])
      outpar <- as.numeric(optim_result[,par.pos])
      names(outpar) <- outnames
      llik <- optim_result$value
    }else if(optimization_method %in% c("DEoptimR","hydroPSO","GenSA")){
      outpar <- optim_result$par
      names(outpar) <- names(init_par$init_par)
      log.likelihood <- optim_result$value
    }else{
      outpar <- optim_result$solution
      names(outpar) <- names(init_par$init_par)
      log.likelihood <- optim_result$objective
    }# if-else method
    
    optim_params <- cxr_retrieve_er_params(optim_params = outpar,
                                           lambda_length = length(init_lambda),
                                           effect_length = length(init_effect),
                                           response_length = length(init_response),
                                           lambda_cov_length = length(init_lambda_cov),
                                           effect_cov_length = length(init_effect_cov),
                                           response_cov_length = length(init_response_cov))
    
  }# if not null
  
  # calculate errors --------------------------------------------------------
  
  # TODO: well...

  if(bootstrap_samples > 0){

    # errors <- cxr_pm_bootstrap(fitness_model = fitness_model,
    #                            optimization_method = optimization_method,
    #                            data = data,
    #                            covariates = covariates,
    #                            init_par = init_par$init_par,
    #                            lower_bounds = init_par$lower_bounds,
    #                            upper_bounds = init_par$upper_bounds,
    #                            fixed_parameters = fixed_parameters,
    #                            bootstrap_samples = bootstrap_samples)
    
    error_params <- cxr_retrieve_er_params(optim_params = errors,
                                        lambda_length = length(init_lambda),
                                        alpha_length = length(init_alpha),
                                        lambda_cov_length = length(init_lambda_cov),
                                        alpha_cov_length = length(init_alpha_cov))
  }else{
    error_params <- list(lambda = NULL,
                         effect = NULL,
                         response = NULL,
                         lambda_cov = NULL,
                         effect_cov = NULL,
                         response_cov = NULL
                         )
  }
  
  # return cxr object ---------------------------------------------------
  
  list_names <- c("model_name",
                  "model",
                  "data",
                  "model_family",
                  "covariates",
                  "sp",
                  "optimization_method",
                  "initial_values",
                  "fixed_terms",
                  "lambda",
                  "effect",
                  "response",
                  "lambda_cov",
                  "effect_cov",
                  "response_cov",
                  "lambda_standard_error",
                  "effect_standard_error",
                  "response_standard_error",
                  "lambda_cov_standard_error",
                  "effect_cov_standard_error",
                  "response_cov_standard_error",
                  "log_likelihood")
  
  fit <- sapply(list_names,function(x) NULL)
  
  fit$model_name <- model_name
  fit$model <- fitness_model
  fit$data <- spdf
  fit$model_family <- model_family
  fit$covariates <- covdf
  fit$sp <- sp.list
  fit$optimization_method <- optimization_method
  fit$initial_values <- initial_values
  
  # for returning explicit NULL values
  if(!is.null(fixed_terms)){
    fit$fixed_terms <- fixed_terms
  }
  if(!is.null(optim_params$lambda)){
    fit$lambda <- optim_params$lambda
  }
  if(!is.null(optim_params$effect)){
    fit$effect <- optim_params$effect
  }
  if(!is.null(optim_params$response)){
    fit$response <- optim_params$response
  }
  if(!is.null(optim_params$lambda_cov)){
    fit$lambda_cov <- optim_params$lambda_cov
  }
  if(!is.null(optim_params$effect_cov)){
    fit$effect_cov <- optim_params$effect_cov
  }
  if(!is.null(optim_params$response_cov)){
    fit$response_cov <- optim_params$response_cov
  }
  if(!is.null(error_params$lambda)){
    fit$lambda <- error_params$lambda
  }
  if(!is.null(error_params$effect)){
    fit$effect <- error_params$effect
  }
  if(!is.null(error_params$response)){
    fit$response <- error_params$response
  }
  if(!is.null(error_params$lambda_cov)){
    fit$lambda_cov <- error_params$lambda_cov
  }
  if(!is.null(error_params$effect_cov)){
    fit$effect_cov <- error_params$effect_cov
  }
  if(!is.null(error_params$response_cov)){
    fit$response_cov <- error_params$response_cov
  }
  
  fit$log_likelihood <- llik
  # define two classes, cxr_pm_fit/cxr_er_fit
  
  class(fit) <- "cxr_er_fit"
  fit
}


# summary method ----------------------------------------------------------

summary.cxr_er_fit <- function(x){
  cat("Effect/Response model '",x$model_name,"' fitted for ",length(x$sp)," species, with ",nrow(x$data)/length(x$sp)," observations per species,", 
      "\nand ",ifelse(is.null(x$covariates),0,ncol(x$covariates))," covariates",
      ", using optimization method '",x$optimization_method,"'",
      # "\n* focal lambda: ",ifelse(is.null(x$lambda)," - not fit - ",x$lambda),
      # # "\n* mean alpha: ",ifelse(is.null(x$alpha)," - not fit - ",mean(x$alpha)),
      # "\n* mean lambda_cov: ",ifelse(is.null(x$lambda_cov),"- not fit - ",mean(x$lambda_cov)),
      # "\n* mean alpha_cov: ",ifelse(is.null(x$alpha_cov),"- not fit - ",mean(x$alpha_cov)),
      "\n* log-likelihood of the fit: ",x$log_likelihood,sep="")
  
}

