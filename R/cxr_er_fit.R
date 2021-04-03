#' General optimization for effect-response models
#' 
#' Estimates parameters of user-specified models of competitive effects and responses.
#'
#' @param data either a list of dataframes or a single dataframe. if 'data' is a list, each element is a dataframe with the following columns:
#' * fitness: fitness metric for each observation
#' * neighbours: named columns giving the number of neighbours of each column
#' the names of the list elements are taken to be the names of the focal species. 
#' 
#' 
#' If 'data' is a dataframe, it also needs a 'focal' column.
#' Regardless of the data structure, all focal species need to have the same number of observations (i.e. same number of rows),
#' and the set of neighbour species needs to be the same as the set of focal species, so that
#' the neighbours columns correspond to the names of the list elements or, if 'data' is a dataframe, 
#' to the values of the 'focal' column. Future versions will relax this requirement.
#' @param covariates a data structure equivalent to 'data', in which each column are the values of a covariate.
#' @param effect_cov_form form of the covariate effects on competitive effects. 
#' Either "none" (no covariate effects) or "global" (one estimate per covariate)
#' @param response_cov_form form of the covariate effects on competitive responses. 
#' Either "none" (no covariate effects) or "global" (one estimate per covariate)
#' @param initial_values list with components "lambda","effect","response", and optionally
#' "lambda_cov", "effect_cov", "response_cov", specifying the initial values
#' for numerical optimization. Single values are allowed.
#' @param lower_bounds optional list with single values for "lambda", "effect","response", 
#' and optionally "lambda_cov", "effect_cov", "response_cov".
#' @param upper_bounds optional list with single values for "lambda", "effect","response", 
#' and optionally "lambda_cov", "effect_cov", "response_cov".
#' @param fixed_terms optional list specifying which model parameters are fixed.
#' @inheritParams cxr_pm_fit
#' @md
#' @return an object of type 'cxr_er_fit' which is a list with the following components:
#' * model_name: string with the name of the fitness model
#' * model: model function
#' * data: data supplied
#' * taxa: names of the taxa fitted 
#' * covariates: covariate data supplied
#' * optimization_method: optimization method used
#' * initial_values: list with initial values
#' * fixed_terms: list with fixed terms
#' * lambda: fitted values for lambdas, or NULL if fixed
#' * effect: fitted values for competitive effects, or NULL if fixed
#' * response: fitted values for competitive responses, or NULL if fixed
#' * lambda_cov: fitted values for effect of covariates on lambdas, or NULL if fixed
#' * effect_cov: fitted values for effect of covariates on competitive effects, or NULL if fixed
#' * response_cov: fitted values for effect of covariates on competitive responses, or NULL if fixed
#' * lambda_standard_error: standard errors for lambdas, if calculated
#' * effect_standard_error: standard errors for competitive effects, if calculated
#' * response_standard_error: standard errors for competitive responses, if calculated
#' * lambda_cov_standard_error: standard errors for effect of covariates on lambdas, if calculated
#' * effect_cov_standard_error: standard errors for effect of covariates on competitive effects, if calculated
#' * response_cov_standard_error: standard errors for effect of covariates on competitive responses, if calculated
#' * log_likelihood: log-likelihood of the fits
#' @export
#'
#' @examples
#' \donttest{
#' # fit three species at once
#' data("neigh_list")
#' # these species all have >250 observations
#' example_sp <- c("BEMA","LEMA","HOMA")
#' sp.pos <- which(names(neigh_list) %in% example_sp)
#' data <- neigh_list[sp.pos]
#' n.obs <- 250
#' # keep only fitness and neighbours columns
#' for(i in 1:length(data)){
#'   data[[i]] <- data[[i]][1:n.obs,c(2,sp.pos+2)]#2:length(data[[i]])]
#' }
#' 
#' # covariates: salinity
#' data("salinity_list")
#' salinity <- salinity_list[example_sp]
#' # keep only salinity column
#' for(i in 1:length(salinity)){
#'   salinity[[i]] <- salinity[[i]][1:n.obs,2:length(salinity[[i]])]
#' }
#' 
#' initial_values = list(lambda = 1, 
#'                      effect = 1, 
#'                      response = 1
#'                      # lambda_cov = 0, 
#'                      # effect_cov = 0, 
#'                      # response_cov = 0
#')
#'lower_bounds = list(lambda = 0, 
#'                    effect = 0, 
#'                    response = 0
#'                    # lambda_cov = 0, 
#'                    # effect_cov = 0, 
#'                    # response_cov = 0
#')
#'upper_bounds = list(lambda = 100, 
#'                     effect = 10, 
#'                     response = 10
#'                    # lambda_cov = 0, 
#'                    # effect_cov = 0, 
#'                    # response_cov = 0
#' )
#' 
#' er_3sp <- cxr_er_fit(data = data,
#'                      model_family = "BH",
#'                      # fit without covariates, 
#'                      # as it may be very computationally expensive
#'                      # covariates = salinity,
#'                      optimization_method = "bobyqa",
#'                      lambda_cov_form = "none",
#'                      effect_cov_form = "none",
#'                      response_cov_form = "none",
#'                      initial_values = initial_values,
#'                      lower_bounds = lower_bounds,
#'                      upper_bounds = upper_bounds,
#'                      # syntaxis for fixed values
#'                      # fixed_terms = list("response"),
#'                      bootstrap_samples = 3)
#' # brief summary
#' summary(er_3sp)
#' }
cxr_er_fit <- function(data, 
                       model_family = c("BH"),
                       covariates = NULL, 
                       optimization_method = c("Nelder-Mead", "BFGS", "CG", 
                                               "ucminf","L-BFGS-B", "nlm", "nlminb", 
                                               "Rcgmin", "Rvmmin", "spg", 
                                               "bobyqa", "nmkb", "hjkb",
                                               "nloptr_CRS2_LM","nloptr_ISRES",
                                               "nloptr_DIRECT_L_RAND","DEoptimR",
                                               "hydroPSO","GenSA"), 
                       lambda_cov_form = c("none","global"),
                       effect_cov_form = c("none","global"),
                       response_cov_form = c("none","global"),
                       initial_values = list(lambda = 1, 
                                             effect = 1, 
                                             response = 1, 
                                             lambda_cov = 0, 
                                             effect_cov = 0, 
                                             response_cov = 0),
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
  data.ok <- cxr_check_input_er(data,covariates)
  if(!data.ok){
    stop("cxr_er_fit ERROR: check the consistency of your input data: 
    1) 'data' is either a list of dataframes or a single dataframe.
    2) if 'data' is a list, each component is a dataframe; all
    dataframes have the same number of observations and same columns.
    3) if 'data' is a dataframe, it has a 'focal' column, and all 
    focal species have the same number of observations.
    4) The set of focal species is given by the names of the list elements
    or by the 'focal' column. These species need to be the same as the
    ones in the neighbour columns, and in the same order.
    4) 'data' and 'covariates' (if present) are of the same class and
    have the same number of observations.
    5) All variables are integer/numeric, with no NAs.")
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
  model_name <- paste(model_family,"_er","_lambdacov_",lambda_cov_form,"_effectcov_",effect_cov_form,"_responsecov_",response_cov_form,sep="")
  
  fitness_model <- try(get(model_name),silent = TRUE)
  if(inherits(fitness_model,"try-error")){
    stop(paste("cxr_er_fit ERROR: model '",model_name,"' could not be retrieved. 
  Make sure it is defined and available in the cxr package or in the global environment.\n",sep=""))
  }
  
  # check that lower/upper bounds are provided if the method requires it
  bound.ok <- cxr_check_method_boundaries(optimization_method,lower_bounds,upper_bounds, type = "er")
  if(!bound.ok){
    stop("cxr_er_fit ERROR: check the optimization method selected and lower/upper bounds.
         The following methods require explicit lower and upper parameter boundaries to be set:
         L-BFGS-B, nlm, nlminb, Rcgmin, Rvmmin, spg, bobyqa, nmkb, hjkb, nloptr_CRS2_LM,
         nloptr_ISRES, nloptr_DIRECT_L_RAND, GenSA, hydroPSO, DEoptimR.")
  }
  
  # warning if initial values are default ones
  if(identical(initial_values,list(lambda = 1, 
                                   effect = 1, 
                                   response = 1, 
                                   lambda_cov = 0, 
                                   effect_cov = 0, 
                                   response_cov = 0))){
    message("cxr_er_fit: Using default initial values. Note that these may not be appropriate for your data/model, or
    for the optimization method selected.")
  }
  
  message("cxr_er_fit: note that computation time grows exponentially with number of species and covariates.")
  
  # prepare data ------------------------------------------------------------
  
  # if list, create a single dataframe
  covdf <- NULL
  if(inherits(data,"list")){
    spdf <- do.call(rbind,data)
    spdf$focal <- rep(names(data),each = nrow(data[[1]]))
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
  
  # species names and number
  sp.list <- as.character(unique(spdf$focal))
  num.sp <- length(sp.list)
  
  for(i.sp in 1:num.sp){
    
    target.my.sp <- integer(nrow(spdf))
    target.my.sp <- ifelse(spdf$focal == sp.list[i.sp],1,0)
    density.my.sp <- spdf[[sp.list[i.sp]]]
    
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
    if(length(initial_values$lambda) == 1){
      fixed_parameters[["lambda"]] <- rep(initial_values$lambda,num.sp)
    }else{
      fixed_parameters[["lambda"]] <- initial_values$lambda
    }
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
    if(length(initial_values$effect) == 1){
      fixed_parameters[["effect"]] <- rep(initial_values$effect,num.sp)
    }else{
      fixed_parameters[["effect"]] <- initial_values$effect
    }
  }else{
    if(length(initial_values$effect) == 1){
      init_effect <- rep(initial_values$effect,num.sp)
    }else{
      init_effect <- initial_values$effect
    }
    names(init_effect) <- paste("effect_",sp.list,sep="")
  }
  
  if("response" %in% fixed_terms){
    if(length(initial_values$response) == 1){
      fixed_parameters[["response"]] <- rep(initial_values$response,num.sp)
    }else{
      fixed_parameters[["response"]] <- initial_values$response
    }
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
  
  if(lambda_cov_form != "none" & !is.null(covdf)){
    
    lc_names <- as.vector(t(outer("lambda_cov", sp.list, paste, sep="_"))) 
    lc_names <- as.vector(t(outer(lc_names, names(covdf), paste, sep="_"))) 
    
    if("lambda_cov" %in% fixed_terms){
      fixed_parameters[["lambda_cov"]] <- cxr_return_init_length(lambda_cov_form,
                                                                 initial_values$lambda_cov,
                                                                 lc_names,"er")
    }else{
      init_lambda_cov <- cxr_return_init_length(lambda_cov_form,
                                                initial_values$lambda_cov,
                                                lc_names,"er")
    }
  }
  
  if(effect_cov_form != "none" & !is.null(covdf)){
    ec_names <- as.vector(t(outer("effect_cov", sp.list, paste, sep="_"))) 
    ec_names <- as.vector(t(outer(ec_names, names(covdf), paste, sep="_"))) 
    if("effect_cov" %in% fixed_terms){
      fixed_parameters[["effect_cov"]] <- cxr_return_init_length(effect_cov_form,
                                                                 initial_values$effect_cov,
                                                                 ec_names,"er")
    }else{
      init_effect_cov <- cxr_return_init_length(effect_cov_form,
                                                initial_values$effect_cov,
                                                ec_names,"er")
    }
  }
  
  if(response_cov_form != "none" & !is.null(covdf)){
    rc_names <- as.vector(t(outer("response_cov", sp.list, paste, sep="_"))) 
    rc_names <- as.vector(t(outer(rc_names, names(covdf), paste, sep="_"))) 
    if("response_cov" %in% fixed_terms){
      fixed_parameters[["response_cov"]] <- cxr_return_init_length(response_cov_form,
                                                                   initial_values$response_cov,
                                                                   rc_names,"er")
    }else{
      init_response_cov <- cxr_return_init_length(response_cov_form,
                                                  initial_values$response_cov,
                                                  rc_names,"er")
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
    lower_sigma <- 1e-5
    upper_sigma <- 1e5
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
  
  # TODO beware in previous versions order was lambda-e-r-lcov-ecov-rcov
  # it is apparently all good now, but it may pop up somewhere?
  
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
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
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
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
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
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
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
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
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
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
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
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
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
      
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
    
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
    }, error=function(e){cat("cxr_er_fit optimization ERROR :",conditionMessage(e), "\n")})
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
  
  if(bootstrap_samples > 0){
    
    errors <- cxr_er_bootstrap(fitness_model = fitness_model,
                               optimization_method = optimization_method,
                               data = spdf,
                               covariates = covdf,
                               init_par = init_par$init_par,
                               lower_bounds = init_par$lower_bounds,
                               upper_bounds = init_par$upper_bounds,
                               fixed_parameters = fixed_parameters,
                               bootstrap_samples = bootstrap_samples)
    
    error_params <- cxr_retrieve_er_params(optim_params = errors,
                                           lambda_length = length(init_lambda),
                                           effect_length = length(init_effect),
                                           response_length = length(init_response),
                                           lambda_cov_length = length(init_lambda_cov),
                                           effect_cov_length = length(init_effect_cov),
                                           response_cov_length = length(init_response_cov))
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
                  "taxa",
                  "covariates",
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
  # fit$model_family <- model_family
  fit$covariates <- covdf
  fit$taxa <- sp.list
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
    fit$lambda_standard_error <- error_params$lambda
  }
  if(!is.null(error_params$effect)){
    fit$effect_standard_error <- error_params$effect
  }
  if(!is.null(error_params$response)){
    fit$response_standard_error <- error_params$response
  }
  if(!is.null(error_params$lambda_cov)){
    fit$lambda_cov_standard_error <- error_params$lambda_cov
  }
  if(!is.null(error_params$effect_cov)){
    fit$effect_cov_standard_error <- error_params$effect_cov
  }
  if(!is.null(error_params$response_cov)){
    fit$response_cov_standard_error <- error_params$response_cov
  }
  
  fit$log_likelihood <- llik
  
  class(fit) <- "cxr_er_fit"
  
  if(!is.null(fit$lambda) & !is.null(lower_lambda) & !is.null(upper_lambda)){
    if(any(fit$lambda == lower_lambda) | any(fit$lambda == upper_lambda)){
      message("cxr_er_fit: A fitted lambda is equal to lower or upper bounds. Consider refitting
              with different boundaries.")
    }
  }
  
  if(!is.null(fit$response) & !is.null(lower_response) & !is.null(upper_response)){
    if(any(fit$response == lower_response) | any(fit$response == upper_response)){
      message("cxr_er_fit: One or more fitted responses are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$effect) & !is.null(lower_effect) & !is.null(upper_effect)){
    if(any(fit$effect == lower_effect) | any(fit$effect == upper_effect)){
      message("cxr_er_fit: One or more fitted effects are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$lambda_cov) & !is.null(lower_lambda_cov) & !is.null(upper_lambda_cov)){
    if(any(fit$lambda_cov == lower_lambda_cov) | any(fit$lambda_cov == upper_lambda_cov)){
      message("cxr_er_fit: A fitted lambda_cov is equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$response_cov) & !is.null(lower_response_cov) & !is.null(upper_response_cov)){
    if(any(fit$response_cov == lower_response_cov) | any(fit$response_cov == upper_response_cov)){
      message("cxr_er_fit: One or more fitted response_covs are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$effect_cov) & !is.null(lower_effect_cov) & !is.null(upper_effect_cov)){
    if(any(fit$effect_cov == lower_effect_cov) | any(fit$effect_cov == upper_effect_cov)){
      message("cxr_er_fit: One or more fitted effects are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  fit
}
