#' General optimization for population models
#' 
#' Estimates parameters of user-specified population dynamics models.
#'
#' @param data dataframe with observations in rows and two sets of columns:
#' * fitness: fitness metric for the focal individual
#' * neighbours: columns with user-defined names, giving number of neighbours for each group
#' @param model_family family of model to use. Available families are BH (Beverton-Holt) as default.
#' Users may define their own families and models (see vignette XXXXX).
#' @param covariates optional named matrix or dataframe with observations (rows) of any number of environmental covariates (columns)
#' @param optimization_method numerical optimization method
#' @param alpha_form what form does the alpha parameter take? one of "none" (no alpha in the model), 
#' "global" (a single alpha for all pairwise interactions), or "pairwise" (one alpha value for every interaction)
#' @param lambda_cov_form form of the covariate effects on lambda. Either "none" (no covariate effects) or "global" (one estimate per covariate)
#' @param alpha_cov_form form of the covariate effects on alpha. One of "none" (no covariate effects), "global" (one estimate per covariate on every alpha),
#' or "pairwise" (one estimate per covariate and pairwise alpha)
#' @param initial_values list with components "lambda","alpha","lambda_cov", "alpha_cov", specifying the initial values
#' for numerical optimization. Single values are allowed.
#' @param lower_bounds optional list with single values for "lambda","alpha","lambda_cov", "alpha_cov".
#' @param upper_bounds optional list with single values for "lambda","alpha","lambda_cov", "alpha_cov".
#' @param fixed_terms optional list specifying which model parameters are fixed, among "lambda","alpha","lambda_cov", and "alpha_cov".
#' @param bootstrap_samples number of bootstrap samples for error calculation. Defaults to 0, i.e. no error is calculated.
#' @return an object of type 'cxr_pm_fit' which is a list with the following components:
#' * model_name: string with the name of the fitness model
#' * model: model function
#' * data: data supplied 
#' * covariates: covariate data supplied
#' * optimization_method: optimization method used
#' * initial_values: list with initial values
#' * fixed_terms: list with fixed terms
#' * lambda: fitted value for lambda, or NULL if fixed
#' * alpha: fitted value(s) for alpha, or NULL if fixed
#' * lambda_cov: fitted value(s) for lambda_cov, or NULL if fixed
#' * alpha_cov: fitted value(s) for alpha_cov, or NULL if fixed
#' * lambda_standard_error: standard error for lambda, if computed
#' * alpha_standard_error: standard error for alpha, if computed
#' * lambda_cov_standard_error: standard error for lambda_cov, if computed
#' * alpha_cov_standard_error: standard error for alpha_cov, if computed
#' * log_likelihood: log-likelihood of the fit
#' @export
#' @md
#' @examples
#' data("neigh_list")
#' # data for a single species, keep only fitness and neighbours columns
#' sp_data <- neigh_list[[1]][2:ncol(neigh_list[[1]])]
#' \dontrun{
#'   sp_fit <- cxr_pm_fit(data = sp_data,
#'                        optimization_method = "bobyqa",
#'                        alpha_form = "pairwise",
#'                        lambda_cov_form = "none",
#'                        alpha_cov_form = "none",
#'                        initial_values = list(lambda = 1,alpha = 0.1),
#'                        lower_bounds = list(lambda = 0,alpha = 0),
#'                        upper_bounds = list(lambda = 100,alpha = 1),
#'                        bootstrap_samples = 3)
#'   summary(sp_fit)
#' }
cxr_pm_fit <- function(data, 
                       model_family = c("BH"),
                       covariates = NULL, 
                       optimization_method = c("BFGS", "CG", "Nelder-Mead", 
                                               "ucminf","L-BFGS-B", "nlm", "nlminb", 
                                               "Rcgmin", "Rvmmin", "spg", 
                                               "bobyqa", "nmkb", "hjkb",
                                               "nloptr_CRS2_LM","nloptr_ISRES",
                                               "nloptr_DIRECT_L_RAND","DEoptimR",
                                               "hydroPSO","GenSA"), 
                       alpha_form = c("none","global","pairwise"), 
                       lambda_cov_form = c("none","global"),
                       alpha_cov_form = c("none","global","pairwise"),
                       initial_values = list(lambda = 0, alpha = 0, lambda_cov = 0, alpha_cov = 0),
                       lower_bounds = NULL,
                       upper_bounds = NULL,
                       fixed_terms = NULL,
                       # errors
                       bootstrap_samples = 0
){
  
  # TODO add cxr:: to the internal functions once they are added
  
  # argument checks ---------------------------------------------------------
  
  optimization_method <- match.arg(optimization_method)
  alpha_form <- match.arg(alpha_form)
  lambda_cov_form <- match.arg(lambda_cov_form)
  alpha_cov_form <- match.arg(alpha_cov_form)
  
  # TODO fixed terms?
  
  # check installed packages for optimization method
  if (optimization_method %in% c("nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND") & !requireNamespace("nloptr", quietly = TRUE)) {
    stop("cxr_pm_fit ERROR: Package \"nloptr\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "GenSA" & !requireNamespace("GenSA", quietly = TRUE)) {
    stop("cxr_pm_fit ERROR: Package \"GenSA\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "hydroPSO" & !requireNamespace("hydroPSO", quietly = TRUE)) {
    stop("cxr_pm_fit ERROR: Package \"hydroPSO\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "DEoptimR" & !requireNamespace("DEoptimR", quietly = TRUE)) {
    stop("cxr_pm_fit ERROR: Package \"DEoptimR\" needed for the method selected to work.",
         call. = FALSE)
  }
  
  # check input data
  data.ok <- cxr_check_input_data(data,covariates)
  if(!data.ok){
    stop("cxr_pm_fit ERROR: check the consistency of your input data: 
    1) All variables are integer/numeric, with no NAs; 
    2) first column in 'data' is named 'fitness'; 
    3) abundances of at least one neighbour species in 'data';
    4) data and covariates (if present) have the same number of observations")
  }
  
  # check covariates if alpha_cov or lambda_cov are to be fit
  if(is.null(covariates) & (alpha_cov_form != "none" | lambda_cov_form != "none")){
    stop("cxr_pm_fit ERROR: need to specify covariates if lambda_cov and/or alpha_cov are to be fit")
  }
  
  # retrieve model ----------------------------------------------------------
  # character string giving the name of the model
  model_name <- paste("pm_",model_family,"_alpha_",alpha_form,"_lambdacov_",lambda_cov_form,"_alphacov_",alpha_cov_form,sep="")
  
  # try to retrieve the function from its name
  # using function "get"
  fitness_model <- try(get(model_name),silent = TRUE)
  if(class(fitness_model) == "try-error"){
    stop(paste("cxr_pm_fit ERROR: model '",model_name,"' could not be retrieved. 
  Make sure it is defined and available in the cxr package or in the global environment.\n",sep=""))
  }
  
  # check that lower/upper bounds are provided if the method requires it
  bound.ok <- cxr_check_method_boundaries(optimization_method,lower_bounds,upper_bounds, type = "pm")
  if(!bound.ok){
    stop("cxr_pm_fit ERROR: check the optimization method selected and lower/upper bounds.
         The following methods require explicit lower and upper parameter boundaries to be set:
         L-BFGS-B, nlm, nlminb, Rcgmin, Rvmmin, spg, bobyqa, nmkb, hjkb, nloptr_CRS2_LM,
         nloptr_ISRES, nloptr_DIRECT_L_RAND, GenSA, hydroPSO, DEoptimR.")
  }
  
  # warning if initial values are not set
  if(identical(initial_values,list(lambda = 0, alpha = 0, lambda_cov = 0, alpha_cov = 0))){
    message("cxr_pm_fit: Using default initial values. Note that these may not be appropriate for your data/model, or
    for the optimization method selected.")
  }
  
  # prepare data ------------------------------------------------------------
  # neighbour matrix
  
  # just to avoid a note in R CMD CHECK
  dropname <- "fitness"
  neigh_matrix <- as.matrix(data[ , !(names(data) %in% dropname)])
  # neigh_matrix <- subset(data, select = -c(fitness))
  # neigh_matrix <- as.matrix(neigh_matrix)
  
  # neighbour species?
  neigh <- colnames(neigh_matrix)
  
  # are covariates named?
  if(!is.null(covariates)){
    if(is.null(names(covariates))){
      names(covariates) <- paste("cov",1:ncol(covariates),sep="")
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
  init_alpha <- NULL
  init_lambda_cov <- NULL
  init_alpha_cov <- NULL
  
  if("lambda" %in% fixed_terms){
    fixed_parameters[["lambda"]] <- initial_values$lambda
  }else{
    init_lambda <- initial_values$lambda
    names(init_lambda) <- "lambda"
  }
  
  if(!is.null(initial_values$sigma)){
    init_sigma <- initial_values$sigma
  }else{
    init_sigma <- 0.1
  }
  names(init_sigma) <- "sigma"
  
  # return_init_length is an auxiliary function
  # in case initial values are not of the same length of the expected parameters
  # e.g. if we want to fit pairwise alphas but only provide a single initial value
  
  if(alpha_form != "none"){
    if("alpha" %in% fixed_terms){
      fixed_parameters[["alpha"]] <- cxr_return_init_length(alpha_form,
                                                            initial_values$alpha,
                                                            neigh,"pm")
    }else{
      init_alpha <- cxr_return_init_length(alpha_form,
                                           initial_values$alpha,
                                           neigh,"pm")
    }
  }
  
  if(lambda_cov_form != "none" & !is.null(covariates)){
    if("lambda_cov" %in% fixed_terms){
      fixed_parameters[["lambda_cov"]] <- cxr_return_init_length(lambda_cov_form,
                                                                 initial_values$lambda_cov,
                                                                 names(covariates),"pm")
    }else{
      init_lambda_cov <- cxr_return_init_length(lambda_cov_form,
                                                initial_values$lambda_cov,
                                                names(covariates),"pm")
      names(init_lambda_cov) <- paste("lambda_cov_",names(covariates),sep="")
    }
  }
  
  if(alpha_cov_form != "none" & !is.null(covariates)){
    if(alpha_cov_form == "global"){
      name.alpha.cov <- paste("alpha_cov_",names(covariates),sep="")
    }else{
      name.alpha.cov <- paste("alpha_cov",rep(names(covariates),each = length(neigh)),rep(neigh,ncol(covariates)),sep="_")
    }
    if("alpha_cov" %in% fixed_terms){
      fixed_parameters[["alpha_cov"]] <- cxr_return_init_length(alpha_cov_form,
                                                                initial_values$alpha_cov,
                                                                name.alpha.cov,"pm")
    }else{
      init_alpha_cov <- cxr_return_init_length(alpha_cov_form,
                                               initial_values$alpha_cov,
                                               name.alpha.cov,"pm")
    }
  }
  
  # retrieve lower and upper bounds, if present
  
  lower_lambda <- NULL
  upper_lambda <- NULL
  lower_sigma <- NULL
  upper_sigma <- NULL
  lower_alpha <- NULL
  upper_alpha <- NULL
  lower_lambda_cov <- NULL
  upper_lambda_cov <- NULL
  lower_alpha_cov <- NULL
  upper_alpha_cov <- NULL
  
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
  
  if(!is.null(lower_bounds$alpha) & 
     !is.null(upper_bounds$alpha) &
     !"alpha" %in% fixed_terms){
    lower_alpha <- lower_bounds$alpha
    upper_alpha <- upper_bounds$alpha
  }
  
  if(!is.null(lower_bounds$lambda_cov) &
     !is.null(upper_bounds$lambda_cov) &
     !"lambda_cov" %in% fixed_terms){
    lower_lambda_cov <- lower_bounds$lambda_cov
    upper_lambda_cov <- upper_bounds$lambda_cov
  }
  
  if(!is.null(lower_bounds$alpha_cov) &
     !is.null(upper_bounds$alpha_cov) &
     !"alpha_cov" %in% fixed_terms){
    lower_alpha_cov <- lower_bounds$alpha_cov
    upper_alpha_cov <- upper_bounds$alpha_cov
  }
  
  # sort parameters for optim routine ---------------------------------------
  init_par <- cxr_init_params(init_lambda = init_lambda,
                              init_sigma = init_sigma,
                              init_alpha = init_alpha,
                              init_lambda_cov = init_lambda_cov,
                              init_alpha_cov = init_alpha_cov,
                              lower_lambda = lower_lambda,
                              upper_lambda = upper_lambda,
                              lower_sigma = lower_sigma,
                              upper_sigma = upper_sigma,
                              lower_alpha = lower_alpha,
                              upper_alpha = upper_alpha,
                              lower_lambda_cov = lower_lambda_cov,
                              upper_lambda_cov = upper_lambda_cov,
                              lower_alpha_cov = lower_alpha_cov,
                              upper_alpha_cov = upper_alpha_cov)
  
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
                                     fitness = log(data$fitness), 
                                     neigh_matrix = neigh_matrix,
                                     covariates = covariates, 
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
                                     fitness = log(data$fitness), 
                                     neigh_matrix = neigh_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_CRS2_LM"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = init_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e4),
                                     lb = init_par$lower_bounds,
                                     ub = init_par$upper_bounds,
                                     fitness = log(data$fitness), 
                                     neigh_matrix = neigh_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_ISRES"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = init_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e4),
                                     lb = init_par$lower_bounds,
                                     ub = init_par$upper_bounds,
                                     fitness = log(data$fitness), 
                                     neigh_matrix = neigh_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_DIRECT_L_RAND"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = init_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e4),
                                     lb = init_par$lower_bounds,
                                     ub = init_par$upper_bounds,
                                     fitness = log(data$fitness), 
                                     neigh_matrix = neigh_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "GenSA"){
    tryCatch({
      optim_result <- GenSA::GenSA(par = init_par$init_par,
                                   fn = fitness_model,
                                   lower = init_par$lower_bounds,
                                   upper = init_par$upper_bounds, 
                                   control = list(maxit = 1e3), 
                                   fitness = log(data$fitness), 
                                   neigh_matrix = neigh_matrix,
                                   covariates = covariates, 
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
                                         fitness = log(data$fitness), 
                                         neigh_matrix = neigh_matrix,
                                         covariates = covariates, 
                                         fixed_parameters = fixed_parameters)

    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
    
  }else if(optimization_method == "DEoptimR"){
    tryCatch({
      optim_result <- DEoptimR::JDEoptim(lower = init_par$lower_bounds,
                                         upper = init_par$upper_bounds,
                                         fn = fitness_model,
                                         fitness = log(data$fitness), 
                                         neigh_matrix = neigh_matrix,
                                         covariates = covariates, 
                                         fixed_parameters = fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }
  
  # gather output -----------------------------------------------------------
  
  # this is cumbersome and boring because different methods 
  # have different output types
  
  if(is.null(optim_result)){
    optim_params <- cxr_retrieve_params(optim_par = rep(NA_real_,length(init_par$init_par)),
                                        lambda_length = length(init_lambda),
                                        alpha_length = length(init_alpha),
                                        lambda_cov_length = length(init_lambda_cov),
                                        alpha_cov_length = length(init_alpha_cov))
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
    
    optim_params <- cxr_retrieve_params(optim_par = outpar,
                                        lambda_length = length(init_lambda),
                                        alpha_length = length(init_alpha),
                                        lambda_cov_length = length(init_lambda_cov),
                                        alpha_cov_length = length(init_alpha_cov))
    
  }# if not null
  
  # calculate errors --------------------------------------------------------

  if(bootstrap_samples > 0){

    errors <- cxr_pm_bootstrap(fitness_model = fitness_model,
                               optimization_method = optimization_method,
                               data = data,
                               covariates = covariates,
                               init_par = init_par$init_par,
                               lower_bounds = init_par$lower_bounds,
                               upper_bounds = init_par$upper_bounds,
                               fixed_parameters = fixed_parameters,
                               bootstrap_samples = bootstrap_samples)
    
    error_params <- cxr_retrieve_params(optim_par = errors,
                                        lambda_length = length(init_lambda),
                                        alpha_length = length(init_alpha),
                                        lambda_cov_length = length(init_lambda_cov),
                                        alpha_cov_length = length(init_alpha_cov))
  }else{
    error_params <- list(lambda = NULL, 
                         alpha = NULL,
                         lambda_cov = NULL,
                         alpha_cov = NULL)
  }
  
  # return cxr object ---------------------------------------------------
  
  list_names <- c("model_name",
                  "model",
                  "data",
                  "covariates",
                  "optimization_method",
                  "initial_values",
                  "fixed_terms",
                  "lambda","alpha","lambda_cov","alpha_cov",
                  "lambda_standard_error","alpha_standard_error",
                  "lambda_cov_standard_error","alpha_cov_standard_error",
                  "log_likelihood")
  
  fit <- sapply(list_names,function(x) NULL)
  
  fit$model_name <- model_name
  fit$model <- fitness_model
  fit$data <- data
  fit$covariates <- covariates
  fit$optimization_method <- optimization_method
  fit$initial_values <- initial_values
  
  # for returning explicit NULL values
  if(!is.null(fixed_terms)){
    fit$fixed_terms <- fixed_terms
  }
  if(!is.null(optim_params$lambda)){
    fit$lambda <- optim_params$lambda
  }
  if(!is.null(optim_params$alpha)){
    fit$alpha <- optim_params$alpha
  }
  if(!is.null(optim_params$lambda_cov)){
    fit$lambda_cov <- optim_params$lambda_cov
  }
  if(!is.null(optim_params$alpha_cov)){
    fit$alpha_cov <- optim_params$alpha_cov
  }
  if(!is.null(error_params$lambda)){
    fit$lambda_standard_error <- error_params$lambda
  }
  if(!is.null(error_params$alpha)){
    fit$alpha_standard_error <- error_params$alpha
  }
  if(!is.null(error_params$lambda_cov)){
    fit$lambda_cov_standard_error <- error_params$lambda_cov
  }
  if(!is.null(error_params$alpha_cov)){
    fit$alpha_cov_standard_error <- error_params$alpha_cov
  }
  
  fit$log_likelihood <- llik
  # define two classes, cxr_pm_fit/cxr_er_fit
  
  class(fit) <- "cxr_pm_fit"
  
  if(!is.null(fit$lambda) & !is.null(lower_lambda) & !is.null(upper_lambda)){
    if(fit$lambda == lower_lambda | fit$lambda == upper_lambda){
      message("cxr_pm_fit: A fitted lambda is equal to lower or upper bounds. Consider refitting
              with different boundaries.")
    }
  }
  
  if(!is.null(fit$alpha) & !is.null(lower_alpha) & !is.null(upper_alpha)){
    if(any(fit$alpha == lower_alpha) | any(fit$alpha == upper_alpha)){
      message("cxr_pm_fit: One or more fitted alphas are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$lambda_cov) & !is.null(lower_lambda_cov) & !is.null(upper_lambda_cov)){
    if(any(fit$lambda_cov == lower_lambda_cov) | any(fit$lambda_cov == upper_lambda_cov)){
      message("cxr_pm_fit: A fitted lambda_cov is equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$alpha_cov) & !is.null(lower_alpha_cov) & !is.null(upper_alpha_cov)){
    if(any(fit$alpha_cov == lower_alpha_cov) | any(fit$alpha_cov == upper_alpha_cov)){
      message("cxr_pm_fit: One or more fitted alpha_covs are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }  
  
  fit
}


# summary method ----------------------------------------------------------

summary.cxr_pm_fit <- function(x){
  cat(
    "\nmodel: '",x$model_name,"'",
    "\noptimization method: '",x$optimization_method,"'",
    "\n----------",
    "\nobservations: ",nrow(x$data),
    "\nneighbours: ",length(x$data)-1,
    "\ncovariates: ",ifelse(is.null(x$covariates),0,ncol(x$covariates)),
    "\n----------",
    "\nfocal lambda: ",ifelse(is.null(x$lambda)," - not fit - ",x$lambda),
    "\nmean alpha: ",ifelse(is.null(x$alpha)," - not fit - ",mean(x$alpha)),
    "\nmean lambda_cov: ",ifelse(is.null(x$lambda_cov),"- not fit - ",mean(x$lambda_cov)),
    "\nmean alpha_cov: ",ifelse(is.null(x$alpha_cov),"- not fit - ",mean(x$alpha_cov)),
    "\nlog-likelihood of the fit: ",x$log_likelihood,
    "\n----------",
    sep="")
}

