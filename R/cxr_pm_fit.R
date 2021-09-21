#' General optimization for population models
#' 
#' Estimates parameters of user-specified population dynamics models.
#'
#' @param data dataframe with observations in rows and two sets of columns:
#' * fitness: fitness metric for the focal individual
#' * neighbours: numeric columns with user-defined names, giving number of neighbours for each group
#' @param focal_column optional integer or character giving the column
#' with neighbours from the same species as the focal one. This field is necessary if "alpha_intra" is specified
#' in \code{initial_values}, \code{lower_bounds}, \code{upper_bounds}, or \code{fixed_terms}.
#' @param model_family family of model to use. Available families are BH (Beverton-Holt), LV (Lotka-Volterra),
#' RK (Ricker), and LW (Law-Watkinson). Users may also define their own families and models (see vignette 4).
#' @param covariates optional named matrix or dataframe with observations (rows) of any number of environmental covariates (columns).
#' @param optimization_method numerical optimization method.
#' @param alpha_form what form does the alpha parameter take? one of "none" (no alpha in the model), 
#' "global" (a single alpha for all pairwise interactions), or "pairwise" (one alpha value for every interaction).
#' @param lambda_cov_form form of the covariate effects on lambda. Either "none" (no covariate effects) or "global" (one estimate per covariate).
#' @param alpha_cov_form form of the covariate effects on alpha. One of "none" (no covariate effects), "global" (one estimate per covariate on every alpha),
#' or "pairwise" (one estimate per covariate and pairwise alpha)
#' @param initial_values list with components "lambda","alpha_intra","alpha_inter","lambda_cov", "alpha_cov", specifying the initial values
#' for numerical optimization. Single values are allowed.
#' @param lower_bounds optional list with single values for "lambda","alpha_intra","alpha_inter","lambda_cov", "alpha_cov".
#' @param upper_bounds optional list with single values for "lambda","alpha_intra","alpha_inter","lambda_cov", "alpha_cov".
#' @param fixed_terms optional list of numeric vectors specifying the value of fixed model parameters, among 
#' "lambda","alpha_intra","alpha_inter","lambda_cov", and "alpha_cov".
#' @param bootstrap_samples number of bootstrap samples for error calculation. Defaults to 0, i.e. no error is calculated.
#' @return an object of type 'cxr_pm_fit' which is a list with the following components:
#' * model_name: string with the name of the fitness model
#' * model: model function
#' * data: data supplied 
#' * focal_ID: name/ID of the focal taxa, if provided in 'focal_column'
#' * covariates: covariate data supplied
#' * optimization_method: optimization method used
#' * initial_values: list with initial values
#' * fixed_terms: list with fixed terms
#' * lambda: fitted value for lambda, or NULL if fixed
#' * alpha_intra: fitted value for intraspecific alpha, or NULL if fixed
#' * alpha_inter: fitted value for interspecific alpha, or NULL if fixed
#' * lambda_cov: fitted value(s) for lambda_cov, or NULL if fixed.
#' * alpha_cov: fitted value(s) for alpha_cov, or NULL if fixed. 
#' These are structured as a list with one element for each covariate.
#' * lambda_standard_error: standard error for lambda, if computed
#' * alpha_intra_standard_error: standard error for intraspecific alpha, if computed
#' * alpha_inter_standard_error: standard error for interspecific alpha, if computed
#' * lambda_cov_standard_error: standard error for lambda_cov, if computed
#' * alpha_cov_standard_error: standard error for alpha_cov, if computed
#' * log_likelihood: log-likelihood of the fit
#' @export
#' @md
#' @examples
#' data("neigh_list")
#' my.sp <- "BEMA"
#' # data for a single species, keep only fitness and neighbours columns
#' sp_data <- neigh_list[[my.sp]][2:ncol(neigh_list[[1]])]
#' \donttest{
#'   sp_fit <- cxr_pm_fit(data = sp_data,
#'                        focal_column = my.sp,
#'                        optimization_method = "bobyqa",
#'                        model_family = "BH",
#'                        alpha_form = "pairwise",
#'                        lambda_cov_form = "none",
#'                        alpha_cov_form = "none",
#'                        initial_values = list(lambda = 1,alpha_intra = 0.1,alpha_inter = 0.1),
#'                        lower_bounds = list(lambda = 0,alpha_intra = 0,alpha_inter = 0),
#'                        upper_bounds = list(lambda = 100,alpha_intra = 1,alpha_inter = 1),
#'                        bootstrap_samples = 3)
#'   summary(sp_fit)
#' }
#' 
cxr_pm_fit <- function(data, 
                       focal_column = NULL,
                       model_family,
                       covariates = NULL, 
                       optimization_method = c("Nelder-Mead","BFGS", "CG", 
                                               "ucminf","L-BFGS-B", "nlm", "nlminb", 
                                               "Rcgmin", "Rvmmin", "spg", 
                                               "bobyqa", "nmkb", "hjkb",
                                               "nloptr_CRS2_LM","nloptr_ISRES",
                                               "nloptr_DIRECT_L_RAND","DEoptimR",
                                               "hydroPSO","GenSA"), 
                       alpha_form = c("none","global","pairwise"), 
                       lambda_cov_form = c("none","global"),
                       alpha_cov_form = c("none","global","pairwise"),
                       initial_values = list(lambda = 0, 
                                             alpha_intra = 0,
                                             alpha_inter = 0, 
                                             lambda_cov = 0, 
                                             alpha_cov = 0),
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
  
  # wrapper for several consistency checks on the input arguments
  input.ok <- cxr_check_pm_input(data,
                                 focal_column,
                                 model_family,
                                 covariates,
                                 optimization_method,
                                 alpha_form,
                                 lambda_cov_form,
                                 alpha_cov_form,
                                 initial_values,
                                 lower_bounds,
                                 upper_bounds,
                                 fixed_terms)
  if(input.ok[[1]] == "error"){
    message(input.ok[[2]])
    return(NULL)
    # stop(input.ok[[2]],call. = FALSE)
  }else if(input.ok[[1]] == "warning"){
    message(input.ok[[2]])
  }
  
  # retrieve model ----------------------------------------------------------
  # character string giving the name of the model
  model_name <- paste(model_family,"_pm",
                      "_alpha_",alpha_form,
                      "_lambdacov_",lambda_cov_form,
                      "_alphacov_",alpha_cov_form,sep="")
  
  # try to retrieve the function from its name
  # using function "get"
  fitness_model <- try(get(model_name),silent = TRUE)
  if(inherits(fitness_model,"try-error")){
    message(paste("cxr_pm_fit ERROR: model '",model_name,"' could not be retrieved. 
  Make sure it is defined and available 
               in the cxr package or in the global environment.\n"
               ,sep=""))
    return(NULL)
  }
  
  # prepare data ------------------------------------------------------------
  
  # neighbour matrix
  
  # just to avoid a note in R CMD CHECK
  dropname <- "fitness"
  neigh_matrix <- as.matrix(data[ , !(names(data) %in% dropname)])
  
  # initial check to remove any neighbour with no presences
  # also ensure there is at least one neighbour with presences
  if(any(colSums(neigh_matrix) == 0)){
    empty_neigh <- colnames(neigh_matrix)[which(colSums(neigh_matrix) == 0)]
    present_neigh <- colnames(neigh_matrix)[which(colSums(neigh_matrix) != 0)]
    
    if(length(present_neigh) == 0){
      message("cxr_pm_fit ERROR: No neighbours with densities > 0 in any observation.")      
      return(NULL)
    }else{
      neigh_matrix <- neigh_matrix[,present_neigh]
    }
    
    # neigh_matrix <- neigh_matrix[,present_neigh]
  }else{
    empty_neigh <- NULL
    present_neigh <- colnames(neigh_matrix)
  }# if-else any neighbour with no presences

  # differentiate alpha_intra/inter
  
  # if alpha_intra is present, intraspecific neighbours
  # should be in a different 1-column matrix
  # there are other options, but this 
  # clearly separates both types of observations
  # and avoids mixing columns in neigh_matrix, etc
  
  if(is.null(focal_column)){
    # no alpha_intra
    neigh_inter_matrix <- neigh_matrix
    neigh_intra_matrix <- NULL
    # set also names
    neigh_inter <- colnames(neigh_inter_matrix)
    neigh_intra <- NULL
  }else{
    
    # which column number
    if(inherits(focal_column,"character")){
      focal_column_num <- which(colnames(neigh_matrix) == focal_column)
    }else{
      focal_column_num <- focal_column -1 # because data has fitness in column 1
    }
    
    if(length(focal_column_num) == 0){
      # no alpha_intra
      neigh_inter_matrix <- neigh_matrix
      neigh_intra_matrix <- NULL
      # set also names
      neigh_inter <- colnames(neigh_inter_matrix)
      neigh_intra <- NULL
      
      message("cxr_pm_fit: a focal column is provided, but it 
            contains no densities > 0. It will be discarded.")
    }else{
      
      # intra and inter observations in different matrices
      neigh_inter_matrix <- as.matrix(neigh_matrix[,-focal_column_num])
      neigh_intra_matrix <- as.matrix(neigh_matrix[,focal_column_num])
      colnames(neigh_intra_matrix) <- colnames(neigh_matrix)[focal_column_num]
      if(ncol(neigh_inter_matrix) == 1){
        colnames(neigh_inter_matrix) <- colnames(neigh_matrix)[-focal_column_num]
      }
      # set also names
      neigh_inter <- colnames(neigh_inter_matrix)
      neigh_intra <- colnames(neigh_matrix)[focal_column_num]
    }
  }
  

  
  # are covariates named?
  if(!is.null(covariates)){
    if(is.null(colnames(covariates))){
      colnames(covariates) <- paste("cov",1:ncol(covariates),sep="")
    }
  }

  # clean initial parameters ------------------------------------------------
  
  # prepare initial parameters, 
  # checking which ones are fixed or not,
  # their length, and names.
  # these functions do not check for errors
  # since that is already done with 'cxr_check_pm_input'
  init_par <- cxr_get_init_params(initial_values,
                                  fixed_terms,
                                  alpha_form,
                                  lambda_cov_form,
                                  alpha_cov_form,
                                  model_type = "pm",
                                  neigh_intra,
                                  neigh_inter,
                                  covariates)
  
  
  # retrieve lower and upper bounds, if present
  bounds <- cxr_get_model_bounds(lower_bounds,upper_bounds,fixed_terms)
  
  # sort parameters for optim routine ---------------------------------------
  
  vector_par <- cxr_sort_params(init_lambda = init_par$init_lambda,
                              init_sigma = init_par$init_sigma,
                              init_alpha_intra = init_par$init_alpha_intra,
                              init_alpha_inter = init_par$init_alpha_inter,
                              init_lambda_cov = init_par$init_lambda_cov,
                              init_alpha_cov = init_par$init_alpha_cov,
                              lower_lambda = bounds$lower_lambda,
                              upper_lambda = bounds$upper_lambda,
                              lower_sigma = bounds$lower_sigma,
                              upper_sigma = bounds$upper_sigma,
                              lower_alpha_intra = bounds$lower_alpha_intra,
                              upper_alpha_intra = bounds$upper_alpha_intra,
                              lower_alpha_inter = bounds$lower_alpha_inter,
                              upper_alpha_inter = bounds$upper_alpha_inter,
                              lower_lambda_cov = bounds$lower_lambda_cov,
                              upper_lambda_cov = bounds$upper_lambda_cov,
                              lower_alpha_cov = bounds$lower_alpha_cov,
                              upper_alpha_cov = bounds$upper_alpha_cov)
  
  # fit parameters ----------------------------------------------------------

  optim_result <- NULL
  error.message <- NULL

  if(optimization_method %in% c("BFGS", "CG", "Nelder-Mead", "ucminf")){
    tryCatch({
      optim_result <- optimx::optimx(par = vector_par$init_par, 
                                     fn = fitness_model, 
                                     gr = NULL, 
                                     method = optimization_method,
                                     # lower = init_par$lower_bounds,
                                     # upper = init_par$upper_bounds,
                                     control = list(), 
                                     hessian = F,
                                     fitness = log(data$fitness), 
                                     neigh_intra_matrix = neigh_intra_matrix,
                                     neigh_inter_matrix = neigh_inter_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = init_par$fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method %in% c("L-BFGS-B", "nlm", "nlminb", 
                               "Rcgmin", "Rvmmin", "spg", 
                               "bobyqa", "nmkb", "hjkb")){
    tryCatch({
      optim_result <- optimx::optimx(par = vector_par$init_par, 
                                     fn = fitness_model, 
                                     gr = NULL, 
                                     method = optimization_method,
                                     lower = vector_par$lower_bounds,
                                     upper = vector_par$upper_bounds,
                                     control = list(), 
                                     hessian = F,
                                     fitness = log(data$fitness), 
                                     neigh_intra_matrix = neigh_intra_matrix,
                                     neigh_inter_matrix = neigh_inter_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = init_par$fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_CRS2_LM"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = vector_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e4),
                                     lb = vector_par$lower_bounds,
                                     ub = vector_par$upper_bounds,
                                     fitness = log(data$fitness), 
                                     neigh_intra_matrix = neigh_intra_matrix,
                                     neigh_inter_matrix = neigh_inter_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = init_par$fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_ISRES"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = vector_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e4),
                                     lb = vector_par$lower_bounds,
                                     ub = vector_par$upper_bounds,
                                     fitness = log(data$fitness), 
                                     neigh_intra_matrix = neigh_intra_matrix,
                                     neigh_inter_matrix = neigh_inter_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = init_par$fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "nloptr_DIRECT_L_RAND"){
    tryCatch({
      optim_result <- nloptr::nloptr(x0 = vector_par$init_par,
                                     eval_f = fitness_model,
                                     opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e4),
                                     lb = vector_par$lower_bounds,
                                     ub = vector_par$upper_bounds,
                                     fitness = log(data$fitness), 
                                     neigh_intra_matrix = neigh_intra_matrix,
                                     neigh_inter_matrix = neigh_inter_matrix,
                                     covariates = covariates, 
                                     fixed_parameters = init_par$fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "GenSA"){
    tryCatch({
      optim_result <- GenSA::GenSA(par = vector_par$init_par,
                                   fn = fitness_model,
                                   lower = vector_par$lower_bounds,
                                   upper = vector_par$upper_bounds, 
                                   control = list(maxit = 1e3), 
                                   fitness = log(data$fitness), 
                                   neigh_intra_matrix = neigh_intra_matrix,
                                   neigh_inter_matrix = neigh_inter_matrix,
                                   covariates = covariates, 
                                   fixed_parameters = init_par$fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }else if(optimization_method == "hydroPSO"){
    tryCatch({
      # suppress annoying output??
      # sink("/dev/null")
      optim_result <- hydroPSO::hydroPSO(par = vector_par$init_par,
                                         fn = fitness_model,
                                         lower = vector_par$lower_bounds,
                                         upper = vector_par$upper_bounds, 
                                         control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
                                         fitness = log(data$fitness), 
                                         neigh_intra_matrix = neigh_intra_matrix,
                                         neigh_inter_matrix = neigh_inter_matrix,
                                         covariates = covariates, 
                                         fixed_parameters = init_par$fixed_parameters)

    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
    
  }else if(optimization_method == "DEoptimR"){
    tryCatch({
      optim_result <- DEoptimR::JDEoptim(lower = vector_par$lower_bounds,
                                         upper = vector_par$upper_bounds,
                                         fn = fitness_model,
                                         fitness = log(data$fitness), 
                                         neigh_intra_matrix = neigh_intra_matrix,
                                         neigh_inter_matrix = neigh_inter_matrix,
                                         covariates = covariates, 
                                         fixed_parameters = init_par$fixed_parameters)
    }, error=function(e){cat("cxr_pm_fit optimization ERROR :",conditionMessage(e), "\n")})
  }
  
  # gather output -----------------------------------------------------------
  
  # this is cumbersome and boring because different methods 
  # have different output types
  
  if(is.null(optim_result)){
    # error has been printed by trycatch method above,
    # so simply stop
    return(NULL)
    
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
      llik <- optim_result$value
    }else{
      outpar <- optim_result$solution
      names(outpar) <- names(init_par$init_par)
      llik <- optim_result$objective
    }# if-else method
    
    if(any(is.na(outpar)) | is.na(llik)){
      message("cxr_pm_fit ERROR: one or more fitted values are NA. This is likely to be because the optimization algorithm failed to find a solution.
              Please check the lower and upper bounds provided and/or the optimization method selected.")
      return(NULL)
    }
    
    optim_params <- cxr_retrieve_params(optim_par = outpar,
                                        lambda_length = length(init_par$init_lambda),
                                        alpha_intra_length = length(init_par$init_alpha_intra),
                                        alpha_inter_length = length(init_par$init_alpha_inter),
                                        lambda_cov_length = length(init_par$init_lambda_cov),
                                        alpha_cov_length = length(init_par$init_alpha_cov),
                                        empty_neigh = empty_neigh,
                                        alpha_cov_form = alpha_cov_form)
    
  }# if not null
  
  # calculate errors --------------------------------------------------------

  # double check 
  # 1) whether errors are calculated
  # 2) whether errors are *correctly* calculated
  if(bootstrap_samples > 0){
    
    errors <- cxr_pm_bootstrap(fitness_model = fitness_model,
                               optimization_method = optimization_method,
                               data = data,
                               focal_column = focal_column,
                               covariates = covariates,
                               init_par = vector_par$init_par,
                               lower_bounds = vector_par$lower_bounds,
                               upper_bounds = vector_par$upper_bounds,
                               fixed_parameters = init_par$fixed_parameters,
                               bootstrap_samples = bootstrap_samples)
    if(!is.null(errors)){
      error_params <- cxr_retrieve_params(optim_par = errors,
                                          lambda_length = length(init_par$init_lambda),
                                          alpha_intra_length = length(init_par$init_alpha_intra),
                                          alpha_inter_length = length(init_par$init_alpha_inter),
                                          lambda_cov_length = length(init_par$init_lambda_cov),
                                          alpha_cov_length = length(init_par$init_alpha_cov),
                                          empty_neigh = empty_neigh,
                                          alpha_cov_form = alpha_cov_form,error_par = TRUE)
    }else{
      error_params <- list(lambda = NULL, 
                           alpha_intra = NULL,
                           alpha_inter = NULL,
                           lambda_cov = NULL,
                           alpha_cov = NULL)
    }
  }else{
    error_params <- list(lambda = NULL, 
                         alpha_intra = NULL,
                         alpha_inter = NULL,
                         lambda_cov = NULL,
                         alpha_cov = NULL)
  }
  
  # prepare and return cxr_pm_fit object --------------------------------

  list_names <- c("model_name",
                  "model",
                  "data",
                  "focal_ID",
                  "covariates",
                  "optimization_method",
                  "initial_values",
                  "fixed_terms",
                  "lambda","alpha_intra","alpha_inter","lambda_cov","alpha_cov",
                  "lambda_standard_error","alpha_intra_standard_error",
                  "alpha_inter_standard_error",
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
  if(!is.null(focal_column)){
    fit$focal_ID <- neigh_intra
  }
  if(!is.null(fixed_terms)){
    fit$fixed_terms <- fixed_terms
  }
  if(!is.null(optim_params$lambda)){
    fit$lambda <- optim_params$lambda
  }
  if(!is.null(optim_params$alpha_intra)){
    fit$alpha_intra <- optim_params$alpha_intra
  }
  if(!is.null(optim_params$alpha_inter)){
    fit$alpha_inter <- optim_params$alpha_inter
  }
  if(!is.null(optim_params$lambda_cov)){
    fit$lambda_cov <- optim_params$lambda_cov
  }
  
  if(!is.null(optim_params$alpha_cov)){
    tidy_ac <- list()
    if(length(optim_params$alpha_cov) == ncol(covariates)){
      acnames <- substr(names(optim_params$alpha_cov),11,(nchar(names(optim_params$alpha_cov))))
    }else{
      # extract covariate names from names(alpha_cov)
      # trickier than i thought
      # 1 - remove species names
      neighnames <- colnames(data)[-1]
      acnames <- names(optim_params$alpha_cov)
      
      # gsub is not vectorized and alternatives are difficult to understand/other packages
      # so go the loop way
      for(i.acn in 1:length(optim_params$alpha_cov)){
        acnames[i.acn] <- gsub("alpha_cov_","",acnames[i.acn])
        for(i.n in 1:length(neighnames)){
          acnames[i.acn] <- gsub(paste("_",neighnames[i.n],sep=""),"",acnames[i.acn])
        }
      }
      
      # anames <- substr(names(optim_params$alpha_cov),11,(nchar(names(optim_params$alpha_cov))-5))
    }
    for(i.cov in 1:ncol(covariates)){
      # grepl does not give exact match
      # my.cov <- which(grepl(colnames(covariates)[i.cov],names(optim_params$alpha_cov)))
      
      my.cov <- which(acnames == colnames(covariates)[i.cov])
      tidy_ac[[i.cov]] <- optim_params$alpha_cov[my.cov]
    }
    names(tidy_ac) <- colnames(covariates)
    fit$alpha_cov <- tidy_ac
  }
  
  if(!is.null(error_params$lambda)){
    fit$lambda_standard_error <- error_params$lambda
  }
  if(!is.null(error_params$alpha_intra)){
    fit$alpha_intra_standard_error <- error_params$alpha_intra
  }
  if(!is.null(error_params$alpha_inter)){
    fit$alpha_inter_standard_error <- error_params$alpha_inter
  }
  if(!is.null(error_params$lambda_cov)){
    fit$lambda_cov_standard_error <- error_params$lambda_cov
  }
  
  if(!is.null(error_params$alpha_cov)){
    tidy_acr <- list()
    
    if(length(error_params$alpha_cov) == ncol(covariates)){
      ernames <- substr(names(error_params$alpha_cov),11,(nchar(names(error_params$alpha_cov))-3))
    }else{
      
      # extract covariate names from names(alpha_cov)
      # trickier than i thought
      # 1 - remove species names
      neighnames <- colnames(data)[-1]
      ernames <- names(error_params$alpha_cov)
      
      # gsub is not vectorized and alternatives are difficult to understand/other packages
      # so go the loop way
      for(i.acn in 1:length(error_params$alpha_cov)){
        ernames[i.acn] <- gsub("alpha_cov_","",ernames[i.acn])
        for(i.n in 1:length(neighnames)){
          ernames[i.acn] <- gsub(paste("_",neighnames[i.n],"_se",sep=""),"",ernames[i.acn])
        }
      }
      
    }
    
    for(i.cov in 1:ncol(covariates)){
      # my.cov <- which(grepl(colnames(covariates)[i.cov],names(error_params$alpha_cov)))
      my.cov <- which(ernames == colnames(covariates)[i.cov])
      tidy_acr[[i.cov]] <- error_params$alpha_cov[my.cov]
    }
    names(tidy_acr) <- colnames(covariates)
    fit$alpha_cov_standard_error <- tidy_acr
  }
  
  fit$log_likelihood <- llik
  # define two classes, cxr_pm_fit/cxr_er_fit
  
  class(fit) <- "cxr_pm_fit"
  
  if(!is.null(fit$lambda) & !is.null(bounds$lower_lambda) & !is.null(bounds$upper_lambda)){
    if(fit$lambda == bounds$lower_lambda | fit$lambda == bounds$upper_lambda){
      message("cxr_pm_fit: A fitted lambda is equal to lower or upper bounds. Consider refitting
              with different boundaries.")
    }
  }
  
  if(!is.null(fit$alpha_intra) & !is.null(bounds$lower_alpha_intra) & !is.null(bounds$upper_alpha_intra)){
    if(fit$alpha_intra == bounds$lower_alpha_intra | fit$alpha_intra == bounds$upper_alpha_intra){
      message("cxr_pm_fit: Fitted Intraspecific alpha is equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$alpha_inter) & !is.null(bounds$lower_alpha_inter) & !is.null(bounds$upper_alpha_inter)){
    if(any(fit$alpha_inter == bounds$lower_alpha_inter,na.rm = TRUE) | 
       any(fit$alpha_inter == bounds$upper_alpha_inter,na.rm = TRUE)){
      message("cxr_pm_fit: One or more fitted interspecific alphas are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$lambda_cov) & !is.null(bounds$lower_lambda_cov) & !is.null(bounds$upper_lambda_cov)){
    if(any(fit$lambda_cov == bounds$lower_lambda_cov,na.rm = TRUE) | 
       any(fit$lambda_cov == bounds$pper_lambda_cov,na.rm = TRUE)){
      message("cxr_pm_fit: A fitted lambda_cov is equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }
  
  if(!is.null(fit$alpha_cov) & !is.null(bounds$lower_alpha_cov) & !is.null(bounds$upper_alpha_cov)){
    if(any(unlist(fit$alpha_cov) == bounds$lower_alpha_cov,na.rm = TRUE) | 
       any(unlist(fit$alpha_cov) == bounds$upper_alpha_cov,na.rm = TRUE)){
      message("cxr_pm_fit: One or more fitted alpha_covs are equal to lower or upper bounds. 
      Consider refitting with different boundaries.")
    }
  }  
  
  fit
}