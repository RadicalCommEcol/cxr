# 
# 
# # load test data
# library(cxr)
# data("competition")
# 
# # spread the data from long to wide format
# competition.data <- tidyr::spread(competition,competitor,number,fill = 0)
# focal.sp <- unique(competition$focal)
# mindata <- subset(competition.data,focal == "CHFU")
# mindata$fitness <- log(mindata$seed)
# mindata <- mindata[,c("fitness",focal.sp)]
# data <- mindata
# 
# initial_values <- list(lambda = 1,alpha = 0,lambda_cov = 0, alpha_cov = 0)
# lower_bounds <- list(lambda = 0,alpha = -1,lambda_cov = 0, alpha_cov = 0)
# upper_bounds <- list(lambda = 10,alpha = 1,lambda_cov = 1, alpha_cov = 1)
# 
# # covariates: rows are observations, columns are different covariates
# # either matrix or dataframe, will be transformed to matrix in the function
# covariates <- data.frame(c1 = rnorm(nrow(mindata),0,1))
# 
# model_family <- "BH"
# optimization_method <- "bobyqa"
# alpha_form <- "pairwise"
# lambda_cov_form <- "none"
# alpha_cov_form <- "none"
# fixed_terms <- NULL
# bootstrap_samples <- 3

# NOTE single species 
# TODO write a wrapper for several sp?

#' General optimization for population models
#' 
#' Estimates parameters of user-specified population dynamics models.
#'
#' @param data dataframe with observations in rows and two sets of columns:
#' * fitness: fitness metric for the focal individual
#' * neighbours: columns with user-defined names with number of neighbours for each group
#' @param model_family family of model to use. Available families are BH (Beverton-Holt) as default.
#' Users may define their own families and models (see vignette XXXXX).
#' @param covariates optional matrix or dataframe with observations (rows) of any number of environmental covariates (columns)
#' @param optimization_method numerical optimization method
#' @param alpha_form what form does the alpha parameter take? one of "none" (no alpha in the model), 
#' "global" (a single alpha for all pairwise interactions), or "pairwise" (one alpha value for every interaction)
#' @param lambda_cov_form form of the covariate effects on lambda. Either "none" (no covariate effects) or "global" (one estimate per covariate)
#' @param alpha_cov_form form of the covariate effects on alpha. One of "none" (no covariate effects), "global" (one estimate per covariate on every alpha),
#' or "pairwise" (one estimate per covariate and pairwise alpha)
#' @param initial_values list with components "lambda","alpha","lambda_cov", "alpha_cov", and "sigma", specifying the initial values
#' for numerical optimization. Single values are allowed.
#' @param lower_bounds optional list with single values for "lambda","alpha","lambda_cov", "alpha_cov".
#' @param upper_bounds optional list with single values for "lambda","alpha","lambda_cov", "alpha_cov".
#' @param fixed_terms optional list specifying which model parameters are fixed, among "lambda","alpha","lambda_cov", and "alpha_cov".
#' @param bootstrap_samples number of bootstrap samples for error calculation. Defaults to 0, i.e. no error is calculated.
#'
#' @return
#' @export
#' @md
#' @examples
cxr_pm_fit <- function(data, 
                    # model, # is it necessary to specify the full model when model_family, alpha_form, etc are given?
                    # --> NO
                    model_family = c("BH"),
                    # add parameter to differentiate population dynamics estimations from effect-response ones 
                    # or code two different functions
                    covariates = NULL, 
                    optimization_method = c(), 
                    alpha_form = c("none","global","pairwise"), 
                    lambda_cov_form = c("none","global"),
                    alpha_cov_form = c("none","global","pairwise"),
                    initial_values = list(lambda = 0, alpha = 0, lambda_cov = 0, alpha_cov = 0, sigma = 0),
                    lower_bounds = NULL,
                    upper_bounds = NULL,
                    fixed_terms = NULL,
                    # errors
                    bootstrap_samples = 0
                    ){
  
  # TODO add cxr:: to the internal functions once they are added
  
  # match arguments ---------------------------------------------------------
  optimization_method <- match.arg(optimization_method)
  alpha_form <- match.arg(alpha_form)
  lambda_cov_form <- match.arg(lambda_cov_form)
  alpha_cov_form <- match.arg(alpha_cov_form)

  # TODO fixed terms?
  
# more checks
  if (optimization_method %in% c("nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND") & !requireNamespace("nloptr", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"nloptr\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "GenSA" & !requireNamespace("GenSA", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"GenSA\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "hydroPSO" & !requireNamespace("hydroPSO", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"hydroPSO\" needed for the method selected to work.",
         call. = FALSE)
  }
  if (optimization_method == "DEoptimR" & !requireNamespace("DEoptimR", quietly = TRUE)) {
    stop("pm_optim ERROR: Package \"DEoptimR\" needed for the method selected to work.",
         call. = FALSE)
  }
  
  # prepare data ------------------------------------------------------------
  # neighbour matrix
  neigh_matrix <- subset(data, select = -c(fitness))
  neigh_matrix <- as.matrix(neigh_matrix)
  
  # how many neighbour species?
  num_neigh <- ncol(neigh_matrix)
  
  # initial parameters
  # check wheter each model parameter is to be fitted or is fixed
  # note how initial_values also function as values 
  # for those parameters that are fixed
  
  fixed_parameters <- list(lambda = NULL, 
                       alpha = NULL,
                       lambda_cov = NULL,
                       alpha_cov = NULL)
  
  init_lambda <- NULL
  init_alpha <- NULL
  init_lambda_cov <- NULL
  init_alpha_cov <- NULL
  
  if("lambda" %in% fixed_terms){
    fixed_parameters[["lambda"]] <- initial_values$lambda
  }else{
    init_lambda <- initial_values$lambda
  }

  if(!is.null(initial_values$sigma)){
    init_sigma <- initial_values$sigma
  }else{
    init_sigma <- 0.1
  }
  
  # return_init_length is an auxiliary function
  # in case initial values are not the same length of the expected parameters
  # e.g. if we want to fit pairwise alphas but only provide a single initial value
  
  if(alpha_form != "none"){
    if("alpha" %in% fixed_terms){
      fixed_parameters[["alpha"]] <- cxr_return_init_length(alpha_form,initial_values$alpha,num_neigh)
    }else{
      init_alpha <- cxr_return_init_length(alpha_form,initial_values$alpha,num_neigh)
    }
  }
  
  if(lambda_cov_form != "none"){
    if("lambda_cov" %in% fixed_terms){
      fixed_parameters[["lambda_cov"]] <- cxr_return_init_length(lambda_cov_form,initial_values$lambda_cov,num_neigh)
    }else{
      init_lambda_cov <- cxr_return_init_length(lambda_cov_form,initial_values$lambda_cov,num_neigh)
    }
  }
    
  if(alpha_cov_form != "none"){
    if("alpha_cov" %in% fixed_terms){
      fixed_parameters[["alpha_cov"]] <- cxr_return_init_length(alpha_cov_form,initial_values$alpha_cov,num_neigh)
    }else{
      init_alpha_cov <- cxr_return_init_length(alpha_cov_form,initial_values$alpha_cov,num_neigh)
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
      lower_sigma <- 1e-10
      upper_sigma <- 1
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
       !"alpha" %in% fixed_terms){
      lower_alpha_cov <- lower_bounds$alpha_cov
      upper_alpha_cov <- upper_bounds$alpha_cov
    }
  
  # retrieve model ----------------------------------------------------------
  fitness_model <- paste("pm_",model_family,"_alpha_",alpha_form,"_lambdacov_",lambda_cov_form,"_alphacov_",alpha_cov_form,sep="")
  
  # TODO this is temporary, while I rename models
  fitness_model <- model_BH3
  
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
                                   upper_alpha_cov = upper_alpha_cov,
                                   num_neigh = num_neigh,
                                   num_covariates = ncol(covariates))

  
  
  # fit parameters ----------------------------------------------------------
  # TODO copy from EBD computer and, in general, rename
  tryCatch({
    optim.result <- optimx::optimx(par = init_par$init_par, 
                                   fn =  fitness_model, 
                                   gr = NULL, 
                                   method = "bobyqa", 
                                   lower = init_par$lower_bounds, 
                                   upper = init_par$upper_bounds,
                                   control = list(parscale = abs(init_par$init_par)), 
                                   hessian = F,
                                   # this will be dropped
                                   param.list = c("lambda","alpha"),#param.list,
                                   # change to "fitness"
                                   log.fitness = data$fitness, 
                                   focal.comp.matrix = neigh_matrix,
                                   num.covariates = ncol(covariates), 
                                   num.competitors = num_neigh, 
                                   focal.covariates = covariates,
                                   # change to fixed_parameters
                                   fixed.terms = fixed_parameters)
  }, error=function(e){cat("pm_optim ERROR :",conditionMessage(e), "\n")})

# gather output -----------------------------------------------------------

  #cxr_retrieve_params

# calculate errors --------------------------------------------------------

  if(bootstrap_samples > 0){
    # TODO check when updated
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
    errors <- rep(NA_real_,length(init_par[[1]]))
  }

  error_params <- cxr_retrieve_params(optim.params = errors,
                                      # param.list = param.list,
                                      alpha.length = length(init.alpha),
                                      alpha.cov.length = length(init.alpha.cov),
                                      num.competitors = num.competitors,
                                      num.covariates = num.covariates)
  
  
# return cxr object ---------------------------------------------------
fit <- list()
  fit$model <- fitness_model
  fit$data <- data
  fit$model_family <- model_family
  fit$covariates <- covariates
  fit$optimization_method <- optimization_method
  fit$initial_values <- initial_values
  fit$fixed_terms <- fixed_terms
  fit$lambda <- NULL
  fit$alpha <- NULL
  fit$lambda_cov <- NULL
  fit$alpha_cov <- NULL
  fit$lambda_lower_bound <- NULL
  fit$lambda_upper_bound <- NULL
  fit$alpha_lower_bound <- NULL
  fit$alpha_upper_bound <- NULL
  fit$lambda_cov_lower_bound <- NULL
  fit$lambda_cov_upper_bound <- NULL
  fit$alpha_cov_lower_bound <- NULL
  fit$alpha_cov_upper_bound <- NULL
  
  # define two classes, cxr_pm_fit/cxr_er_fit
  
  class(fit) <- "cxr_pm_fit"
  fit
}

# seedbank parameters? in principle, not necessary





