#' Multi-species parameter optimization
#' 
#' This function is a wrapper for estimating parameters for several
#' focal species, instead of making separate calls to \code{cxr_pm_fit}.
#'
#' @param data named list in which each component is 
#' a dataframe with a fitness column and a number of columns representing neigbhours
#' @param focal_column character vector with the same length as data,
#' giving the names of the columns representing
#' intraspecific observations for each species, 
#' or numeric vector giving the position of such columns.
#' @inheritParams cxr_pm_fit
#' @param covariates optional named list in which each component is
#' a dataframe with values of each covariate for each observation. The ith component
#' of \code{covariates} are the covariate values that correspond to 
#' the ith component of \code{data}, so they must have the same number of observations.
#'
#' @return an object of type 'cxr_pm_multifit' which is a list with the following components:
#' * model_name: string with the name of the fitness model
#' * model: model function
#' * data: data supplied 
#' * covariates: covariate data supplied
#' * optimization_method: optimization method used
#' * initial_values: list with initial values
#' * fixed_terms: list with fixed terms
#' * lambda: fitted values for lambda, or NULL if fixed
#' * alpha_intra: fitted values for alpha_intra, or NULL if fixed
#' * alpha_inter: fitted values for alpha_inter, or NULL if fixed
#' * lambda_cov: fitted values for lambda_cov, or NULL if fixed
#' * alpha_cov: fitted values for alpha_cov, or NULL if fixed
#' * lambda_standard_error: standard errors for lambda, if computed
#' * alpha_standard_error: standard errors for alpha, if computed
#' * lambda_cov_standard_error: standard errors for lambda_cov, if computed
#' * alpha_cov_standard_error: standard errors for alpha_cov, if computed
#' * log_likelihood: log-likelihood of the fits
#' @export
#' @md
#' @examples
#' # fit three species at once
#' data("neigh_list")
#' data <- neigh_list[1:3]
#' # keep only fitness and neighbours columns
#' for(i in 1:length(data)){
#'   data[[i]] <- data[[i]][,2:length(data[[i]])]
#' }
#' # covariates: salinity
#' data("salinity_list")
#' salinity <- salinity_list[1:3]
#' # keep only salinity column
#' for(i in 1:length(salinity)){
#'   salinity[[i]] <- salinity[[i]][,2:length(salinity[[i]])]
#' }
#' \donttest{
#'   fit_3sp <- cxr_pm_multifit(data = data,
#'                              optimization_method = "bobyqa",
#'                              covariates = salinity,
#'                              alpha_form = "pairwise",
#'                              lambda_cov_form = "global",
#'                              alpha_cov_form = "global",
#'                              initial_values = list(lambda = 1,
#'                                                    alpha = 0.1,
#'                                                    lambda_cov = 0.1, 
#'                                                    alpha_cov = 0.1),
#'                              lower_bounds = list(lambda = 0.01,
#'                                                  alpha = 0,
#'                                                  lambda_cov = 0, 
#'                                                  alpha_cov = 0),
#'                              upper_bounds = list(lambda = 100,
#'                                                  alpha = 1,
#'                                                  lambda_cov = 1, 
#'                                                  alpha_cov = 1),
#'                              bootstrap_samples = 3)
#'   # brief summary
#'   summary(fit_3sp)
#'   # interaction matrix
#'   fit_3sp$alpha
#' }
cxr_pm_multifit <- function(data, 
                            model_family = c("BH"),
                            focal_column = NULL,
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
                            initial_values = NULL,
                            lower_bounds = NULL,
                            upper_bounds = NULL,
                            fixed_terms = NULL,
                            # errors
                            bootstrap_samples = 0
){
  
  # argument checks ---------------------------------------------------------
  
  optimization_method <- match.arg(optimization_method)
  alpha_form <- match.arg(alpha_form)
  lambda_cov_form <- match.arg(lambda_cov_form)
  alpha_cov_form <- match.arg(alpha_cov_form)
  
  # TODO fixed terms?
  
  # check input data
  if(class(data) != "list"){
    stop("cxr_pm_multifit ERROR: check the consistency of your input data: 
    1) data is a named list containing dataframes with observations for each focal species;
    2) No NAs; 
    3) first column in 'data' is named 'fitness'; 
    4) abundances of at least one neighbour species in 'data';
    5) data and covariates (if present) have the same number of observations")
  }else{
    for(i.sp in 1:length(data)){
      temp <- cxr_check_pm_input(data = data[[i.sp]],
                                 focal_column = focal_column[i.sp],
                                 model_family = model_family,
                                 covariates = covariates[[i.sp]],
                                 optimization_method = optimization_method,
                                 alpha_form = alpha_form,
                                 lambda_cov_form = lambda_cov_form,
                                 alpha_cov_form = alpha_cov_form,
                                 initial_values = initial_values,
                                 lower_bounds = lower_bounds,
                                 upper_bounds = upper_bounds,
                                 fixed_terms = fixed_terms,
                                 bootstrap_samples = bootstrap_samples)
      if(temp[[1]] == "error"){
        message(paste("cxr_pm_multifit ERROR: check data format for sp ",names(data)[i.sp]," and/or input parameters.\n",sep=""))
        message(paste("more info on the error:\n",temp[[2]],sep=""))
        return(NULL)
      }
    }
}
  if(!data.ok){
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
  
  
# prepare multisp data ----------------------------------------------------
spnames <- names(data)

# fit every sp ------------------------------------------------------------
spfits <- list()
for(i.sp in 1:length(data)){
  
  spfits[[i.sp]] <- try(cxr_pm_fit(data = data[[i.sp]],
                                   focal_column = focal_column[i.sp],
                                   model_family = model_family,
                                   covariates = covariates[[i.sp]],
                                   optimization_method = optimization_method,
                                   alpha_form = alpha_form,
                                   lambda_cov_form = lambda_cov_form,
                                   alpha_cov_form = alpha_cov_form,
                                   initial_values = initial_values,
                                   lower_bounds = lower_bounds,
                                   upper_bounds = upper_bounds,
                                   fixed_terms = fixed_terms,
                                   bootstrap_samples = bootstrap_samples
  ))
}

# output ------------------------------------------------------------------

  # return a cxr_pm_multifit object
  # which is basically the same as the base object
  # but with info on more sp. e.g. lambda is a 1d vector,
  # alpha_inter is a n x n matrix

splambda <- NULL
spalpha_intra <- NULL
spalpha_inter <- NULL
splambda_cov <- NULL
spalpha_cov <- NULL

er_splambda <- NULL
er_spalpha_intra <- NULL
er_spalpha_inter <- NULL
er_splambda_cov <- NULL
er_spalpha_cov <- NULL

spllik <- NULL

for(i.sp in 1:length(spnames)){
  # lambda
  if(!is.null(spfits[[i.sp]]$lambda)){
    mylambda <- spfits[[i.sp]]$lambda
    names(mylambda) <- spnames[i.sp]
    splambda <- c(splambda,mylambda)
  }
  # alpha_intra
  if(!is.null(spfits[[i.sp]]$alpha_intra)){
    myalpha_intra <- spfits[[i.sp]]$alpha_intra
    names(myalpha_intra) <- spnames[i.sp]
    spalpha_intra <- c(spalpha_intra,myalpha_intra)
  }
  # alpha
  if(!is.null(spfits[[i.sp]]$alpha_inter)){
    myalpha_inter <- spfits[[i.sp]]$alpha_inter
    spalpha_inter <- rbind(spalpha_inter,myalpha_inter)
    rownames(spalpha_inter)[i.sp] <- spnames[i.sp]
  }
  # lambda_cov
  if(!is.null(spfits[[i.sp]]$lambda_cov)){
    mylambda_cov <- spfits[[i.sp]]$lambda_cov
    names(mylambda_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$lambda_cov),sep="")
    splambda_cov <- c(splambda_cov,mylambda_cov)
  }
  # alpha_cov
  if(!is.null(spfits[[i.sp]]$alpha_cov)){
    myalpha_cov <- spfits[[i.sp]]$alpha_cov
    names(myalpha_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$alpha_cov),sep="")
    spalpha_cov <- c(spalpha_cov,myalpha_cov)
  }
  
  # errors
  
  # lambda
  if(!is.null(spfits[[i.sp]]$lambda_standard_error)){
    erlambda <- spfits[[i.sp]]$lambda_standard_error
    if(!is.null(names(erlambda))){
      names(erlambda) <- paste(spnames[i.sp],"_",names(erlambda),sep="")
    }else{
      names(erlambda) <- paste(spnames[i.sp],"_lambda_se",sep="")
    }
    er_splambda <- c(er_splambda,erlambda)
  }
  # alpha_intra
  if(!is.null(spfits[[i.sp]]$alpha_intra_standard_error)){
    eralpha_intra <- spfits[[i.sp]]$alpha_intra_standard_error
    er_spalpha_intra <- c(er_spalpha_intra,eralpha_intra)
  }
  # alpha_inter
  if(!is.null(spfits[[i.sp]]$alpha_inter_standard_error)){
    eralpha_inter <- spfits[[i.sp]]$alpha_inter_standard_error
    er_spalpha_inter <- rbind(er_spalpha_inter,eralpha_inter)
    rownames(er_spalpha_inter)[i.sp] <- spnames[i.sp]
  }
  # lambda_cov
  if(!is.null(spfits[[i.sp]]$lambda_cov_standard_error)){
    erlambda_cov <- spfits[[i.sp]]$lambda_cov_standard_error
    names(erlambda_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$lambda_cov_standard_error),sep="")
    er_splambda_cov <- c(er_splambda_cov,erlambda_cov)
  }
  if(!is.null(spfits[[i.sp]]$alpha_cov_standard_error)){
    eralpha_cov <- spfits[[i.sp]]$alpha_cov_standard_error
    names(eralpha_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$alpha_cov_standard_error),sep="")
    er_spalpha_cov <- c(er_spalpha_cov,eralpha_cov)
  }
  
  # log-likelihood
  myllik <- spfits[[i.sp]]$log_likelihood
  names(myllik) <- spnames[i.sp]
  spllik <- c(spllik,myllik)
}

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

fit$model_name <- spfits[[1]]$model_name
fit$model <- spfits[[1]]$fitness_model
fit$data <- data
fit$covariates <- covariates
fit$optimization_method <- optimization_method
fit$initial_values <- initial_values

# for returning explicit NULL values
if(!is.null(fixed_terms)){
  fit$fixed_terms <- fixed_terms
}
if(!is.null(splambda)){
  fit$lambda <- splambda
}
if(!is.null(spalpha)){
  fit$alpha <- spalpha
}
if(!is.null(splambda_cov)){
  fit$lambda_cov <- splambda_cov
}
if(!is.null(spalpha_cov)){
  fit$alpha_cov <- spalpha_cov
}
if(!is.null(er_splambda)){
  fit$lambda_standard_error <- er_splambda
}
if(!is.null(er_spalpha)){
  fit$alpha_standard_error <- er_spalpha
}
if(!is.null(er_splambda_cov)){
  fit$lambda_cov_standard_error <- er_splambda_cov
}
if(!is.null(er_spalpha_cov)){
  fit$alpha_cov_standard_error <- er_spalpha_cov
}

fit$log_likelihood <- spllik

class(fit) <- "cxr_pm_multifit"
fit

}