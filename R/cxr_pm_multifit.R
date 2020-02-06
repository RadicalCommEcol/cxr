# load test data
# library(cxr)
# data("competition")

# TEMP
# source("R/cxr_return_init_length.R")
# source("R/cxr_init_params.R")
# source("R/cxr_retrieve_params.R")
# source("R/pm_BH_alpha_pairwise_lambdacov_none_alphacov_none.R")
# source("R/pm_BH_alpha_pairwise_lambdacov_global_alphacov_global.R")
# source("R/cxr_pm_bootstrap.R")
# source("R/cxr_pm_fit.R")
# source("R/cxr_check_input_data.R")

#' Multi-species parameter optimization
#' 
#' This function is a wrapper for estimating parameters for several
#' focal species, instead of making separate calls to \code{cxr_pm_fit}.
#'
#' @param data named list in which each component is 
#' a dataframe with a fitness column and a number of columns representing neigbhours
#' @inheritParams model_family cxr_pm_fit
#' @param covariates optional named list in which each component is
#' a dataframe with values of each covariate for each observation. The ith component
#' of \code{covariates} are the covariate values that correspond to 
#' the ith component of \code{data}, so they must have the same number of observations.
#' @inheritParams optimization_method cxr_pm_fit
#' @inheritParams alpha_form cxr_pm_fit
#' @inheritParams lambda_cov_form cxr_pm_fit
#' @inheritParams alpha_cov_form cxr_pm_fit
#' @inheritParams initial_values cxr_pm_fit
#' @inheritParams lower_bounds cxr_pm_fit
#' @inheritParams upper_bounds cxr_pm_fit
#' @inheritParams fixed_terms cxr_pm_fit
#' @inheritParams bootstrap_samples cxr_pm_fit
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
#' * alpha: fitted values for alpha, or NULL if fixed
#' * lambda_cov: fitted values for lambda_cov, or NULL if fixed
#' * alpha_cov: fitted values for alpha_cov, or NULL if fixed
#' * lambda_standard_error: standard errors for lambda, if computed
#' * alpha_standard_error: standard errors for alpha, if computed
#' * lambda_cov_standard_error: standard errors for lambda_cov, if computed
#' * alpha_cov_standard_error: standard errors for alpha_cov, if computed
#' * log_likelihood: log-likelihood of the fits
#' @export
#'
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
#' \dontrun{
#'   fit_3sp <- cxr_pm_multifit(data = data,
#'                              optimization_method = "bobyqa",
#'                              covariates = salinity,
#'                              alpha_form = "pairwise",
#'                              lambda_cov_form = "global",
#'                              alpha_cov_form = "global",
#'                              initial_values = list(lambda = 1,alpha = 0.1,lambda_cov = 0.1, alpha_cov = 0.1),
#'                              lower_bounds = list(lambda = 0.01,alpha = 0,lambda_cov = 0, alpha_cov = 0),
#'                              upper_bounds = list(lambda = 100,alpha = 1,lambda_cov = 1, alpha_cov = 1),
#'                              bootstrap_samples = 3)
#'   # brief summary
#'   summary(fit_3sp)
#'   # interaction matrix
#'   fit_3sp$alpha
#' }
cxr_pm_multifit <- function(data, 
                            model_family = c("BH"),
                            covariates = NULL, 
                            optimization_method = c(), 
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
  

# prepare multisp data ----------------------------------------------------
spnames <- names(data)

# fit every sp ------------------------------------------------------------
spfits <- list()
for(i.sp in 1:length(data)){
  spfits[[i.sp]] <- try(cxr_pm_fit(data = data[[i.sp]],
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
  # alpha is a n x n matrix

splambda <- NULL
spalpha <- NULL
splambda_cov <- NULL
spalpha_cov <- NULL

er_splambda <- NULL
er_spalpha <- NULL
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
  # alpha
  if(!is.null(spfits[[i.sp]]$alpha)){
    myalpha <- spfits[[i.sp]]$alpha
    spalpha <- rbind(spalpha,myalpha)
    rownames(spalpha)[i.sp] <- spnames[i.sp]
  }
  # lambda_cov
  if(!is.null(spfits[[i.sp]]$lambda_cov)){
    mylambda_cov <- spfits[[i.sp]]$lambda_cov
    names(mylambda_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$lambda_cov),sep="")
    splambda_cov <- c(splambda_cov,mylambda_cov)
  }
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
  # alpha
  if(!is.null(spfits[[i.sp]]$alpha_standard_error)){
    eralpha <- spfits[[i.sp]]$alpha_standard_error
    er_spalpha <- rbind(er_spalpha,eralpha)
    rownames(er_spalpha)[i.sp] <- spnames[i.sp]
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

summary.cxr_pm_multifit <- function(x){
  cat("model: '",x$model_name,"'",
      "\noptimization method: '",x$optimization_method,"'",
      "\n----------",sep="")
  # for printing null or valid values
  # ifelse returns single values over single conditions
  if(is.null(x$lambda)){
    sl <- rep("- not fit -",nrow(x$data))
  }else{
    sl <- x$lambda
  }
  if(is.null(x$alpha)){
    sa <- rep("- not fit -",nrow(x$data))
  }else{
    sa <- rowMeans(x$alpha)
  }
  if(is.null(x$lambda_cov)){
    slc <- rep("- not fit -",nrow(x$data))
  }else{
    slc <- x$lambda_cov
  }
  if(is.null(x$alpha_cov)){
    sac <- rep("- not fit -",nrow(x$data))
  }else{
    sac <- x$alpha_cov
  }
  summary_table <- data.frame(sp = names(x$data),
                              observations = unlist(lapply(x$data, nrow)),
                              neighbours = unlist(lapply(x$data,length))-1,
                              covariates = ifelse(is.null(x$covariates[[1]]),0,
                                                  ncol(x$covariates[[1]])),
                              # lambda = NULL,
                              lambda = sl,
                              mean_alpha = sa,
                              # mean_lambda_cov = ifelse(is.null(x$lambda_cov),
                              #                          "- not fit - ",
                              #                          mean(x$lambda_cov)),
                              # mean_alpha_cov = ifelse(is.null(x$alpha_cov),
                              #                         "- not fit - ",
                              #                         mean(x$alpha_cov)),
                              row.names = NULL)
  cat("\n")
  summary_table
  # for(i.sp in 1:length(x$data)){
  #   cat("\n",names(x$data)[i.sp],":",
  #       "\nobservations: ",nrow(x$data[[i.sp]]),
  #       "\nneighbours: ",length(x$data[[i.sp]])-1,
  #       "\ncovariates: ",ifelse(is.null(x$covariates[[i.sp]]),0,
  #                               ncol(x$covariates[[i.sp]])),
  #       "\n----------",
  #       "\nfocal lambda: ",ifelse(is.null(x$lambda)," - not fit - ",x$lambda[i.sp]),
  #       "\nmean alpha: ",ifelse(is.null(x$alpha)," - not fit - ",mean(x$alpha[i.sp])),
  #       "\nmean lambda_cov: ",ifelse(is.null(x$lambda_cov),"- not fit - ",mean(x$lambda_cov[i.sp])),
  #       "\nmean alpha_cov: ",ifelse(is.null(x$alpha_cov),"- not fit - ",mean(x$alpha_cov[i.sp])),
  #       "\nlog-likelihood of the fit: ",x$log_likelihood[i.sp],
  #       "\n",
  #       # "\n------------------------",
  #       sep="")
  # }
}
# summary(fit_3sp)

