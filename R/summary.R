#' CXR summary method for population model fits
#'
#' @param object a \code{cxr_pm_fit} object, from the function with the same name
#' @param ... other arguments, not used
#' @return console output
#' @export
summary.cxr_pm_fit <- function(object,...){
  cat(
    "\nmodel: '",object$model_name,"'",
    "\noptimization method: '",object$optimization_method,"'",
    "\n----------",
    "\nobservations: ",nrow(object$data),
    "\nneighbours: ",length(object$data)-1,
    "\ncovariates: ",ifelse(is.null(object$covariates),0,ncol(object$covariates)),
    "\n----------",
    "\nfocal lambda: ",ifelse(is.null(object$lambda)," - not fit - ",object$lambda),
    "\nmean alpha: ",ifelse(is.null(object$alpha)," - not fit - ",mean(object$alpha)),
    "\nmean lambda_cov: ",ifelse(is.null(object$lambda_cov),"- not fit - ",mean(object$lambda_cov)),
    "\nmean alpha_cov: ",ifelse(is.null(object$alpha_cov),"- not fit - ",mean(object$alpha_cov)),
    "\nlog-likelihood of the fit: ",object$log_likelihood,
    "\n----------",
    sep="")
}

#' CXR summary method for effect response model fits
#'
#' @param object a \code{cxr_er_fit} object, from the function with the same name
#' @param ... other arguments, not used
#' @return console output
#' @export
summary.cxr_er_fit <- function(object,...){
  covar <- 0
  if("list" %in% class(object$data)){
    obs <- nrow(object$data[[1]])
    if(!is.null(object$covariates)){
      covar <- length(object$covariates)
    }
  }else{
    obs <- unique(table(object$data$focal))
    if(!is.null(object$covariates)){
      covar <- ncol(object$covariates)
    }
  }

  cat("model:",object$model_name,"",
      "\noptimization method:",object$optimization_method,"",
      "\nspecies:", object$sp,
      "\ncovariates:", covar,
      "\nobservations:", obs,
      "\n----------",sep=" ")

  # for printing null or valid values
  # ifelse returns single values over single conditions
  if(is.null(object$lambda)){
    sl <- rep("-not fit-",length(object$sp))
  }else{
    sl <- object$lambda
  }
  if(is.null(object$effect)){
    se <- rep("-not fit-",length(object$sp))
  }else{
    se <- object$effect
  }
  if(is.null(object$response)){
    sr <- rep("-not fit-",length(object$sp))
  }else{
    sr <- object$response
  }
  if(is.null(object$lambda_cov)){
    slc <- rep("-not fit-",length(object$sp))
  }else{
    slc <- object$lambda_cov
  }
  if(is.null(object$effect_cov)){
    sec <- rep("-not fit-",length(object$sp))
  }else{
    sec <- object$effect_cov
  }
  if(is.null(object$response_cov)){
    src <- rep("-not fit-",length(object$sp))
  }else{
    src <- object$response_cov
  }
  summary_table <- data.frame(lambda = sl,
                              effect = se,
                              response = sr,
                              lambda_cov = slc,
                              effect_cov = sec,
                              response_cov = src,
                              row.names = object$sp)
  cat("\n")
  summary_table
}

#' CXR summary method for multispecies fits
#'
#' @param object a \code{cxr_pm_multifit} object, from the function with the same name
#' @param ... other arguments, not used
#' @return console output
#' @export
summary.cxr_pm_multifit <- function(object,...){
  cat("model: '",object$model_name,"'",
      "\noptimization method: '",object$optimization_method,"'",
      "\n----------",sep="")
  # for printing null or valid values
  # ifelse returns single values over single conditions
  if(is.null(object$lambda)){
    sl <- rep("- not fit -",nrow(object$data))
  }else{
    sl <- object$lambda
  }
  if(is.null(object$alpha)){
    sa <- rep("- not fit -",nrow(object$data))
  }else{
    sa <- rowMeans(object$alpha)
  }
  if(is.null(object$lambda_cov)){
    slc <- rep("- not fit -",nrow(object$data))
  }else{
    slc <- object$lambda_cov
  }
  if(is.null(object$alpha_cov)){
    sac <- rep("- not fit -",nrow(object$data))
  }else{
    sac <- object$alpha_cov
  }
  summary_table <- data.frame(sp = names(object$data),
                              observations = unlist(lapply(object$data, nrow)),
                              neighbours = unlist(lapply(object$data,length))-1,
                              covariates = ifelse(is.null(object$covariates[[1]]),0,
                                                  ncol(object$covariates[[1]])),
                              # lambda = NULL,
                              lambda = sl,
                              mean_alpha = sa,
                              # mean_lambda_cov = ifelse(is.null(object$lambda_cov),
                              #                          "- not fit - ",
                              #                          mean(object$lambda_cov)),
                              # mean_alpha_cov = ifelse(is.null(object$alpha_cov),
                              #                         "- not fit - ",
                              #                         mean(object$alpha_cov)),
                              row.names = NULL)
  cat("\n")
  summary_table
  
}


