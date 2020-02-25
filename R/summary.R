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
    "\nfocal taxa ID: ",ifelse(is.null(object$focal_ID),"-",object$focal_ID),
    "\nobservations: ",nrow(object$data),
    "\nneighbours: ",length(object$data)-1,
    "\ncovariates: ",ifelse(is.null(object$covariates),0,ncol(object$covariates)),
    "\n----------",
    "\nfocal lambda: ",ifelse(is.null(object$lambda)," - not fit - ",object$lambda),
    "\nalpha_intra: ",ifelse(is.null(object$alpha_intra)," - not fit - ",object$alpha_intra),
    "\nmean alpha_inter: ",ifelse(is.null(object$alpha_inter)," - not fit - ",mean(object$alpha_inter)),
    "\nmean lambda_cov: ",ifelse(is.null(object$lambda_cov),"- not fit - ",mean(object$lambda_cov)),
    "\nmean alpha_cov: ",ifelse(is.null(object$alpha_cov),"- not fit - ",mean(unlist(object$alpha_cov))),
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
  # for printing null or valid values
  # ifelse returns single values over single conditions
  if(is.null(object$lambda)){
    sl <- rep("- not fitted -",length(object$data))
  }else{
    sl <- object$lambda
  }
  if(is.null(object$alpha_matrix)){
    sa <- "- not fit -"
  }else{
    sa <- object$alpha_matrix
  }
  if(is.null(object$lambda_cov)){
    slc <- as.matrix(rep("- not fit -",length(object$data)))
    colnames(slc) <- "lambda_cov"
  }else{
    slc <- object$lambda_cov
    colnames(slc) <- paste("lambda_cov_",colnames(slc),sep="")
  }
  if(is.null(object$alpha_cov)){
    sac <- as.matrix(rep("- not fit -",length(object$data)))
    colnames(sac) <- "alpha_cov"
  }else{
    # average over covariates and alpha values, 
    # hence the double rowmeans
    sac <- sapply(object$alpha_cov,function(x)rowMeans(x))
    colnames(sac) <- paste("mean_alpha_cov_",colnames(sac),sep="")
    
  }
  summary_table <- data.frame(sp = names(object$data),
                              observations = unlist(lapply(object$data, nrow)),
                              neighbours = unlist(lapply(object$data,length))-1,
                              covariates = ifelse(is.null(object$covariates[[1]]),0,
                                                  ncol(object$covariates[[1]])),
                              # lambda = NULL,
                              lambda = sl,
                              # alpha_intra = sa,
                              # mean_alpha_inter = sai,
                              # mean_lambda_cov = slc,
                              row.names = NULL)
  summary_table <- cbind(summary_table,slc,sac)
  row.names(summary_table) <- NULL
  
  cat("model: '",object$model_name,"'",
      "\noptimization method: '",object$optimization_method,"'",
      "\n----------",sep="")
  cat("\n")
  print(summary_table)
  cat("\n----------\nalpha matrix:\n")
  object$alpha_matrix
}


