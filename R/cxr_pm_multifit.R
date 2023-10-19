#' Multi-species parameter optimization
#' 
#' This function is a wrapper for estimating parameters for several
#' focal species, instead of making separate calls to \code{cxr_pm_fit}.
#'
#' @param data named list in which each component is 
#' a dataframe with a fitness column and a number of columns representing neighbours
#' @param focal_column character vector with the same length as data,
#' giving the names of the columns representing
#' intraspecific observations for each species, 
#' or numeric vector giving the position of such columns.
#' @inheritParams cxr_pm_fit
#' @param covariates optional named list in which each component is
#' a dataframe with values of each covariate for each observation. The ith component
#' of \code{covariates} are the covariate values that correspond to 
#' the ith component of \code{data}, so they must have the same number of observations.
#' @param fixed_terms optional named list in which each component is 
#' itself a list containing fixed terms for each focal species.
#'
#' @return an object of class 'cxr_pm_multifit' which is a list with the following components:
#' * model_name: string with the name of the fitness model
#' * model: model function
#' * data: data supplied 
#' * taxa: names of the taxa fitted
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
#' * log_likelihood: log-likelihoods of the fits
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
#' # be explicit about the focal species
#' focal.sp <- names(data)
#' # covariates: salinity
#' data("salinity_list")
#' salinity <- salinity_list[1:3]
#' # keep only salinity column
#' for(i in 1:length(salinity)){
#'   salinity[[i]] <- data.frame(salinity = salinity[[i]][,2:length(salinity[[i]])])
#' }
#' \donttest{
#'   fit_3sp <- cxr_pm_multifit(data = data,
#'                              optimization_method = "bobyqa",
#'                              model_family = "BH",
#'                              focal_column = focal.sp,
#'                              covariates = salinity,
#'                              alpha_form = "pairwise",
#'                              lambda_cov_form = "global",
#'                              alpha_cov_form = "global",
#'                              initial_values = list(lambda = 1,
#'                                                    alpha_intra = 0.1,
#'                                                    alpha_inter = 0.1,
#'                                                    lambda_cov = 0.1, 
#'                                                    alpha_cov = 0.1),
#'                              lower_bounds = list(lambda = 0.01,
#'                                                  alpha_intra = 0,
#'                                                  alpha_inter = 0,
#'                                                  lambda_cov = 0, 
#'                                                  alpha_cov = 0),
#'                              upper_bounds = list(lambda = 100,
#'                                                  alpha_intra = 1,
#'                                                  alpha_inter = 1,
#'                                                  lambda_cov = 1, 
#'                                                  alpha_cov = 1),
#'                              bootstrap_samples = 3)
#'   # brief summary
#'   summary(fit_3sp)
#'   # interaction matrix
#'   fit_3sp$alpha_matrix
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
                                                    "GenSA"), 
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
  
  multifit_data_check <- TRUE
  
  # check input data
  if(!inherits(data,"list")){
    multifit_data_check <- FALSE
  }else{
    
    if(!is.null(covariates)){
      if(!inherits(covariates,"list")){
        multifit_data_check <- FALSE
      }else if(length(covariates) != length(data)){
        multifit_data_check <- FALSE
      }
    }
    
    if(!is.null(fixed_terms)){
      if(!inherits(fixed_terms,"list")){
        multifit_data_check <- FALSE
      }else if(length(fixed_terms) != length(data)){
        multifit_data_check <- FALSE
      }
    }
    
    if(!multifit_data_check){
      message("cxr_pm_multifit ERROR: Check your input data
              1) 'data' is a list;
              2) if present, 'covariates' and/or 'fixed_terms' are also lists;
              3) these lists have the same number of elements, one for each focal species.")
      return(NULL)
    }
    
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
                                 fixed_terms = fixed_terms[[i.sp]])
      if(temp[[1]] == "error"){
        message(paste("cxr_pm_multifit ERROR: check data format for sp ",names(data)[i.sp]," and/or input parameters.\n",sep=""))
        message(paste("more info on the error:\n",temp[[2]],sep=""))
        return(NULL)
      }
    }
  }# if-else
  
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
    Make sure it is defined and available in the cxr package or in the global environment.\n",sep=""))
    return(NULL)
  }
  
  
  # prepare multisp data ----------------------------------------------------
  spnames <- names(data)
  
  # ugly check to see all covariates are named
  if(!is.null(covariates)){
    covnames <- unlist(sapply(covariates,colnames))
    if(any(is.null(covnames))){
      message("at least one covariate column is not named. Please fix this to avoid inconsistencies
              in the fitting process.")
      return(NULL)
    }
  } 
  
  # fit every sp ------------------------------------------------------------
  spfits <- list()
  for(i.sp in 1:length(data)){
    
    spfits[[i.sp]] <- cxr_pm_fit(data = data[[i.sp]],
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
                                 fixed_terms = fixed_terms[[i.sp]],
                                 bootstrap_samples = bootstrap_samples)
    # summary(spfits[[i.sp]])
    # spfits[[i.sp]]$alpha_inter
    # spfits[[i.sp]]$alpha_cov_standard_error
  }
  
  if(length(spfits) == length(spnames)){
    names(spfits) <- spnames
  }  
  
  # output ------------------------------------------------------------------
  
  # in case some or all fits fail
  if(length(spfits) < length(spnames)){
    if(length(spfits) == 0){
      message("cxr_pm_multifit ERROR: Parameter fitting failed for all focal species.
            Please check the data, initial values, and bounds provided.")
    }else{
      validfits <- which(!sapply(spfits,is.null))
      message(paste("cxr_pm_multifit ERROR: Parameter fitting failed for the following taxa:\n",
                    paste(!which(sapply(spfits,is.null)),collapse = ","),
            "\nPlease check the data, initial values, and bounds provided.",sep=""))
    }
    return(NULL)
  }
  
  # if fits are ok,
  # integrate output from the different sp
  
  splambda <- NULL
  alpha_matrix <- NULL
  splambda_cov <- NULL
  spalpha_cov <- NULL
  spllik <- NULL

  for(i.sp in 1:length(spnames)){
    # lambda
    if(!is.null(spfits[[i.sp]]$lambda)){
      mylambda <- spfits[[i.sp]]$lambda
      names(mylambda) <- spnames[i.sp]
      splambda <- c(splambda,mylambda)
    }
  }
  
  # a single alpha matrix, 
  # including all focal and neighbours
  
  # get the names of all sp and sort them
  # these will be the columns of the alpha matrix
  matrix.names <- NULL
  for(i.sp in 1:length(spfits)){
    if(!is.null(spfits[[i.sp]]$alpha_intra)){
      matrix.names <- c(matrix.names,names(spfits[[i.sp]]$alpha_intra))
    }
    if(!is.null(spfits[[i.sp]]$alpha_inter)){
      matrix.names <- c(matrix.names,names(spfits[[i.sp]]$alpha_inter))
    }
  }
  matrix.names <- sort(unique(matrix.names))
  
  # build the matrix
  if(is.null(matrix.names)){
    alpha_matrix <- NULL
  }else{
    alpha_matrix <- matrix(nrow = length(spfits),ncol = length(matrix.names),dimnames = list(names(spfits),matrix.names))
    
    # go sp by sp filling the matrix
    for(i.sp in 1:nrow(alpha_matrix)){
      # intraspecific term
      if(!is.null(spfits[[i.sp]]$alpha_intra)){
        intra.col <- which(colnames(alpha_matrix) == names(spfits[[i.sp]]$alpha_intra))
        alpha_matrix[i.sp,intra.col] <- spfits[[i.sp]]$alpha_intra
      }
      # interspecific terms
      if(!is.null(spfits[[i.sp]]$alpha_inter)){
        inter.col <- match(names(spfits[[i.sp]]$alpha_inter),colnames(alpha_matrix))
        # which(colnames(alpha_matrix) == names(spfits[[i.sp]]$alpha_inter))
        alpha_matrix[i.sp,inter.col] <- spfits[[i.sp]]$alpha_inter
      }
    }
  }# if-else
  
  # lambda_cov also as a matrix, with covariates in columns
  # and focal species in rows
  # splambda_cov <- NULL
  if(!is.null(covariates) & 
     !all(sapply(spfits,function(x){is.null(x$lambda_cov)}))){
    # covariate names
    # in case different sp have different covariates
    cov.names <- NULL
    for(i.cov in 1:length(covariates)){
      if(!is.null(colnames(covariates[[i.cov]]))){
      cov.names <- c(cov.names,colnames(covariates[[i.cov]]))
      }else{
        cov.names <- c(cov.names,paste("cov",i.cov,sep=""))
      }
    }
    cov.names <- sort(unique(cov.names))
    
    splambda_cov <- matrix(nrow = length(spfits),ncol = length(cov.names),dimnames = list(names(spfits),cov.names))
    for(i.sp in 1:length(spnames)){
      # lambda_cov
      if(!is.null(spfits[[i.sp]]$lambda_cov)){
        for(i.cov in 1:length(cov.names)){
          # look for i.cov covariate in the vector of lambda_covs affecting i.sp
          # and place it in the matrix
          spfitnames <- substr(names(spfits[[i.sp]]$lambda_cov),12,nchar(names(spfits[[i.sp]]$lambda_cov)))
          splambda_cov[i.sp,i.cov] <- spfits[[i.sp]]$lambda_cov[which(spfitnames == cov.names[i.cov])]
          # splambda_cov[i.sp,i.cov] <- spfits[[i.sp]]$lambda_cov[which(grepl(cov.names[i.cov],names(spfits[[i.sp]]$lambda_cov)))]
        }
      }
    }# for i.sp
  }# if covariates
  
  # alpha_cov
  
  # this is a list where, for each cov, a focalsp x neighsp matrix gives the alpha_covs
  # if a global alpha_cov is fitted, does not matter
  
  # spalpha_cov <- NULL
  if(!is.null(covariates) & 
     !all(sapply(spfits,function(x){is.null(x$alpha_cov)}))){
    
    # covariate names
    # in case different sp have different covariates
    cov.names <- NULL
    for(i.cov in 1:length(covariates)){
      cov.names <- c(cov.names,colnames(covariates[[i.cov]]))
    }
    cov.names <- sort(unique(cov.names))
    
    # get the names of all sp and sort them
    # these will be the columns of the alpha_cov matrix
    matrix.names <- NULL
    for(i.sp in 1:length(spfits)){
      if(!is.null(spfits[[i.sp]])){
        matrix.names <- c(matrix.names,spnames[i.sp],names(data[[i.sp]]))
      }
    }
    matrix.names <- sort(unique(matrix.names[which(!matrix.names == "fitness")]))
    
    # template
    # I could simply copy the alpha one, but not sure if its robust enough
    # for the cases in which alpha is fixed. anyway.
    ac_matrix_template <- matrix(nrow = length(spfits),ncol = length(matrix.names),dimnames = list(spnames,matrix.names))
    
    spalpha_cov <- list()
    for(i.cov in 1:length(cov.names)){
      spalpha_cov[[i.cov]] <- ac_matrix_template
      
      # fill up matrix
      for(i.sp in 1:length(spnames)){

        # spfitnames <- substr(names(spfits[[i.sp]]$alpha_cov),12,nchar(names(spfits[[i.sp]]$alpha_cov)))
        # splambda_cov[i.sp,i.cov] <- spfits[[i.sp]]$lambda_cov[which(spfitnames == cov.names[i.cov])]

        listpos <- which(names(spfits[[i.sp]]$alpha_cov) == cov.names[i.cov])
        mycov <- spfits[[i.sp]]$alpha_cov[[listpos]]

        if(length(mycov)==1){
          # same alpha_cov for all interactions
          spalpha_cov[[i.cov]][i.sp,] <- mycov
        }else{
          # specific alpha_cov
          # sort sp just in case
          sp.pos <- sapply(matrix.names,function(x)which(grepl(x,names(mycov))))
          spalpha_cov[[i.cov]][i.sp,] <- mycov[sp.pos]
        }
      }# for each focal sp
    }# for each covariate
    names(spalpha_cov) <- cov.names
  }
  
  # log-likelihood
  # spllik <- NULL
  for(i.sp in 1:length(spnames)){
    myllik <- spfits[[i.sp]]$log_likelihood
    if(!is.null(myllik)){
      names(myllik) <- spnames[i.sp]
      spllik <- c(spllik,myllik)
    }
  }
  
  # error output ------------------------------------------------------------
  
  # lambda
  er_splambda <- NULL
  for(i.sp in 1:length(spnames)){
    # lambda
    if(!is.null(spfits[[i.sp]]$lambda_standard_error)){
      myer_lambda <- spfits[[i.sp]]$lambda_standard_error
      names(myer_lambda) <- spnames[i.sp]
      er_splambda <- c(er_splambda,myer_lambda)
    }
  }
  
  # a single alpha matrix, 
  # including all focal and neighbours
  
  # get the names of all sp and sort them
  # these will be the columns of the alpha matrix
  er_matrix.names <- NULL
  
  for(i.sp in 1:length(spfits)){
    if(!is.null(spfits[[i.sp]]$alpha_intra_standard_error)){
      er_matrix.names <- c(er_matrix.names,names(spfits[[i.sp]]$alpha_intra_standard_error))
    }
    if(!is.null(spfits[[i.sp]]$alpha_inter_standard_error)){
      er_matrix.names <- c(er_matrix.names,names(spfits[[i.sp]]$alpha_inter_standard_error))
    }
  }
  er_matrix.names <- sort(unique(er_matrix.names))
  
  # build the matrix
  if(is.null(er_matrix.names)){
    er_alpha_matrix <- NULL
  }else{
    er_alpha_matrix <- matrix(nrow = length(spfits),
                              ncol = length(er_matrix.names),
                              dimnames = list(names(spfits),er_matrix.names))
    
    # go sp by sp filling the matrix
    for(i.sp in 1:nrow(er_alpha_matrix)){
      # intraspecific term
      if(!is.null(spfits[[i.sp]]$alpha_intra_standard_error)){
        intra.col <- which(colnames(er_alpha_matrix) == 
                             names(spfits[[i.sp]]$alpha_intra_standard_error))
        er_alpha_matrix[i.sp,intra.col] <- spfits[[i.sp]]$alpha_intra_standard_error
      }
      # interspecific terms
      if(!is.null(spfits[[i.sp]]$alpha_inter_standard_error)){
        inter.col <- match(names(spfits[[i.sp]]$alpha_inter_standard_error),
                           colnames(er_alpha_matrix))
        # which(colnames(alpha_matrix) == names(spfits[[i.sp]]$alpha_inter))
        er_alpha_matrix[i.sp,inter.col] <- spfits[[i.sp]]$alpha_inter_standard_error
      }
    }
  }# if-else
  
  # lambda_cov also as a matrix, with covariates in columns
  # and focal species in rows
  
  er_splambda_cov <- NULL
  if(!is.null(covariates) & 
     !all(sapply(spfits,function(x){is.null(x$lambda_cov_standard_error)}))){
    # covariate names
    # in case different sp have different covariates
    cov.names <- NULL
    for(i.cov in 1:length(covariates)){
      cov.names <- c(cov.names,colnames(covariates[[i.cov]]))
    }
    cov.names <- sort(unique(cov.names))
    
    er_splambda_cov <- matrix(nrow = length(spfits),
                              ncol = length(cov.names),
                              dimnames = list(names(spfits),cov.names))
    for(i.sp in 1:length(spnames)){
      # lambda_cov
      if(!is.null(spfits[[i.sp]]$lambda_cov_standard_error)){
        for(i.cov in 1:length(cov.names)){
          # look for i.cov covariate in the vector of lambda_covs affecting i.sp
          # and place it in the matrix
          spfitnames <- substr(names(spfits[[i.sp]]$lambda_cov_standard_error),12,(nchar(names(spfits[[i.sp]]$lambda_cov_standard_error))-3))
          er_splambda_cov[i.sp,i.cov] <- spfits[[i.sp]]$lambda_cov_standard_error[which(spfitnames == cov.names[i.cov])]
          
          # grepl does not return exact match
          # er_splambda_cov[i.sp,i.cov] <- 
          #   spfits[[i.sp]]$lambda_cov_standard_error[which(
          #     grepl(cov.names[i.cov],
          #           names(spfits[[i.sp]]$lambda_cov_standard_error)))]
        }
      }
    }# for i.sp
  }# if covariates
  
  # alpha_cov
  
  # this is a list where, for each cov, a focalsp x neighsp matrix gives the alpha_covs
  # if a global alpha_cov is fitted, does not matter
  
  er_spalpha_cov <- NULL
  if(!is.null(covariates) & 
     !all(sapply(spfits,function(x){is.null(x$alpha_cov_standard_error)}))){
    
    # covariate names
    # in case different sp have different covariates
    cov.names <- NULL
    for(i.cov in 1:length(covariates)){
      cov.names <- c(cov.names,colnames(covariates[[i.cov]]))
    }
    cov.names <- sort(unique(cov.names))
    
    # get the names of all sp and sort them
    # these will be the columns of the alpha_cov matrix
    er_matrix.names <- NULL
    for(i.sp in 1:length(spfits)){
      if(!is.null(spfits[[i.sp]])){
        er_matrix.names <- c(er_matrix.names,spnames[i.sp],names(data[[i.sp]]))
      }
    }
    er_matrix.names <- sort(unique(matrix.names[which(!matrix.names == "fitness")]))
    er_matrix.names <- paste(er_matrix.names,"_se",sep="")
    
    # template
    # I could simply copy the alpha one, but not sure if its robust enough
    # for the cases in which alpha is fixed. anyway.
    er_ac_matrix_template <- matrix(nrow = length(spfits),
                                    ncol = length(er_matrix.names),
                                    dimnames = list(spnames,er_matrix.names))
    
    er_spalpha_cov <- list()
    for(i.cov in 1:length(cov.names)){
      er_spalpha_cov[[i.cov]] <- er_ac_matrix_template
      
      # fill up matrix
      for(i.sp in 1:length(spnames)){
        
        # mycov <- spfits[[i.sp]]$alpha_cov_standard_error[[cov.names[i.cov]]]
        if(!is.null(spfits[[i.sp]]$alpha_cov_standard_error)){
          
          listpos <- which(names(spfits[[i.sp]]$alpha_cov_standard_error) == cov.names[i.cov])
          mycov <- spfits[[i.sp]]$alpha_cov_standard_error[[listpos]]
          
          if(length(mycov)==1){
            # same alpha_cov for all interactions
            er_spalpha_cov[[i.cov]][i.sp,] <- mycov
          }else{
            # specific alpha_cov
            # sort sp just in case
            sp.pos <- sapply(er_matrix.names,function(x)which(grepl(x,names(mycov))))
            er_spalpha_cov[[i.cov]][i.sp,] <- mycov[sp.pos]
          }
        }
      }# for each focal sp
    }# for each covariate
    names(er_spalpha_cov) <- cov.names
  }
  
  
  # prepare cxr_pm_multifit object ------------------------------------------
  
  # return a cxr_pm_multifit object
  # which is basically the same as the base object
  # but with info on more sp. e.g. lambda is a 1d vector,
  # alpha is a matrix including the intra and interspecific terms
  
  list_names <- c("model_name",
                  "model",
                  "data",
                  "taxa",
                  "covariates",
                  "optimization_method",
                  "initial_values",
                  "fixed_terms",
                  "lambda","alpha_matrix","lambda_cov","alpha_cov",
                  "lambda_standard_error","alpha_standard_error",
                  "lambda_cov_standard_error","alpha_cov_standard_error",
                  "log_likelihood")
  
  fit <- sapply(list_names,function(x) NULL)
  
  fit$model_name <- spfits[[1]]$model_name
  fit$model <- spfits[[1]]$fitness_model
  fit$data <- data
  fit$taxa <- spnames
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
  if(!is.null(alpha_matrix)){
    fit$alpha_matrix <- alpha_matrix
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
  if(!is.null(er_alpha_matrix)){
    fit$alpha_matrix_standard_error <- er_alpha_matrix
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
