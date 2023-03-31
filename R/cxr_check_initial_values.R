#' Internal, check the initial values and bounds provided to cxr_pm_fit
#'
#' @inheritParams cxr_pm_fit
#'
#' @return boolean, whether the values are consistent
#' @noRd
cxr_check_initial_values <- function(initial_values,
                                     focal_column,
                                     lower_bounds,
                                     upper_bounds,
                                     covariates,
                                     fixed_terms){
  iv.ok <- TRUE
  
  # double check, this should already be in cxr_pm_fit as match.arg
  valid.names <- c("lambda","alpha_intra","alpha_inter","lambda_cov","alpha_cov")
  if(!all(names(initial_values) %in% valid.names)){
    iv.ok <- FALSE
  }
  
  # focal_column and alpha_intra come together
  # either both specified or none
  if((is.null(focal_column) & "alpha_intra" %in% names(initial_values)) | 
     (!is.null(focal_column) & !("alpha_intra" %in% names(initial_values)))){
    iv.ok <- FALSE
  }
  
  # either both or none lower/upper
  if(any(is.null(lower_bounds),is.null(upper_bounds)) & 
     !all(is.null(lower_bounds),is.null(upper_bounds))){
    iv.ok <- FALSE
  }
  
  # identical elements in the three lists
  if(!is.null(upper_bounds) & !is.null(lower_bounds)){
    if(!identical(names(initial_values),names(lower_bounds)) |
       !identical(names(initial_values),names(upper_bounds)) |
       !identical(names(upper_bounds),names(lower_bounds))){
      iv.ok <- FALSE
    }
  }
  
  # check that the number of initial values is equal to either:
  # 1, so that all parameters have same initial values
  # number of covariates, so that each covariate effect has different starting point
  if(!is.null(covariates)){
    
    if("lambda_cov" %in% names(initial_values) |
       "alpha_cov" %in% names(initial_values)){
      
      if("lambda_cov" %in% names(initial_values)){
        il <- length(initial_values$lambda_cov)
        if(!il %in% c(1,ncol(covariates))){
          iv.ok <- FALSE
        }
      }
      
      if("alpha_cov" %in% names(initial_values)){
        ia <- length(initial_values$alpha_cov)
        if(!ia %in% c(1,ncol(covariates))){
          iv.ok <- FALSE
        }
      }
    }
  }
  
  # check fixed terms is a list
  # and if a param is here, it is not in
  # initial values or bounds
  if(!is.null(fixed_terms)){
    if(!is.list(fixed_terms)){
      iv.ok <- FALSE
    }else{
      if("lambda" %in% names(fixed_terms)){
        if(!is.numeric(fixed_terms[["lambda"]])){
          iv.ok <- FALSE
        }
        if("lambda" %in% c(names(initial_values), 
                           names(lower_bounds), 
                           names(upper_bounds))){
          iv.ok <- FALSE
        }
      }
      
      if("alpha_intra" %in% names(fixed_terms)){
        if(!is.numeric(fixed_terms[["alpha_intra"]])){
          iv.ok <- FALSE
        }
        if("alpha_intra" %in% c(names(initial_values), 
                                names(lower_bounds), 
                                names(upper_bounds))){
          iv.ok <- FALSE
        }
      }
      
      if("alpha_inter" %in% names(fixed_terms)){
        if(!is.numeric(fixed_terms[["alpha_inter"]])){
          iv.ok <- FALSE
        }
        if("alpha_inter" %in% c(names(initial_values), 
                                names(lower_bounds), 
                                names(upper_bounds))){
          iv.ok <- FALSE
        }
      }
      
      if("lambda_cov" %in% names(fixed_terms)){
        if(!is.numeric(fixed_terms[["lambda_cov"]])){
          iv.ok <- FALSE
        }
        if("lambda_cov" %in% c(names(initial_values), 
                               names(lower_bounds), 
                               names(upper_bounds))){
          iv.ok <- FALSE
        }
      }
      
      if("alpha_cov" %in% names(fixed_terms)){
        if(!is.numeric(fixed_terms[["alpha_cov"]])){
          iv.ok <- FALSE
        }
        if("alpha_cov" %in% c(names(initial_values), 
                              names(lower_bounds), 
                              names(upper_bounds))){
          iv.ok <- FALSE
        }
      }
      
    }
  }
  
  # check that the number of lower/upper bounds is equal to either
  # 1, so that all parameters have same bounds
  # number of covariates, so that each covariate has different boundaries
  
  if(!is.null(upper_bounds) & !is.null(lower_bounds)){
    if(!is.null(covariates)){
      if("lambda_cov" %in% names(upper_bounds) & "lambda_cov" %in% names(lower_bounds)){
        lu <- length(upper_bounds$lambda_cov)
        ll <- length(lower_bounds$lambda_cov)
        if(lu != ll | !lu %in% c(1,ncol(covariates))){
          iv.ok <- FALSE
        }
      }
      
      if("alpha_cov" %in% names(upper_bounds) & "alpha_cov" %in% names(lower_bounds)){
        au <- length(upper_bounds$alpha_cov)
        al <- length(lower_bounds$alpha_cov)
        if(au != al | !au %in% c(1,ncol(covariates))){
          iv.ok <- FALSE
        }
      }
      
    }else{
      # if no covariates
      # check all list elements have length 1
      l1 <- all(sapply(lower_bounds,length) == 1)
      a1 <- all(sapply(upper_bounds,length) == 1)
      
      if(!l1 | !a1){
        iv.ok <- FALSE
      }
    }
    
  }
  
  iv.ok
  
}
