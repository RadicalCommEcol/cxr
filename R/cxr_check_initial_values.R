cxr_check_initial_values <- function(initial_values,
                                     focal_column,
                                     lower_bounds,
                                     upper_bounds,
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
  
  iv.ok
  
}