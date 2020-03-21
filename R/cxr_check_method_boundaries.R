#' Internal, check whether a combination of optimization method and bounds is valid
#' 
#' Some methods need explicit lower/upper bounds, but cxr_pm/er fit
#' have default NULL values. This function checks that, for such methods,
#' lower and upper bounds are given as lists with appropriate components.
#' Note that it does not check the values of the bounds themselves.
#' 
#' @param optimization_method character
#' @param lower_bounds either NULL or a list
#' @param upper_bounds either NULL or a list
#'
#' @return boolean, whether appropriate lower/upper bounds are provided.
#' @noRd
cxr_check_method_boundaries <- function(optimization_method,
                                        lower_bounds = NULL,
                                        upper_bounds = NULL,
                                        type = c("pm","er")){
  method.ok <- TRUE
  if(optimization_method %in% c("L-BFGS-B", "nlm", "nlminb", 
                                "Rcgmin", "Rvmmin", "spg", 
                                "bobyqa", "nmkb", "hjkb",
                                "nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND",
                                "GenSA","hydroPSO","DEoptimR") & 
     (is.null(lower_bounds) | is.null(upper_bounds))){
    method.ok <- FALSE
  }else if(optimization_method %in% c("BFGS", "CG", "Nelder-Mead", "ucminf") & 
           (!is.null(lower_bounds) | !is.null(upper_bounds))){
    method.ok <- FALSE
  }else if(!is.null(lower_bounds) & !is.null(upper_bounds)){
    if(type == "pm"){
      bnames <- c("lambda", "alpha_intra","alpha_inter", "lambda_cov", "alpha_cov")
    }else{
      bnames <- c("lambda", "effect", "response", "lambda_cov", "effect_cov", "response_cov")
    }
    my.names <- unique(c(names(lower_bounds),names(upper_bounds)))#,"yeah"))
    bad.names <- sum(!my.names %in% bnames)
    if(bad.names > 0){
      method.ok <- FALSE
    }
  }
  method.ok
}