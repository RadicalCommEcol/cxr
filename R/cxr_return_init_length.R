#' internal, return a parameter with its appropriate length
#'
#' @param par_type global or pairwise
#' @param par_value numeric value
#' @param par_names character vector with the names of the parameter
#' @param fit_type "global" params have different lengths in
#' er and pm fits 
#'
#' @return numeric vector of appropriate length
#' @noRd
#'
cxr_return_init_length <- function(par_type, 
                                   par_value, 
                                   par_names, 
                                   fit_type = c("er","pm")){
  
  fit_type <- match.arg(fit_type)
  
  return_par <- NA
  
  if(fit_type == "er"){
    if(par_type == "global"){
      if(length(par_value) == 1){
        return_par <- rep(par_value,length(par_names))
      }else{
        if(length(par_value) == length(par_names)){
          return_par <- par_value
        }else{
          return_par <- rep(par_value[1],length(par_names))
        }
      }
    }
  }else{
  # global means single param
    if(par_type == "global"){
    if(length(par_value) != 1){
      if(length(par_value) == length(par_names)){
        return_par <- par_value
      }else{
        return_par <- rep(par_value[1],length(par_names))
      }
    }else{
      return_par <- rep(par_value,length(par_names))
    }
      # else, pairwise means a single param for each pair
      # so extend if needed
    }else if(par_type == "pairwise"){
    if(length(par_value) == 1){
      return_par <- rep(par_value,length(par_names))
    }else{
      if(length(par_value) == length(par_names)){
        return_par <- par_value
      }else{
        num.names <- length(par_names)
        num.values <- length(par_value)
        if(num.names %% num.values == 0){
          num.rep <- num.names/num.values
          return_par <- rep(par_value,each = num.rep)
        }else{
          return_par <- rep(par_value[1],num.names)
        }
      }
    }
  }# par_type
  }
  names(return_par) <- par_names
  return_par
}