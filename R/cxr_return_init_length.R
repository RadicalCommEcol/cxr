cxr_return_init_length <- function(par_type, par_value, par_names, fit_type = c("er","pm")){
  
  fit_type <- match.arg(fit_type)
  
  return_par <- NA
  
  if(fit_type == "er"){
    if(par_type == "global"){
      if(length(par_value) == 1){
        return_par <- rep(par_value,length(par_names))
      }else{
        return_par <- par_value
      }
    }
  }else{
  # global means single param
    if(par_type == "global"){
    if(length(par_value) != 1){
      return_par <- par_value[1]
    }else{
      return_par <- par_value
    }
      # else, pairwise means a single param for each pair
      # so extend if needed
    }else if(par_type == "pairwise"){
    if(length(par_value) == 1){
      return_par <- rep(par_value,length(par_names))
    }else{
      return_par <- par_value
    }
  }# par_type
  }
  names(return_par) <- par_names
  return_par
}