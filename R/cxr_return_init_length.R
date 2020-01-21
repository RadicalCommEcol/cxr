cxr_return_init_length <- function(par_type, par_value, par_names){
  return_par <- NA
  # global means single param
    if(par_type == "global"){
    if(length(par_value) != 1){
      return_par <- par_value[1]
    }
      # else, pairwise means a single param for each pair
      # so extend if needed
    }else if(par_type == "pairwise"){
    if(length(par_value) == 1){
      return_par <- rep(par_value,length(par_names))
    }
  }# par_type
  names(return_par) <- par_names
  return_par
}