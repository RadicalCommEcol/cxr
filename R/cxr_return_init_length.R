cxr_return_init_length <- function(par_type, par_value, pairwise_length){
  
  # global means single param
    if(par_type == "global"){
    if(length(par_value) != 1){
      return_par <- par_value[1]
    }
      # else, pairwise means a single param for each pair
      # so extend if needed
    }else if(par_type == "pairwise"){
    if(length(par_value) == 1){
      return_par <- rep(par_value,pairwise_length)
    }
  }# par_type
  
}