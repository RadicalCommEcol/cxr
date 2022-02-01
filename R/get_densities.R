
#' sort densities
#' 
#' internal function, to pass a vector in a proper form to the vital_rate function.
#' In particular, the first element must be the focal species, and the rest in
#' order (e.g. 2,1,3 if 2 is focal, 3,1,2 if 3 is focal, etc.).
#' It also sums the densities from the three life stages.
#'
#' @param focal.sp integer
#' @param site integer
#' @param current.densities list of length sp, each element is a matrix site*stages
#'
#' @return densities vector
#' @noRd
get_densities <- function(focal.sp,site,current.densities){
  
  if(is.null(focal.sp) | is.null(site) | is.null(current.densities)){
    message(paste("get_densities ERROR: one or more argumets are not specified",sep=""))
    return(NULL)
  }
  
  num.sp <- length(current.densities)
  non.focal.sp <- c(1:num.sp)[-focal.sp]
  
  dens <- rep(NA_real_,num.sp)
  dens[1] <- sum(current.densities[[focal.sp]][site,])
  
  for(i.nf in 2:length(dens)){
    dens[i.nf] <- sum(current.densities[[non.focal.sp[i.nf-1]]][site,])
  }
  
  return(dens)
}
