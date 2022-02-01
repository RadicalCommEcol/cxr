
#' Obtain species densities from transition matrices
#' 
#' Using the vec-permutation approach as defined in: 
#' Hunter and Caswell 2005, doi:10.1016/j.ecolmodel.2005.05.002, 
#' Ozgul et al. 2009, doi: 10.1086/597225
#' In particular, it uses the arrangement by patches, and calculates
#' first demography, then dispersal (Table 1 of Hunter and Caswell 2005).
#'
#' @param focal.sp integer, focal species
#' @param vpm data structure holding all vector-permutation matrices; see 
#' `vec_permutation_matrices`. If not in an appropriate format, it is likely
#' to fail without warning.
#' @param current.densities list of length sp, each element is a matrix sites*stages. If
#' not in that format, it is likely to fail without warning.
#' 
#' @return matrix of sites x stages, each element is the density of a given life stage
#' (juvenile, non-reproductive adult, reproductive adult) at a given site.
#' @export
calculate_densities <- function(focal.sp,vpm,current.densities){
  
  if(is.null(focal.sp) | is.null(vpm) | is.null(current.densities)){
    message(paste("calculate_densities ERROR: one or more argumets are not specified",sep=""))
    return(NULL)
  }
  
  # transform to numeric vector
  dv <- c(t(current.densities[[focal.sp]]))
  d2 <- t(vpm[["permutation"]][[focal.sp]])%*%vpm[["dispersal"]][[focal.sp]]%*%vpm[["permutation"]][[focal.sp]]%*%vpm[["demography"]][[focal.sp]]%*%dv
  d3 <- matrix(d2,ncol = 3,byrow = T)
  return(d3)
}
