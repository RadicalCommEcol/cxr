#' Fill the vec-permutation demography matrix
#' 
#' Fill for a given species, across all sites.
#'
#' @param focal.sp integer, focal species.
#' @param vpm data structure holding all vector-permutation matrices; see 
#' `vec_permutation_matrices`. If not in an appropriate format, it is likely
#' to fail without warning.
#' @param transition_matrices nested list species x sites, in which each element
#' holds a 3x3 transition matrix. If not in that format, it is likely to fail
#' without warning.
#'
#' @return vec-permutation demography matrix for a given species across sites.
#' @export
fill_demography_matrix <- function(focal.sp,vpm,transition_matrices){
  
  if(is.null(focal.sp) | is.null(vpm) | is.null(transition_matrices)){
    message(paste("fill_demography_matrix ERROR: one or more argumets are not specified",sep=""))
    return(NULL)
  }
  
  temp.dem.orig <- vpm[["demography"]][[focal.sp]]
  
  temp.dem <- as.matrix(Matrix::bdiag(transition_matrices[[focal.sp]]))
  
  if(nrow(temp.dem.orig) != nrow(temp.dem) | 
     ncol(temp.dem.orig) != ncol(temp.dem)){
    stop("fill_demography_matrix ERROR: vpm matrix and the resulting demography matrix have different dimensions.")
  }
  
  return(temp.dem)
}
