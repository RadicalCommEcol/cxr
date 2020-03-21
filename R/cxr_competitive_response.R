#' competitive response ratio
#' 
#' according to eq. 4 of Godoy et al. (2014)
#' 
#' @param pair_matrix 2x2 matrix, competition coefficients between the two species and intraspecific terms
#'
#' @return single numeric value or NA if any coefficient is negative (facilitation)
#' @noRd
#'
cxr_competitive_response <- function(pair_matrix){
  if(all(pair_matrix>=0)){
    return(sqrt((pair_matrix[2,1]/pair_matrix[1,1]) * (pair_matrix[2,2]/pair_matrix[1,2])))
  }else{
    return(NA_real_)
  }
}