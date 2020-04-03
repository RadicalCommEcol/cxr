#' niche overlap function from a structural perspective
#'
#' formulation from Saavedra et al. (2017) Ecological Monographs
#'
#' @param pair_matrix pairwise interaction matrix
#'
#' @return numeric value or NA if error
#' @noRd
#'
#' @examples niche_overlap_SA(matrix(c(.1,.2,.03,.3),nrow = 2))
niche_overlap_SA <- function(pair_matrix){
  if(any(is.na(pair_matrix))){
    message("niche_overlap: one or more interaction coefficients are NA.")
    return(NA_real_)
  }
  n <- nrow(pair_matrix)
  Sigma <- try(solve(t(pair_matrix) %*% pair_matrix, tol = 1e-40),silent = TRUE)
  if(class(Sigma) != "try-error"){
    d <- mvtnorm::pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
    out <- log10(d[1]) + n * log10(2)
    return(1-10^out) 
  }else{
    message("niche_overlap: Structural niche overlap could not be numerically calculated,
            most likely because the pair matrix provided is singular/has determinant 0.")
    return(NA_real_)
  }
}