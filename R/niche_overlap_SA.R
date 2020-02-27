#' niche overlap function from a structural perspective
#'
#' formulation from Saavedra et al. (2017) Ecological Monographs
#'
#' @param pair_matrix pairwise interaction matrix
#'
#' @return numeric value
#' @export
#'
#' @examples niche_overlap_MCT(matrix(c(.1,.2,.03,.3),nrow = 2))
niche_overlap_SA <- function(pair_matrix){
  n <- nrow(pair_matrix)
  Sigma <-solve(t(pair_matrix) %*% pair_matrix, tol = 1e-40)
  d <- mvtnorm::pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(10^out) 
}