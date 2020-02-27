#' niche overlap function for a suite of models (LV,RK,LW,BH)
#'
#' formulation from the modern coexistence theory (e.g. Chesson 2012)
#'
#' @param pair_matrix pairwise interaction matrix
#'
#' @return numeric value
#' @export
#'
#' @examples niche_overlap_MCT(matrix(c(.1,.2,.03,.3),nrow = 2))
niche_overlap_MCT <- function(pair_matrix){
  sqrt((pair_matrix[1,2]/pair_matrix[2,2])*(pair_matrix[2,1]/pair_matrix[1,1]))
}