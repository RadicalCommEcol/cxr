#' niche overlap function for Lqw-Watkinson (LW) models
#'
#' @param pair_matrix pairwise interaction matrix
#'
#' @return numeric value
#' @export
#'
#' @examples LW_niche_overlap(matrix(c(.1,.2,.03,.3),nrow = 2))
LW_niche_overlap <- function(pair_matrix){
  sqrt((pair_matrix[1,2]/pair_matrix[2,2])*(pair_matrix[2,1]/pair_matrix[1,1]))
}