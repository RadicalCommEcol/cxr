#' Niche overlap between two species
#'
#' quoting Godoy et al. (2014):
#' reflects the average degree to which species limit individuals of their own species relative to competitors. 
#' Low niche overlap causes species to have greater per capita growth rates when rare than when common. 
#' If species limit individuals of their own species and their competitors equally, then niche overlap is 1, 
#' and coexistence is not possible unless species are otherwise identical. 
#' At the other extreme, if species have no interspecific effects, then niche overlap is 0.
#'
#' @param pair.matrix 2x2 matrix with competition coefficients between the two species, and intraspecific terms
#'
#' @return niche overlap value, in the range 0-1.
#' @export
NicheOverlap <- function(pair.matrix){
  sqrt((pair.matrix[1,2]/pair.matrix[2,2])*(pair.matrix[2,1]/pair.matrix[1,1]))
}