
#' Average fitness differences between a pair of species
#' Calculates the product of (1) the demographic ratio, and (2) the competitive response ratio between two species
#' according to their vital rates and competition coefficients. First species in the parameters is numerator 
#' (j in eq. 4 of Godoy et al. 2014).
#'
#' @param lambda vector of length 2, per germinant fecundity of the species in the absence of competition
#' @param germ.rate vector of length 2, germination rate of the two species
#' @param survival.rate vector of length 2, annual survival of ungerminated seed in the soil
#' @param pair.matrix 2x2 matrix, competition coefficients between the two species and intraspecific terms
#'
#' @return single value giving the average fitness ratio between the two species
#' @export
#'
#' @examples
AvgFitnessRatio <- function(lambda, germ.rate, survival.rate, pair.matrix){
  nu <- (lambda*germ.rate)/(1-(1-germ.rate)*survival.rate)
  demographic.ratio <- (nu[1]-1)/(nu[2]-1)
  comp.response.ratio <- sqrt((pair.matrix[1,2]/pair.matrix[2,2])*(pair.matrix[1,1]/pair.matrix[2,1]))
  demographic.ratio*comp.response.ratio
}
