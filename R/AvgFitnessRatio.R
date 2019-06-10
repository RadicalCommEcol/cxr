
#' Average fitness differences between a pair of species
#' Calculates the product of (1) the demographic ratio, and (2) the competitive response ratio between two species
#' according to their vital rates and competition coefficients. First species in the parameters is numerator 
#' (j in eq. 4 of Godoy et al. 2014). If germ.rate and survival.rate are provided, it calculates the ratio
#' according to eq. 4 of Godoy et al. (2014), for annual plants. Otherwise, if only lambda is provided, it
#' returns the general version of the demographic ratio, where nu = lambda
#'
#' @param lambda vector of length 2, per capita fecundity of the species in the absence of competition
#' @param germ.rate optional vector of length 2, germination rate of the two species
#' @param survival.rate optional vector of length 2, annual survival of ungerminated seed in the soil
#' @param pair.matrix 2x2 matrix, competition coefficients between the two species and intraspecific terms
#'
#' @return list with three numeric values, giving 1) the demographic ratio, 2) the competitive response ratio, and 
#' 3) the average fitness ratio between the two species
#' @export
#'
AvgFitnessRatio <- function(lambda, germ.rate = NULL, survival.rate = NULL, pair.matrix){
  if(!is.null(germ.rate) & !is.null(survival.rate)){
    nu <- (lambda*germ.rate)/(1-(1-germ.rate)*survival.rate)
  }else{
    nu <- lambda
  }
  demographic.ratio <- (nu[1]-1)/(nu[2]-1)
  comp.response.ratio <- sqrt((pair.matrix[1,2]/pair.matrix[2,2])*(pair.matrix[1,1]/pair.matrix[2,1]))
  names(demographic.ratio) <- NULL
  names(comp.response.ratio) <- NULL
  
  # return the two components and the resulting ratio
  list(demographic.ratio = demographic.ratio, 
       comp.response.ratio = comp.response.ratio, 
       avg.fitness.ratio = demographic.ratio*comp.response.ratio)
}
