#' Fitness of a species
#' 
#' Calculates the fitness of a species following eq. S2 of Godoy et al. (2014). If germ.rate and survival.rate are provided, 
#' it calculates nu according to eq. 4 of Godoy et al. (2014), for annual plants. Otherwise, if only lambda is provided, nu = lambda
#'
#' @param lambda per capita fecundity of the species in the absence of competition
#' @param germ.rate optional, germination rate of the species
#' @param survival.rate optional, annual survival of ungerminated seed in the soil
#' @param competitive.response parameter reflecting the species' sensitivity to competition
#'
#' @return single numeric value, species fitness
#' @export
#'
SpeciesFitness <- function(lambda, germ.rate = NULL, survival.rate = NULL,competitive.response){
  
  if(!is.null(germ.rate) & !is.null(survival.rate)){
    nu <- (lambda*germ.rate)/(1-(1-germ.rate)*survival.rate)
  }else{
    nu <- lambda
  }
  
  (nu-1)/competitive.response
  
}