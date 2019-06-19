#' Competitive ability of a species
#' 
#' Calculates the competitive ability of a species following eq. S2 of Godoy et al. (2014). If germ.rate and survival.rate are provided, 
#' it calculates nu according to eq. 4 of Godoy et al. (2014), for annual plants. Otherwise, if only lambda is provided, nu = lambda
#'
#' @param lambda per capita fecundity of the species in the absence of competition
#' @param germ.rate germination rate of the species
#' @param survival.rate annual survival of ungerminated seed in the soil
#' @param competitive.response parameter reflecting the species' sensitivity to competition
#'
#' @return single numeric value, competitive ability of the species
#' @export
#'
CompetitiveAbility <- function(lambda, germ.rate = NULL, survival.rate = NULL,competitive.response){
  
  if(!is.null(germ.rate) & !is.null(survival.rate)){
    nu <- (lambda*germ.rate)/(1-(1-germ.rate)*survival.rate)
  }else{
    nu <- lambda
  }
  
  (nu-1)/competitive.response
  
}