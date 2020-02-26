#' Fitness of a species
#' 
#' Calculates the fitness of a species, a.k.a. its competitive ability (Godoy et al. 2014, Hart et al. 2018 Journal of Ecology).
#' Note that its definition is model-specific, i.e. it depends on the model family from which interaction coefficients were estimated.
#' The function given here assumes a community of n-species, so that species fitness is calculated according to a 
#' general competitive response (r) substituting the 2-sp denominator terms of table A1 of Hart et al. 2018. 
#' This competitive response can be calculated with the function 'cxr_er_fit'.
#' 
#' Thus, the function accepts two sets of parameters. First, a 'cxr_er_fit' object returned from that function. In this case,
#' species fitness will be calculated for all focal taxa included in the 'cxr_er_fit' object.
#' 
#' Otherwise, users may enter a specification of the model to follow as well as lambda and competitive response parameters.
#' Note that there is no 'default' way of calculating species fitness without specifying the underlying model.
#' 
#' @param effect_response_fit cxr_er_fit object with valid lambda and response terms.
#' @param model_family model with which to calculate species fitness.
#' @param lambda per capita fecundity of the species in the absence of competition.
#' @param competitive_response parameter reflecting the species' sensitivity to competition.
#'
#' @return single numeric value, species fitness
#' @export
#'
species_fitness <- function(effect_response_fit = NULL, 
                            model_family = NULL, 
                            lambda = NULL, 
                            competitive_response = NULL){
  res <- NULL
  if(!is.null(effect_response_fit)){
    
    sf_f <- substr(effect_response_fit$model_name,4,5)
    
    sf_fun <- paste(sf_f,"_species_fitness",sep="")
    sf_model <- try(get(sf_fun),silent = TRUE)
    
    if(class(sf_model) == "try-error"){
      message(paste("species_fitness ERROR: function '",sf_fun,"' could not be retrieved. 
      Make sure it is defined and available in the cxr package or in the global environment.\n"
                    ,sep=""))
      return(NULL)
    }else{
      res <- mapply(sf_model,effect_response_fit$lambda,effect_response_fit$response)
      names(res) <- paste("fitness_",effect_response_fit$sp,sep = "")
    }
  }else{
    sf_fun <- paste(model_family,"_species_fitness")
    sf_model <- try(get(sf_fun),silent = TRUE)
    
    if(class(sf_model) == "try-error"){
      message(paste("species_fitness ERROR: function '",sf_fun,"' could not be retrieved. 
      Make sure it is defined and available in the cxr package or in the global environment.\n"
                    ,sep=""))
      return(NULL)
    }else{
      res <- sf_model(lambda,competitive_response)
    }
  }
  
  res
  
}