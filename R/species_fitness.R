#' Fitness of a species
#' 
#' Calculates the fitness of a species sensu Godoy et al. (2014).Note that its definition is model-specific, 
#' i.e. it depends on the model family from which interaction coefficients were estimated.
#' The function given here assumes a community of n-species, so that species fitness is calculated according to a 
#' general competitive response (r) substituting the 2-sp denominator terms of table A1 of Hart et al. 2018. 
#' This competitive response can be calculated for a series of species with the function 'cxr_er_fit'.
#' 
#' Thus, the function accepts two sets of parameters. First, a 'cxr_er_fit' object returned from that function. 
#' In this case, species fitness will be calculated for all focal taxa included in the 'cxr_er_fit' object.
#' 
#' Otherwise, users may enter a specification of the model to use, 
#' as well as lambda and competitive response parameters of a single species.
#' 
#' If no model family is provided, or a model family for which there is no associated 'XX_species_fitness'
#' function, the function resorts to the standard Lotka-Volterra formulation (Hart et al. 2018). 
#' Overall, we strongly suggest that you use the standard formulation ONLY if you are completely confident 
#' that the model from which you obtained your parameters is consistent with it. 
#' Otherwise, you should include your own formulation of species fitness (see vignette 4).
#' 
#' @param effect_response_fit cxr_er_fit object with valid lambda and response terms.
#' @param lambda per capita fecundity of the species in the absence of competition.
#' @param competitive_response parameter reflecting the species' sensitivity to competition.
#' @param model_family model family for which to calculate species fitness.
#'
#' @return single numeric value/vector, species fitness of one or several taxa
#' @export
#'
species_fitness <- function(effect_response_fit = NULL, 
                            lambda = NULL, 
                            competitive_response = NULL,
                            model_family = NULL){
  res <- NULL
  if(!is.null(effect_response_fit)){
    
    if(!is.null(model_family) | !is.null(lambda) | !is.null(competitive_response)){
      message("cxr species_fitness: the 'cxr_er_fit' object will be used, other
              arguments will be discarded.")
    }
    
    sf_f <- substr(effect_response_fit$model_name,1,2)
    
    sf_fun <- paste(sf_f,"_species_fitness",sep="")
    sf_model <- try(get(sf_fun),silent = TRUE)
    
    if(class(sf_model) == "try-error"){
      message(paste("cxr species_fitness ERROR: function '",sf_fun,"' could not be retrieved. 
      Make sure it is defined and available in the cxr package or in the global environment.\n"
                    ,sep=""))
      return(NULL)
    }else{
      res <- mapply(sf_model,effect_response_fit$lambda,effect_response_fit$response)
      names(res) <- paste("fitness_",effect_response_fit$taxa,sep = "")
    }
  }else{
    sf_fun <- paste(model_family,"_species_fitness")
    sf_model <- try(get(sf_fun),silent = TRUE)
    
    if(class(sf_model) == "try-error"){
      
      sf_model <- BH_species_fitness
      
      message(paste("cxr species_fitness: function '",sf_fun,"' could not be retrieved. 
      The default formulation for Lotka-Volterra models (lambda - 1 in numerator term) will be used. 
      Be aware that this may yield incorrect results for your model family.\n"
                    ,sep=""))
    }
     res <- sf_model(lambda,competitive_response)
  }# if-else parameters
  
  res
  
}