#' Fitness ratio among two or more species
#'
#' @param effect_response_fit cxr_er_fit object
#' @param fitness_sp1 numeric value representing the fitness (a.k.a. competitive ability) of the first taxa
#' @param fitness_sp2 numeric value representing the fitness (a.k.a. competitive ability) of the second taxa
#'
#' @return either a matrix with fitness ratios for all pairs of fitted species, or a single numeric value.
#' The matrix elements represent the ratios of species in columns over species in rows, and conversely, 
#' the numeric value represents the ratio of sp1 over sp2.
#' 
#' @export
#'
#' @examples fitness_ratio(fitness_sp1 = 0.6, fitness_sp2 = 0.3)
fitness_ratio <- function(effect_response_fit = NULL,fitness_sp1 = NULL,fitness_sp2 = NULL){
  res <- NULL
  if(!is.null(effect_response_fit)){
    
    spfitness <- try(cxr::species_fitness(effect_response_fit = effect_response_fit),silent = TRUE)
    
    if(!is.null(spfitness)){
      spnames <- effect_response_fit$taxa
      res <- matrix(nrow = length(spnames),ncol = length(spnames),dimnames = list(spnames,spnames))
      for(i.sp in 1:nrow(res)){
        for(j.sp in 1:ncol(res)){
          res[i.sp,j.sp] <- spfitness[j.sp]/spfitness[i.sp]
        }
      }

    }else{
      warning("fitness_ratio: species fitness could not be calculated from the 'cxr_er_fit' object supplied.")
      return(NULL)
    }
    
  }else{
    if(!is.null(fitness_sp1) & !is.null(fitness_sp2)){
      res <- fitness_sp1/fitness_sp2
    }
  }#if-else
  res
}