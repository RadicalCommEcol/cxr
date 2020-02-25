LW_species_fitness <- function(lambda, competitive_response){
  if(lambda > 1){
    log((lambda-1))/competitive_response
  }else{
    NA
  }
}