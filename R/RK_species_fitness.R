RK_species_fitness <- function(lambda, competitive_response){
  if(lambda > 0){
    log(lambda)/competitive_response
  }else{
    NA
  }
}