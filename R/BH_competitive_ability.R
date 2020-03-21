BH_competitive_ability <- function(lambda, pair_matrix){
  if(all(pair_matrix >= 0)){
    (lambda - 1)/sqrt(pair_matrix[1,1] * pair_matrix[1,2])
  }else{
    NA_real_
  }
}