RK_competitive_ability <- function(lambda, pair_matrix){
  if(lambda > 0 & all(pair_matrix >= 0)){
    log(lambda)/sqrt(pair_matrix[1,1] * pair_matrix[1,2])
  }else{
    NA_real_
  }
}