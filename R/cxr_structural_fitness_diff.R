#' structural fitness difference (in degree)
#'
#' @param alpha interaction matrix
#' @param r 1d vector of growth rates
#' @param mf model familiy
#'
#' @return numeric value, structural fitness difference between two species
#' @noRd
cxr_structural_fitness_diff <- function(alpha,r,mf){
  
  #vector defining the centroid of the feasibility domain
  r_centroid <- function(alpha){
    n <- nrow(alpha)
    D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
    alpha_n <- alpha %*% D
    r_c <- rowSums(alpha_n) /n 
    r_c <- t(t(r_c))
    return(r_c)
  }
  
  r_c <- r_centroid(alpha)
  
  # transform r to LV-like growth rates
  if(mf == "BH"){
    my.r <- r - 1
  }else if(mf == "RK"){
    my.r <- log(r)
  }else if(mf == "LW"){
    my.r <- log(r - 1)
  }
  else{
    
    if(mf != "LV"){
      message(paste("structural_fitness: Calculating structural fitness differences for a custom model family,
      using lambdas as growth rates without transforming them. 
      Be aware that this may be incorrect depending on your underlying model.\n"
                    ,sep=""))
    }
    
    my.r <- r
  }
  
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}