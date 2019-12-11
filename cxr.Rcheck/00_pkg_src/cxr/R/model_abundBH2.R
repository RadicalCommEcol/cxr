#' Project abundance of individuals according to the Beverton-Holt second model
#'
#' @param sp.par dataframe with species in rows, and the following columns:
#' lambda: fecundity term
#' germ.rate: seed germination rate
#' survival.rate: annual survival of ungerminated seed
#' @param init.abund number of individuals at time t
#' @param cov.values Not used in model_abundBH2
#' @param alpha.matrix competition value, same for all interactions
#' @param lambda.cov.matrix Not used in model_abundBH2
#' @param alpha.cov.matrix Not used in model_abundBH2
#' @param return.seeds boolean flag, whether the prediction should return 
#' number of seeds (i.e. \eqn{N_{i,t+1}}, eq. 1 of Lanuza et al. 2018), or number of
#' adult individuals, (i.e. \eqn{N_{i,t+1} * g} )
#'
#' @return 1d vector with number of individuals of each species at time t+1
#' @export
model_abundBH2 <- function(sp.par,init.abund,cov.values,alpha.matrix,lambda.cov.matrix,alpha.cov.matrix,return.seeds = TRUE){
  expected.abund <- rep(0,nrow(sp.par))
  for(i.sp in 1:nrow(sp.par)){
    # numerator
    num <- sp.par$lambda[i.sp]
    # denominator
    den <- 0
    
    if(return.seeds){
      # if init.abund are seeds
      for(j.sp in 1:nrow(sp.par)){
        den <- den + alpha.matrix*(init.abund[j.sp]*sp.par$germ.rate[j.sp])
      }
    }else{
      # if init.abund are adult individuals
      for(j.sp in 1:nrow(sp.par)){
        den <- den + alpha.matrix*init.abund[j.sp]
      }
    }# if-else return.seeds
    
    den <- 1+den
    # overall fitness metric
    fitness <- num/den
    expected.abund[i.sp] <- ((1-sp.par$germ.rate[i.sp])*sp.par$survival.rate[i.sp]) + sp.par$germ.rate[i.sp]*fitness
  }
  # just in case
  expected.abund[expected.abund < 0] <- 0
  # the above formulation gives number of seeds. If necessary, multiply by the germination rate to get number of individuals
  if(return.seeds){
    return(expected.abund)
  }else{
    return(expected.abund*sp.par$germ.rate)
  }
}

