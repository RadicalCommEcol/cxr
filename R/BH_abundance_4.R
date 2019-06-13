#' Project abundance of individuals according to the Beverton-Holt fourth model
#'
#' @param sp.par dataframe with species in rows, and the following columns:
#' lambda: fecundity term
#' germ.rate: seed germination rate
#' survival.rate: annual survival of ungerminated seed
#' @param init.abund number of individuals at time t
#' @param cov.values 1d vector with values of each covariate
#' @param alpha.matrix competition matrix
#' @param lambda.cov.matrix matrix of  num.sp x num.cov
#' giving the effect of each covariate over the fecundity (lambda) of each species
#' @param alpha.cov.matrix list of dimension number of covariates. 
#' Each component of the list is a single value,
#' giving the effect of the covariate in question over the interaction matrix. 
#' @param return.seeds boolean flag, whether the prediction should return 
#' number of seeds (i.e. $N_{i,t+1}$ eq. 1 of Lanuza et al. 2018), or number of
#' adult individuals, (i.e. $N_{i,t+1} * g$ )
#' 
#' @return 1d vector with number of individuals of each species at time t+1
#' @export
BH_abundance_4 <- function(sp.par,init.abund,cov.values,alpha.matrix,lambda.cov.matrix,alpha.cov.matrix,return.seeds = TRUE){
  expected.abund <- rep(0,nrow(sp.par))
  for(i.sp in 1:nrow(sp.par)){
    # numerator
    lambda.cov <- 0
    for(i.cov in 1:ncol(lambda.cov.matrix)){
      lambda.cov <- lambda.cov + cov.values[i.cov]*lambda.cov.matrix[i.sp,i.cov]
    }
    num <- sp.par$lambda[i.sp] * lambda.cov 
    # denominator
    den <- 0
    for(j.sp in 1:nrow(sp.par)){
      alpha.term <- alpha.matrix[i.sp,j.sp]
      for(i.cov in 1:ncol(lambda.cov.matrix)){
        alpha.term <- alpha.term + alpha.cov.matrix[[i.cov]]*cov.values[i.cov]
      }# for i.cov
      den <- den + alpha.term*init.abund[j.sp]
    }# for j.sp
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