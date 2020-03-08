
BH_project_alpha_none_lambdacov_none_alphacov_none <- function(lambda,
                                                               alpha,
                                                               lambda_cov,
                                                               alpha_cov,
                                                               abundances,
                                                               covariates){
  
  numsp <- length(lambda)
  expected_abund <- rep(0,numsp)
  
  for(i.sp in 1:numsp){
    expected_abund[i.sp] <- lambda[i.sp] * abundances[i.sp]
  }# for each sp
  
  expected_abund
}

