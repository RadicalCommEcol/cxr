# Function that takes the par of model BH_5 and return the a list of TRUE if we keep the coefficient and false if we set it to 0, 
# depending on the criterion.
Choose_coefficients <-function(par, focal.comp.matrix, num.covariates, num.competitors, focal.covariates, criterion){

# First, we set a numeric vector, lambda.cov.comparable, containig the coeffictients multiplied by the mean of the covariate
  lambda.cov.comparable <- vector("numeric",num.covariates)
  for(v in 1:num.covariates){
  lambda.cov.comparable [v] <- lambda.cov[v]*mean(focal.cov.matrix[,v])
    }

# We get a booleen vector containig TRUE if the coefficient is kept and FALSE if it is set to 0
  lambda.cov.comparison  <- lambda.cov.comparable > criterion*max(lambda.cov.comparable)
  
  
## We compare the alpha covariates 
# We obtain a vector containing the coefficients of the alpha covariates 
#We weigh it with the effect mean 
  alpha_coeff.comparable <- vector("numeric",length(alpha_coeff))
  
  # For the alphas, we multiply by the mean of the number of competitors :
  for (v in 1:length(num.competitors)){
  alpha_coeff.comparable[v] <- alpha_coeff[v]*mean(focal.comp.matrix[,v])
    }
 # For the other covariates we multiply by the mean of the covariate times the number of competitors
  focal.cov.matrix <- as.matrix(focal.covariates)
   for (u in 1:num.covariates){
    for (v in (num.competitors*u+1) : (num.competitors*(u+1))){
	
      alpha_coeff.comparable[v] <- alpha_coeff[v]*mean(focal.comp.matrix[,v-num.competitors*(u)]*focal.cov.matrix[,u])
      }
    }
 # We get a booleen vector containig TRUE if the coefficient is kept and FALSE if it is set to 0
  alpha.coeff.comparison  <- alpha_coeff.comparable > (criterion*max(alpha_coeff.comparable))
	
new_par<-function(par, focal.comp.matrix, num.covariates, num.competitors, focal.covariates, criterion){
