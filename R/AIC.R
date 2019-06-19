# METHOD 1 (Burnham and Anderson, 2002)
#AIC=-2log(max likl)+2k #k=number of parameters
# METHOD 2 
#AIC=2log(max lik)-2np #np=number of parameters for a given model
AIC <- function(negative_llik, num.covariates, num.competitors,model.number,param.list){
  if(model.number==1){
    num.parameters <- 1 #lambda
    return(2*negative_llik + 2*num.parameters)
  }
  
  else if(model.number==2){
    num.parameters <- 2 #lambda and alpha
    return(2*negative_llik + 2*num.parameters)
  }
  
  else if(model.number==3){
    num.parameters <- 1 + num.competitors #lambda and one alpha per competitor
    return(2*negative_llik + 2*num.parameters)
  }
  
  else if(model.number==4){
    num.parameters <- 1 + num.competitors + 2*num.covariates #lambda, lambda.cov, one alpha per competitor, one alpha.cov per covariate
    if("lambda.cov_NL" %in% param.list){num.parameters<-num.parameters + num.covariates}
    if("alpha.cov_NL" %in% param.list){num.parameters<-num.parameters + num.covariates}
    return(2*negative_llik + 2*num.parameters)
  }
  
  else if(model.number==5){
    num.parameters <- 1 + num.competitors + num.covariates*(num.competitors+1) #lambda, lambda.cov, one alpha per competitor, one alpha.cov per covariate and per competirors
    if("lambda.cov_NL" %in% param.list){num.parameters <- num.parameters + num.covariates*num.competitors}
    if("alpha.cov_NL" %in% param.list){num.parameters <- num.parameters + num.covariates*num.competitors}
    return(2*negative_llik + 2*num.parameters)
  }
}