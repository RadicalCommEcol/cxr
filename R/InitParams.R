#' Title gather initial parameters for the optimization procedures
#' 
#' the list of initial parameters varies depending on the model to optimize (1-5, see nested_models.R). This function
#' does the dirty work of sorting the parameters appropriately for feeding the optim functions.
#' Basically it needs to sort the parameters in a 1-dim vector, with the following order:
#' lambda-lambda.covariates-alphas-alphas.covariates-sigma.
#' Depending on the model, some of these will not be necessary, or have different lengths.
#' This function returns this 1d vector of initial parameters alongside the lower and upper bounds for the optimization
#' procedures.
#'
#' @param model.index the number of the model to use
#' @param log.fitness in case model = 1, it uses as a starting value the mean of this vector
#' @param lambda.results in case model > 1, take initial lambdas and sigmas from the previous model
#' @param num.competitors for initializing the alpha matrices
#' @param num.covariates idem
#' @param init.par.method in case it needs to gather initial values from previous results and there are different optimization methods
#' calculated, from which should it take the starting values?
#' @param focal.sp which focal species
#' @param param.matrices in case model > 2, it gathers initial estimates from previous values of alpha, etc.
#' @param lower.lambda lower bounds
#' @param upper.lambda upper bounds
#' @param lower.sigma lower bounds
#' @param upper.sigma upper bounds
#' @param init.alpha initial estimates if no previous model
#' @param lower.alpha lower bounds
#' @param upper.alpha upper bounds
#' @param init.lambda.cov initial estimates if no previous model
#' @param lower.lambda.cov lower bounds
#' @param upper.lambda.cov upper bounds
#' @param init.alpha.cov initial estimates if no previous model
#' @param lower.alpha.cov lower bounds
#' @param upper.alpha.cov upper bounds
#'
#' @return a list with three components: init.par with initial values, lower.bounds, and upper.bounds
#' @export
#'
#' @examples
InitParams <- function(model.index,
                       log.fitness,
                       lambda.results,
                       num.competitors,
                       num.covariates,
                       init.par.method,
                       focal.sp,
                       param.matrices,
                       lower.lambda = 1,
                       upper.lambda = 1e5,
                       lower.sigma = 1e-10,
                       upper.sigma = 1e5,
                       init.alpha = 1e-4,
                       lower.alpha = 0,
                       upper.alpha = 1e5,
                       init.lambda.cov = 1e-3,
                       lower.lambda.cov = 1e-4,
                       upper.lambda.cov = 1e5,
                       init.alpha.cov = 1e-3,
                       lower.alpha.cov = 1e-4,
                       upper.alpha.cov = 1e5){
  
  
  if(model.index == 1){
    # fitness.model <- model1
    
    init.lambda <- mean(log.fitness)
    init.sigma <- sd(log.fitness)
    
    init.par <- c(init.lambda, init.sigma)
    lower.bounds <- c(lower.lambda, lower.sigma) 
    upper.bounds <- c(upper.lambda, upper.sigma) 
    
  }else if(model.index == 2){
    # fitness.model <- model2
    
    init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
                                           lambda.results$model == 1 & 
                                           lambda.results$optim.method == init.par.method]
    
    init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
                                         lambda.results$model == 1 & 
                                         lambda.results$optim.method == init.par.method]
    
    init.par <- c(init.lambda,init.alpha,init.sigma)
    lower.bounds <- c(lower.lambda, lower.alpha, lower.sigma)
    upper.bounds <- c(upper.lambda, upper.alpha, upper.sigma)
    
  }else if(model.index == 3){
    # fitness.model <- model3
    
    init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
                                           lambda.results$model == 2 & 
                                           lambda.results$optim.method == init.par.method]
    init.alphas <- rep(param.matrices[[focal.sp]][[2]][[init.method.num]]$alpha.matrix, times=num.competitors)
    init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
                                         lambda.results$model == 2 & 
                                         lambda.results$optim.method == init.par.method]
    
    init.par <- c(init.lambda,init.alphas,init.sigma)
    
    lower.bounds <- c(lower.lambda, rep(lower.alpha, times=num.competitors),lower.sigma)
    upper.bounds <- c(upper.lambda, rep(upper.alpha, times=num.competitors),upper.sigma)
    
  }else if(model.index == 4){
    # fitness.model <- model4
    
    init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
                                           lambda.results$model == 3 & 
                                           lambda.results$optim.method == init.par.method]
    init.alphas <- param.matrices[[focal.sp]][[3]][[init.method.num]]$alpha.matrix
    init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
                                         lambda.results$model == 3 & 
                                         lambda.results$optim.method == init.par.method]
    
    init.l_cov <- rep(init.lambda.cov, times = num.covariates)
    init.a_cov <- rep(init.alpha.cov, times = num.covariates)
    
    init.par <- c(init.lambda,init.l_cov,init.alphas,init.a_cov,init.sigma)
    
    lower.bounds <- c(lower.lambda, 
                      rep(lower.lambda.cov, times=ncol(covariates)), 
                      rep(lower.alpha, times=num.competitors), 
                      rep(lower.alpha.cov, times=ncol(covariates)), 
                      lower.sigma)
    upper.bounds <- c(upper.lambda, 
                      rep(upper.lambda.cov, times=ncol(covariates)), 
                      rep(upper.alpha, times=num.competitors), 
                      rep(upper.alpha.cov, times=ncol(covariates)), 
                      upper.sigma)
    
  }else if(model.index == 5){
    # fitness.model <- model5
    
    init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
                                           lambda.results$model == 4 & 
                                           lambda.results$optim.method == init.par.method]
    init.alphas <- param.matrices[[focal.sp]][[3]][[init.method.num]]$alpha.matrix
    init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
                                         lambda.results$model == 4 & 
                                         lambda.results$optim.method == init.par.method]
    
    init.l_cov <- param.matrices[[focal.sp]][[4]][[init.method.num]]$lambda.cov.matrix
    a_cov4 <- param.matrices[[focal.sp]][[4]][[init.method.num]]$alpha.cov.matrix
    
    # rep(each)?
    init.a_cov <- c()
    for(w in 1:length(a_cov4)){
      init.a_cov <- c(init.a_cov, rep(a_cov4[w], times = num.competitors))
    }
    
    init.par <- c(init.lambda,init.l_cov,init.alphas,init.a_cov,init.sigma)
    lower.bounds <- c(lower.lambda, 
                      rep(lower.lambda.cov, times=ncol(covariates)), 
                      rep(lower.alpha, times=num.competitors), 
                      rep(lower.alpha.cov, times=(ncol(covariates)*num.competitors)), 
                      lower.sigma)
    upper.bounds <- c(upper.lambda, 
                      rep(upper.lambda.cov, times=ncol(covariates)), 
                      rep(upper.alpha, times=num.competitors), 
                      rep(upper.alpha.cov, times=(ncol(covariates)*num.competitors)), 
                      upper.sigma)
    
  }# if-else model
  return(list(init.par = init.par, lower.bounds = lower.bounds, upper.bounds = upper.bounds))
}



