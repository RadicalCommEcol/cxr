
InitParams_2 <- function(init.lambda = NULL,
                         init.sigma = 0,
                         init.alpha = NULL,
                         init.lambda.cov = NULL,
                         init.alpha.cov = NULL,
                       lower.lambda = 1,
                       upper.lambda = 1e5,
                       lower.sigma = 1e-10,
                       upper.sigma = 1e5,
                       lower.alpha = 0,
                       upper.alpha = 1e5,
                       lower.lambda.cov = 1e-4,
                       upper.lambda.cov = 1e5,
                       lower.alpha.cov = 1e-4,
                       upper.alpha.cov = 1e5,
                       num.competitors,
                       num.covariates
                       ){
  init.par <- NULL
  lower.bounds <- NULL
  upper.bounds <- NULL
  
  # if lambda is not null, it goes first
  if(!is.null(init.lambda)){
    init.par <- init.lambda
    lower.bounds <- rep(lower.lambda,length(init.lambda))
    upper.bounds <- rep(upper.lambda,length(init.lambda))
  }
  
  # effect of covariates on lambda
  if(!is.null(init.lambda.cov)){
    init.par <- c(init.par,init.lambda.cov)
    lower.bounds <- c(lower.bounds,rep(lower.lambda.cov,length(init.lambda.cov)))
    upper.bounds <- c(upper.bounds,rep(upper.lambda.cov,length(init.lambda.cov)))
  }
  
  # alpha value/matrix
  if(!is.null(init.alpha)){
    init.par <- c(init.par,init.alpha)
    lower.bounds <- c(lower.bounds,rep(lower.alpha,length(init.alpha)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha,length(init.alpha)))                  
  }
  
  # effect of covariates on alpha
  if(!is.null(init.alpha.cov)){
    init.par <- c(init.par,init.alpha.cov)
    lower.bounds <- c(lower.bounds,rep(lower.alpha.cov,length(init.alpha.cov)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha.cov,length(init.alpha.cov)))
  }
  
  # sigma goes at the end
  init.par <- c(init.par,init.sigma)
  lower.bounds <- c(lower.bounds,lower.sigma)
  upper.bounds <- c(upper.bounds,upper.sigma)
  
  return(list(init.par = init.par, lower.bounds = lower.bounds, upper.bounds = upper.bounds))
  
  # if(model.index == 1){
  #   # fitness.model <- model1
  #   
  #   init.lambda <- mean(log.fitness)
  #   init.sigma <- sd(log.fitness)
  #   
  #   init.par <- c(init.lambda, init.sigma)
  #   lower.bounds <- c(lower.lambda, lower.sigma) 
  #   upper.bounds <- c(upper.lambda, upper.sigma) 
  #   
  # }else if(model.index == 2){
  #   # fitness.model <- model2
  #   
  #   init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
  #                                          lambda.results$model == 1 & 
  #                                          lambda.results$optim.method == init.par.method]
  #   
  #   init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
  #                                        lambda.results$model == 1 & 
  #                                        lambda.results$optim.method == init.par.method]
  #   if(include.lambda){
  #     init.par <- c(init.lambda,init.alpha,init.sigma)
  #     lower.bounds <- c(lower.lambda, lower.alpha, lower.sigma)
  #     upper.bounds <- c(upper.lambda, upper.alpha, upper.sigma)
  #   }else{
  #     init.par <- c(init.alpha,init.sigma)
  #     lower.bounds <- c(lower.alpha, lower.sigma)
  #     upper.bounds <- c(upper.alpha, upper.sigma)
  #   }
  #   
  # }else if(model.index == 3){
  #   # fitness.model <- model3
  #   
  #   init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
  #                                          lambda.results$model == 2 & 
  #                                          lambda.results$optim.method == init.par.method]
  #   init.alphas <- rep(param.matrices[[focal.sp]][[2]][[init.method.num]]$alpha.matrix, times=num.competitors)
  #   init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
  #                                        lambda.results$model == 2 & 
  #                                        lambda.results$optim.method == init.par.method]
  #   
  #   if(include.lambda){
  #     init.par <- c(init.lambda,init.alphas,init.sigma)
  #     lower.bounds <- c(lower.lambda, rep(lower.alpha, times=num.competitors),lower.sigma)
  #     upper.bounds <- c(upper.lambda, rep(upper.alpha, times=num.competitors),upper.sigma)
  #   }else{
  #     init.par <- c(init.alphas,init.sigma)
  #     lower.bounds <- c(rep(lower.alpha, times=num.competitors),lower.sigma)
  #     upper.bounds <- c(rep(upper.alpha, times=num.competitors),upper.sigma)
  #   }    
  # }else if(model.index == 4){
  #   # fitness.model <- model4
  #   
  #   init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
  #                                          lambda.results$model == 3 & 
  #                                          lambda.results$optim.method == init.par.method]
  #   init.alphas <- param.matrices[[focal.sp]][[3]][[init.method.num]]$alpha.matrix
  #   init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
  #                                        lambda.results$model == 3 & 
  #                                        lambda.results$optim.method == init.par.method]
  #   
  #   init.l_cov <- rep(init.lambda.cov, times = num.covariates)
  #   init.a_cov <- rep(init.alpha.cov, times = num.covariates)
  #   
  #   if(include.lambda){
  #     init.par <- c(init.lambda,init.l_cov,init.alphas,init.a_cov,init.sigma)
  #     
  #     lower.bounds <- c(lower.lambda, 
  #                       rep(lower.lambda.cov, times=ncol(covariates)), 
  #                       rep(lower.alpha, times=num.competitors), 
  #                       rep(lower.alpha.cov, times=ncol(covariates)), 
  #                       lower.sigma)
  #     upper.bounds <- c(upper.lambda, 
  #                       rep(upper.lambda.cov, times=ncol(covariates)), 
  #                       rep(upper.alpha, times=num.competitors), 
  #                       rep(upper.alpha.cov, times=ncol(covariates)), 
  #                       upper.sigma)
  #   }else{
  #     init.par <- c(init.l_cov,init.alphas,init.a_cov,init.sigma)
  #     
  #     lower.bounds <- c(#lower.lambda, 
  #       rep(lower.lambda.cov, times=ncol(covariates)), 
  #       rep(lower.alpha, times=num.competitors), 
  #       rep(lower.alpha.cov, times=ncol(covariates)), 
  #       lower.sigma)
  #     upper.bounds <- c(#upper.lambda, 
  #       rep(upper.lambda.cov, times=ncol(covariates)), 
  #       rep(upper.alpha, times=num.competitors), 
  #       rep(upper.alpha.cov, times=ncol(covariates)), 
  #       upper.sigma)
  #   }
  # }else if(model.index == 5){
  #   # fitness.model <- model5
  #   
  #   init.lambda <- lambda.results$lambda[lambda.results$focal.sp == focal.sp & 
  #                                          lambda.results$model == 4 & 
  #                                          lambda.results$optim.method == init.par.method]
  #   init.alphas <- param.matrices[[focal.sp]][[3]][[init.method.num]]$alpha.matrix
  #   init.sigma <- lambda.results$sigma[lambda.results$focal.sp == focal.sp & 
  #                                        lambda.results$model == 4 & 
  #                                        lambda.results$optim.method == init.par.method]
  #   
  #   init.l_cov <- param.matrices[[focal.sp]][[4]][[init.method.num]]$lambda.cov.matrix
  #   a_cov4 <- param.matrices[[focal.sp]][[4]][[init.method.num]]$alpha.cov.matrix
  #   
  #   # rep(each)?
  #   init.a_cov <- c()
  #   for(w in 1:length(a_cov4)){
  #     init.a_cov <- c(init.a_cov, rep(a_cov4[w], times = num.competitors))
  #   }
  #   
  #   if(include.lambda){
  #     init.par <- c(init.lambda,init.l_cov,init.alphas,init.a_cov,init.sigma)
  #     lower.bounds <- c(lower.lambda, 
  #                       rep(lower.lambda.cov, times=ncol(covariates)), 
  #                       rep(lower.alpha, times=num.competitors), 
  #                       rep(lower.alpha.cov, times=(ncol(covariates)*num.competitors)), 
  #                       lower.sigma)
  #     upper.bounds <- c(upper.lambda, 
  #                       rep(upper.lambda.cov, times=ncol(covariates)), 
  #                       rep(upper.alpha, times=num.competitors), 
  #                       rep(upper.alpha.cov, times=(ncol(covariates)*num.competitors)), 
  #                       upper.sigma)
  #   }else{
  #     init.par <- c(init.l_cov,init.alphas,init.a_cov,init.sigma)
  #     lower.bounds <- c(#lower.lambda, 
  #       rep(lower.lambda.cov, times=ncol(covariates)), 
  #       rep(lower.alpha, times=num.competitors), 
  #       rep(lower.alpha.cov, times=(ncol(covariates)*num.competitors)), 
  #       lower.sigma)
  #     upper.bounds <- c(#upper.lambda, 
  #       rep(upper.lambda.cov, times=ncol(covariates)), 
  #       rep(upper.alpha, times=num.competitors), 
  #       rep(upper.alpha.cov, times=(ncol(covariates)*num.competitors)), 
  #       upper.sigma)
  #   }
  #   
  # }# if-else model
  # return(list(init.par = init.par, lower.bounds = lower.bounds, upper.bounds = upper.bounds))
}



