InitParams <- function(init.lambda = NULL,
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
                       init.alpha_NL = NULL,
                       init.lambda.cov_NL = NULL,
                       init.alpha.cov_NL = NULL,
                       lower.alpha_NL = 0,
                       upper.alpha_NL = 1e5,
                       lower.lambda.cov_NL = 1e-4,
                       upper.lambda.cov_NL = 1e5,
                       lower.alpha.cov_NL = 1e-4,
                       upper.alpha.cov_NL = 1e5,
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
  
  # effect of non linearity of covariates on lambda
  if(!is.null(init.lambda.cov_NL)){
    init.par <- c(init.par,init.lambda.cov_NL)
    lower.bounds <- c(lower.bounds,rep(lower.lambda.cov_NL,length(init.lambda.cov_NL)))
    upper.bounds <- c(upper.bounds,rep(upper.lambda.cov_NL,length(init.lambda.cov_NL)))
  }
  
  # alpha value/matrix
  if(!is.null(init.alpha)){
    init.par <- c(init.par,init.alpha)
    lower.bounds <- c(lower.bounds,rep(lower.alpha,length(init.alpha)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha,length(init.alpha)))                  
  }
  
  # alpha value/matrix non linear
  if(!is.null(init.alpha_NL)){
    init.par <- c(init.par,init.alpha_NL)
    lower.bounds <- c(lower.bounds,rep(lower.alpha,length(init.alpha_NL)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha,length(init.alpha_NL)))                  
  }
  
  # effect of covariates on alpha
  if(!is.null(init.alpha.cov)){
    init.par <- c(init.par,init.alpha.cov)
    lower.bounds <- c(lower.bounds,rep(lower.alpha.cov,length(init.alpha.cov)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha.cov,length(init.alpha.cov)))
  }
  
   # effect of non linearity of covariates on alpha
  if(!is.null(init.alpha.cov_NL)){
    init.par <- c(init.par,init.alpha.cov_NL)
    lower.bounds <- c(lower.bounds,rep(lower.alpha.cov_NL,length(init.alpha.cov_NL)))
    upper.bounds <- c(upper.bounds,rep(upper.alpha.cov_NL,length(init.alpha.cov_NL)))
  }
  
  # sigma goes at the end
  init.par <- c(init.par,init.sigma)
  lower.bounds <- c(lower.bounds,lower.sigma)
  upper.bounds <- c(upper.bounds,upper.sigma)
  
  return(list(init.par = init.par, lower.bounds = lower.bounds, upper.bounds = upper.bounds))
}
