#' Vital rate calculation
#' 
#' Calculates vital rates from their effect sizes and terms. This
#' is equivalent to predicting from a binomial glm with
#' given coefficients.
#' In this version, the user needs to ensure that `param` and `env` match,
#' i.e. that if the `param` list is defined with environmental forcing, it is
#' passed here, and viceversa. In future versions I may implement checks for that
#' here, but for now, be aware that it will fail.
#' 
#' @param vr integer or char, vital rate to obtain, from the ones defined in `param`.
#' So far, valid names are "Sj","Sn","Sr","Rn","Rr","D","Ds,"O".
#' @param sp integer or char, species
#' @param site intger or char, site
#' @param param param nested list (see `build_param`)
#' @param env optional numeric, environmental forcing
#' @param densities densities of all sp in the site, including individuals from
#' all three life stages
#'
#' @return numeric value
#' @export
vital_rate <- function(vr,sp,site,param,env = NULL,densities){
  
  if(is.null(vr) | is.null(sp) | is.null(site) | is.null(param) | is.null(densities)){
    message(paste("vital_rate ERROR: one or more argumets are not specified",sep=""))
    return(NULL)
  }
  
  rates <- rownames(param[[1]][[1]])
  # rates <- c("Sj","Sn","Sr","Rn","Rr","D","Ds","O")
  # if(is.character(vr)){
    if(vr %in% rates){
      vr.num <- which(rates == vr)
    }else{
      if(vr == "Ds"){
        # this is for a specific situation in which Ds is not included in "param",
        # as it was the case for the first iterations of the model,
        # but we want the model to run anyway, so assume full dispersal survival.
        return(1)
      }else{
        message("vital_rate ERROR: please provide a valid vr param, among Sj,Sn,Sr,Rn,Rr,D,Ds,O.")
        return(NULL)
      }
    }
  # }else if(is.numeric(vr)){
  #   vr.num <- vr
  # }
  
  coefs <- param[[sp]][[site]][vr.num,]
  
  # TODO beware of this, I guess it's ok, but may crash at some point?
  # values <- length(coefs)
  if(!is.null(env)){
    # values <- c(1,env,densities,densities[1]*env)
    # consider all interactions dens:env, not only focal sp
    values <- c(1,env,densities,densities*env)
    # this is a hack for the situation of not having all interactions
    # so that the length of the coefs vector should give the number of 
    # interactions considered. I assume that these are always in the same order
    # first the focal.sp:env, then the second sp:env, etc.
    if(length(values)>length(coefs)){
      values <- values[1:length(coefs)]
    }
  }else{
    values <- c(1,densities)
  }
  
  mean <- exp(sum(coefs*values))
  
  if(vr!= "O"){
    output <- mean/(1+mean)
  }else{ # for reproductive output
    output <- mean
  }
  
  return(output)
}
