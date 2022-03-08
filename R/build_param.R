#' Build param structure
#' 
#' Builds a nested list for the parameters of a given metapopulation
#'
#' @param sp character vector with species names
#' @param sites character vector with site names
#' @param rates character vector, vital rate names
#' @param env boolean, whether environment is accounted for
#' @param num.params optional, integer giving the number of parameters to account for.
#' If not specified, it will include environment interactions with all species densities.
#' E.g. if 3 sp and env = TRUE, there will be 7 params (intercept + 6 betas)
#'
#' @return nested list of the form `list[[sp]][[site]]`. Each of these elements
#' is a NA matrix with vital rates in rows and expected parameters in columns.
#' @export
#'
#' @examples
#' sp <- c("s1","s2","s3")
#' sites <- c("sa","sb")
#' rates <- c("Sj","Sn","Sr","Rn","Rr","D","O")
#' env <- TRUE
#' param <- build_param(sp = sp,sites = sites,rates = rates,env = env)
build_param <- function(sp,sites,rates,env,num.params = NULL){
  
  if(is.null(sp) | is.null(sites) | is.null(rates) | is.null(env)){
    message(paste("build_param ERROR: one or more argumets are not specified",sep=""))
    return(NULL)
  }
  
  ########################
  # auxiliary function
  # number of params depend on the number of species
  # and on whether there is environmental forcing
  get_nparams <- function(nsp,env = FALSE){
    num <- 1 + nsp
    if(env){
      num <- num + 1 + nsp # all sp interactions with env
      # num <- num + 2 # only interaction between focal sp and env
    }
    return(num)
  }
  #######################
  
  nparams <- ifelse(is.null(num.params),get_nparams(length(sp),env),num.params)
  nameparams <- c("alpha",paste("beta",1:(nparams-1),sep=""))
  
  testd <- list(length(sp))
  
  for(i.sp in 1:length(sp)){
    testd[[i.sp]] <- list(length(sites))
    
    for(i.site in 1:length(sites)){
      testd[[i.sp]][[i.site]] <- matrix(NA_real_,ncol = nparams,
                                        nrow = length(rates),
                                        dimnames = list(rates,nameparams))
    }# for i.site
    names(testd[[i.sp]]) <- sites
  }# for i.sp
  names(testd) <- sp
  return(testd)
}
