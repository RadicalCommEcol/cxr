#' Fill a transition matrix
#' 
#' Calculates the elements of a site-specific transition matrix for a given sp.
#' Note that here, and through all functions, we fix three life stages.
#' Also note that `param` and `env` must match, as for the `vital_rate` function.
#'
#' @param focal.sp integer, species
#' @param site integer, site
#' @param param param structure (see `build_param` function)
#' @param env optional numeric, environmental forcing for a given timestep
#' @param current.densities list of length sp, each element is a matrix site*stages
#'
#' @return 3x3 transition matrix
#' @export
fill_transition_matrix <- function(focal.sp,site,param,env = NULL,current.densities){
  
  if(is.null(focal.sp) | is.null(site) | is.null(param) | is.null(current.densities)){
    message(paste("fill_transition_matrix ERROR: one or more argumets are not specified",sep=""))
    return(NULL)
  }
  
  # get a vector of densities for this site
  dens <- get_densities(focal.sp,site,current.densities)
  
  # each transition matrix is n.stages x n.stages
  # for now, assume three life stages
  tm <- matrix(NA,3,3)
  
  # fill the different transitions between stages
  tm[1,1] <- 0
  tm[2,1] <- vital_rate(vr = 1,
                        sp = focal.sp,
                        site = site,
                        param = param,
                        env = env,
                        densities = dens)
  tm[3,1] <- 0
  
  tm[1,2] <- 0
  tm[2,2] <- vital_rate(vr = 2,
                        sp = focal.sp, 
                        site = site, 
                        param = param,
                        env = env,
                        densities = dens) * (1 - vital_rate(vr = 4,
                                                            sp = focal.sp,
                                                            site = site,
                                                            param = param,
                                                            env = env,
                                                            densities = dens))
  tm[3,2] <- vital_rate(vr = 2,
                        sp = focal.sp, 
                        site = site, 
                        param = param,
                        env = env,
                        densities = dens) * vital_rate(vr = 4,
                                                       sp = focal.sp,
                                                       site = site,
                                                       param = param,
                                                       env = env,
                                                       densities = dens)
  
  tm[1,3] <- vital_rate(vr = 3,
                        sp = focal.sp, 
                        site = site, 
                        param = param,
                        env = env,
                        densities = dens) * vital_rate(vr = 7,
                                                       sp = focal.sp,
                                                       site = site,
                                                       param = param,
                                                       env = env,
                                                       densities = dens)
  tm[2,3] <- vital_rate(vr = 3,
                        sp = focal.sp, 
                        site = site, 
                        param = param,
                        env = env,
                        densities = dens) * (1 - vital_rate(vr = 5,
                                                            sp = focal.sp,
                                                            site = site,
                                                            param = param,
                                                            env = env,
                                                            densities = dens))
  tm[3,3] <- vital_rate(vr = 3,
                        sp = focal.sp, 
                        site = site, 
                        param = param,
                        env = env,
                        densities = dens) * vital_rate(vr = 5,
                                                       sp = focal.sp,
                                                       site = site,
                                                       param = param,
                                                       env = env,
                                                       densities = dens)
  
  return(tm)
  
}