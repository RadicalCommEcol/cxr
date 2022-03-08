
#' Generate templates for dispersal, demography, and permutation matrices
#'
#' this follows the vec-permutation approach
#' as defined in: 
#' Hunter and Caswell 2005, doi:10.1016/j.ecolmodel.2005.05.002, 
#' Ozgul et al. 2009, doi: 10.1086/597225
#'
#' @param num.sp integer, number of species
#' @param num.sites integer, number of sites
#' @param num.stages integer, number of stages
#'
#' @return nested list, of the form `list[[type]][[sp]]`, where `type` is 
#' demography, dispersal, or permutation.
#' @export
#'
#' @examples
#' # number of demographic stages - this should be always fixed to 3 for 
#' # compatibility with other functions
#' num.stages <- 3
#' num.sp <- 4
#' num.sites <- 5
#' vpm <- vec_permutation_matrices(num.sp,num.sites,num.stages)
vec_permutation_matrices <- function(num.sp, num.sites, num.stages){
  
  if(is.null(num.sp) | is.null(num.sites) | is.null(num.stages)){
    message(paste("vec_permutation_matrices ERROR: one or more argumets are not specified",sep=""))
    return(NULL)
  }
  
  if(num.sp <= 0 | num.sites <= 0 | num.stages <= 0){
    message(paste("vec_permutation_matrices ERROR: all arguments must be > 0",sep=""))
    return(NULL)
  }
  
  vpm <- list()
  vpm[[1]] <- list()
  vpm[[2]] <- list()
  vpm[[3]] <- list()
  names(vpm) <- c("demography","dispersal","permutation")
  
  dispersal.template <- matrix(0,num.sites*num.stages,num.sites*num.stages)
  diag(dispersal.template) <- 1
  
  demo.template <- matrix(0,num.sites*num.stages,num.sites*num.stages)
  
  perm.template <- matrix(0,num.sites*num.stages,num.sites*num.stages)

  for(i.stage in 1:num.stages){
    for(i.site in 1:num.sites){
      E <- matrix(0,num.stages,num.sites)
      E[i.stage,i.site] <- 1
      perm.template <- perm.template+(E%x%t(E))
    }
  }
  
  for(i.sp in 1:num.sp){
    vpm[["demography"]][[i.sp]] <- demo.template
    vpm[["dispersal"]][[i.sp]] <- dispersal.template
    vpm[["permutation"]][[i.sp]] <- perm.template
  }
  
  return(vpm)
}
