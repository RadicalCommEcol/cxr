#' Fill the vec-permutation dispersal matrix
#' 
#' Fill for a given species, all sites
#'
#' @param focal.sp integer, focal species
#' @param num.sites integer, how many sites
#' @param param param nested list,see `build_param` function
#' @param vpm data structure holding all vector-permutation matrices; see 
#' `vec_permutation_matrices`
#' @param env optional numeric, environmental forcing for a given timestep
#' @param current.densities list of length sp, each element is a matrix sites*stages
#'
#' @return dispersal matrix, stages*sites
#' @export
fill_dispersal_matrix <- function(focal.sp,num.sites,param,vpm,env = NULL,current.densities){
  
  if(is.null(focal.sp) | is.null(vpm) | is.null(num.sites) | is.null(current.densities)){
    message(paste("fill_dispersal_matrix ERROR: one or more arguments are not specified",sep=""))
    return(NULL)
  }
  
  temp.disp <- vpm[["dispersal"]][[focal.sp]] 
  
  # if this has dimensions of the M matrix, I create another one:
  temp.mat <- matrix(0,num.sites,num.sites)
  
  # fill up the diagonal of temp.mat, the "staying"
  for(i.site in 1:length(diag(temp.mat))){
    
    # get a vector of densities for this site
    dens <- get_densities(focal.sp,i.site,current.densities)
    
    diag(temp.mat)[i.site] <- 1 - vital_rate(vr = "D",
                                             sp = focal.sp,
                                             site = i.site,
                                             param = param,
                                             env = env,
                                             densities = dens) 
  }
  
  # Assuming temp.disp is the M matrix and the non-reproductives disperse:
  temp.disp[(num.sites+1):(2*num.sites),(num.sites+1):(2*num.sites)] <- temp.mat
  
  # reset temp.mat for using it for the non-diagonal
  diag(temp.mat) <- 0
  
  # moving ( = all the off diagonal)
  # these are the non-reproductive adults that disperse,
  # becoming reproductive adults in the new site
  # since the dispersal matrix is stages*sites, the 
  # non-rep adults are "in the middle" of the matrix,
  # and as they disperse, they go "downwards", 
  # i.e. to the "reproductive adult" part.
  # this loop, then, gives the transition terms
  # from each nr in each site to each r in each other site
  
  # effective dispersal is mediated by "disperser survival", 
  # the vital rate "Ds", so that it is D*Ds.
  
  # i.row is the arriving site
  # i.col is the source
  for(i.row in 1:num.sites){
    for(i.col in 1:num.sites){
      # we don't want the diagonal here,
      # it makes no sense
      if(i.row != i.col){
        
        # get a vector of densities for this site
        dens <- get_densities(focal.sp,i.col,current.densities)
        arriving.dens  <- get_densities(focal.sp,i.row,current.densities)
        
        temp.mat[i.row,i.col] <- (vital_rate(vr = "D",
                                             sp = focal.sp,
                                             site = i.col,
                                             param = param,
                                             env = env,
                                             densities = dens) *
                                    # 1
                                    vital_rate(vr = "Ds",
                                               sp = focal.sp,
                                               site = i.col, # TODO this or i.row? should be i.col, because it refers to the Ds of the sp at the source
                                               param = param,
                                               env = env,
                                               densities = arriving.dens)
        ) # but the densities are those of the receiving site
        
        # TODO ask Maria about /num.sites:
        
        # Also note that I am dividing the dispersal by 
        # how many sites can be dispersed to 
        # to comply with mass conservation
        
      }# not diag
    }# for i.col
  }# for i.row
  
  # Assuming temp.disp is the M matrix and the non-reproductives disperse:
  temp.disp[(2*num.sites+1):(3*num.sites),(num.sites+1):(2*num.sites)] <- temp.mat
  return(temp.disp)
}
