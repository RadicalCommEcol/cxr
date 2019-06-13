#' Project abundances according to specified models and parameters
#'
#' @param par list with the following components: 
#' 1 - dataframe sp.par, with the parameters to be passed to the predictive model
#' 2 - dataframe initial.values, with fields "site","species","abundance"
#' 3 - covariates, either 0 if there are no covariates, or a dataframe with fields "site","timestep","covariate","value"
#' 4 - list other.par, other parameters to \code{abundance model}, such as lambda.cov.matrix, alpha.cov.matrix.
#' species should be a unique identifier, character or numeric. idem for site and covariate.
#' @param timesteps number of timesteps to project
#' @param abundance.model a function that accepts parameters from sp.par, a set of initial abundances, and optionally other parameters.
#' The function returns the projected abundances at t+1
#' @param return.seeds boolean flag, whether the prediction should return 
#' number of seeds (i.e. $N_{i,t+1}$ eq. 1 of Lanuza et al. 2018), or number of
#' adult individuals, (i.e. $N_{i,t+1} * g$ )
#'
#' @return dataframe with fields "timestep","site","sp","abundance", giving the expected abundance for each species, timestep, and site.
#' @export
PredictAbundances <- function(par,timesteps,abundance.model, return.seeds = TRUE){
  sites <- unique(par$initial.values$site)
  sp <- unique(par$initial.values$species)
  num.sp <- nrow(par$sp.par)
  
  predicted.abundances <- expand.grid(1:timesteps,sites,sp)
  names(predicted.abundances) <- c("timestep","site","sp")
  predicted.abundances$abundance <- 0
  predicted.abundances$abundance[predicted.abundances$timestep == 1] <- as.numeric(par$initial.values$abundance)
  
  for(i.timestep in 2:timesteps){
    for(i.site in 1:length(sites)){
      init.abund <- predicted.abundances$abundance[predicted.abundances$timestep == (i.timestep-1) & 
                                                     predicted.abundances$site == sites[i.site]]
      if(is.data.frame(par$covariates)){
        cov.values <- par$covariates$value[par$covariates$site == sites[i.site] & 
                                             par$covariates$timestep == i.timestep]
      }else{
        cov.values <- 0
      }
      predicted.abundances$abundance[predicted.abundances$timestep == i.timestep & 
                                       predicted.abundances$site == sites[i.site]] <- abundance.model(sp.par = par$sp.par, 
                                                                                               init.abund = init.abund, 
                                                                                               cov.values = cov.values, 
                                                                                               alpha.matrix = par$other.par$alpha.matrix,
                                                                                               lambda.cov.matrix = par$other.par$lambda.cov,
                                                                                               alpha.cov.matrix = par$other.par$alpha.cov,
                                                                                               return.seeds = return.seeds)
    }# i.site
  }# i.timestep
  return(predicted.abundances)
}