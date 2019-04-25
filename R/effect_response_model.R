
# function to optimize lambda, r, and e parameters (competition responses/effects)
# see Godoy et al. 2014.
# here, all sp have to be parameterized in the same function, because e
# is not pair-specific.
# Therefore, the vector of parameters includes all lambdas, all r, all e, and the sigma term

# species are sorted alphabetically

effect_response_model <- function(par, sp.data){
  
  # sort data just in case
  sp.data <- dplyr::arrange(sp.data,site,focal,competitor)
  
  # set of focal sp
  focal.sp <- sort(unique(sp.data$focal))
  num.focal <- length(focal.sp)
  
  # set of competitor sp need not be = set of focal sp
  comp.sp <- unique(sp.data$competitor)
  num.comp <- length(comp.sp)
  
  sites <- unique(sp.data$site)
  
  log.fitness <- log(sp.data$fitness)
  
  lambda.vector <- par[1:num.focal]
  r.vector <- par[(num.focal+1):(num.focal+num.focal)]
  e.vector <- par[(num.focal+1+num.focal):(length(par)-1)]
  sigma <- par[length(par)]
  
  pred <- NULL
  
  for(i.site in 1:length(sites)){
    
    pred.site <- rep(0,num.focal)
    comp.abundances <- matrix(0,nrow = num.focal,ncol = num.comp)
    
    my.site.data <- subset(sp.data, site == sites[i.site])
    site.focal <- unique(my.site.data$focal)
    # fill temporary competition matrix, 
    # giving how many competitors for each focal sp
    # in this specific site
    # with zero entries if no focal/no competitors
    for(i.focal in 1:length(site.focal)){
      focal.index <- which(focal.sp == site.focal[i.focal])
      comp.abundances[focal.index,] <- my.site.data$number[my.site.data$focal == site.focal[i.focal]]
    }
    
    for(i.sp in 1:num.focal){ 
      
      num <- lambda.vector[i.sp]
      
      term <- 0
      for(j.sp in 1:num.comp){
        term <- term + e.vector[j.sp]*comp.abundances[i.sp,j.sp]
      }
      den <- 1+(r.vector[i.sp]*term)
      pred.site[i.sp] <- num/den
    }# for i.sp
    
    pred <- c(pred,pred.site)
  }# for i.site
  
  llik <- dnorm(log.fitness, mean = as.numeric(log(pred)), sd = (sigma), log = TRUE)
  # return sum of negative log likelihoods
  return(sum(-1*llik)) 
}


