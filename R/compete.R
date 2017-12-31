#' Calculates pairwise competitive effects between species
#'
#' Calculates how reproduction decays with number of intra and inter-specific
#'  compettors. The idea is to fit a series of nested models. First model fits
#'  an intercept. Second model fits a common response to all competitors, 
#'  Third model fits a pairwise effct of competition. Parameters obtained in each 
#'  model serves as priors for the next model. 
#' 
#' @param focal vector with the focal species name of each observation
#' @param reprod Vector of the reproductive success of each observation
#' @param comp_matrix data.frame with as many columns as competitor species. Each cell is number of competitors per observtion
#' @param log Should reproductive success be loed. Forced to be TRUE
#' @param op number of iterations to ptimize the models.
#'
#' @note This script follows previous methodology developed by Godoy et al 2014,
#'  and Kraft et al. 2015 papers
#'  
#' @return a data frame?
#' @examples
#'  

compete <- function(focal, reprod, comp_matrix, treatment, #need to implement model 4
                    drop0 = FALSE, log = TRUE, op = 25){
  #check focal is c() and match reprod and nrow(comp_matrix)
  #check names(focal) match focal colnames(comp_matrix)
  
  d <- data.frame(focal, reprod, comp_matrix)
  #calculate total competitors
  d$background <- rowSums(comp_matrix)
  
  #subset lamda events (response in absence of competition) and competition events
  lam <- subset(d, d$background=="0") #!
  comp <- subset(d, d$background>"0") #!
  if(drop0 == TRUE){
    ## drop lambda rows without seed production:
    lam <- subset(lam, reprod!=0)
    comp <- subset(comp, reprod!=0)
  }
  #test and warn if any is = 0
  
  ## get a list of target species to work through sequentially:
  splist <- unique(focal) #this assumes all focal and background are equal.
  #What when you have bg sp not as focal?
  ##alpha order splist:
  splist<-splist[order(splist)]
  
  ## objects to hold the final parameter estimates from model 3: 
  alpha_matrix <- matrix(0, nrow=length(splist), ncol=length(splist)) #better with NA?
  row.names(alpha_matrix) <- splist
  colnames(alpha_matrix) <- splist
  lambda_est <- rep(NA, length(splist))
  sigma_est <- rep(NA, length(splist))
  convergence_code <- rep(NA, length(splist))
  
  ## for each species in turn as a target:
  for(i in 1:length(splist)){
    ## subset out the rows that are needed from the competition df
    comp_points <- subset(comp, comp$focal==splist[i], drop=TRUE)
    ## and the correct lambda plants
    lam_points <- subset(lam, lam$focal==splist[i])
    
    ## now need to build up a vector of nonzero seed production
    ## and corresponding density vectors (held in a matrix) for background species 
    ## to use in model fitting
    
    ## start with the lambda seeds- will add on to this for each background:
    reprod <- lam_points$reprod 
    ##build density matrix (each row will be density of a species)
    dens <- matrix(0, nrow=length(splist), ncol=(nrow(lam_points)+ nrow(comp_points)))
    row.names(dens) <- splist 
    ##for each background species the target competes against:
    background_list <- row.names(dens)
    
    ## use this counter to keep track of which column is next to have data added to
    ## set it to begin after the lambda points:
    start <- nrow(lam_points) +1
    for(j in 1:length(background_list)){
      ## LOOP creates a seeds vector and a corresponding dens ("density") matrix with 
      ## background species as rows, and columns corresponding to the seed production
      ## vector.
      ## beginning columns correspond to lambda plants
      
      #take just the rows pertaining to a specific background sp:
      bg_points <- subset(comp_points, select=c(background_list[j], "reprod"))
      ## which row of the density matrix corresponds?
      rownum <- which(row.names(dens)==background_list[j])
      #column to end with:
      end <- start + nrow(bg_points)-1
      ## drop in density values into matrix:
      dens[rownum, start:end] <- bg_points[,1] #the first column is background_list[j]
      ## add seed numbers into the seeds vector
      reprod<-c(reprod, bg_points$reprod)
    } #close j

    ## should now have "reprod" vector and corresponding "dens"ity matrix
    ## can test ncol(dens)==length(reprod) #BUT IT DO NOT!
    ## we'll be working with log reprod (for giving a lognormal error structure):
    log_reprod <- log(reprod) 
    if(log == FALSE){message("only log = TRUE implmented.log needed for giving a lognormal error structure")}
    #model fitting using optim and earlier likelihood functions
    # model 1, no competition
    #recall parameters are lambda and sigma- initialize these with estimates from the data:
    par1 <- c(mean(log_reprod), sd(log_reprod)) 
    ##repeat optimization until we get convergence (or we try 25 times)
    for(k in 1:op){
      testcomp1<-optim(par = par1, compmodel1)
      ##update start parameters to final estimate to use in next run in case of nonconvergence
      par1<-testcomp1$par
      if(testcomp1$convergence==0){
        print(paste(splist[i], "model 1 converged on rep", k, sep=" "))
        break
      }
    }
    # model 2, one common alpha
    ## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates
    par2<-c(testcomp1$par[1],0.001,testcomp1$par[2])
    ##as before:
    for(k in 1:op){
      ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
      testcomp2<-optim(par2,compmodel2, method="L-BFGS-B", lower=c(1,0,0.0000000001)) #outsource those params
      par2<-testcomp2$par
      if(testcomp2$convergence==0){
        print(paste(splist[i],  "model 2 converged on rep", k, sep=" "))
        break
      }
    }
    # model 3, unique alphas
    par3<-c(testcomp2$par[1], rep(testcomp2$par[2], times=3), testcomp2$par[3])
    ##as before
    for(k in 1:op){
      ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
      testcomp3<-optim(par3,compmodel3, method="L-BFGS-B", lower=c(1, rep(0, times=3),0.0000000001), control=list(maxit=1000))
      par3<-testcomp3$par
      if(testcomp3$convergence==0){
        print(paste(splist[i], "model 3 converged on rep", k, sep=" "))
        break
      }
    }
    # save estimates from model 3 
    lambda_est[i] <- par3[1]
    sigma_est[i] <- par3[length(par3)]
    convergence_code[i] <- testcomp3$convergence
    ## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
    ## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
    ## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
    alphas <- par3[2:(length(par3)-1)]
    alpha_matrix[i,] <- alphas
    ##note that in cases where there is no data for a particular species the 
    #alpha estimate for that species ends up as the starting value- we need to 
    #be careful of these as they are basically gargbage numbers.  Keeping them 
    #in up to now to keep the structure of the data constant, but will set them 
    #to NA here:
    #identify which species have no data in this fit:
    no_data <- which(apply(dens, MARGIN=1, FUN=mean)==0)
    ## set their alphas to NA in the matrix:
    alpha_matrix[i,no_data] <- NA
    ## some diagnostics
    ## print an error to the console if any one of the three models failed to converge:
    if(testcomp1$convergence + testcomp2$convergence + testcomp3$convergence !=0){
      warning(paste("at least one model did not converge for", splist[i], sep=" "))
    }
  }
  #store results
  results <- data.frame(splist, lambda_est, sigma_est, convergence_code)
  out <- list(results, alpha_matrix)
  out
}
  
  #basic models to fit
  #CRAP, I need to call objects from outside the function!

  #model 1 - no effect of density (no competitive effects)
  compmodel1 <- function(par, log_reprod = log_reprod){
    # lambda and sigma parameters for the normal distribution
    #(assuming lognormal error- seed data are logged) #!
    lambda <- par[1]
    sigma <- par[2]
    #this is the predictive model- here is just fitting a horizontal
    #line through the data:
    pred <- rep(lambda, times=length(log_reprod)) #will this work in a function? 
    #these are the log likelihoods of the data given the model + parameters
    llik <- dnorm(log_reprod,log(pred), sd=sigma, log=TRUE) #WTF??
    #return the sum of negative log likelihoods - what optim minimizes
    return(sum(-1*llik)) 
  }
  
  #model 2 - competition, but no difference between species
  compmodel2<-function(par, log_reprod = log_reprod, background = background){
    lambda<-par[1] ## same as model 1
    alpha<-par[2]  ## new parameter introduced in model 2
    sigma<-par[3] ## same as model 1
    ## apply is used to get a single vector of densities rather than a matrix
    ## there is only one competitor at a time here- just adding each density to a 
    ## column of other zeros:
    ## predictive model:
    pred <- lambda/(1+alpha*(background)) #is this the right background?? Oscar uses this one, but we should be in a  
    ## log likelihoods of data given the model + parameters:
    llik<-dnorm(log_reprod,log(pred), sd=sigma, log=TRUE)
    ## return sum of negative log likelihoods:
    return(sum(-1*llik)) 
  }
  
  #model 3 - all species have different competitive effects
  compmodel3 <- function(par, comp_matrix = comp_matrix, log_reprod = log_reprod){
    lambda <- par[1] #same as model 2
    a_comp <- par[2:(length(par)-1)]
    ## new parameters- use alpha estimate from model 2 as start 
    #value for fitting
    sigma <- par[length(par)] ## same as model 2
    ## predictive model:
    term = 1 #create the denominator term for the model
    for(i in 1:ncol(comp_matrix)){
      term <- term + a_comp[i]*comp_matrix[,i] #idem, why the not subseted one?
    }
    pred<- lambda/ term
    # likelihood as before:
    llik<-dnorm(log_reprod,log(pred), sd=sigma, log=TRUE)
    # return sum of negative log likelihoods
    return(sum(-1*llik)) #sum of negative log likelihoods
  }
  