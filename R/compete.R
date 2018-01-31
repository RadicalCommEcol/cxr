#' Calculates pairwise competitive effects between species
#'
#' Calculates how reproduction decays with number of intra and inter-specific
#'  compettors. The idea is to fit a series of nested models. First model fits
#'  an intercept (lambda). Second model fits a common response to all competitors (alpha), 
#'  Third model fits a pairwise effct of competition (alpha_ij). And a fourth and fiveth models 
#'  fit covariates on lambda's and alpha's (l_cov and a_cov). Parameters obtained in each model serves as priors for the
#'  next model. 
#' 
#' @param focal vector with the focal species name of each observation
#' @param fitness Vector of the reproductive success or biomass or any other fitness measure of each observation
#' @param comp_matrix data.frame with as many columns as competitor species. Each cell is the number of competitors per observtion
#' @param covariates data.frame with as many columns as covariates measured.
#' @param model you can provide a function with a model to test if you don't like the default. This is not implemented yet. 
#' In any case models need to have the same initial par's as the ones implemented to be easy to plug.
#' @param drop0 logical. Should observations with zero fitness be discarded?
#' @param log logical. Should fitnes be loged. Forced to be TRUE for now. 
#' For accepting non log fitness, the dnorm has to be modified in the models.
#' @param op number of iterations to ptimize the models.
#' @param method Options passed to `optim()` Default "L-BFGS-B". 
#' @param lower1 Options passed to `optim()` One value per par. For model one it needs 2 values.
#' @param upper1 Options passed to `optim()` One value per par. For model one it needs 2 values.
#' @param control1 Options passed to `optim()`. It should be a list.
#' @param hessian1 Options passed to `optim()`. Default = TRUE. This provides SE for all params
#' @param gr1 Options passed to `optim()`. Usually not used by us, but for completness.
#' @param lower2 Options passed to `optim()` One value per par. For model two it needs 3 values.
#' @param upper2 Options passed to `optim()` One value per par. For model two it needs 3 values.
#' @param control2 Options passed to `optim()`. It should be a list.
#' @param hessian2 Options passed to `optim()`. Default = TRUE. This provides SE for all params
#' @param gr2 Options passed to `optim()`. Usually not used by us, but for completness.
#' @param lower3 Options passed to `optim()` One value per par. For model three it needs 2+ncol(comp_matrix) values.
#' @param upper3 Options passed to `optim()` One value per par. For model three it needs 2+ncol(comp_matrix) values.
#' @param control3 Options passed to `optim()`. It should be a list.
#' @param hessian3 Options passed to `optim()`. Default = TRUE. This provides SE.
#' @param gr3 Options passed to `optim()`. Usually not used by us, but for completness. 
#' @param lower4 Options passed to `optim()` One value per par. For model four it needs 2+ncol(comp_matrix)+(2*ncol(covariates)) values.
#' @param upper4 Options passed to `optim()` One value per par. For model four it needs 2+ncol(comp_matrix)+(2*ncol(covariates)) values.
#' @param control4 Options passed to `optim()`. It should be a list.
#' @param hessian4 Options passed to `optim()`. Default = TRUE. This provides SE.
#' @param gr4 Options passed to `optim()`. Usually not used by us, but for completness.
#' @param lower5 Options passed to `optim()` One value per par. For model five it needs 
#' 2+ncol(comp_matrix)+ncol(covariates)+(ncol(comp_matrix)*ncol(covariates)) values.
#' @param upper5 Options passed to `optim()` One value per par. For model five it needs 
#' 2+ncol(comp_matrix)+ncol(covariates)+(ncol(comp_matrix)*ncol(covariates)) values.
#' @param control5 Options passed to `optim()`. It should be a list.
#' @param hessian5 Options passed to `optim()`. Default = TRUE. This provides SE.
#' @param gr5 Options passed to `optim()`. Usually not used by us, but for completness. 
#' 
#' @note This script follows previous methodology developed by Godoy et al 2014,
#'  Kraft et al. 2015 and Lanuza et al 2018 papers
#'  
#' @details Useful control params are maxit=1000, or par scale (with as much elements as pars has the model)  
#'  
#' @return a list of estimated parameters for each model.
#' @examples
#'  
compete <- function(focal, fitness, comp_matrix, #basic data needed
                    covariates = NULL, 
                    model = NULL, 
                    drop0 = FALSE, log = TRUE, op = 25, 
                    method="L-BFGS-B", 
                    #m1
                    lower1 = -Inf, upper1 = Inf, 
                    control1 = list(), gr1 = NULL, hessian1 = TRUE,
                    #m2
                    lower2 = -Inf, upper2 = Inf, 
                    control2 = list(), gr2 = NULL, hessian2 = TRUE,
                    #m3
                    lower3 = -Inf, upper3 = Inf, 
                    control3 = list(), gr3 = NULL, hessian3 = TRUE,
                    #m4
                    lower4 = -Inf, upper4 = Inf, 
                    control4 = list(), gr4 = NULL, hessian4 = TRUE,
                    #m5
                    lower5 = -Inf, upper5 = Inf, 
                    control5 = list(), gr5 = NULL, hessian5 = TRUE){ 
  #We can check focal is c() and length match fitness and nrow(comp_matrix), etc...
  if(is.null(covariates)){
    d <- data.frame(focal, fitness, comp_matrix)
    n_cov <- 0
  } else{
    d <- data.frame(focal, fitness, comp_matrix, covariates)
    n_cov <- ncol(covariates)
  }
  #DELETE  
  #calculate total competitors
  #d$background <- rowSums(comp_matrix) #not used in d Delete?
  if(drop0 == TRUE){
    ## drop lambda rows without seed production:
    d <- subset(d, fitness!=0)
  }
  #DELETE
  #subset lamda events (response in absence of competition) and competition events
  #lam <- subset(d, d$background=="0") #!
  #comp <- subset(d, d$background>"0") #!
  #if(drop0 == TRUE){
    ## drop lambda rows without seed production:
   # lam <- subset(lam, reprod!=0)
  #  comp <- subset(comp, reprod!=0)
  #}
  #test and warn if any is = 0
  
  ## get a list of target species to work through sequentially:
  splist <- unique(focal) 
  bglist <- colnames(comp_matrix)
  ##alpha order splist:
  splist <- splist[order(splist)]
  bglist <- bglist[order(bglist)]
  
  ## objects to hold the final parameter estimates from all model: 
  #m3
  alpha_matrix3 <- matrix(NA, nrow=length(splist), ncol=length(bglist)) 
  row.names(alpha_matrix3) <- splist
  colnames(alpha_matrix3) <- bglist
  lambda_est3 <- rep(NA, length(splist))
  sigma_est3 <- rep(NA, length(splist))
  convergence_code3 <- rep(NA, length(splist))
  loglike3 <- rep(NA, length(splist))
  alpha_lower_error3 <- matrix(NA, nrow = length(splist), ncol = length(bglist))
  alpha_upper_error3 <- matrix(NA, nrow = length(splist), ncol = length(bglist))
  row.names(alpha_lower_error3) <- splist
  colnames(alpha_lower_error3) <- bglist
  row.names(alpha_upper_error3) <- splist
  colnames(alpha_upper_error3) <- bglist
  lambda_error3 <- matrix(NA, nrow = length(splist), ncol = 2)
  row.names(lambda_error3) <- splist
  colnames(lambda_error3) <- c("lower.error", "upper.error")
  #m1
  lambda_est1 <- lambda_est3
  sigma_est1 <- sigma_est3
  convergence_code1 <- convergence_code3
  loglike1 <- loglike3
  lambda_error1 <- lambda_error3
  #m2
  alpha_est2 <- rep(NA, length(splist))
  lambda_est2 <- rep(NA, length(splist))
  sigma_est2 <- sigma_est3
  convergence_code2 <- convergence_code3
  loglike2 <- loglike3
  alpha_error2 <- lambda_error3
  lambda_error2 <- lambda_error3
  #m4
  alpha_matrix4 <- alpha_matrix3
  lambda_est4 <- lambda_est3
  sigma_est4 <- lambda_est3
  l_cov_est4 <- rep(NA, n_cov)
  a_cov_est4 <- rep(NA, n_cov)
  convergence_code4 <- convergence_code3
  loglike4 <- loglike3
  alpha_lower_error4 <- alpha_lower_error3
  alpha_upper_error4 <- alpha_upper_error3
  lambda_error4 <- lambda_error3
  l_cov_error4 <- matrix(NA, nrow = length(splist), ncol = n_cov*2)
  row.names(l_cov_error4) <- splist
  colnames(l_cov_error4) <- c(rep("lower.error", n_cov), rep("upper.error", n_cov))
  a_cov_error4 <- l_cov_error4
  #m5
  alpha_matrix5 <- alpha_matrix3
  lambda_est5 <- lambda_est3
  sigma_est5 <- lambda_est3
  l_cov_est5 <- l_cov_est4
  a_cov_est5 <- matrix(NA, nrow = length(splist), ncol = n_cov*ncol(comp_matrix))
  row.names(a_cov_est5) <- splist
  colnames(a_cov_est5) <- rep("a_cov", n_cov*ncol(comp_matrix))
  convergence_code5 <- convergence_code3
  loglike5 <- loglike3
  alpha_lower_error5 <- alpha_lower_error3
  alpha_upper_error5 <- alpha_upper_error3
  lambda_error5 <- lambda_error3
  l_cov_error5 <- l_cov_error4
  a_cov_error5 <- matrix(NA, nrow = length(splist), ncol = n_cov*ncol(comp_matrix)*2)
  row.names(a_cov_error5) <- splist
  colnames(a_cov_error5) <- c(rep("lower.error", n_cov*ncol(comp_matrix)), rep("upper.error", n_cov*ncol(comp_matrix)))

  ## for each species in turn as a target:
  for(i in 1:length(splist)){
    #DELETE
    ## subset out the rows that are needed from the competition df
    #comp_points <- subset(comp, comp$focal==splist[i], drop=TRUE)
    ## and the correct lambda plants
    #lam_points <- subset(lam, lam$focal==splist[i])
    
    ## now need to build up a vector of nonzero seed production
    ## and corresponding density vectors (held in a matrix) for background species 
    ## to use in model fitting
    
    ## start with the lambda seeds- will add on to this for each background:
    #fitness <- lam_points$fitness 
    ##build density matrix (each row will be density of a species)
    #dens <- matrix(NA, nrow=length(splist), ncol=(nrow(lam_points) + nrow(comp_points)))
    #row.names(dens) <- splist 
    ##for each background species the target competes against:
    #background_list <- row.names(dens)
    
    ## use this counter to keep track of which column is next to have data added to
    ## set it to begin after the lambda points:
    #start <- nrow(lam_points) +1
    #for(j in 1:length(background_list)){
      ## LOOP creates a seeds vector and a corresponding dens ("density") matrix with 
      ## background species as rows, and columns corresponding to the seed production
      ## vector.
      ## beginning columns correspond to lambda plants
      
      #take just the rows pertaining to a specific background sp:
      #bg_points <- subset(comp_points, select=c(background_list[j], "reprod"))
      ## which row of the density matrix corresponds?
      #rownum <- which(row.names(dens)==background_list[j])
      #column to end with:
      #end <- start + nrow(bg_points)-1
      ## drop in density values into matrix:
      #dens[rownum, start:end] <- bg_points[,1] #the first column is background_list[j]
      ## add seed numbers into the seeds vector
      #reprod<-c(reprod, bg_points$reprod)
    #} #close j
    
    ######THIS LOOP can be changes to a t(subset(d, focal = i)[comp matrix]) for dens
    ## and reprod is just reprod subseted for sp i. 
    
    ## subset out the rows that are needed from the competition df
    #comp <- subset(d, focal == splist[i]) #DELETE if nothing is broken.
    comp_matrix_i <- comp_matrix[which(d$focal == splist[i]),]
    assign("comp_matrix_i", comp_matrix_i, envir = .GlobalEnv) #needed to run the models within the function.
    n_bg <- dim(comp_matrix_i)[2]
    assign("n_bg", n_bg, envir = .GlobalEnv)
    background <- rowSums(comp_matrix_i)
    assign("background", background, envir = .GlobalEnv)
    ## we'll be working with log fitness (for giving a lognormal error structure)
    log_fitness <- log(fitness[which(d$focal == splist[i])]) 
    assign("log_fitness", log_fitness, envir = .GlobalEnv)
    
    if(log == FALSE){message("only log = TRUE implmented. log needed for giving a lognormal error structure")}
    #model fitting using optim and earlier likelihood functions
    # model 1, no competition
    #recall parameters are lambda and sigma- initialize these with estimates from the data:
    par1 <- c(mean(log_fitness), #lambda
              sd(log_fitness))  #sigma
    ##repeat optimization until we get convergence (or we try 25 times)
    for(k in 1:op){
      testcomp1 <- optim(par = par1, fn = compmodel1,
                         method = method,
                         gr = gr1, lower = lower1, upper = upper1,
                         control = control1, hessian = hessian1)
      ##update start parameters to final estimate to use in next run in case of nonconvergence
      par1 <- testcomp1$par
      if(testcomp1$convergence == 0){
        message(paste(splist[i], "model 1 converged on rep", k, sep = " "))
        break
      }
    }
    if(hessian1 == TRUE){
      #calculate Confidence intervals at 95% of the estimates, both lambda and alphas
      inverse <- solve(testcomp1$hessian)
      errors <- sqrt(diag(inverse))  
      # save estimates errors
      lambda_error1[i,1] <- testcomp1$par[1]-1.96*errors[1]
      lambda_error1[i,2] <- testcomp1$par[1]+1.96*errors[1]    
    }
    # model 2, one common alpha
    ## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates
    par2 <- c(testcomp1$par[1], #lambda
              0.0001,            # alfa 
              testcomp1$par[2]) # sigma
    ##as before:
    for(k in 1:op){
      ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
      testcomp2 <- optim(par2, compmodel2, method = method,
                       gr = gr2, lower = lower2, upper = upper2,
                       control = control2, hessian = hessian2)
      par2 <- testcomp2$par
      if(testcomp2$convergence == 0){
        message(paste(splist[i],  "model 2 converged on rep", k, sep = " "))
        break
      }
    }
    if(hessian2 == TRUE){
      #calculate Confidence intervals at 95% of the estimates, both lambda and alphas
      inverse <- solve(testcomp2$hessian)
      errors <- sqrt(diag(inverse))  
      # save estimates errors
      alpha_error2[i,] <- c(testcomp2$par[2]-1.96*errors[2], testcomp2$par[2]+1.96*errors[2])
      lambda_error2[i,1] <- testcomp2$par[1]-1.96*errors[1]
      lambda_error2[i,2] <- testcomp2$par[1]+1.96*errors[1]    
    }
    # model 3, unique alphas
    par3 <- c(testcomp2$par[1], #lambda
              rep(testcomp2$par[2], times = n_bg),  #paired alfas
              testcomp2$par[3]) #sigma
    ##as before
    for(k in 1:op){
      ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
      testcomp3 <- optim(par3, compmodel3, 
                       gr = gr3, method = method, lower = lower3, upper = upper3,
                       control = control3, hessian = hessian3) #check hessian can be done with all methods.
      par3 <- testcomp3$par
      if(testcomp3$convergence == 0){
        message(paste(splist[i], "model 3 converged on rep", k, sep = " "))
        break
      }
    }
    if(hessian3 == TRUE){
      #calculate Confidence intervals at 95% of the estimates, both lambda and alphas
      inverse <- solve(testcomp3$hessian)
      errors <- sqrt(diag(inverse))  
      # save estimates errors
      alpha_lower_error3[i,] <- testcomp3$par[2:(n_bg+1)]-1.96*errors[2:(n_bg+1)] 
      alpha_upper_error3[i,] <- testcomp3$par[2:(n_bg+1)]+1.96*errors[2:(n_bg+1)]
      lambda_error3[i,1] <- testcomp3$par[1]-1.96*errors[1]
      lambda_error3[i,2] <- testcomp3$par[1]+1.96*errors[1]    
    }
    
    #When covariates are in:
    if(!is.null(covariates)){
      n_cov <- ncol(covariates) 
      assign("n_cov", n_cov, envir = .GlobalEnv)
      covariates_i <- covariates[which(focal == splist[i]), ,drop = FALSE]
      assign("covariates_i", covariates_i, envir = .GlobalEnv)
      # model 4!
      par4 <- c(testcomp3$par[1],  #lambda 1
                rep(0.0001, times = n_cov), #l_cov n_cov
                testcomp3$par[2:(n_bg+1)], #alfas n_bg
                rep(0.0001, times = n_cov), #a_cov n_cov
                testcomp3$par[length(par3)]) # sigma 1
      ##as before
      for(k in 1:op){
        ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
        testcomp4 <- optim(par4, compmodel4, 
                           gr = gr4, method = method, lower = lower4, upper = upper4,
                           control = control4, hessian = hessian4) #check hessian can be done with all methods.
        par4 <- testcomp4$par
        if(testcomp4$convergence == 0){
          message(paste(splist[i], "model 4 converged on rep", k, sep = " "))
          break
        }
      }
      if(hessian4 == TRUE){
        #calculate Confidence intervals at 95% of the estimates, both lambda and alphas
        inverse <- solve(testcomp4$hessian)
        errors <- sqrt(diag(inverse))  
        # save estimates errors
        lambda_error4[i,1] <- testcomp4$par[1]-1.96*errors[1]
        lambda_error4[i,2] <- testcomp4$par[1]+1.96*errors[1]   
        l_cov_error4[i,] <- c(testcomp4$par[2:(1+n_cov)]-1.96*errors[2:(1+n_cov)], 
                              testcomp4$par[2:(1+n_cov)]+1.96*errors[2:(1+n_cov)]) 
        alpha_lower_error4[i,] <- testcomp4$par[(1+n_cov+1):(1+n_cov+n_bg)]-1.96*errors[(1+n_cov+1):(1+n_cov+n_bg)]
        alpha_upper_error4[i,] <- testcomp4$par[(1+n_cov+1):(1+n_cov+n_bg)]+1.96*errors[(1+n_cov+1):(1+n_cov+n_bg)]
        a_cov_error4[i,] <- c(testcomp4$par[(2+n_cov+n_bg):(n_cov+n_bg+1+n_cov)]-1.96*errors[(2+n_cov+n_bg):(n_cov+n_bg+1+n_cov)], 
                              testcomp4$par[(2+n_cov+n_bg):(n_cov+n_bg+1+n_cov)]+1.96*errors[(2+n_cov+n_bg):(n_cov+n_bg+1+n_cov)])
      }
      #model 5
      #subset a_cov's
      par4_a_cov <- testcomp4$par[(1+n_cov+n_bg+1):(1+n_cov+n_bg+n_cov)]
      par5_a_cov <- c()
      for(w in 1:length(par4_a_cov)){
        par5_a_cov <- c(par5_a_cov, rep(par4_a_cov[w], times = n_bg))
      }
      par5 <- c(testcomp4$par[1],  #lambda
                testcomp4$par[(1+1):(1+n_cov)], #l_cov
                testcomp4$par[(1+n_cov+1):(1+n_cov+n_bg)], #alfas
                par5_a_cov, #a_cov pairwise
                testcomp4$par[3]) # sigma
      ##as before
      for(k in 1:op){
        ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
        testcomp5 <- optim(par5, compmodel5, 
                           gr = gr5, method = method, lower = lower5, upper = upper5,
                           control = control5, hessian = hessian5) #check hessian can be done with all methods.
        par5 <- testcomp5$par
        if(testcomp5$convergence == 0){
          message(paste(splist[i], "model 5 converged on rep", k, sep = " "))
          break
        }
      }
      if(hessian5 == TRUE){
        #calculate Confidence intervals at 95% of the estimates, both lambda and alphas
        inverse <- solve(testcomp5$hessian)
        errors <- sqrt(diag(inverse))  
        # save estimates errors
        lambda_error5[i,1] <- testcomp5$par[1]-1.96*errors[1]
        lambda_error5[i,2] <- testcomp5$par[1]+1.96*errors[1]   
        l_cov_error5[i,] <- c(testcomp5$par[2:(1+n_cov)]-1.96*errors[2:(1+n_cov)], 
                              testcomp5$par[2:(1+n_cov)]+1.96*errors[2:(1+n_cov)]) 
        alpha_lower_error5[i,] <- testcomp5$par[(1+n_cov+1):(1+n_cov+n_bg)]-1.96*errors[(1+n_cov+1):(1+n_cov+n_bg)]
        alpha_upper_error5[i,] <- testcomp5$par[(1+n_cov+1):(1+n_cov+n_bg)]+1.96*errors[(1+n_cov+1):(1+n_cov+n_bg)]
        a_cov_error5[i,] <- c(testcomp5$par[(1+n_cov+n_bg+1):(1+n_cov+n_bg+(n_bg*n_cov))]-1.96*errors[(1+n_cov+n_bg+1):(1+n_cov+n_bg+(n_bg*n_cov))], 
                              testcomp5$par[(1+n_cov+n_bg+1):(1+n_cov+n_bg+(n_bg*n_cov))]+1.96*errors[(1+n_cov+n_bg+1):(1+n_cov+n_bg+(n_bg*n_cov))])
      }
      #delete objects
      rm(n_cov,
         covariates_i,
         inherits = TRUE)
    }
    # delete objects from environment
    rm(log_fitness,
       background,
       comp_matrix_i, inherits = TRUE) #THIS is not really cleaning it...
    # save estimates from all model 
    #m1
    lambda_est1[i] <- par1[1]
    sigma_est1[i] <- par1[length(par1)]
    convergence_code1[i] <- testcomp1$convergence
    loglike1[i] <- -1*testcomp1$value
    #m2
    lambda_est2[i] <- par2[1]
    sigma_est2[i] <- par3[length(par2)]
    convergence_code2[i] <- testcomp2$convergence
    loglike2[i] <- -1*testcomp2$value
    #m3
    lambda_est3[i] <- par3[1]
    sigma_est3[i] <- par3[length(par3)]
    convergence_code3[i] <- testcomp3$convergence
    loglike3[i] <- -1*testcomp3$value
    if(!is.null(covariates)){
      #m4
      lambda_est4[i] <- par4[1]
      l_cov_est4 <- par4[2:(1+n_cov)] 
      a_cov_est4 <- par4[(2+n_cov+n_bg):(n_cov+n_bg+1+n_cov)] 
      sigma_est4[i] <- par4[length(par4)]
      convergence_code4[i] <- testcomp4$convergence
      loglike4[i] <- -1*testcomp4$value
      #m5
      lambda_est5[i] <- par5[1]
      l_cov_est5 <- par5[2:(1+n_cov)] 
      a_cov_est5 <- par5[(1+n_cov+n_bg+1):(1+n_cov+n_bg+(n_bg*n_cov))] 
      sigma_est5[i] <- par5[length(par5)]
      convergence_code5[i] <- testcomp5$convergence
      loglike5[i] <- -1*testcomp5$value
      }
    
    ## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
    ## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
    ## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
    alpha_est2[i] <- par2[2]
    alpha_matrix3[i,] <- par3[2:(length(par3)-1)]
    if(!is.null(covariates)){
      alpha_matrix4[i,] <- par4[(1+n_cov+1):(1+n_cov+n_bg)]
      alpha_matrix5[i,] <- par5[(1+n_cov+1):(1+n_cov+n_bg)]
    }
  } #close i
    #DELETE
    ##note that in cases where there is no data for a particular species the 
    #alpha estimate for that species ends up as the starting value- we need to 
    #be careful of these as they are basically gargbage numbers.  Keeping them 
    #in up to now to keep the structure of the data constant, but will set them 
    #to NA here:
    #identify which species have no data in this fit:
    #no_data <- which(apply(dens, MARGIN=1, FUN=mean)==0)
    ## set their alphas to NA in the matrix:
    #alpha_matrix[i,no_data] <- NA
    ### I THINK THIS IS DONE.
    
    ## some diagnostics
    ## print an error to the console if any one of the three models failed to converge:
    if(!is.null(covariates)){
      if(testcomp1$convergence + testcomp2$convergence + testcomp3$convergence
         + testcomp4$convergence + testcomp5$convergence !=0){
        warning(paste("Sorry, at least one model did not converge for", splist[i], sep=" "))
        } 
      } else {
      if(testcomp1$convergence + testcomp2$convergence + testcomp3$convergence !=0){
        warning(paste("Sorry, at least one model did not converge for", splist[i], sep=" "))
        }
      }
  #store results
  results <- data.frame(splist, lambda_est1, lambda_error1, sigma_est1, convergence_code1, loglike1,
                        alpha_est2, alpha_error2, lambda_est2, lambda_error2, sigma_est2, convergence_code2, loglike2,
                        lambda_est3, lambda_error3, sigma_est3, convergence_code3, loglike3,
                        lambda_est4, lambda_error4, sigma_est4, convergence_code4, loglike4,
                        lambda_est5, lambda_error5, sigma_est5, convergence_code5, loglike5)
  colnames(results)[c(3,4,9,10)] <- c("lower.error.1", "upper.error.1", "alpha.lower.error.2", "alpha.upper.error.2")
  out <- list(results, alpha_matrix3, alpha_lower_error3, alpha_upper_error3,
              alpha_matrix4, alpha_lower_error4, alpha_upper_error4,
              l_cov_est4, l_cov_error4, a_cov_est4, a_cov_error4,
              alpha_matrix5, alpha_lower_error5, alpha_upper_error5,
              l_cov_est5, l_cov_error5, a_cov_est5, a_cov_error5 
              )
  names(out) <- c("lambdas$co", "alpha_matrix3", "alpha_lower_error3", "alpha_upper_error3",
                  "alpha_matrix4", "alpha_lower_error4", "alpha_upper_error4",
                  "l_cov_est4", "l_cov_error4", "a_cov_est4", "a_cov_error4",
                  "alpha_matrix5", "alpha_lower_error5", "alpha_upper_error5",
                  "l_cov_est5", "l_cov_error5", "a_cov_est5", "a_cov_error5")
  out
  }

  #basic models to fit

  #model 1 - no effect of density (no competitive effects)
  compmodel1 <- function(par){ 
    #lambda and sigma parameters for the normal distribution
    #(assuming lognormal error- seed data are logged) 
    lambda <- par[1]
    sigma <- par[2]
    #this is the predictive model- here is just fitting a horizontal
    #line through the data:
    pred <- rep(lambda, times=length(log_fitness)) 
    #these are the log likelihoods of the data given the model + parameters
    llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
    #return the sum of negative log likelihoods - what optim minimizes
    return(sum(-1*llik)) 
  }
  
  #model 2 - competition, but no difference between species
  compmodel2 <- function(par){ 
    lambda <- par[1] ## same as model 1
    alpha <- par[2]  ## new parameter introduced in model 2
    sigma <- par[3] ## same as model 1
    ## apply is used to get a single vector of densities rather than a matrix
    ## there is only one competitor at a time here- just adding each density to a 
    ## column of other zeros:
    ## predictive model:
    pred <- lambda/(1+alpha*(background))  
    ## log likelihoods of data given the model + parameters:
    llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
    ## return sum of negative log likelihoods:
    return(sum(-1*llik)) 
  }
  
  #model 3 - all species have different competitive effects
  compmodel3 <- function(par){
    lambda <- par[1] #same as model 2
    a_comp <- par[2:(length(par)-1)] # new parameters- use alpha estimate from model 2 as start 
    #value for fitting
    sigma <- par[length(par)] ## same as model 2
    ## predictive model:
    term = 1 #create the denominator term for the model
    for(z in 1:ncol(comp_matrix_i)){
      term <- term + a_comp[z]*comp_matrix_i[,z] 
      #THIS IS AS DONE IN LINCX, CHECK WITH OSCAR
    }
    pred <- lambda/ term
    # likelihood as before:
    llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log = TRUE)
    # return sum of negative log likelihoods
    return(sum(-1*llik)) #sum of negative log likelihoods
  }
  

  #model 4 - 
  compmodel4 <- function(par){
    lambda <- par[1] 
    l_cov <- par[(1+1):(1+n_cov)] #effect of cov 1, 2, ... on lambda
    a_comp <- par[(1+n_cov+1):(1+n_cov+n_bg)] # alfas_ij
    a_cov <- par[(1+n_cov+n_bg+1):(1+n_cov+n_bg+n_cov)] # common effect of cov 1, 2... on alphas
    sigma <- par[length(par)]
    num = 1
    for(z in 1:n_cov){
      num <- num + l_cov[z]*covariates_i[,z] 
    }
    cov_term <- 0 
      for(v in 1:n_cov){
        cov_term <- cov_term + a_cov[v] * covariates_i[,v]
      }
    term <- 1 #create the denominator term for the model
    for(z in 1:ncol(comp_matrix_i)){
        term <- term + (a_comp[z] + cov_term) * comp_matrix_i[,z] 
    }
    pred <- lambda * (num) / term 
    # likelihood as before:
    llik<-dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
    # return sum of negative log likelihoods
    return(sum(-1*llik)) #sum of negative log likelihoods
  }
  
  #model 5 - 
  compmodel5 <- function(par){
    lambda <- par[1] 
    l_cov <- par[(1+1):(1+n_cov)] #effect of cov 1, 2... on lambda
    a_comp <- par[(1+n_cov+1):(1+n_cov+n_bg)] #alfas_ij
    a_cov <- par[(1+n_cov+n_bg+1):(1+n_cov+n_bg+(n_cov*n_bg))] #effects of cov 1, 2... on alpha_i
    sigma <- par[length(par)]
    num = 1
    for(v in 1:n_cov){
      num <- num + l_cov[v]*covariates_i[,v] 
    }
    cov_term_x <- list()
    for(v in 1:n_cov){
      cov_temp <- covariates_i[,v]
      for(z in 1:n_bg){
          cov_term_x[[z+(n_bg*(v-1))]] <- a_cov[z+(n_bg*(v-1))] * cov_temp  #create  a_cov_i*cov_i vector
      }
    }
    cov_term <- list()
    #here I need to reformat cov_term_x to sumatories of the form a_cov_i* cov_i + a_covj* cov_j + ...
    for(z in 0:(n_bg-1)){
      cov_term_x_sum <- cov_term_x[[z+1]] 
      for(v in 1:n_cov){ #IT WAS 2: bfbfbf...
        cov_term_x_sum <- cov_term_x_sum + cov_term_x[[v+n_bg]]
      } 
      cov_term[[z+1]] <- cov_term_x_sum
    }
    term <- 1 #create the denominator term for the model
    for(z in 1:n_bg){
      term <- term + (a_comp[z] + cov_term[[z]]) * comp_matrix_i[,z]  
    }
    pred <- lambda * (num) / term 
    # likelihood as before:
    llik <- dnorm(log_fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
    # return sum of negative log likelihoods
    return(sum(-1*llik)) #sum of negative log likelihoods
    }

  
  
  
  
  