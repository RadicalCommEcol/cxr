library(tidyverse)

#' Title
#'
#' @param focal.sp number of focal species
#' @param num.sp total number of species, including focal ones
#' @param num.cov number of covariates
#' @param num.obs number of observations/sites
#' @param fitness.model model to generate data from, from 1 to 5 in increasing levels of complexity. TODO: document explaining the different models
#' @param focal.lambda vector with lambdas of the focal sp
#' @param alpha.matrix num.sp x num.sp interaction matrix
#' @param alpha.cov list containing num.sp matrices of size num.sp * num.cov, giving the effect of a given covariate -column- on the interaction alpha i -list index- ,j -matrix row-
#' TODO: make it more understandable, maybe alpha.cov[[num.cov]][i.sp,j.sp]
#' @param lambda.cov list containing num.sp vectors giving the effect of each covariate on the lambda of species i
#'
#' @return dataset with a fitness metric calculated for each focal species and observation, according to the fitness model selected
#' @export
#'
#' @examples
GenerateTestData <- function(focal.sp = 1,num.sp = 2,num.cov = 2,num.obs = 10, fitness.model = 1,focal.lambda,alpha.matrix,alpha.cov,lambda.cov){
  
  full.data <- NULL
  
  # focal.sp <- paste("sp",as.character(focal.sp),sep="")
  
  for(i.focal in 1:length(focal.sp)){
    
    temp.data <- data.frame(focal = rep(focal.sp[i.focal],num.obs))
    
    # neighbour occurrences
    neigh.data <- matrix(data = round(rpois(num.sp*num.obs,0.5)),nrow = num.obs)
    colnames(neigh.data) <- paste("sp",as.character(1:num.sp),sep="")
    temp.data <- cbind(temp.data,neigh.data)
    
    # covariates
    cov.data <- matrix(data = round(runif(num.cov*num.obs,0,10)),nrow = num.obs)
    colnames(cov.data) <- paste("cov",as.character(1:num.cov),sep="")
    temp.data <- cbind(temp.data,cov.data)
    
    full.data <- rbind(full.data,temp.data)
  }# for i.focal

  #############
  # fitness model
  full.data$focal.lambda <- rep(focal.lambda, each=num.obs)
  full.data$fitness <- 0
  
  if(fitness.model == 1){
    
    full.data$fitness <- full.data$focal.lambda
    
  }else if(fitness.model == 2){
    
    for(i.obs in 1:nrow(full.data)){
      my.focal.sp <- which(focal.sp == full.data$focal[i.obs])
      den <- 1+(sum(alpha.matrix[my.focal.sp,])*sum(full.data[i.obs,2:(num.sp+1)]))
      full.data$fitness[i.obs] <- full.data$focal.lambda[i.obs]/den
    }
    
  }else if(fitness.model == 3){
    
    for(i.obs in 1:nrow(full.data)){
      my.focal.sp <- which(focal.sp == full.data$focal[i.obs])
      den <- 1
      focal.comp <- full.data[i.obs,2:(num.sp+1)]
      for(i.comp in 1:length(focal.comp)){
        den <- den + alpha.matrix[my.focal.sp,i.comp]*focal.comp[i.comp]
      }
      full.data$fitness[i.obs] <- full.data$focal.lambda[i.obs]/den
    }
    
  }else if(fitness.model == 4){
    
    for(i.obs in 1:nrow(full.data)){
      my.focal.sp <- which(focal.sp == full.data$focal[i.obs])
      focal.comp <- full.data[i.obs,2:(num.sp+1)]
      
    num = 1
    for(i.cov in 1:num.cov){
      num <- num + lambda.cov[[my.focal.sp]][i.cov]*full.data[i.obs,1+num.sp+i.cov] 
    }
    
    cov_term <- 0 
    for(j.cov in 1:num.cov){
      # sum(alpha.cov...) in case a full alpha.cov matrix is specified (i.e. with different effects for different alphas), 
      # which could be the case for generating with model 5
      cov_term <- cov_term + sum(alpha.cov[[my.focal.sp]][,j.cov]) * full.data[i.obs,1+num.sp+i.cov]
    }
    
    term <- 1 #create the denominator term for the model
    for(i.comp in 1:length(focal.comp)){
      term <- term + (alpha.matrix[my.focal.sp,i.comp] + cov_term) * focal.comp[i.comp] 
    }
    full.data$fitness[i.obs] <- full.data$focal.lambda[i.obs] * (num) / term 
    }# for i.obs
    
  }else if(fitness.model == 5){
    
    for(i.obs in 1:nrow(full.data)){
      my.focal.sp <- which(focal.sp == full.data$focal[i.obs])
      focal.comp <- full.data[i.obs,2:(num.sp+1)]
      
      num = 1
      for(i.cov in 1:num.cov){
        num <- num + lambda.cov[[my.focal.sp]][i.cov]*full.data[i.obs,1+num.sp+i.cov] 
      }
      
      # for each alpha term and each covariate, compute
      # alpha_ij + (covariate_k * covariate effect_k_ij)
      # then, each summed alpha is multiplied by the density of the competitor i.comp
      # and the overall sum (+1) is the denominator
      term <- 1
      for(i.comp in 1:num.sp){
        alpha.cov.term <- alpha.matrix[my.focal.sp,i.comp]
        for(i.cov in 1:num.cov){
          alpha.cov.term <- (alpha.cov.term + (full.data[i.obs,1+num.sp+i.cov] * alpha.cov[[my.focal.sp]][i.comp,i.cov]))
        }
        term <- term + alpha.cov.term*focal.comp[i.comp]
      }
      
      full.data$fitness[i.obs] <- as.numeric(full.data$focal.lambda[i.obs] * (num) / term) 

    }# for i.obs
  }else{ # if no generative model
    full.data$fitness <- 0
  }
  
  full.data[,!names(full.data) == "focal.lambda"]
  
}
