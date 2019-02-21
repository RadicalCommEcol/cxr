#true parameters I want to recover
lambda = 300 #lambda
l_cov1 = 1 #l_cov
a_sp1 = 0.2 #alphas
a_sp2 = 1
a_sp3 = 2
a_cov1 = 0 #a_cov's
a_cov2 = 0
a_cov3 = 0

fitness_calc <- function(lambda, sp1, sp2, sp3, cov1){
  (lambda * (1+ l_cov1 * cov1)) / 
    (1 + ((a_sp1 + (a_cov1*cov1))*sp1) + 
       ((a_sp2 + (a_cov2*cov1))*sp2) +
       ((a_sp3 + (a_cov3*cov1))*sp3))  
}

# community data

n <- 500*3 #multiple of 3

# one or three focal species
focal.sp = c(rep("sp1", n/3))#, 
          # rep("sp2", n/3),
          # rep("sp3", n/3))

# data with one or two covariates
test.data <- data.frame(
  #focal = focal.sp,
  focal = rep(focal.sp,3),
  sp1 = c(round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5)))), 
  sp2 = c(round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5)))), 
  sp3 = c(round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5)))),
  cov1 = c(round((runif(n/3, 0, 6))), round((runif(n/3, 0, 6))), round((runif(n/3, 0, 6)))) 
  #               , cov2 = c(round(sort(runif(n/3, 0, 6))), round((runif(n/3, 0, 6))), round((runif(n/3, 0, 6)))) 
) 

# number of competitors at each observation
comp_matrix <- test.data[,2:4]

# covariates matrix
covariates <- test.data[,5, drop = FALSE]

# fitness metric
test.data$fitness <- fitness_calc(lambda, test.data$sp1, test.data$sp2, test.data$sp3, test.data$cov1)

#####################

# test1 - no effect of covariates on lambdas and alphas (models 1-3)

test1 <- compete(focal = test.data$focal, 
                fitness = test.data$fitness, 
                comp_matrix = comp_matrix
                , lower1 = c(1,0.0001)
                , lower2 = c(1,0,0.0001)
                , lower3 = c(1, rep(0, times=ncol(comp_matrix)),0.0000000001)
                , hessian3 = TRUE) #works

# test2 - full set of models (1-5)

test2 <- compete(focal, test.data$fitness, comp_matrix, covariates
                     , lower1 = c(1,0.0001)
                     , lower2 = c(1,0,0.0001)
                     , lower3 = c(1, rep(0, times=ncol(comp_matrix)),0.0000000001)
                     , lower4 = c(1, rep(0.001, times=ncol(covariates)), #n_cov
                                  rep(0, times=ncol(comp_matrix)), #alfas
                                  rep(0.001, times=ncol(covariates)), #n_cov
                                  0.0000000001)
                     , lower5 = c(1, rep(0.001, times=ncol(covariates)), #n_cov
                                  rep(0, times=ncol(comp_matrix)), #alfas
                                  rep(0.001, times=(ncol(covariates)*ncol(comp_matrix))), #n_cov*n_bg
                                  0.0000000001)
                     , hessian5 = TRUE)

#####################
# diagnostic plots
library(tidyverse)

predicted.values <- test1$`lambdas$co`

# this value cannot be compared to a real one, as it merges all effects into one alpha

alpha2.plot <- ggplot(predicted.values, aes(x = splist, y = alpha_est2)) + 
  geom_errorbar(aes(ymin=alpha.lower.error.2, ymax = alpha.upper.error.2), width = 0.05, size = 0.3) +
  geom_point(size = 2) + 
  ggtitle("alpha from model 2") + 
  NULL
alpha2.plot

##########
# models 3,4,5
# convert to tidy format
alpha.values <- NULL
for(i.model in c(3,4,5)){
  
  if(i.model == 3){
    my.predicted.values <- test2$alpha_matrix3
    my.upper.errors <- test2$alpha_upper_error3
    my.lower.errors <- test2$alpha_lower_error3
  }else if(i.model == 4){
    my.predicted.values <- test2$alpha_matrix4
    my.upper.errors <- test2$alpha_upper_error4
    my.lower.errors <- test2$alpha_lower_error4
  }else if(i.model == 5){
    my.predicted.values <- test2$alpha_matrix5
    my.upper.errors <- test2$alpha_upper_error5
    my.lower.errors <- test2$alpha_lower_error5
  }
  
  predicted.values <- as.data.frame(my.predicted.values)
  predicted.values <- gather(predicted.values,key = "sp",value = "alpha")
  
  upper.errors <- as.data.frame(my.upper.errors)
  upper.errors <- gather(upper.errors,key = "sp",value = "alpha.upper.error")
  
  lower.errors <- as.data.frame(my.lower.errors)
  lower.errors <- gather(lower.errors,key = "sp",value = "alpha.lower.error")
  
  temp.alpha.values <- left_join(predicted.values,upper.errors)
  temp.alpha.values <- left_join(temp.alpha.values,lower.errors)
  
  temp.alpha.values$model <- i.model
  
  alpha.values <- rbind(alpha.values, temp.alpha.values)
}

alpha.values$focal.sp <- 1
alpha.values$model <- as.factor(alpha.values$model)

# append true values
true.values <- data.frame(sp = c("sp1","sp2","sp3"), alpha = c(a_sp1,a_sp2,a_sp3), alpha.upper.error = 0, alpha.lower.error = 0,
                          model = "truth", focal.sp = 1)

alpha.values <- rbind(alpha.values, true.values)

dodge <- position_dodge(.3)
alpha.plot <- ggplot(alpha.values, aes(x = sp, y = alpha, group = model)) + 
  geom_errorbar(aes(ymin=alpha.lower.error, ymax = alpha.upper.error), width = 0.05, size = 0.3, position = dodge) +
  geom_point(aes(color = model), size = 2, position = dodge) + 
  facet_grid(focal.sp~.) + 
  NULL
alpha.plot


