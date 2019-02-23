
set.seed(123)

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
