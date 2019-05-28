#generate fake data

#NOTE: if lambda of m1 is always better, make m3 to use m1 as starting point?

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

n <- 500*3 #multiple of 3
focal = c(rep("sp1", n/3)#, 
          #rep("sp2", n/3),
          #rep("sp3", n/3)
          )
d <- data.frame(focal = rep(focal,3),
                sp1 = c(round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5)))), 
                sp2 = c(round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5)))), 
                sp3 = c(round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5))), round((rpois(n/3, 0.5)))),
                cov1 = c(round((runif(n/3, 0, 6))), round((runif(n/3, 0, 6))), round((runif(n/3, 0, 6)))) 
#               , cov2 = c(round(sort(runif(n/3, 0, 6))), round((runif(n/3, 0, 6))), round((runif(n/3, 0, 6)))) 
) 
head(d) 
comp_matrix <- d[,2:4]
covariates <- d[,5, drop = FALSE]
d$fitness <- fitness_calc(lambda, d$sp1, d$sp2, d$sp3, d$cov1)
hist(d$fitness)
scatter.smooth(d$fitness ~ (d$sp1+ d$sp2+ d$sp3), col = d$focal) #model 1 and 2
scatter.smooth(d$fitness ~ d$sp1, col = d$focal) #model 3 alphas 
scatter.smooth(d$fitness ~ d$sp2, col = d$focal) #model 3 alphas 
scatter.smooth(d$fitness ~ d$sp3, col = d$focal) #model 3 alphas  
scatter.smooth(d$fitness ~ d$cov1, col = d$focal) #l_cov1 
#scatter.smooth(d$fitness ~ d$cov2, col = d$focal) #l_cov1 0
tapply(d$fitness, d$focal, max)
tapply(d$fitness, d$focal, mean)

#simple model, no cov's
test <- compete(focal = d$focal, 
        fitness = d$fitness, 
        comp_matrix = comp_matrix
        , lower1 = c(1,0.0001)
        , lower2 = c(1,0,0.0001)
        , lower3 = c(1, rep(0, times=ncol(comp_matrix)),0.0000000001)
        , hessian3 = TRUE) #works

compete(focal, d$fitness, comp_matrix, covariates
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
        , hessian5 = TRUE) #Lapack routine dgesv: system is exactly singular: U[10,10] = 0 


#lambda's are estimated too low and errors go up with model number, which is odd.
#alpha2 & 3 are not detected!
# cov effects on lambda

scatter.smooth(d$fitness ~ (d$sp1+ d$sp2+ d$sp3), col = d$focal) #model 1 and 2
curve(expr = test[[1]]$lambda_est2
      / (1+ (test[[1]]$alpha_est2)
         *x), add=T, lwd=3, col="red")



