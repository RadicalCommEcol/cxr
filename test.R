#generate fake data

n <- 500*3 #multiple of 3
focal = c(rep("sp1", n/3), 
          rep("sp2", n/3),
          rep("sp3", n/3))
fitness = c(sort(runif(n/3, 3,10), decreasing = TRUE), sort(runif(n/3, 5,13), decreasing = TRUE), sort(runif(n/3, 0,8), decreasing = TRUE))
d <- data.frame(focal = focal,
                fitness = fitness,
                sp1 = c(round(sort(rpois(n/3, 3))), round(sort(rpois(n/3, 3))), round(sort(rpois(n/3, 3)))), #high compe with all
                sp2 = c(round(sort(rpois(n/3, 3))), round((rpois(n/3, 3))), round((rpois(n/3, 3)))), #high intra, low inter
                sp3 = c(round((rpois(n/3, 3))), round(sort(rpois(n/3, 3))), round(sort(rpois(n/3, 3)))), #low intra, high compe
                cov1 = c(round((runif(n/3, 0, 6))), round(sort(runif(n/3, 0, 6))), round(sort(runif(n/3, 0, 6)))), #high impact on sp 2 $ 3
                cov2 = c(round(sort(runif(n/3, 0, 6))), round((runif(n/3, 0, 6))), round((runif(n/3, 0, 6))))) #only on 1
d
comp_matrix <- d[,3:5]
covariates <- d[,6:7]
scatter.smooth(d$fitness ~ (d$sp1+ d$sp2+ d$sp3), col = d$focal) #model 1 and 2
scatter.smooth(d$fitness ~ d$sp1, col = d$focal) #model 3 alphas - - - 
scatter.smooth(d$fitness ~ d$sp2, col = d$focal) #model 3 alphas - 0 0
scatter.smooth(d$fitness ~ d$sp3, col = d$focal) #model 3 alphas 0 - - 
scatter.smooth(d$fitness ~ d$cov1, col = d$focal) #l_cov1 - (at least moderatelly)
scatter.smooth(d$fitness ~ d$cov2, col = d$focal) #l_cov1 0
tapply(d$fitness, d$focal, max)
tapply(d$fitness, d$focal, mean)

#compete(focal, fitness, comp_matrix
#       , lower1 = c(1,0.0001)
#        , lower2 = c(1,0,0.0001)
#        , lower3 = c(1, rep(0, times=length(unique(focal))),0.0000000001)
#        , hessian3 = TRUE) #works

compete(focal, fitness, comp_matrix, covariates
        , lower1 = c(1,0.0001)
        , lower2 = c(1,0,0.0001)
        , lower3 = c(1, rep(0, times=length(unique(focal))),0.0000000001)
        , lower4 = c(1, rep(0, times=ncol(covariates)), #n_cov
                     rep(0, times=length(unique(focal))), #alfas
                     rep(0, times=ncol(covariates)), #n_cov
                     0.0000000001)
        , lower5 = c(1, rep(0, times=ncol(covariates)), #n_cov
                     rep(0.001, times=length(unique(focal))), #alfas
                     rep(0.001, times=(ncol(covariates)*length(unique(focal)))), #n_cov*n_bg
                     0.0000000001)
        , hessian5 = FALSE) #Lapack routine dgesv: system is exactly singular: U[10,10] = 0 


#lambda's are estimated too low and errors go up with model number, which is odd.
#alpha2 & 3 are not detected!
# cov effects on lambda
