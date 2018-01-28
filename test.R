#generate fake data

n <- 500*3 #multiple of 3
focal = c(rep("sp1", n/3), 
          rep("sp2", n/3),
          rep("sp3", n/3))
fitness = sort(runif(n, 0,10), decreasing = TRUE)
d <- data.frame(focal = focal,
                fitness = fitness,
                sp1 = round(sort(runif(n, 0,6))),
                sp2 = round(sort(runif(n, 0,6))),
                sp3 = round(sort(runif(n, 0,6))),
                cov1 = round((runif(n, 0,10))),
                cov2 = round((runif(n, 0,10))))
d
comp_matrix <- d[,3:5]
covariates <- d[,6:7]
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



