#generate fake data

n <- 100*3 #multiple of 3
focal = c(rep("sp1", n/3), 
          rep("sp2", n/3),
          rep("sp3", n/3))
fitness = sort(runif(n, 0,10), decreasing = TRUE)
d <- data.frame(focal = focal,
                fitness = fitness,
                sp1 = round(sort(runif(n, 0,3))),
                sp2 = round(sort(runif(n, 0,3))),
                sp3 = round(sort(runif(n, 0,3))))
d
comp_matrix <- d[,3:5]

compete(focal, fitness, comp_matrix
        , lower2 = c(1,0,0.0001)
        , lower3 = c(1, rep(0, times=length(unique(focal))),0.0000000001) 
        , hessian3 = TRUE)
