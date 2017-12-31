#generate fake data

n <- 21 #multiple of 3
d <- data.frame(focal = c(rep("sp1", n/3), 
                          rep("sp2", n/3),
                          rep("sp3", n/3)),
                reprod = runif(n, 0,10),
                sp1 = round(runif(n, 0,3)),
                sp2 = round(runif(n, 0,3)),
                sp3 = round(runif(n, 0,3)))
d
comp_matrix <- d[,3:5]

compete(focal, reprod, comp_matrix)
