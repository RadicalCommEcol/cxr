# test of pairwise niche overlap functions
# premise: MCT formulation does not allow facilitation coefficients
# in that case, the structural analogue may be used

# for competitive interactions, MCT and structural approaches 
# should be equivalent

# MCT function
niche_overlap_MCT <- function(pair_matrix){
  sqrt((pair_matrix[1,2] * pair_matrix[2,1])/(pair_matrix[1,1] * pair_matrix[2,2]))
}

# saavedra et al. 2017 eq. 7
niche_diff_SA1 <- function(pair_matrix){
  (2/pi) * asin((pair_matrix[1,1] * pair_matrix[2,2] - pair_matrix[2,1]*pair_matrix[1,2])/
                  (sqrt(pair_matrix[1,1]^2 + pair_matrix[2,1]^2)*sqrt(pair_matrix[1,2]^2 + pair_matrix[2,2]^2))) 
}

# saavedra et al. 2017 - supp.mat code
niche_diff_SA2 <- function(pair_matrix){
  n <- nrow(pair_matrix)
  Sigma <-solve(t(pair_matrix) %*% pair_matrix, tol = 1e-40)
  d <- mvtnorm::pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

# song et al. 2018 A guideline... - supp.mat code
niche_diff_SA3 <- function(pair_matrix) {
  S <- nrow(pair_matrix)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- mvtnorm::pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    return(out)
  }
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (f(pair_matrix) == FALSE) {
    return(0)
  }
  else {
    Sigma <- solve(t(pair_matrix) %*% pair_matrix)
    return(omega(S, Sigma))
  }
}

niche_overlap_SA1 <- function(pair_matrix){1 - niche_diff_SA1(pair_matrix)}
niche_overlap_SA2 <- function(pair_matrix){1 - 10^niche_diff_SA2(pair_matrix)}
niche_overlap_SA3 <- function(pair_matrix){1 - niche_diff_SA3(pair_matrix)}

# test --------------------------------------------------------------------

set.seed(42)
# make sure intra>inter
posmatrix <- matrix(c(runif(1,0.5,1),runif(2,0,0.5),runif(1,0.5,1)),nrow = 2)

# MCT formulation
niche_overlap_MCT(posmatrix)
# saavedra et al. 2017 eq. 7
niche_overlap_SA1(posmatrix)
# saavedra et al. 2017 supp.mat
niche_overlap_SA2(posmatrix)
# song et al. 2018
niche_overlap_SA3(posmatrix)

