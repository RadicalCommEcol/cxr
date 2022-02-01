context("metapopulation functions")
# test metapopulation functions

# define species names
sp <- c("s1","s2","s3")
# number of species
num.sp <- length(sp)

# define site names
sites <- c("sa","sb")
# number of sites
num.sites <- length(sites)

# number of demographic stages - this should be always fixed
num.stages <- 3

# vital rates to account for - these names should be fixed
rates <- c("Sj","Sn","Sr","Rn","Rr","D","Ds","O")

# environmental forcing?
env <- TRUE

# this builds the data structure - still empty
param <- build_param(sp = sp,
                     sites = sites,
                     rates = rates,
                     env = env, 
                     num.params = 6)

years <- 100

# -------------------------------------------------------------------------

param <- generate_vital_rate_coefs(param = param,
                                   sp = "s1",
                                   # sites = "sa",
                                   vital.rate = c("Sj"),#,"Sn","Sr","Rn","Rr","O"),
                                   vr.coef = "alpha",
                                   mean.coef = .1,sd.coef = 0)

p1 <- param[["s1"]][["sa"]]["Sj","alpha"]

test_that("vital rate coefficients are correctly generated", {
  
  expect_is(p1, "numeric")
  
})

# -------------------------------------------------------------------------

vpm <- vec_permutation_matrices(num.sp,num.sites,num.stages)

env <- rnorm(years, mean=0, sd=1)

initial.densities <- list()
# sp1
initial.densities[[1]] <- matrix(c(20,20,15,20,20,13),
                                 nrow = num.sites,
                                 ncol = num.stages,
                                 byrow = TRUE)
# sp2
initial.densities[[2]] <- matrix(c(25,15,3,25,15,3),
                                 nrow = num.sites,
                                 ncol = num.stages,
                                 byrow = TRUE)
# sp2
initial.densities[[3]] <- matrix(c(5,5,2,3,4,2),
                                 nrow = num.sites,
                                 ncol = num.stages,
                                 byrow = TRUE)
current.densities <- initial.densities

# for(i.year in 1:years){
  
  # this is a list per species and site
  transition_matrices <- list()
  
  for(i.sp in 1:length(sp)){
    transition_matrices[[i.sp]] <- list()
    for(i.site in 1:length(sites)){
      
      # store the transition matrix for this sp and site
      transition_matrices[[i.sp]][[i.site]] <- fill_transition_matrix(focal.sp = i.sp,
                                                                      site = i.site,
                                                                      param = param,
                                                                      env = env[1],
                                                                      current.densities = current.densities)
      
    }# for each site
    names(transition_matrices[[i.sp]]) <- sites
  }# for each species
  names(transition_matrices) <- sp
  
  # update demography and dispersal matrices --------------------------------
  
  for(i.sp in 1:length(sp)){
    
    # demography
    vpm[["demography"]][[i.sp]] <- fill_demography_matrix(focal.sp = i.sp,
                                                          vpm = vpm,
                                                          transition_matrices = transition_matrices)
    # dispersal
    vpm[["dispersal"]][[i.sp]] <- fill_dispersal_matrix(focal.sp = i.sp,
                                                        num.sites = num.sites,
                                                        param = param,
                                                        vpm = vpm,
                                                        env = env[1],
                                                        current.densities = current.densities) 
  }# for i.sp
  
# }# for i.year

test_that("transition matrices are correctly generated", {
  
  expect_is(transition_matrices[["s1"]][["sa"]], "matrix")
  
})

test_that("demography matrices are correctly generated", {
  
  expect_is(vpm[["demography"]][[1]], "matrix")
  expect_is(vpm[["dispersal"]][[1]], "matrix")
})

# -------------------------------------------------------------------------

dens <- calculate_densities(focal.sp = 1,
                            vpm,
                            current.densities)

test_that("densities are correctly generated", {
  
  expect_is(dens, "matrix")
  
})

