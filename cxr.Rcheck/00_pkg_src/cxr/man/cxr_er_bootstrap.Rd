% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cxr_er_bootstrap.R
\name{cxr_er_bootstrap}
\alias{cxr_er_bootstrap}
\title{standard error estimates for effect and response parameters}
\usage{
cxr_er_bootstrap(effect.response.model, optim.method, sp.data, init.par,
  lower.bounds, upper.bounds, covariates, optimize.lambda, lambda.vector,
  nsamples)
}
\arguments{
\item{effect.response.model}{effect/response function}

\item{optim.method}{optimization method to use. One of the following: "optim_NM","optim_L-BFGS-B","nloptr_CRS2_LM", 
"nloptr_ISRES","nloptr_DIRECT_L_RAND","GenSA","hydroPSO","DEoptimR".}

\item{sp.data}{dataframe with all the necessary information in long format. It should have the following columns:
- site: character ID
- focal: character ID of the focal species. Any number of focal species is allowed, but the number of focal species
must match the number of initial parameters (one lambda, e, and r per species).
- fitness: numeric, a fitness metric
- competitor: character, ID of a competitor for that observation. The set of competitors must be, for now, the same
as the set of focal species.
- number: number of neighbouring/competitor individuals from the associated species. Observations without competitors of a given species
must be explicit, i.e. setting number to zero.}

\item{init.par}{1d vector of initial parameters}

\item{lower.bounds}{1d vector of lower bounds}

\item{upper.bounds}{1d vector of upper bounds}

\item{covariates}{dataframe/matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
from an observation.}

\item{optimize.lambda}{boolean, whether to optimize the values of lambda or not.}

\item{lambda.vector}{in case lambda is not to be optimized, fixed values for it.}

\item{nsamples}{how many bootstrap samples to compute.}
}
\value{
1d vector, the standard error of each parameter in init.par
}
\description{
Computes bootstrap standard errors for a given effect/response function
}
