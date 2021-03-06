% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_BH5.R
\name{model_BH5}
\alias{model_BH5}
\title{Title Beverton-Holt fecundity, fifth model}
\usage{
model_BH5(par, param.list, log.fitness, focal.comp.matrix, num.covariates,
  num.competitors, focal.covariates, fixed.terms)
}
\arguments{
\item{par}{vector of variable length, with the following order: first, lambda of focal sp; 
second lambda.cov, effects of every covariate on lambda; 
third alpha, interaction coefficients with every species; 
fourth alpha.cov, effects of every covariate on alpha values (varying effect of every covariate over every interaction coefficient); 
last, sigma value. If any element is not to be optimized, it must not be present in this vector, but rather in the \code{fixed.terms} list}

\item{param.list}{string listing parameters to optimize. Possible elements are \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}.}

\item{log.fitness}{log of fitness value}

\item{focal.comp.matrix}{dataframe with as many rows as observations, and one column for each competitor sp. 
Values of the dataframe are number of competitors of each sp per observation.}

\item{num.covariates}{number of covariates}

\item{num.competitors}{number of competitor species}

\item{focal.covariates}{dataframe with as many rows as observationes, and one column for each covariate.
Values of the dataframe are covariate values for every observation.}

\item{fixed.terms}{list with elements \code{lambda}, \code{lambda.cov}, \code{alpha}, \code{alpha.cov}. It contains parameters not to be optimized.
Each element of the list must be of its appropriate length. Note that adding an element in "param.list" will force the function
to look for it in \code{par}, and will not consider it here. In this model, \code{lambda.cov} and \code{alpha.cov} are not considered.}
}
\value{
log-likelihood value
}
\description{
These functions return the negative log-likelihood of the data
given the model and parameters. model_BH5 is 
\eqn{F_i = \frac{\lambda_i + \sum_k b_k \Theta_{i,k}}{1+\sum_j(\alpha_{ij} +\sum_k b_k \phi_{i,j,k})N_{j}}}{(lambda_i + sum(b_k*theta_k))/(1 + sum(alpha_ij + sum(b_k*phi_ijk))N_j)}
}
