% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_abundBH3.R
\name{model_abundBH3}
\alias{model_abundBH3}
\title{Project abundance of individuals according to the Beverton-Holt third model}
\usage{
model_abundBH3(sp.par, init.abund, cov.values, alpha.matrix,
  lambda.cov.matrix, alpha.cov.matrix, return.seeds = TRUE)
}
\arguments{
\item{sp.par}{dataframe with species in rows, and the following columns:
lambda: fecundity term
germ.rate: seed germination rate
survival.rate: annual survival of ungerminated seed}

\item{init.abund}{number of individuals at time t}

\item{cov.values}{Not used in model_abundBH3}

\item{alpha.matrix}{competition matrix}

\item{lambda.cov.matrix}{Not used in model_abundBH3}

\item{alpha.cov.matrix}{Not used in model_abundBH3}

\item{return.seeds}{boolean flag, whether the prediction should return 
number of seeds (i.e. \eqn{N_{i,t+1}}, eq. 1 of Lanuza et al. 2018), or number of
adult individuals, (i.e. \eqn{N_{i,t+1} * g} )}
}
\value{
1d vector with number of individuals of each species at time t+1
}
\description{
Project abundance of individuals according to the Beverton-Holt third model
}
