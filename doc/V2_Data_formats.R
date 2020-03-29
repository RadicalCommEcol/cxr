## ----setup,echo=FALSE----------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## ------------------------------------------------------------------------
library(cxr)
data("neigh_list", package = "cxr")

## ------------------------------------------------------------------------
?neigh_list
names(neigh_list)
head(neigh_list[[1]])

## ------------------------------------------------------------------------
data("abundance", package = "cxr")
head(abundance)

## ------------------------------------------------------------------------
data("spatial_sampling")
names(spatial_sampling)
head(spatial_sampling[["BEMA"]])

## ------------------------------------------------------------------------
data("species_rates", package = "cxr")
species_rates 

## ------------------------------------------------------------------------
data("salinity_list", package = "cxr")
names(salinity_list)
head(salinity_list[[1]])

