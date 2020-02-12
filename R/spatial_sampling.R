#' spatial arrangement of the observations
#' 
#' A dataset giving the spatial arrangement of observations.
#' The dataset is a list of 16 elements following the structure of
#' 'neigh_list'. Each list component is a dataframe
#' with columns:
#' 
#' \itemize{          
#'   \item obs_ID: unique identifier for each observation
#'   \item plot: one of 9 plots of 8.5 x 8.5 m
#'   \item subplot: one of 36 subplots of 1x1 m within each plot
#'   }
#'   
#' @note For details, see Lanuza et al. 2018 Ecology Letters. 
#' @docType data
#' @keywords datasets
#' @name spatial_sampling
#' @usage data(spatial_sampling)
#' @format A list with 16 elements, each of which a dataframe
#' of variable number of rows and 18 columns
"spatial_sampling"