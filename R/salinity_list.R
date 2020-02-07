#' Salinity measurements
#' 
#' A list containing salinity values associated to the data from 'neigh_list'.
#' The list has 16 elements, one for each focal species considered.
#' Each element of the list is a dataframe with 2 columns:
#' 
#' \itemize{          
#'   \item obs_ID: unique identifier of each observation
#'   \item salinity: salinity measurement for that observation,
#'   in accumulated microsiemens/m2
#'   }
#'   
#' @note For details, see Lanuza et al. 2018 Ecology Letters. 
#' @docType data
#' @keywords datasets
#' @name salinity_list
#' @usage data(salinity_list)
#' @format A list with 16 elements, each of which a dataframe
#' of variable number of rows and 2 numeric columns
"salinity_list"