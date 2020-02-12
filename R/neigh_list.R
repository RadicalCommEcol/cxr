#' neighbours and fitness observations
#' 
#' A dataset containing fitness and neighbours 
#' for plant individuals of 16 species. 
#' The dataset is a named list with 16 elements,
#' each of which is a dataframe with the following columns:
#' * obs_ID: unique identifier for each observation
#' * fitness: number of viable seeds of the focal individual
#' * 16 columns indicating the number of neighbours from each plant sp.
#'   in a radius of 7.5 cm from the focal individual
#'   
#' @note For details, see Lanuza et al. 2018 Ecology Letters. 
#' @docType data
#' @keywords datasets
#' @name neigh_list
#' @usage data(neigh_list)
#' @format A list with 16 elements, each of which a dataframe
#' of variable number of rows and 18 columns
#' @md 
"neigh_list"