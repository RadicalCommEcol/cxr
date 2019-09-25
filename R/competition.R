#' Competition measurements
#' 
#' A dataset containing fitness and neighbours for each plant individual
#' 
#' \itemize{          
#'   \item year: year
#'   \item month: month
#'   \item day: day
#'   \item plot: plot
#'   \item subplot: subplot code
#'   \item focal: focal plant species identity
#'   \item individual_ID: unique identifier for each individual focal plant
#'   \item fruit: total fruits produced by the focal individual
#'   \item seed: total seeds produced by the focal individual
#'   \item competitor: competitor identity
#'   \item number: number of competitors
#'   }
#'   
#' @note For details, see Lanuza et al. 2018 Ecology Letters. 
#' @docType data
#' @keywords datasets
#' @name competition
#' @usage data(competition)
#' @format A data frame with 57432 rows and 11 variables
"competition"