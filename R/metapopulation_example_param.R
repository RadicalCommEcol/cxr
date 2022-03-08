#' Metapopulation dynamics coefficients
#' 
#' A nested list containing vital rate coefficients for projecting
#' metapopulation dynamics.
#' The first level of the list has 3 elements, one for each species modelled.
#' The second level of the list has 2 elements, one for each site modelled.
#' For each combination species-site, there is a data.frame of eight rows - 
#' one per each vital rate, and eight columns - one per coefficient, that 
#' correspond to the coefficients of a GLM. These are named as `alpha`,`beta1`, 
#' etc, in the data.frame, and correspond to the intercept, environmental effect, 
#' effects of each of the three species' density, and environment:density interactions
#'   
#' @docType data
#' @keywords datasets
#' @name metapopulation_example_param
#' @usage data(metapopulation_example_param)
#' @format A nested list with 3x2 elements, each of which a dataframe
#' of 8 rows and 8 numeric columns
#' @md
"metapopulation_example_param"