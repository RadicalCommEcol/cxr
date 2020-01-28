#' internal checks for data consistency
#'
#' @param data dataframe with at least a column "fitness"
#' @param covariates optional dataframe/matrix with covariates in columns
#'
#' @return boolean value, whether input data are valid or not
#' @export
#'
#' @examples
cxr_check_input_data <- function(data,covariates = NULL){
  data.ok <- TRUE
  if(names(data)[1] != "fitness"){
    data.ok <- FALSE
  } 
  if(ncol(data)<2){
    data.ok <- FALSE
  }
  if(!is.null(covariates)){
    if(nrow(covariates) != nrow(data)){
      data.ok <- FALSE
    }
  }
  if(sum(is.na(data)) > 0 | sum(is.na(covariates))>0){
    data.ok <- FALSE
  }
  data.ok
}