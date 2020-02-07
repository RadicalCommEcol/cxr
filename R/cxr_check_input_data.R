#' Internal checks for data consistency
#'
#' checks whether data have a column 'fitness' and at least one more column.
#' All columns must be numeric or integer.
#'
#' @param data dataframe with at least a column "fitness"
#' @param covariates optional dataframe/matrix with covariates in columns
#'
#' @return boolean value, whether input data are valid or not
#'
#' @noRd
cxr_check_input_data <- function(data,covariates = NULL){
  data.ok <- TRUE
  if(names(data)[1] != "fitness"){
    data.ok <- FALSE
  } 
  if(ncol(data)<2){
    data.ok <- FALSE
  }else{
    classes <- sapply(data,class)
    if(any(!classes %in% c("integer","numeric"))){
      data.ok <- FALSE
    }
  }
  if(!is.null(covariates)){
    if(nrow(covariates) != nrow(data)){
      data.ok <- FALSE
    }
    cclasses <- sapply(covariates,class)
    if(any(!cclasses %in% c("integer","numeric"))){
      data.ok <- FALSE
    }
  }
  if(sum(is.na(data)) > 0 | sum(is.na(covariates))>0){
    data.ok <- FALSE
  }

  data.ok
}