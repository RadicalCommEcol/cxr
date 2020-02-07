#' Check the input data for cxr_er_fit
#'
#' Several requirements:
#' 1-equal nrows for all focal sp
#' 2-same for covariates
#' 3-all dataframes same columns, all numeric
#' 4-same set of focals and neighbour sp
#'
#' @param data as in cxr_er_fit
#' @param covariates as in cxr_er_fit
#'
#' @return logical, whether conditions are met
#'
#' @noRd
cxr_check_input_er <- function(data,covariates = NULL){
  data.ok <- TRUE
  if("list" %in% class(data)){
    # check nrow
    nrows <- unlist(lapply(data,nrow))
    nna <- sum(unlist(lapply(data,function(x){sum(is.na(x))})))
    for(i.sp in 1:length(data)){
      classes <- sapply(data[[i.sp]],class)
      if(any(!classes %in% c("integer","numeric"))){
        data.ok <- FALSE
      }
    }    
    if(nna > 0 | length(unique(nrows))>1){
      data.ok <- FALSE
    }else{
      # check cols
      mynames <- names(data[[1]])
      namesok <- all(unlist(lapply(data,function(x){identical(names(x),mynames)})))
      focalneigh <- identical(names(nrows),mynames[which(mynames != "fitness")])
      if(!namesok | length(mynames)<3 | !focalneigh){
        data.ok <- FALSE
      }else{
        if(!is.null(covariates)){
          if(!class(covariates) == "list"){
            data.ok <- FALSE
          }else{
            nrcov <- unlist(lapply(covariates,nrow))
            if(nrcov[1] != nrows[1] | length(unique(nrcov))>1){
              data.ok <- FALSE
            }
          }# if-else covariates list
        }# covariates
      }# if-else cols ok
    }# if-else length ok
  }else if(class(data) == "data.frame"){
    if(sum(is.na(data))>0){
      data.ok <- FALSE
    }
    classes <- sapply(data,class)
    if(any(!classes %in% c("integer","numeric"))){
      data.ok <- FALSE
    }
    if(!c("focal") %in% names(data) | !c("fitness") %in% names(data)){
      data.ok <- FALSE
    }else{
      obsp <- table(data$focal)
      if(length(unique(obsp))>1){
        data.ok <- FALSE
      }else{
        focalneigh <- identical(unique(data$focal),names(data)[which(!names(data) %in% c("fitness","focal"))])
        if(ncol(data)<3 | !focalneigh){
          data.ok <- FALSE
        }else{
          if(!is.null(covariates)){
            if(!is.data.frame(covariates)){
              data.ok <- FALSE
            }else{
              if(nrow(covariates) != nrow(data)){
                data.ok <- FALSE
              }# if nrows data and cov
            }# if-else covariates
          }# if covariates
        }# if-else colnums ok
      }# if-else num obs ok
      
    }# if-else names ok
  }else{
    data.ok <- FALSE
  }# if-else list or df
  
  data.ok
  
}