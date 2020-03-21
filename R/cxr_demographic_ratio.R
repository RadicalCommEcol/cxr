#' competitive response ratio
#' 
#' according to eq. 4 of Godoy et al. (2014)
#' 
#' @param model_family acronym 
#' @param pair_lambdas vector of length 2
#'
#' @return single numeric value
#' @noRd
#'
cxr_demographic_ratio <- function(model_family, pair_lambdas){
  dr_fun <- paste(model_family,"_demographic_ratio",sep="")
  
  dr_model <- try(get(dr_fun),silent = TRUE)
  
  if(class(dr_model) == "try-error"){
    message(paste("demographic_ratio: function '",dr_fun,"' could not be retrieved, which means
    that you have not defined a specific demographic ratio formulation for your custom model family. 
    Demographic ratio will be calculated as the ratio between the lambda values provided. Be aware
    that this may be incorrect depending on your underlying model.\n"
                  ,sep=""))
    return(pair_lambdas[1]/pair_lambdas[2])
  }else{
    dr_model(pair_lambdas)
  }
  
}
