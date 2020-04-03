#' Internal, retrieve parameters from the vector returned by the optimization procedures
#'
#' @param optim_par 1d vector, the result of an optimization method
#' @param lambda_length either 0 (lambda not fit) or 1
#' @param alpha_intra_length either 0 (alpha_intra not fit) or 1
#' @param alpha_inter_length either 0 (alpha_inter not fit) or a positive number
#' @param lambda_cov_length either 0 (lambda_cov not fit) or a positive number
#' @param alpha_cov_length either 0 (alpha_cov not fit) or a positive number
#' @param empty_neigh optional character vector with the names of columns without
#' observations. these should be added with coefficients "NA".
#' @param alpha_cov_form whether and how to include covariates on alpha, for the case
#' in which empty neighbours are added with NA values.
#' @param error_par flag to know if the function is to retrieve error terms or normal fits
#' @return list with elements "lambda", "alpha", "lambda_cov", "alpha_cov", "sigma". 
#' If one of these elements is not present, returns NULL.
#' @noRd
cxr_retrieve_params <- function(optim_par, 
                                lambda_length = 0, 
                                alpha_intra_length = 0, 
                                alpha_inter_length = 0, 
                                lambda_cov_length = 0, 
                                alpha_cov_length = 0,
                                empty_neigh = NULL,
                                alpha_cov_form, error_par = FALSE){
  
  lambda <- NULL
  alpha_intra <- NULL
  alpha_inter <- NULL
  lambda_cov <- NULL
  alpha_cov <- NULL
  sigma <- optim_par[length(optim_par)]
  
  pos <- 1
  
  if(lambda_length == 1){
    lambda <- optim_par[pos]
    pos <- pos + 1
  }
  
  if(lambda_cov_length > 0){
    lambda_cov <- optim_par[pos:(pos+lambda_cov_length-1)]
    pos <- pos + lambda_cov_length
  }
  
  if(alpha_intra_length == 1){
    alpha_intra <- optim_par[pos]
    pos <- pos + 1
  }
  
  if(alpha_inter_length > 0){
    alpha_inter <- optim_par[pos:(pos+alpha_inter_length-1)]
    pos <- pos + alpha_inter_length
    
    # if alpha is pairwise, and there are empty neighbours
    # fill their slots with NA
    if(alpha_inter_length>1){
      if(!is.null(empty_neigh)){
        emptyv <- rep(NA_real_,length(empty_neigh))
        names(emptyv) <- empty_neigh
        if(error_par){
          names(emptyv) <- paste(names(emptyv),"_se",sep="")
        }
        alpha_inter <- c(alpha_inter,emptyv)
        alpha_inter <- alpha_inter[sort(names(alpha_inter))]
      }
    }
  }
  
  if(alpha_cov_length > 0){
    alpha_cov <- optim_par[pos:(pos+alpha_cov_length-1)]
    
    if(!is.null(empty_neigh)){
      if(alpha_cov_form == "pairwise"){
        # the only case in which this is needed
        
        # retrieve covariate names
        # i could also pass them as an argument
        if(error_par){
          cnames <- gsub("lambda_cov_","",names(lambda_cov))
          cnames <- gsub("_se","",cnames)
            # substr(names(lambda_cov),12,nchar(names(lambda_cov)))
        }else{
          cnames <- gsub("lambda_cov_","",names(lambda_cov))
            #substr(names(lambda_cov),12,nchar(names(lambda_cov)))
        }
        
        # beware a single paste command wont work when cnames > 1 AND
        # empty_neigh > 1
        
        # first, the covariates
        acname <- paste("alpha_cov_",cnames,sep="")
        # second, the species
        acname <- as.character(sapply(acname,function(x){paste(x,"_",empty_neigh,sep="")}))
        
        # if retrieving error terms, add it
        if(error_par){
          acname <- paste(acname,"_se",sep="")
        }
        # append them
        emptyac <- rep(NA_real_,length(acname))
        names(emptyac) <- acname
        
        alpha_cov <- c(alpha_cov,emptyac)
        alpha_cov <- alpha_cov[sort(names(alpha_cov))]
        
      }
    }
    
  }
  
  return(list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter, lambda_cov = lambda_cov, alpha_cov = alpha_cov, sigma = sigma))
  
}

