#' Niche overlap between two species
#'
#' quoting Godoy et al. (2014):
#' reflects the average degree to which species limit individuals of their own species relative to competitors. 
#' Low niche overlap causes species to have greater per capita growth rates when rare than when common. 
#' If species limit individuals of their own species and their competitors equally, then niche overlap is 1, 
#' and coexistence is not possible unless species are otherwise identical. 
#' At the other extreme, if species have no interspecific effects, then niche overlap is 0.
#' 
#' Niche overlap has a common functional form for a series of models, including those specified in table A1 of
#' Hart et al. (2018) How to quantify competitive ability. Journal of Ecology 106, 1902-1909. 
#' Other model families may not adhere to the general definition. 
#' 
#' This function will calculate niche overlap, using either the general definition or user-provided formulae.
#' It accepts either two 'cxr_pm_fit' objects calculated with the 'cxr_pm_fit' function,
#' or a 2x2 interaction matrix. In the first case, it checks 
#' whether the model family from which the parameters were calculated has an associated niche overlap formula. 
#' In the second case, it uses the general formula from Godoy et al. (2014).
#' Again, beware that this formula may not be appropriate if the interaction matrix provided is inferred
#' from specific models.
#'  
#' 
#' @param cxr_sp1 cxr_pm_fit object giving the parameters from the first species.
#' @param cxr_sp2 cxr_pm_fit object giving the parameters from the second species. 
#' @param pair_matrix optional 2x2 matrix with intra and interspecific interaction 
#' coefficients between the two species.
#'
#' @return niche overlap value, in the range 0-1.
#' @export
niche_overlap <- function(cxr_sp1 = NULL, cxr_sp2 = NULL, pair_matrix = NULL){
  
  res <- NULL
  
  if(!is.null(cxr_sp1) & !is.null(cxr_sp2)){
    if(!is.null(pair_matrix)){
      message("cxr niche_overlap:both cxr objects and a pairwise matrix were specified. 
              Pairwise matrix will be discarded")
    }
    
    if(cxr_sp1$model_family == cxr_sp2$model_family)
    
    nov_fun <- paste(cxr_sp1$model_family,"_niche_overlap")
    nov_model <- try(get(nov_fun),silent = TRUE)
    
    if(class(nov_model) == "try-error"){
      message(paste("niche_overlap ERROR: function '",nov_fun,"' could not be retrieved. 
      Make sure it is defined and available in the cxr package or in the global environment.\n"
                    ,sep=""))
      return(NULL)
    }else{
      
      sp1 <- cxr_sp1$focal_ID
      sp2 <- cxr_sp2$focal_ID
      
      if(is.null(sp1) | is.null(sp2)){
        message(paste("niche_overlap ERROR: a 'cxr_pm_fit' object passed does not contain information for identifying the 
                      focal taxa, i.e. does not have a valid 'focal_ID' field. 
                      This is probably because the fit was done without specifying 'alpha_intra' and 
                      'focal_column' in the function 'cxr_pm_fit'. Without this information on intra and inter-specific
                      interactions, niche overlap cannot be computed."))
        return(NULL)
      }else{
        intra_sp1 <- cxr_sp1$alpha_intra
        intra_sp2 <- cxr_sp2$alpha_intra
        inter_sp1_sp2 <- cxr_sp1$alpha_inter[which(names(cxr_sp1$alpha_inter) == sp2)]
        inter_sp2_sp1 <- cxr_sp1$alpha_inter[which(names(cxr_sp2$alpha_inter) == sp1)]
        
        if(is.null(intra_sp1) | is.null(intra_sp2) | is.null(inter_sp2_sp1) | is.null(inter_sp1_sp2)){
          message(paste("niche_overlap ERROR: pairwise interactions could not be retrieved from the 'cxr_pm_fit' objects. This could
                        be because terms in the 'alpha_inter' fields do not match the names of the other taxa, 
                        or because the 'alpha_intra' fields are NULL."))
          return(NULL)
        }else{
          nov_matrix <- matrix(intra_sp1,inter_sp2_sp1,inter_sp1_sp2,intra_sp2,nrow = 2)
          res <- nov_model(nov_matrix)
        }# if-else null
        
      }# if-else null

    }# if-else model found
  }else if(!is.null(pair_matrix)){
    
    message("Niche overlap calculated with the standard formula. You should be sure that
            this formula is applicable to the model with which the matrix coefficients
            were calculated.")
    
    res <- sqrt((pair_matrix[1,2]/pair_matrix[2,2])*(pair_matrix[2,1]/pair_matrix[1,1]))
    
  }
  
  # get model family
  
  
  
  res
  
}