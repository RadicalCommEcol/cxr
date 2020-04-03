#' Niche overlap between two species
#'
#' quoting Godoy et al. (2014):
#' reflects the average degree to which species limit individuals of their own species relative to competitors. 
#' Low niche overlap causes species to have greater per capita growth rates when rare than when common. 
#' If species limit individuals of their own species and their competitors equally, then niche overlap is 1, 
#' and coexistence is not possible unless species are otherwise identical. 
#' At the other extreme, if species have no interspecific effects, then niche overlap is 0.
#' 
#' Niche overlap has a common functional form, in the context of Modern Coexistence Theory (MCT), 
#' for a series of models, including those specified in table A1 of
#' Hart et al. (2018) Journal of Ecology 106, 1902-1909. 
#' Other model families may not adhere to the general definition. 
#' 
#' Furthermore, the MCT definition only accounts for competitive interactions (i.e. positive alpha coefficients
#' in these models). An alternative definition is given in Saavedra et al. (2017) Ecological Monographs 87,470-486. 
#' In this 'structural approach', positive interactions are allowed. Incidentally, both approaches yield
#' qualitatively similar, but not equivalent, results for purely competitive matrices. 
#' 
#' In all cases, these definitions only apply to models whose feasible equilibrium point can be described by a 
#' linear equation (see Saavedra et al. 2017, Hart et al. 2018 for details).
#' 
#' This function calculates niche overlap among two or more taxa, using both the MCT and the structural formulation.
#' The function, as in `avg_fitness_diff` and `competitive_ability`, accepts three different parameterizations:
#' * A cxr_pm_multifit object, from which niche overlap will be computed across all species pairs.
#' * two cxr_pm_fit objects, one for each species.
#' * explicit lambda and alpha values, as well as the model family from which these parameters were obtained.
#' 
#' If negative interactions are present, the MCT niche overlap will be NA.
#' The cxr objects may be calculated with user-defined model families. If this is the case, or
#' if simply a 2x2 matrix is provided, the niche overlap metrics will be calculated and 
#' a warning will be raised.
#' 
#' @param cxr_multifit cxr_pm_multifit object, with parameters for a series of species.
#' @param cxr_sp1 cxr_pm_fit object giving the parameters from the first species.
#' @param cxr_sp2 cxr_pm_fit object giving the parameters from the second species. 
#' @param pair_matrix 2x2 matrix with intra and interspecific interaction 
#' coefficients between the two species.
#'
#' @return either a dataframe with as many rows as species, or a single named numeric vector,
#' containing niche overlap values for the MCT (modern coexistence theory) and SA (structural approach)
#' formulations.
#' @importFrom utils combn
#' @export
#' @md
#' @examples 
#' niche_overlap(pair_matrix = matrix(c(0.33,0.12,0.2,0.4),nrow = 2))
#' 
niche_overlap <- function(cxr_multifit = NULL,cxr_sp1 = NULL, cxr_sp2 = NULL, pair_matrix = NULL){
  
  res <- NULL
  
  if(!is.null(cxr_multifit)){
    if(!is.null(cxr_sp1) | !is.null(cxr_sp2) | !is.null(pair_matrix)){
      message("cxr niche_overlap: the 'cxr_pm_multifit' object will be used, other
              arguments will be discarded.")
    }
    
    mf <- substr(cxr_multifit$model_name,1,2)
    if(!mf %in% c("LV","BH","RK","LW")){
      warning("niche_overlap: calculating niche overlap for coefficients estimated from a
              non-standard model family. Be aware that this may yield incorrect results.",call.=FALSE)
    }
      res <- as.data.frame(t(combn(cxr_multifit$taxa,2)),stringsAsFactors = FALSE)
      names(res) <- c("sp1","sp2")
      res$niche_overlap_MCT <- NA_real_
      res$niche_overlap_SA <- NA_real_
      
      for(ic in 1:nrow(res)){
        nov_matrix <- cxr_multifit$alpha_matrix[c(res$sp1[ic],res$sp2[ic]),c(res$sp1[ic],res$sp2[ic])]
        if(!any(is.na(nov_matrix))){
          res$niche_overlap_SA[ic] <- niche_overlap_SA(nov_matrix)
          if(all(nov_matrix>=0)){
            res$niche_overlap_MCT[ic] <- niche_overlap_MCT(nov_matrix)
          }# if all >=0
        }
      }# for each pair
    
  }# multispecies fit
  
  if(!is.null(cxr_sp1) & !is.null(cxr_sp2)){
    if(!is.null(pair_matrix)){
      message("cxr niche_overlap: both cxr objects and a pairwise matrix were specified. 
              Pairwise matrix will be discarded.")
    }
    
    sp1_model <- substr(cxr_sp1$model_name,4,5)
    sp2_model <- substr(cxr_sp2$model_name,4,5)
    
    if(sp1_model == sp2_model){
    
      if(!sp1_model %in% c("LV","BH","RK","LW")){
        warning("niche_overlap: calculating niche overlap for coefficients estimated from a
              non-standard model family. Be aware that this may yield incorrect results.",call.=FALSE)
      }

      sp1 <- cxr_sp1$focal_ID
      sp2 <- cxr_sp2$focal_ID
      
      if(is.null(sp1) | is.null(sp2)){
        message(paste("niche_overlap ERROR: a 'cxr_pm_fit' object passed does not contain information for identifying the 
                      focal taxa, i.e. does not have a valid 'focal_ID' field. 
                      This is probably because the fit was done without specifying 'alpha_intra' and 
                      'focal_column' in the function 'cxr_pm_fit'. Without this information,
                      intra and inter-specific interactions cannot be identified, 
                      and niche overlap cannot be computed."))
        return(NULL)
      }else{
        intra_sp1 <- cxr_sp1$alpha_intra
        intra_sp2 <- cxr_sp2$alpha_intra
        inter_sp1_sp2 <- cxr_sp1$alpha_inter[which(names(cxr_sp1$alpha_inter) == sp2)]
        inter_sp2_sp1 <- cxr_sp1$alpha_inter[which(names(cxr_sp2$alpha_inter) == sp1)]
        
        if(is.null(intra_sp1) | is.null(intra_sp2) | is.null(inter_sp2_sp1) | is.null(inter_sp1_sp2)){
          message(paste("niche_overlap ERROR: pairwise interactions could not be retrieved from the 'cxr_pm_fit' objects. 
          This could be because terms in the 'alpha_inter' fields do not match the names of the other taxa, 
          or because the 'alpha_intra' fields are NULL."))
          return(NULL)
        }else{
          nov_matrix <- matrix(c(intra_sp1,inter_sp2_sp1,inter_sp1_sp2,intra_sp2),nrow = 2)
          if(!any(is.na(nov_matrix))){
            res <- c(niche_overlap_MCT = niche_overlap_MCT(nov_matrix),
                   niche_overlap_SA = niche_overlap_SA(nov_matrix))
          }else{
            res <- c(niche_overlap_MCT = NA_real_,
                     niche_overlap_SA = NA_real_)
          }# if not NA
        }# if-else null
        
      }# if-else null

    # }# if-else model found
    }else{
      message("niche_overlap ERROR: 'cxr_sp1' and 'cxr_sp2' were fitted for different model families.")
      return(NULL)
    }
  }else if(!is.null(pair_matrix)){
    
    warning("niche_overlap: calculating niche overlap for coefficients estimated from an
              unknown model family. Be aware that this may yield incorrect results.",call.=FALSE)
    if(!any(is.na(pair_matrix))){
      res <- c(niche_overlap_MCT = niche_overlap_MCT(pair_matrix),
               niche_overlap_SA = niche_overlap_SA(pair_matrix))
    }else{
      res <- c(niche_overlap_MCT = NA_real_,
               niche_overlap_SA = NA_real_)
    }# if not NA
    
  }
  res
}
