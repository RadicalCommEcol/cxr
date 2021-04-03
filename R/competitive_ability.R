#' Competitive ability among pairs of species
#'
#' Computes the competitive ability among two species, as defined by Hart et al. (2018). This metric, as others in MCT, is
#' model-specific; the formulation for a series of Lotka-Volterra-like models is given in table A1 of Hart et al. (2018). 
#' We include in `cxr` by default the formulation for Beverton-Holt, Ricker, Law-Watkinson, and Lotka-Volterra families.
#'
#' The function, as in `avg_fitness_diff` and `niche_overlap`, accepts three different parameterizations:
#' * A cxr_pm_multifit object, from which competitive ability of a focal species relative to a given competitor 
#' will be computed across all species pairs.
#' * two cxr_pm_fit objects, one for a focal species and one for a competitor.
#' * explicit lambda and alpha values, as well as the model family from which these parameters were obtained.
#'
#' If the third parameterization is used, the function will try to find a model-specific function
#' for obtaining the competitive ability, by looking at the 'model_family' parameter.
#' If this specific function is not found, it will resort to the standard Lotka-Volterra 
#' formulation (lambda - 1 in the numerator term, Hart et al. 2018). 
#' Overall, we strongly suggest that you use the standard formulation ONLY if you are completely confident 
#' that the model from which you obtained your parameters is consistent with it. 
#' Otherwise, you should include your own formulation of competitive ability (see vignette 4).
#'
#' @param cxr_multifit cxr_pm_multifit object, with parameters for a series of species.
#' @param cxr_sp1 cxr_pm_fit object giving the parameters from the first species.
#' @param cxr_sp2 cxr_pm_fit object giving the parameters from the second species. 
#' @param lambda numeric lambda value of the focal species.
#' @param pair_matrix 2x2 matrix with intra and interspecific interaction 
#' coefficients between the focal and competitor species.
#' @param model_family model family for which to calculate competitive ability.
#'
#' @return data frame with variable number of rows and three columns, specifying taxa identity and the competitive ability
#' of focal species (sp1) relative to the competitor (sp2).
#' @export
#' @md
#'
#' @examples competitive_ability(lambda = runif(1,1,10),
#'                               pair_matrix = matrix(runif(4,0,1),nrow = 2),
#'                               model_family = "BH")
competitive_ability <- function(cxr_multifit = NULL,
                             cxr_sp1 = NULL, 
                             cxr_sp2 = NULL, 
                             lambda = NULL, 
                             pair_matrix = NULL,
                             model_family = NULL){
  
  res <- NULL
  
  if(!is.null(cxr_multifit)){
    if(!is.null(cxr_sp1) | !is.null(cxr_sp2) | !is.null(pair_matrix)){
      message("cxr competitive ability: the 'cxr_pm_multifit' object will be used, other
              arguments will be discarded.")
    }
    
    mf <- substr(cxr_multifit$model_name,1,2)
    if(!mf %in% c("LV","BH","RK","LW")){
      message("cxr competitive ability: calculating average niche differences for coefficients estimated from a
              custom model family.\n")
    }
    
    ca_fun <- paste(mf,"_competitive_ability",sep="")
    ca_model <- try(get(ca_fun),silent = TRUE)
    
    if(inherits(ca_model,"try-error")){
      warning(paste("cxr competitive ability: function '",ca_fun,"' could not be retrieved.
                    The default formulation for Lotka-Volterra models (lambda - 1 in numerator term) will be used. 
                    Be aware that this may yield incorrect results for your model family.\n"
                    ,sep=""),call. = FALSE)
      ca_model <- BH_competitive_ability
    }
    
    res <- expand.grid(sp1 = cxr_multifit$taxa,sp2 = cxr_multifit$taxa,stringsAsFactors = FALSE)#as.data.frame(t(combn(cxr_multifit$taxa,2)),stringsAsFactors = FALSE)
    res <- subset(res,sp1 != sp2)
    res$competitive_ability_sp1 <- NULL
    
    for(ic in 1:nrow(res)){
      my_matrix <- cxr_multifit$alpha_matrix[c(res$sp1[ic],res$sp2[ic]),c(res$sp1[ic],res$sp2[ic])]
      my_lambda <- cxr_multifit$lambda[res$sp1[ic]]
      if(!any(is.na(my_matrix)) & !is.na(my_lambda)){
        res$competitive_ability_sp1[ic] <- ca_model(my_lambda,my_matrix)
      }else{
        res$competitive_ability_sp1[ic] <- NA_real_
      }
      
      }# for each pair
    
  }# multispecies fit
  
  if(!is.null(cxr_sp1) & !is.null(cxr_sp2)){
    if(!is.null(pair_matrix)){
      message("cxr competitive ability: both cxr objects and a pairwise matrix were specified. 
              Pairwise matrix will be discarded.")
    }
    
    sp1_model <- substr(cxr_sp1$model_name,4,5)
    sp2_model <- substr(cxr_sp2$model_name,4,5)
    
    if(sp1_model == sp2_model){
      
      if(!sp1_model %in% c("LV","BH","RK","LW")){
        message("cxr competitive ability: calculating competitive ability for coefficients estimated from a
              custom model family.\n")
      }
        ca_fun <- paste(sp1_model,"_competitive_ability",sep="")
        ca_model <- try(get(ca_fun),silent = TRUE)
        
        if(inherits(ca_model,"try-error")){
          warning(paste("cxr competitive ability: function '",ca_fun,"' could not be retrieved.
                    The default formulation (lambda - 1 in numerator term) will be used. 
                    Be aware that this may yield incorrect results for your model family.\n"
                        ,sep=""),call. = FALSE)
          ca_model <- BH_competitive_ability
        }
      
      sp1 <- cxr_sp1$focal_ID
      sp2 <- cxr_sp2$focal_ID
      
      if(is.null(sp1) | is.null(sp2)){
        message(paste("cxr competitive ability ERROR: a 'cxr_pm_fit' object passed does not contain information for identifying the 
                      focal taxa, i.e. does not have a valid 'focal_ID' field. 
                      This is probably because the fit was done without specifying 'alpha_intra' and 
                      'focal_column' in the function 'cxr_pm_fit'. Without this information,
                      intra and inter-specific interactions cannot be identified, 
                      and competitive ability cannot be computed."))
        return(NULL)
      }else{
        lambda_sp1 <- cxr_sp1$lambda
        # lambda_sp2 <- cxr_sp2$lambda
        
        intra_sp1 <- cxr_sp1$alpha_intra
        intra_sp2 <- cxr_sp2$alpha_intra
        inter_sp1_sp2 <- cxr_sp1$alpha_inter[which(names(cxr_sp1$alpha_inter) == sp2)]
        inter_sp2_sp1 <- cxr_sp1$alpha_inter[which(names(cxr_sp2$alpha_inter) == sp1)]
        
        if(is.null(intra_sp1) | is.null(intra_sp2) | is.null(inter_sp2_sp1) | is.null(inter_sp1_sp2)){
          message(paste("cxr competitive ability ERROR: pairwise interactions could not be retrieved from the 'cxr_pm_fit' objects. 
          This could be because terms in the 'alpha_inter' fields do not match the names of the other taxa, 
          or because the 'alpha_intra' fields are NULL."))
          return(NULL)
        }else{
          
          my_lambda <- lambda_sp1
          my_matrix <- matrix(c(intra_sp1,inter_sp2_sp1,inter_sp1_sp2,intra_sp2),nrow = 2)
          res <- data.frame(sp1 = sp1, 
                            sp2 = sp2)
          
          if(!any(is.na(my_matrix)) & !is.na(my_lambda)){
            res$competitive_ability_sp1 <- ca_model(my_lambda,my_matrix)
          }else{
            res$competitive_ability_sp1 <- NA_real_
          }

        }# if-else null
        
      }# if-else null
      
      # }# if-else model found
    }else{
      message("cxr competitive ability ERROR: 'cxr_sp1' and 'cxr_sp2' were fitted for different model families.")
      return(NULL)
    }
  }else if(!is.null(pair_matrix)){
    
    ca_fun <- paste(model_family,"_competitive_ability",sep="")
    ca_model <- try(get(ca_fun),silent = TRUE)
    
    if(inherits(ca_model,"try-error")){
      
      ca_model <- BH_competitive_ability
      
      warning(paste("cxr competitive ability: function '",ca_fun,"' could not be retrieved. 
      The default formulation for Lotka-Volterra models (lambda - 1 in numerator term) will be used. 
      Be aware that this may yield incorrect results for your model family.\n"
                    ,sep=""),call. = FALSE)
    }

    res <- data.frame(sp1 = NA_character_,sp2 = NA_character_)
    if(!is.null(rownames(pair_matrix))){
      res$sp1 <- rownames(pair_matrix)[1]
      res$sp2 <- rownames(pair_matrix)[2]
    }
    
    if(!any(is.na(pair_matrix)) & !is.na(lambda)){
      res$competitive_ability_sp1 <- ca_model(lambda,pair_matrix)    
    }else{
      res$competitive_ability_sp1 <- NA_real_
    }
    
  }# if-else parameters
  res
}
