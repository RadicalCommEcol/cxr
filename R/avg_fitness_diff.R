#' Average fitness differences
#'
#' computes the average fitness differences among two or more species according to the formulation
#' of the MCT (Chesson 2012, Godoy and Levine 2014), and according to the structural approach (Saavedra et al. 2017).
#' For the MCT version, the average fitness ratio is decomposed in a 'demographic ratio' and a 'competitive response ratio',
#' the product of which is the average fitness ratio (Godoy and Levine 2014). This formulation is only valid for competitive
#' interaction coefficients (i.e. positive alpha values in the interaction matrix). The structural analog can be computed
#' for any interaction matrix, on the other hand. Note that the 'demographic ratio' is model-specific (Hart et al. 2018).
#' 
#' This function, as in `niche_overlap` and `competitive_ability`, accepts three different parameterizations:
#' * A cxr_pm_multifit object, from which average fitness differences will be computed across all species pairs.
#' * two cxr_pm_fit objects, one for each species.
#' * explicit lambda and alpha values, as well as the model family from which these parameters were obtained.
#' 
#' If using the third parameterization, the function will try to find a model-specific function
#' for obtaining the demographic ratio, by looking at the 'model_family' parameter. 
#' If this specific function is not found, it will resort to the standard Lotka-Volterra 
#' formulation (lambda in the numerator term). 
#' Overall, we strongly suggest that you use the standard formulation ONLY if you are completely confident 
#' that your custom model is consistent with it. 
#' Otherwise, you should include your own formulation of the demographic ratio (see vignette 4).
#'
#' @param cxr_multifit cxr_pm_multifit object, with parameters for a series of species.
#' @param cxr_sp1 cxr_pm_fit object giving the parameters from the first species.
#' @param cxr_sp2 cxr_pm_fit object giving the parameters from the second species. 
#' @param pair_lambdas numeric vector of length 2 giving lambda values for the two species.
#' @param pair_matrix 2x2 matrix with intra and interspecific interaction 
#' coefficients between the two species.
#' @param model_family model family for which to calculate fitness differences.
#'
#' @return data frame with variable number of rows, and columns specifying the different components of the MCT average fitness ratio, 
#' as well as its structural analog. The average fitness ratio informs quantitatively about the better competitor. 
#' If the ratio is < 1, sp2 is the better competitor; if = 1, both species are equivalent competitors, if > 1, sp1 is the better competitor.
#' 
#' @export
#' @md
#'
#' @examples 
#' avg_fitness_diff(pair_lambdas = runif(2,1,10),
#'                  pair_matrix = matrix(runif(4,0,1),nrow = 2),
#'                  model_family = "BH")
avg_fitness_diff <- function(cxr_multifit = NULL,
                             cxr_sp1 = NULL, 
                             cxr_sp2 = NULL, 
                             pair_lambdas = NULL, 
                             pair_matrix = NULL,
                             model_family = NULL){
  
  res <- NULL
  
  if(!is.null(cxr_multifit)){
    if(!is.null(cxr_sp1) | !is.null(cxr_sp2) | !is.null(pair_matrix)){
      message("cxr avg_fitness_diff: the 'cxr_pm_multifit' object will be used, other
              arguments will be discarded.")
    }
    
    mf <- substr(cxr_multifit$model_name,1,2)
    if(!mf %in% c("LV","BH","RK","LW")){
      warning("avg_fitness_diff: calculating average niche differences for coefficients estimated from a
              custom model family. Be aware that this may yield incorrect results.",call.=FALSE)
    }
    res <- expand.grid(sp1 = cxr_multifit$taxa,sp2 = cxr_multifit$taxa,stringsAsFactors = FALSE)#as.data.frame(t(combn(cxr_multifit$taxa,2)),stringsAsFactors = FALSE)
    # names(res) <- c("sp1","sp2")
    res$demographic_ratio <- NA_real_
    res$competitive_response_ratio <- NA_real_
    res$average_fitness_ratio <- NA_real_
    res$structural_fitness_diff <- NA_real_
    
    for(ic in 1:nrow(res)){
      my_matrix <- cxr_multifit$alpha_matrix[c(res$sp1[ic],res$sp2[ic]),c(res$sp1[ic],res$sp2[ic])]
      my_lambdas <- cxr_multifit$lambda[c(res$sp1[ic],res$sp2[ic])]
      
      res$competitive_response_ratio[ic] <- cxr_competitive_response(my_matrix)
      res$demographic_ratio[ic] <- cxr_demographic_ratio(mf,my_lambdas)
      res$average_fitness_ratio[ic] <- res$competitive_response_ratio[ic] * res$demographic_ratio[ic]
      
      res$structural_fitness_diff[ic] <- cxr_structural_fitness_diff(my_matrix,my_lambdas,mf)
      
      }# for each pair
    
  }# multispecies fit
  
  if(!is.null(cxr_sp1) & !is.null(cxr_sp2)){
    if(!is.null(pair_matrix)){
      message("cxr avg_fitness_diff: both cxr objects and a pairwise matrix were specified. 
              Pairwise matrix will be discarded.")
    }
    
    sp1_model <- substr(cxr_sp1$model_name,4,5)
    sp2_model <- substr(cxr_sp2$model_name,4,5)
    
    if(sp1_model == sp2_model){
      
      if(!sp1_model %in% c("LV","BH","RK","LW")){
        warning("avg_fitness_diff: calculating niche overlap for coefficients estimated from a
              non-standard model family. Be aware that this may yield incorrect results.",call.=FALSE)
      }
      
      sp1 <- cxr_sp1$focal_ID
      sp2 <- cxr_sp2$focal_ID
      
      if(is.null(sp1) | is.null(sp2)){
        message(paste("avg_fitness_diff ERROR: a 'cxr_pm_fit' object passed does not contain information for identifying the 
                      focal taxa, i.e. does not have a valid 'focal_ID' field. 
                      This is probably because the fit was done without specifying 'alpha_intra' and 
                      'focal_column' in the function 'cxr_pm_fit'. Without this information,
                      intra and inter-specific interactions cannot be identified, 
                      and avg_fitness_diff cannot be computed."))
        return(NULL)
      }else{
        lambda_sp1 <- cxr_sp1$lambda
        lambda_sp2 <- cxr_sp2$lambda
        
        intra_sp1 <- cxr_sp1$alpha_intra
        intra_sp2 <- cxr_sp2$alpha_intra
        inter_sp1_sp2 <- cxr_sp1$alpha_inter[which(names(cxr_sp1$alpha_inter) == sp2)]
        inter_sp2_sp1 <- cxr_sp2$alpha_inter[which(names(cxr_sp2$alpha_inter) == sp1)]
        
        if(is.null(intra_sp1) | is.null(intra_sp2) | is.null(inter_sp2_sp1) | is.null(inter_sp1_sp2)){
          message(paste("avg_fitness_diff ERROR: pairwise interactions could not be retrieved from the 'cxr_pm_fit' objects. 
          This could be because terms in the 'alpha_inter' fields do not match the names of the other taxa, 
          or because the 'alpha_intra' fields are NULL."))
          return(NULL)
        }else{
          my_lambdas <- c(lambda_sp1,lambda_sp2)
          my_matrix <- matrix(c(intra_sp1,inter_sp2_sp1,inter_sp1_sp2,intra_sp2),nrow = 2)
          res <- data.frame(sp1 = sp1, 
                            sp2 = sp2)
          res$demographic_ratio <- cxr_demographic_ratio(sp1_model,my_lambdas)
          res$competitive_response_ratio <- cxr_competitive_response(my_matrix)
          res$average_fitness_ratio <- res$demographic_ratio * res$competitive_response_ratio 
          res$structural_fitness_diff <- cxr_structural_fitness_diff(my_matrix,my_lambdas,sp1_model)
          
        }# if-else null
        
      }# if-else null
      
      # }# if-else model found
    }else{
      message("avf_fitness_diff ERROR: 'cxr_sp1' and 'cxr_sp2' were fitted for different model families.")
      return(NULL)
    }
  }else if(!is.null(pair_matrix)){
    
    # message("avg_fitness_diff: calculating fitness differences for coefficients estimated from an
    #           unknown model family. Be aware that this may yield incorrect results.\n")
    
    res <- data.frame(sp1 = NA_character_,sp2 = NA_character_)
    if(!is.null(rownames(pair_matrix))){
      res$sp1 <- rownames(pair_matrix)[1]
      res$sp2 <- rownames(pair_matrix)[2]
    }
    
    res$demographic_ratio <- cxr_demographic_ratio(model_family,pair_lambdas)
    res$competitive_response_ratio <- cxr_competitive_response(pair_matrix)
    res$average_fitness_ratio <- res$demographic_ratio * res$competitive_response_ratio 
    res$structural_fitness_diff <- cxr_structural_fitness_diff(pair_matrix,pair_lambdas,model_family)
    
  }
  res
}
