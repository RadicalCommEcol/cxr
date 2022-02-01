

#' Generate coefficients for obtaining vital rates
#' 
#' Any vital rate is a function of several parameters, potentially including
#' interactions or environmental effects. This function generates the coefficients
#' for these parameters, so that users do not have to introduce them all manually
#' in a `param` list. Coefficients can be generated from a random sampling of
#' a normal distribution with specified mean and standard deviation, 
#' or they can be retrieved from a model object that accepts a `tidy` function
#' from the broom/broom.mixed packages. This is because coefficients for 
#' vital rates can be understood as coefficients from statistical regressions.
#' 
#' In the current version, we assume that the model coefficients come from a
#' logistic regression with binomial family. Otherwise, the function will probably not fail,
#' but the coefficients will not be interpretable and the results in terms of obtaining the actual
#' vital rates from these will be meaningless.
#' 
#' Also note that you need to take care manually of the signs of the coefficients,
#' if entered through mean/sd pairs.
#'
#' @param param the original list with the structure of species, sites, 
#' vital rates to calculate, and parameters affecting them. See the function `build_param`
#' @param sp number or character of the species to calculate coefficients for. If empty, all species are assumed.
#' @param sites number or character of the sites to calculate coefficients for. If empty, all sites are assumed.
#' @param vital.rate character giving the vital rate to calculate coefficients for. If empty, all vital rates are assumed.
#' @param vr.coef character giving a specific coefficient to calculate. If empty, all coefficients are assumed.
#' @param mean.coef optional numeric value, mean for sampling coefficient values
#' @param sd.coef optional numeric value, standard deviation for sampling coefficient values
#' @param glm.object optional model object/coef table
#' @param glm.coef.equivalence if a glm table is provided and its names differ from the `param` data structure, 
#' you can include a named list in which names are the names from `param` and its elements are the equivalent names from the glm table
#'
#' @return the updated parameter list
#' @export
generate_vital_rate_coefs <- function(param,
                                      sp = NULL,
                                      sites = NULL,
                                      vital.rate = NULL,
                                      vr.coef = NULL,
                                      mean.coef = NULL,
                                      sd.coef = NULL,
                                      glm.object = NULL,
                                      glm.coef.equivalence = NULL){
  
  if(is.null(mean.coef) & is.null(sd.coef) & is.null(glm.object)){
    message("function generate_vital_rate_coefs ERROR: provide either mean.coef/sd.coef or glm.object")
    return(NULL)
  }
  
  # the object to return
  param.result <- param
  
  # coefficient names
  if(is.null(vr.coef)){
    my.coefs <- colnames(param[[1]][[1]])
  }else{
    my.coefs <- vr.coef
  }
  
  # which species? if null, all of them
  if(is.null(sp)){
    my.sp <- names(param)
  }else{
    my.sp <- sp
  }
  
  # which sites? if null, all of them
  if(is.null(sites)){
    my.sites <- names(param[[1]])
  }else{
    my.sites <- sites
  }
  
  # which vital rate to generate? if null, all of them
  if(is.null(vital.rate)){
    my.vital.rates <- rownames(param[[1]][[1]])
  }else{
    my.vital.rates <- vital.rate
  }  
  
  # check whether a glm table is provided, and the equivalence between
  # its coefficients and the names from param
  if(!is.null(glm.object)){
    glm.names <- rownames(glm.object)
    if(is.null(glm.coef.equivalence)){
      if(!my.coefs %in% glm.names){
        stop("names are not equivalent among glm.object and param")
      }# if
    }# if null
  }
  
  for(i.sp in 1:length(my.sp)){
    for(i.site in 1:length(my.sites)){
      for(i.rate in 1:length(my.vital.rates)){
        if(!is.null(mean.coef) & !is.null(sd.coef)){
          
            param.result[[my.sp[i.sp]]][[my.sites[i.site]]][my.vital.rates[i.rate],my.coefs] <- 
              rnorm(length(my.coefs),mean.coef,sd.coef)

        }else{

          for(i.coef in 1:length(my.coefs)){
            
            # identify the coefficient
            this.glm.coef <- glm.coef.equivalence[[which(names(glm.coef.equivalence) == my.coefs[i.coef])]]
            
            if(this.glm.coef %in% rownames(glm.object)){
              param.result[[my.sp[i.sp]]][[my.sites[i.site]]][my.vital.rates[i.rate],my.coefs[i.coef]] <- glm.object[this.glm.coef,1]
            }else{
              param.result[[my.sp[i.sp]]][[my.sites[i.site]]][my.vital.rates[i.rate],my.coefs[i.coef]] <- 0
            }
            
          }# for i.coef

          
        }# if-else mean/sd or glm
        
      }# for i.rate
    }# for i.site
  }# for i.sp
  
  return(param.result)
}
