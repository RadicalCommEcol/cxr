#' Converts a densities list to a tidy dataframe
#'
#' @param densities list, species (optionally x year) with each element holding 
#' a sites x stages matrix. This function assumes three life stages.
#'
#' @return dataframe with columns species-stage-site(-year)-density
#' @export
densities_to_df <- function(densities){
  
  # to consider a densities list with year slots
  include.year <- FALSE
  
  num.sp <- length(densities)
  
  if(class(densities[[1]]) == "list"){
    include.year <- TRUE
    num.stages <- ncol(densities[[1]][[1]])
    num.sites <- nrow(densities[[1]][[1]])
    num.years <- length(densities[[1]])
    
  }else{
    num.stages <- ncol(densities[[1]])
    num.sites <- nrow(densities[[1]])

  }
  
  ddf <- list()
  
  # fill dataframe
  if(include.year){
    for(i.sp in 1:num.sp){
      for(i.year in 1:num.years){
        tdf <- as.data.frame(as.table(densities[[i.sp]][[i.year]]),
                             stringsAsFactors = FALSE)
        names(tdf) <- c("site","stage","density")
        tdf$species <- i.sp
        tdf$year <- i.year
        ddf[[length(ddf) + 1]] <- tdf
      }# for i.year
    }# for i.sp
    
    ddff <- dplyr::bind_rows(ddf)
    ddff <- ddff[,c("species","stage","site","year","density")]
    
  }else{
    for(i.sp in 1:num.sp){
        tdf <- as.data.frame(as.table(densities[[i.sp]]),
                             stringsAsFactors = FALSE)
        names(tdf) <- c("site","stage","density")
        tdf$species <- i.sp
        ddf[[length(ddf) + 1]] <- tdf
    }# for i.sp
    
    ddff <- dplyr::bind_rows(ddf)
    ddff <- ddff[,c("species","stage","site","density")]
  }
  
  ddff$stage <- dplyr::recode(ddff$stage, 
                        A = "juvenile", 
                        B = "non-reproductive-adult", 
                        C = "reproductive-adult")
  ddff$site <- match(ddff$site,LETTERS)
  
  return(ddff)
}