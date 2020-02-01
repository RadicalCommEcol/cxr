# load test data
library(cxr)
# data("competition")

# TEMP
source("R/cxr_return_init_length.R")
source("R/cxr_init_params.R")
source("R/cxr_retrieve_params.R")
source("R/pm_BH_alpha_pairwise_lambdacov_none_alphacov_none.R")
source("R/pm_BH_alpha_pairwise_lambdacov_global_alphacov_global.R")
source("R/cxr_pm_bootstrap.R")
source("R/cxr_pm_fit.R")
source("R/cxr_check_input_data.R")

# spread the data from long to wide format
# competition.data <- tidyr::spread(competition,competitor,number,fill = 0)
# focal.sp <- unique(competition$focal)
# mindata <- subset(competition.data,focal == "LEMA")
# mindata$fitness <- log(mindata$seed)
# mindata <- mindata[,c("fitness",as.character(focal.sp))] #changes as.character,
# mind2 <- subset(competition.data,focal == "HOMA")
# mind2$fitness <- log(mind2$seed)
# mind2 <- mind2[,c("fitness", as.character(focal.sp))] #idem

# test the new competition data, sorted in a list
load("../Caracoles/data/competition.RData")
data <- neigh_list[1:3]

initial_values <- list(lambda = 1,alpha = 0.1,lambda_cov = 0.1, alpha_cov = 0.1)
lower_bounds <- list(lambda = 0.01,alpha = 0.01,lambda_cov = 0.01, alpha_cov = 0.01)
upper_bounds <- list(lambda = 10,alpha = 1,lambda_cov = 1, alpha_cov = 1)

# covariates: rows are observations, columns are different covariates
# either matrix or dataframe, will be transformed to matrix in the function
c1 <- data.frame(c1 = rnorm(nrow(mindata),1,.1))
c2 <- data.frame(c2 = rnorm(nrow(mind2),1,.1))
covariates <- list(c1 = c1, c2 = c2)

model_family <- "BH"
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "global"
fixed_terms <- NULL
bootstrap_samples <- 3

cxr_pm_multifit <- function(data, 
                            model_family = c("BH"),
                            covariates = NULL, 
                            optimization_method = c(), 
                            alpha_form = c("none","global","pairwise"), 
                            lambda_cov_form = c("none","global"),
                            alpha_cov_form = c("none","global","pairwise"),
                            initial_values = NULL,
                            lower_bounds = NULL,
                            upper_bounds = NULL,
                            fixed_terms = NULL,
                            # errors
                            bootstrap_samples = 0
){
  

# prepare multisp data ----------------------------------------------------
spnames <- names(data)

# fit every sp ------------------------------------------------------------
spfits <- list()
for(i.sp in 1:length(data)){
  spfits[[i.sp]] <- cxr_pm_fit(data = data[[i.sp]],
                               model_family = model_family,
                               covariates = covariates[[i.sp]],
                               optimization_method = optimization_method,
                               alpha_form = alpha_form,
                               lambda_cov_form = lambda_cov_form,
                               alpha_cov_form = alpha_cov_form,
                               initial_values = initial_values,
                               lower_bounds = lower_bounds,
                               upper_bounds = upper_bounds,
                               fixed_terms = fixed_terms,
                               bootstrap_samples = bootstrap_samples
                               )
}

# output ------------------------------------------------------------------

  # return a cxr_pm_multifit object
  # which is basically the same as the base object
  # but with info on more sp. e.g. lambda is a 1d vector,
  # alpha is a n x n matrix

splambda <- NULL
spalpha <- NULL
splambda_cov <- NULL
spalpha_cov <- NULL

er_splambda <- NULL
er_spalpha <- NULL
er_splambda_cov <- NULL
er_spalpha_cov <- NULL

spllik <- NULL

for(i.sp in 1:length(spnames)){
  # lambda
  if(!is.null(spfits[[i.sp]]$lambda)){
    mylambda <- spfits[[i.sp]]$lambda
    names(mylambda) <- spnames[i.sp]
    splambda <- c(splambda,mylambda)
  }
  # alpha
  if(!is.null(spfits[[i.sp]]$alpha)){
    myalpha <- spfits[[i.sp]]$alpha
    spalpha <- rbind(spalpha,myalpha)
    rownames(spalpha)[i.sp] <- spnames[i.sp]
  }
  # lambda_cov
  if(!is.null(spfits[[i.sp]]$lambda_cov)){
    mylambda_cov <- spfits[[i.sp]]$lambda_cov
    names(mylambda_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$lambda_cov),sep="")
    splambda_cov <- c(splambda_cov,mylambda_cov)
  }
  if(!is.null(spfits[[i.sp]]$alpha_cov)){
    myalpha_cov <- spfits[[i.sp]]$alpha_cov
    names(myalpha_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$alpha_cov),sep="")
    spalpha_cov <- c(spalpha_cov,myalpha_cov)
  }
  
  # errors
  
  # lambda
  if(!is.null(spfits[[i.sp]]$lambda_standard_error)){
    erlambda <- spfits[[i.sp]]$lambda_standard_error
    if(!is.null(names(erlambda))){
      names(erlambda) <- paste(spnames[i.sp],"_",names(erlambda),sep="")
    }else{
      names(erlambda) <- paste(spnames[i.sp],"_lambda_se",sep="")
    }
    er_splambda <- c(er_splambda,erlambda)
  }
  # alpha
  if(!is.null(spfits[[i.sp]]$alpha_standard_error)){
    eralpha <- spfits[[i.sp]]$alpha_standard_error
    er_spalpha <- rbind(er_spalpha,eralpha)
    rownames(er_spalpha)[i.sp] <- spnames[i.sp]
  }
  # lambda_cov
  if(!is.null(spfits[[i.sp]]$lambda_cov_standard_error)){
    erlambda_cov <- spfits[[i.sp]]$lambda_cov_standard_error
    names(erlambda_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$lambda_cov_standard_error),sep="")
    er_splambda_cov <- c(er_splambda_cov,erlambda_cov)
  }
  if(!is.null(spfits[[i.sp]]$alpha_cov_standard_error)){
    eralpha_cov <- spfits[[i.sp]]$alpha_cov_standard_error
    names(eralpha_cov) <- paste(spnames[i.sp],"_",names(spfits[[i.sp]]$alpha_cov_standard_error),sep="")
    er_spalpha_cov <- c(er_spalpha_cov,eralpha_cov)
  }
  
  # log-likelihood
  myllik <- spfits[[i.sp]]$log_likelihood
  names(myllik) <- spnames[i.sp]
  spllik <- c(spllik,myllik)
}

list_names <- c("model_name",
                "model",
                "data",
                "model_family",
                "covariates",
                "optimization_method",
                "initial_values",
                "fixed_terms",
                "lambda","alpha","lambda_cov","alpha_cov",
                "lambda_standard_error","alpha_standard_error",
                "lambda_cov_standard_error","alpha_cov_standard_error",
                "log_likelihood")

fit <- sapply(list_names,function(x) NULL)

fit$model_name <- spfits[[1]]$model_name
fit$model <- spfits[[1]]$fitness_model
fit$data <- data
fit$model_family <- model_family
fit$covariates <- covariates
fit$optimization_method <- optimization_method
fit$initial_values <- initial_values

# for returning explicit NULL values
if(!is.null(fixed_terms)){
  fit$fixed_terms <- fixed_terms
}
if(!is.null(splambda)){
  fit$lambda <- splambda
}
if(!is.null(spalpha)){
  fit$alpha <- spalpha
}
if(!is.null(splambda_cov)){
  fit$lambda_cov <- splambda_cov
}
if(!is.null(spalpha_cov)){
  fit$alpha_cov <- spalpha_cov
}
if(!is.null(er_splambda)){
  fit$lambda_standard_error <- er_splambda
}
if(!is.null(er_spalpha)){
  fit$alpha_standard_error <- er_spalpha
}
if(!is.null(er_splambda_cov)){
  fit$lambda_cov_standard_error <- er_splambda_cov
}
if(!is.null(er_spalpha_cov)){
  fit$alpha_cov_standard_error <- er_spalpha_cov
}

fit$log_likelihood <- spllik

class(fit) <- "cxr_pm_multifit"
fit

}

