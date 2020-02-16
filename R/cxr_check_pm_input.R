#' internal to check input consistency to cxr_pm_fit
#'
#' this function draws on several other internal functions,
#' performing checks on data consistency, covariates, bounds,
#' initial values, packages installed.
#'
#' @inheritParams cxr_pm_fit
#'
#' @return list with two components. 'input.ok' is a character,
#' either 'ok','warning',or 'error'. 'input.message' is either NULL
#' or the message returned by the errors/warnings
#' @noRd
cxr_check_pm_input <- function(data, 
                               focal_column = NULL,
                               model_family = c("BH"),
                               covariates = NULL, 
                               optimization_method = c("BFGS", "CG", "Nelder-Mead", 
                                                       "ucminf","L-BFGS-B", "nlm", "nlminb", 
                                                       "Rcgmin", "Rvmmin", "spg", 
                                                       "bobyqa", "nmkb", "hjkb",
                                                       "nloptr_CRS2_LM","nloptr_ISRES",
                                                       "nloptr_DIRECT_L_RAND","DEoptimR",
                                                       "hydroPSO","GenSA"), 
                               alpha_form = c("none","global","pairwise"), 
                               lambda_cov_form = c("none","global"),
                               alpha_cov_form = c("none","global","pairwise"),
                               initial_values = list(lambda = 0, alpha = 0, lambda_cov = 0, alpha_cov = 0),
                               lower_bounds = NULL,
                               upper_bounds = NULL,
                               fixed_terms = NULL){
  
  input.ok <- "ok"
  input.message <- NULL
  
  # basic data consistency
  t1 <- cxr_check_input_data(data,covariates)
  if(!t1){
    input.message <- ("cxr_pm_fit ERROR: check the consistency of your input data: 
    1) All variables are integer/numeric, with no NAs; 
    2) first column in 'data' is named 'fitness'; 
    3) abundances of at least one neighbour species in 'data';
    4) data and covariates (if present) have the same number of observations")
  }
  
  # check covariates if alpha_cov or lambda_cov are to be fit
  t2 <- !(is.null(covariates) & (alpha_cov_form != "none" | lambda_cov_form != "none"))
  if(!t2){
    input.message <- ("cxr_pm_fit ERROR: need to specify covariates if lambda_cov and/or alpha_cov are to be fit")
  }
  
  # check that lower/upper bounds are provided if the method requires it
  t3 <- cxr_check_method_boundaries(optimization_method,lower_bounds,upper_bounds, type = "pm")
  if(!t3){
    input.message <- ("cxr_pm_fit ERROR: check the optimization method selected and lower/upper bounds.
         The following methods require explicit lower and upper parameter boundaries to be set:
         L-BFGS-B, nlm, nlminb, Rcgmin, Rvmmin, spg, bobyqa, nmkb, hjkb, nloptr_CRS2_LM,
         nloptr_ISRES, nloptr_DIRECT_L_RAND, GenSA, hydroPSO, DEoptimR.")
  }
  
  t4 <- cxr_check_initial_values(initial_values,
                                 focal_column,
                                 lower_bounds,
                                 upper_bounds,
                                 fixed_terms)
  if(!t4){
    input.message <- ("cxr_pm_fit ERROR: please check the specified initial values/bounds.
                      1) only valid names are allowed, among 'lambda', 'alpha_intra',
                      'alpha_inter','lambda_cov','alpha_cov'.
                      2) if 'focal_column' is provided, you need to specify
                      initial values for 'alpha_intra', and viceversa.
                      3) elements must be the same in the three lists.
                      4) if bounds are provided, you need to specify both lower and upper ones."
                      )
  }
  
  # check installed packages for optimization method
  t5 <- TRUE
  if (optimization_method %in% c("nloptr_CRS2_LM","nloptr_ISRES","nloptr_DIRECT_L_RAND") & !requireNamespace("nloptr", quietly = TRUE)) {
    t5 <- FALSE
    input.message <- ("cxr_pm_fit ERROR: Package \"nloptr\" needed for the method selected to work.")
  }
  if (optimization_method == "GenSA" & !requireNamespace("GenSA", quietly = TRUE)) {
    t5 <- FALSE
    input.message <- ("cxr_pm_fit ERROR: Package \"GenSA\" needed for the method selected to work.")
  }
  if (optimization_method == "hydroPSO" & !requireNamespace("hydroPSO", quietly = TRUE)) {
    t5 <- FALSE
    input.message <- ("cxr_pm_fit ERROR: Package \"hydroPSO\" needed for the method selected to work.")
  }
  if (optimization_method == "DEoptimR" & !requireNamespace("DEoptimR", quietly = TRUE)) {
    t5 <- FALSE
    input.message <- ("cxr_pm_fit ERROR: Package \"DEoptimR\" needed for the method selected to work.")
  }
  
  w1 <- identical(initial_values,list(lambda = 0, alpha_intra = 0, alpha_inter = 0, lambda_cov = 0, alpha_cov = 0))
  if(w1){
    warning.message <- "cxr_pm_fit: Using default initial values. Note that these may not be appropriate for your data/model, or
    for the optimization method selected."
  }  
  w2 <- alpha_form != "pairwise" & (!is.null(focal_column))
  if(w2){
    warning.message <- "cxr_pm_fit: the specified 'alpha_form' does not support differentiating focal and non-focal observations.
  'focal_column' will be discarded, and initial values and bounds, if used, will be taken from 'alpha_inter'."
  }
  
  if(!all(c(t1,t2,t3,t4,t5))){
    input.ok <- "error"
  }
  
  if(any(c(w1,w2))){
    input.ok <- "warning"
    input.message <- warning.message
  }
  
  list(input.ok,input.message)
  
}