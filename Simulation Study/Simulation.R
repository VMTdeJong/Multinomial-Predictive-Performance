# #######################################################################################################################

# Author of code: Valentijn M.T. de Jong.
# File last modified: December 2018
# Code last modified: < April 2018

# #######################################################################################################################

# This is code for a simulation study presented in a manuscript accepted for publication in 
# Statistics in medicine:

# Title: 	Predictive performance of multinomial logistic prediction models

# Authors:
#   Valentijn M. T. de Jong
#   Marinus J. C. Eijkemans
#   Ben van Calster
#   Dirk Timmerman
#   Karel G. M. Moons
#   Ewout W. Steyerberg
#   Maarten van Smeden

# For questions, email: 
#   V.M.T.deJong-2@umcutrecht.nl
#   Valentijn.M.T.de.Jong@gmail.com

# #######################################################################################################################

library(glmnet)
library(mlogit)

Simulation <- function(max_i = 10 ** 3, n, N, coefs, n_bi_var = 0, bi_prob = c(.5), id = NA, n_cv = 1, 
                       pop_seed = matrix(nrow = max_i, ncol = length(.Random.seed)), 
                       sample_seed = matrix(nrow = max_i, ncol = length(.Random.seed)),
                       seed = rep(NA, length(.Random.seed)),
                       cores = 1, computePDI = TRUE, compute_cal = TRUE, repetition = 0, 
                       lasso.which.lambda = "min", lasso.cv.lambda = "mean",
                       ridge.which.lambda = "min", ridge.cv.lambda = "mean",
                       lasso.cv.measure = "deviance", ridge.cv.measure = "deviance",
                       nfolds = 10,
                       max.time = 60 * 60 * 24 * 7 * .98,
                       n_lambda_Lasso = 200, n_lambda_Ridge = 200, proj_loc = getwd(),
                       min.cor = 0, max.cor = 0)
{
  begin             <- proc.time()
  rep_pop_seed      <- rep(FALSE, max_i)        # vector for tracking whether replicated seed was used.
  rep_sample_seed   <- rep(FALSE, max_i)
  if (!is.na(seed[1])) { pop_seed[1, ] <- seed}

  ##### Some important numbers
  n_cat             <- length(n)
  outcomes          <- 0:(n_cat - 1)
  mfor              <- CreateMlogitFormulaInt(ncol(as.matrix(coefs)))
  n_total           <- sum(n)
  N_total           <- sum(N)
  DGM_cm_probs      <- matrix(nrow = max_i, ncol = n_cat)
  dev_mean_cor      <- rep(NA, max_i)
  val_mean_cor      <- rep(NA, max_i)
  
  
  # Shrinkage parameters
  s_lambda_Lasso    <- rep(NA, max_i)
  s_lambda_Ridge    <- rep(NA, max_i)
  s_l_seq_Lasso     <- matrix(nrow = max_i, ncol = max(100, n_lambda_Lasso))
  s_l_seq_Ridge     <- matrix(nrow = max_i, ncol = max(100, n_lambda_Ridge))
  print(paste(lasso.cv.measure, " chosen for Lasso.")) 
  print(paste(ridge.cv.measure, " chosen for Ridge."))
  
  lambda_Lasso_rank <- rep(NA, max_i)
  lambda_Lasso_prop <- rep(NA, max_i)
  lambda_Lasso_max  <- rep(F,  max_i)
  lambda_Lasso_min  <- rep(F,  max_i)
  
  lambda_Ridge_rank <- rep(NA, max_i) 
  lambda_Ridge_prop <- rep(NA, max_i)
  lambda_Ridge_max  <- rep(F,  max_i)
  lambda_Ridge_min  <- rep(F,  max_i)
  user_l_lambda     <- list()
  user_r_lambda     <- list()

  ##### Declare Performance measures.
  # Pairwise One vs Rest, using my own function:
  OneVsRestC_ws_ML     <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  OneVsRestC_ws_Lasso  <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  OneVsRestC_ws_Ridge  <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  
  OneVsRestC_oos_ML    <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  OneVsRestC_oos_Lasso <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  OneVsRestC_oos_Ridge <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  OneVsRestC_oos_DGM   <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  
  # Pairwise C statistic computed with M-index
  M_ws_ML     <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  M_ws_Lasso  <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  M_ws_Ridge  <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  
  M_oos_ML    <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  M_oos_Lasso <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  M_oos_Ridge <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  M_oos_DGM   <- array(dim = list(n_cat, n_cat, max_i), dimnames = list("outcome" = outcomes, "other" = outcomes, "iteration" = 1:max_i))
  
  # PDI
  pdi_names            <- c("pdi1", "pdi2", "pdi3", "pdi", "pdiw")
  PDI_ws_ML            <- data.frame(matrix(NA, ncol = 5, nrow = max_i))
  PDI_ws_Lasso         <- data.frame(matrix(NA, ncol = 5, nrow = max_i))
  PDI_ws_Ridge         <- data.frame(matrix(NA, ncol = 5, nrow = max_i))
  colnames(PDI_ws_ML)  <- colnames(PDI_ws_Lasso) <- colnames(PDI_ws_Ridge) <- pdi_names
  
  PDI_oos_ML           <- data.frame(matrix(NA, ncol = 5, nrow = max_i))
  PDI_oos_Lasso        <- data.frame(matrix(NA, ncol = 5, nrow = max_i))
  PDI_oos_Ridge        <- data.frame(matrix(NA, ncol = 5, nrow = max_i))
  colnames(PDI_oos_ML) <- colnames(PDI_oos_Lasso) <- colnames(PDI_oos_Ridge) <- pdi_names
  
  PDI_oos_DGM          <- data.frame(matrix(NA, ncol = 5, nrow = max_i))
  names(PDI_oos_DGM)   <- pdi_names
  
  # Brier Score
  Brier_ws_ML      <- matrix(NA, ncol = n_cat, nrow = max_i)
  Brier_ws_Lasso   <- matrix(NA, ncol = n_cat, nrow = max_i)
  Brier_ws_Ridge   <- matrix(NA, ncol = n_cat, nrow = max_i)
  
  Brier_oos_ML     <- matrix(NA, ncol = n_cat, nrow = max_i)  
  Brier_oos_Lasso  <- matrix(NA, ncol = n_cat, nrow = max_i) 
  Brier_oos_Ridge  <- matrix(NA, ncol = n_cat, nrow = max_i) 
  Brier_oos_DGM    <- matrix(NA, ncol = n_cat, nrow = max_i) # DGM
  
  
  # Accuracy (mse) of predicted probabilites
  mse_probs_oos_ML      <- matrix(NA, ncol = n_cat, nrow = max_i)
  mse_probs_oos_Lasso   <- matrix(NA, ncol = n_cat, nrow = max_i)
  mse_probs_oos_Ridge   <- matrix(NA, ncol = n_cat, nrow = max_i)
  
  # Pseudo- R^2
  R2_ws_ML     <- data.frame(matrix(NA, ncol = 3, nrow = max_i))
  R2_ws_Lasso  <- data.frame(matrix(NA, ncol = 3, nrow = max_i))
  R2_ws_Ridge  <- data.frame(matrix(NA, ncol = 3, nrow = max_i))
  
  R2_oos_ML    <- data.frame(matrix(NA, ncol = 3, nrow = max_i))
  R2_oos_Lasso <- data.frame(matrix(NA, ncol = 3, nrow = max_i))
  R2_oos_Ridge <- data.frame(matrix(NA, ncol = 3, nrow = max_i))
  R2_oos_DGM   <- data.frame(matrix(NA, ncol = 3, nrow = max_i))    # DGM
  
  
  names(R2_ws_ML)     <- c("Cox", "Nagelkerke", "McFadden")
  names(R2_ws_Lasso)  <- c("Cox", "Nagelkerke", "McFadden")
  names(R2_ws_Ridge)  <- c("Cox", "Nagelkerke", "McFadden")
  
  names(R2_oos_ML)    <- c("Cox", "Nagelkerke", "McFadden")
  names(R2_oos_Lasso) <- c("Cox", "Nagelkerke", "McFadden")
  names(R2_oos_Ridge) <- c("Cox", "Nagelkerke", "McFadden")
  names(R2_oos_DGM)   <- c("Cox", "Nagelkerke", "McFadden")
  
  ### Calibration in the large
  cm_probs_oos_Lasso        <- data.frame(matrix(nrow = max_i, ncol = n_cat))
  cm_probs_oos_Ridge        <- data.frame(matrix(nrow = max_i, ncol = n_cat))
  cm_probs_oos_ML           <- data.frame(matrix(nrow = max_i, ncol = n_cat))
  cm_probs_oos_DGM          <- data.frame(matrix(nrow = max_i, ncol = n_cat))

  names(cm_probs_oos_Lasso) <- outcomes
  names(cm_probs_oos_Ridge) <- outcomes
  names(cm_probs_oos_ML)    <- outcomes
  names(cm_probs_oos_DGM)   <- outcomes
  
  ### Warning and Error handlers
  warns_gen        <- list()
  warns_ML         <- list()
  warns_Lasso      <- list()
  warns_Ridge      <- list()
  warns_cal_ML     <- list()
  warns_cal_Lasso  <- list()
  warns_cal_Ridge  <- list()
  
  errors_gen       <- list()
  errors_ML        <- list()
  errors_Lasso     <- list()
  errors_Ridge     <- list()
  errors_cal_ML    <- list()
  errors_cal_Lasso <- list()
  errors_cal_Ridge <- list()
  
    ## Warning handlers
     # For running the models
  WarningHandlerGeneric <- function(w){ 
    l <- list(w, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    warns_gen[[length(warns_gen) + 1]]     <<- l
    print(w)
    print("Warning in running mlogit on validation sample data.")
    invokeRestart("muffleWarning")
  }    
  
  WarningHandlerMlogit  <- function(w){ 
    print("warning of mlogit");
    l <- list(w, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    warns_ML[[length(warns_ML) + 1]]       <<- l
    print(w)
    print("Warning in running mlogit on development sample data.")
    invokeRestart("muffleWarning")
  }
  
  WarningHandlerLasso   <- function(w){ 
    l <- list(w, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    warns_Lasso[[length(warns_Lasso) + 1]] <<- l
    print(w)
    print("Warning in running lasso/glmnet on development sample data.")
    invokeRestart("muffleWarning")
  }
  
  WarningHandlerRidge   <- function(w){ 
    l <- list(w, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    warns_Ridge[[length(warns_Ridge) + 1]] <<- l
    print(w)
    print("Warning in running ridge/glmnet on development sample data.")
    invokeRestart("muffleWarning")
  }
  
    # For computing the calibration slopes.
  WarningHandlerCalML <- function(w){ 
    l <- list(w, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    warns_cal_ML[[length(warns_cal_ML) + 1]]       <<- l
    print(w)
    print("Warning in computing the calibration slope for ML (on validation sample data).")
    invokeRestart("muffleWarning")
  }  
  
  WarningHandlerCalLasso <- function(w){ 
    l <- list(w, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    warns_cal_Lasso[[length(warns_cal_Lasso) + 1]] <<- l
    print(w)
    print("Warning in computing the calibration slope for Lasso (on validation sample data).")
    invokeRestart("muffleWarning")
  }  
  
  WarningHandlerCalRidge <- function(w){ 
    l <- list(w, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    warns_cal_Ridge[[length(warns_cal_Ridge) + 1]] <<- l
    print(w)
    print("Warning in computing the calibration slope for Ridge (on validation sample data).")
    invokeRestart("muffleWarning")
  } 
    
    ## Error handlers
     # For running the models
  ErrorHandlerGeneric   <- function(e){ 
    errors_gen   <<- c(errors_gen, e,   i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    print(e)
    print("Error in running mlogit on validation sample data.")
    error_data <- list(sample_data = sample_data, sample_seed = sample_seed[i,], pop_data = pop_data, pop_seed = pop_seed[i, ])
    save(error_data, file = paste(proj_loc, paste(paste(paste("id_", paste(id, "_error_data_ML_val_", i, sep = ""), sep = ""), sep = ""), ".txt", sep = "")), sep = "/")
    NA
  }
  
  ErrorHandlerMlogit    <- function(e){ 
    errors_ML    <<- c(errors_ML, e,    i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    print(e)
    print("Error in running mlogit on development sample data.")
    error_data <- list(sample_data = sample_data, sample_seed = sample_seed[i,])
    save(error_data, file = paste(proj_loc, paste(paste(paste("id_", paste(id, "_error_data_ML_dev_", i, sep = ""), sep = ""), sep = ""), ".txt", sep = "")), sep = "/")
    NA
  }
  
  ErrorHandlerLasso     <- function(e){ 
    errors_Lasso <<- c(errors_Lasso, e, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    print(e)
    print("Error in running lasso/glmnet on development sample data.") 
    error_data <- list(sample_data = sample_data, sample_seed = sample_seed[i,])
    print(error_data)
    NA
  }
  
  ErrorHandlerRidge     <- function(e){ 
    errors_Ridge <<- c(errors_Ridge, e, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    print(e)
    print("Error in running ridge/glmnet on development sample data.") 
    error_data <- list(sample_data = sample_data, sample_seed = sample_seed[i,])
    save(error_data, file = paste(proj_loc, paste(paste(paste("id_", paste(id, "_error_data_Ridge_dev_", i, sep = ""), sep = ""), sep = ""), ".txt", sep = "")), sep = "/")
    NA
  }
  
     # For Computing the calibration slopes.
  errorHandlerCalML <- function(e){
    errors_cal_ML[[length(errors_cal_ML) + 1]]       <<- list(e, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    print(e)
    print("Error in computing the calibration slope for ML (on validation sample data).")
    NA
  }
  
  errorHandlerCalLasso <- function(e){
    errors_cal_Lasso[[length(errors_cal_Lasso) + 1]] <<- list(e, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    print(e)
    print("Error in computing the calibration slope for Lasso (on validation sample data).")
    NA
  }
  
  errorHandlerCalRidge <- function(e){
    errors_cal_Ridge[[length(errors_cal_Ridge) + 1]] <<- list(e, i = i, sample_seed = sample_seed[i,], pop_seed = pop_seed[i, ])
    print(e)
    print("Error in computing the calibration slope for Ridge (on validation sample data).")
    NA
  }
  
  # And Model coefficients.
  coef_names            <- CreateNames(ncol(as.matrix(coefs)))
  coefs_ML              <- array(data = NA, dim = c(length(outcomes[-1]), ncol(coefs), max_i), dimnames = list(outcome = outcomes[-1], coef = coef_names, iteration = 1:max_i)) # outcome assumes 0 is the reference category!
  coefs_Lasso           <- array(data = NA, dim = c(length(outcomes[-1]), ncol(coefs), max_i), dimnames = list(outcome = outcomes[-1], coef = coef_names, iteration = 1:max_i)) # Note that glmnet doesnt use a reference category!
  coefs_Ridge           <- array(data = NA, dim = c(length(outcomes[-1]), ncol(coefs), max_i), dimnames = list(outcome = outcomes[-1], coef = coef_names, iteration = 1:max_i)) # But I reparametrize the coefficients, so it does.
  coefs_pop             <- array(data = NA, dim = c(length(outcomes[-1]), ncol(coefs), max_i), dimnames = list(outcome = outcomes[-1], coef = coef_names, iteration = 1:max_i)) 
  se_ML                 <- array(data = NA, dim = c(length(outcomes[-1]), ncol(coefs), max_i), dimnames = list(outcome = outcomes[-1], coef = coef_names, iteration = 1:max_i)) # For ML the SE are kept, so that convergence can be inspected. As (extremely) high SE indicate non-convergence
  
  # Calibration intercepts and slopes
  cal_ML_int    <- array(data = NA, dim = c(n_cat, nrow(coefs) * 3, max_i), dimnames = list(outcome = outcomes, statistic = PolCalIntNames(n_cat), iteration = 1:max_i)) 
  cal_Lasso_int <- array(data = NA, dim = c(n_cat, nrow(coefs) * 3, max_i), dimnames = list(outcome = outcomes, statistic = PolCalIntNames(n_cat), iteration = 1:max_i)) 
  cal_Ridge_int <- array(data = NA, dim = c(n_cat, nrow(coefs) * 3, max_i), dimnames = list(outcome = outcomes, statistic = PolCalIntNames(n_cat), iteration = 1:max_i)) 
  
  cal_ML_slopes    <- array(data = NA, dim = c(n_cat, nrow(coefs) * 3, max_i), dimnames = list(outcome = outcomes, statistic = PolCalSlopesNames(n_cat), iteration = 1:max_i)) 
  cal_Lasso_slopes <- array(data = NA, dim = c(n_cat, nrow(coefs) * 3, max_i), dimnames = list(outcome = outcomes, statistic = PolCalSlopesNames(n_cat), iteration = 1:max_i)) 
  cal_Ridge_slopes <- array(data = NA, dim = c(n_cat, nrow(coefs) * 3, max_i), dimnames = list(outcome = outcomes, statistic = PolCalSlopesNames(n_cat), iteration = 1:max_i)) 
  
  if (cores >= 3) {
    library(doParallel)
    library(foreach)
    mc <- TRUE
    cl <- makePSOCKcluster(cores)
    registerDoParallel(cl)
    } else mc <- FALSE
  

  print("Starting simulation.")

  for (i in 1:max_i)
  { 
    
    ##### Create a population for this iteration.
    # set correlations among predictors for this iteration
    if (min.cor == max.cor) {
      cors <- min.cor
    } else {
      stop("Implement RandCor first")
      cors <- Randcor() #### Change this function if varying correlations are desired!
    }
    
      
    if (is.na(pop_seed[i, 1])) {
      pop_seed[i, ]   <- .Random.seed
    } else {
      .Random.seed     <<- pop_seed[i, ] 
      rep_pop_seed[i]  <- TRUE}
    
    if (n_bi_var == 0) { # Note that creating correlated binary predictors is not possible with this code!
      pop_data_list      <- MakeRetroCorData(n = N,outcome = outcomes, coefs = coefs, Brier = T, c = cors)
    } else {
      pop_data_list      <- MakeRetroData(n = N,outcome = outcomes, n_bi_var = n_bi_var, bi_prob = bi_prob, coefs = coefs, Brier = T)
    }
        
    pop_data           <- data.frame(pop_data_list$dat)
    pop_disease        <- pop_data_list$pop_disease
    Brier_oos_DGM[i, ] <- pop_data_list$pop_Brier
    DGM_probs          <- pop_data_list$y_probs
    
    # And select covariates and outcome:
    pop_x             <- as.matrix(pop_data[ , -1])
    pop_y             <- as.matrix(pop_data[ ,  1])
    
    val_mean_cor[i]   <- mean(cor(pop_x)[upper.tri(cor(pop_x))])
    
    # Convert population to long format:
    pop_data_long     <- mlogit.data(pop_data, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
    
    
    # And ML coefficients in population
    pop_model         <- withCallingHandlers(tryCatch(do.call("mlogit", args = list(formula = mfor, data = pop_data_long)),
                                                      error = ErrorHandlerGeneric), warning = WarningHandlerGeneric) 
    coefs_pop[ , , i] <- matrix(pop_model$coefficients[1:length(coefs)], nrow = nrow(coefs), ncol = ncol(coefs), byrow = F)
    
    
    ##### Draw development sample
    if (is.na(sample_seed[i, 1])) {
      sample_seed[i, ]   <- .Random.seed
    } else {
      .Random.seed       <<- sample_seed[i, ]
      rep_sample_seed[i] <- TRUE}
    
    if (n_bi_var == 0) {
      sample_data          <- data.frame(MakeRetroCorData(n = n, outcome = outcomes, coefs = coefs, Brier = F, c = cors)$dat ) # To sample from DGM
    } else { 
      sample_data          <- data.frame(MakeRetroData(n = n, outcome = outcomes, coefs = coefs, Brier = F, n_bi_var = n_bi_var, bi_prob = bi_prob)$dat ) # To sample from DGM
    }
    
    
    sample_x             <- as.matrix(sample_data[ ,-1])
    sample_y             <- sample_data[ , 1]
    
    dev_mean_cor[i]   <- mean(cor(sample_x)[upper.tri(cor(sample_x))])
    
    # And convert outcome variable to wide format.
    sample_disease <- matrix(0,nrow = sum(n), ncol = n_cat)
    for (cat in 1:n_cat) { sample_disease[sample_y == outcomes[cat], cat] <- 1} # convert sample y # Changed from 0:n_cat to 1:n_cat. correct?
    DGM_cm_probs[i, ] <- colMeans(DGM_probs)
    
    ##### Run models
      ### Run mlogit
        # Convert sample to long format:
      sample_data_long     <- mlogit.data(sample_data, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
      
      mlogit_model         <- withCallingHandlers(tryCatch(do.call("mlogit", args = list(formula = mfor, data = sample_data_long)), 
                                                           error = ErrorHandlerMlogit), warning = WarningHandlerMlogit)  
      coefs_ML[, , i]      <- matrix(mlogit_model$coefficients[1:length(coefs)], nrow = nrow(coefs), ncol = ncol(coefs), byrow = F)
      se_ML[ , , i]        <- matrix(sqrt(diag(solve(-mlogit_model$hessian)))  , nrow = nrow(coefs), ncol = ncol(coefs), byrow = F) # This warning by RStudio makes no sense. 
      sample_mlogit_fitted <- mlogit_model$probabilities
      pop_mlogit_preds     <- predict(mlogit_model, newdata = pop_data_long,    type = "probs")
      
      ### Run Lasso
      foldids              <- matrix(nrow = n_cv, ncol = sum(n))
      
      
      ### If wanted, add lambdas to Lasso's lambda sequence
      if (n_lambda_Lasso >= 100) { 
        lasso_lambdas <- glmnet(x = sample_x, y = sample_y, alpha = 1, family = "multinomial")$lambda
        llrl <- log(lasso_lambdas)
        lindexes <- 1:length(lasso_lambdas)
        lldf <- data.frame(llrl = llrl, lindexes = lindexes)
        
        ld <- lm(llrl ~ lindexes, data = lldf)$coef[2]
        
        l_add_seq <- seq(from = min(llrl), by = (ld), length.out = n_lambda_Lasso - length(lasso_lambdas) + 1)[-1]
        
        l_lambda <- sort(c(lasso_lambdas, exp(l_add_seq)), decreasing = T)
        
      } else {
        l_lambda <- NULL
      }
      
      user_l_lambda[[i]]   <- l_lambda
      lambdas_lasso        <- rep(0, n_cv)
      for (cv in 1:n_cv)
      {
        foldids[cv, ]      <- FoldidLeaveSetOut(n, nfolds)
        lasso_cv           <- withCallingHandlers(tryCatch(cv.glmnet(sample_x, sample_y, alpha = 1, foldid = foldids[cv, ], nfolds = nfolds, lambda = l_lambda,
                                                                     family = "multinomial", parallel = mc, type.measure = lasso.cv.measure), 
                                         error = ErrorHandlerLasso), warning = WarningHandlerLasso)  # alpha = 1, for lasso

        if (is.na(lasso_cv[1])) next() # Together with the following for is.na(...) next(), skips this iteration.
        
        if (lasso.which.lambda == "min")
        {if (i == 1) { print("lambda min chosen for Lasso.")}
          lambdas_lasso[cv]  <- lasso_cv$lambda.min
        } else if (lasso.which.lambda == "1se") {
          lambdas_lasso[cv]  <- lasso_cv$lambda.1se
        } else stop("Invalid lasso.which.lambda selected.")
      }
      
      if (is.na(lasso_cv[1])) next() # Together with the previous for is.na(...) next(), skips this iteration.
      
      if (!is.character(lasso.cv.lambda)) {
        lambda_lasso       <- quantile(lambdas_lasso, prob = lasso.cv.lambda)
      } else if (lasso.cv.lambda == "mean") {
        if (i == 1) { print("Mean chosen for lasso.")}
        lambda_lasso       <- mean(lambdas_lasso)    
      } else {stop("Invalid lasso.cv.lambda chosen.") }
      
      s_lambda_Lasso[i]    <- lambda_lasso

      lasso_model          <- lasso_cv$glmnet.fit
      lasso_l_seq          <- lasso_model$lambda
      
      if (n_cv == 1) {
        lambda_Lasso_rank[i]    <- which(lambda_lasso == lasso_l_seq)
        lambda_Lasso_prop[i]    <- lambda_Lasso_rank[i]/length(lasso_l_seq)
        if (lambda_lasso == lasso_l_seq[1])                   lambda_Lasso_max[i] <- T
        if (lambda_lasso == lasso_l_seq[length(lasso_l_seq)]) lambda_Lasso_min[i] <- T
      }
      
      s_l_seq_Lasso[i,]    <- c(lasso_l_seq, rep(NA, max(0, ncol(s_l_seq_Lasso) - length(lasso_l_seq))))
  
      coefs_Lasso[ , , i]  <- GlmnetAsRefCat(t(do.call(cbind, lapply(coef(lasso_model, s = lambda_lasso), matrix))))[2:length(outcomes), ] # Makes it a baseline-category model
      
      sample_lasso_fitted  <- data.frame(predict(lasso_model, newx = sample_x, s = lambda_lasso, type = "response"))
      pop_lasso_preds      <- data.frame(predict(lasso_model, newx = pop_x,    s = lambda_lasso, type = "response"))
      
      ### If wanted, add lambdas to Ridge's lambda sequence
      if (n_lambda_Ridge >= 100) { 
        ridge_lambdas <- glmnet(x = sample_x, y = sample_y, alpha = 0, family = "multinomial")$lambda
        lrl <- log(ridge_lambdas)
        indexes <- 1:length(ridge_lambdas)
        ldf <- data.frame(lrl = lrl, indexes = indexes)
        
        #min_val <- min(log(lasso_lambdas))
        d <- lm(lrl ~ indexes, data = ldf)$coef[2]
        
        #add_seq <- seq(from = min(lrl), to = min_val, by = (d))
        add_seq <- seq(from = min(lrl), by = (d), length.out = n_lambda_Ridge - length(ridge_lambdas) + 1)[-1]
        
        rlambda <- sort(c(ridge_lambdas, exp(add_seq)), decreasing = T)
        
      } else {
        rlambda <- NULL
      }
      user_r_lambda[[i]]   <- rlambda
      
      ### Run Ridge
      lambdas_ridge        <- rep(0, n_cv)
      for (cv in 1:n_cv)
      {
        ridge_cv           <- withCallingHandlers(tryCatch(cv.glmnet(sample_x, sample_y, alpha = 0, lambda = rlambda, foldid = foldids[cv, ], nfolds = nfolds, family = "multinomial", parallel = mc, type.measure = ridge.cv.measure),
                                                           error = ErrorHandlerRidge), warning = WarningHandlerRidge)                                    # alpha = 0, for ridge
        
        if (ridge.which.lambda == "min")
        { if (i == 1) { print("lambda min chosen for Ridge.")}
          lambdas_ridge[cv]  <- ridge_cv$lambda.min
        } else if (ridge.which.lambda == "1se") {
          lambdas_ridge[cv]  <- ridge_cv$lambda.1se
        } else stop("Invalid which.lambda selected.")
  }
  
  if (!is.character(ridge.cv.lambda )) {
    lambda_ridge      <- quantile(lambdas_ridge, prob = ridge.cv.lambda)
  } else if (ridge.cv.lambda == "mean") {
    if (i == 1) { print("Mean chosen for Ridge")}
    lambda_ridge       <- mean(lambdas_ridge)    
  } else {stop("Invalid cv.lambda chosen.") }
  
      s_lambda_Ridge[i]     <- lambda_ridge
      
      
      ridge_model           <- ridge_cv$glmnet.fit
      ridge_l_seq           <- ridge_model$lambda
      
      if (n_cv == 1) {
        lambda_Ridge_rank[i]    <- which(lambda_ridge == ridge_l_seq)
        lambda_Ridge_prop[i]    <- lambda_Ridge_rank[i]/length(ridge_l_seq)
        
        if (lambda_ridge == ridge_l_seq[1])                   lambda_Ridge_max[i] <- T
        if (lambda_ridge == ridge_l_seq[length(ridge_l_seq)]) lambda_Ridge_min[i] <- T
      }
      
      
      
      s_l_seq_Ridge[i, ]    <-  c(ridge_l_seq, rep(NA, max(0, ncol(s_l_seq_Ridge) - length(ridge_l_seq))))  
      
      coefs_Ridge[ , , i]   <- GlmnetAsRefCat(t(do.call(cbind, lapply(coef(ridge_model, s = lambda_ridge), matrix)))) [2:length(outcomes), ]
      sample_ridge_fitted   <- data.frame(predict(ridge_model, newx = sample_x, s = lambda_ridge, type = "response"))
      pop_ridge_preds       <- data.frame(predict(ridge_model, newx = pop_x,    s = lambda_ridge, type = "response"))

    ##### Compute Performance Measures
  
      ### Overall: Brier Score
        # Overall: Apparent Brier Score
        Brier_ws_ML[i, ]      <- colMeans((sample_mlogit_fitted - sample_disease) ^ 2)
        Brier_ws_Lasso[i, ]   <- colMeans((sample_lasso_fitted  - sample_disease) ^ 2)
        Brier_ws_Ridge[i, ]   <- colMeans((sample_ridge_fitted  - sample_disease) ^ 2)
                  
        # Overall: Actual Brier Score
        Brier_oos_ML[i, ]     <- colMeans((pop_mlogit_preds     - pop_disease) ^ 2)
        Brier_oos_Lasso[i, ]  <- colMeans((pop_lasso_preds      - pop_disease) ^ 2)
        Brier_oos_Ridge[i, ]  <- colMeans((pop_ridge_preds      - pop_disease) ^ 2)
        
      ### Overall: Accuracy of predicted probabilities using DGM's probabilities
        mse_probs_oos_ML[i, ]      <- colMeans((pop_mlogit_preds     - DGM_probs) ^ 2)
        mse_probs_oos_Lasso[i, ]   <- colMeans((pop_lasso_preds      - DGM_probs) ^ 2)
        mse_probs_oos_Ridge[i, ]   <- colMeans((pop_ridge_preds      - DGM_probs) ^ 2) 
    
    
      ### Overall: R^2 : Cox, Nagelkerke and McFadden
    
        # Overall: Apparent R^2
        L0_sample             <- mult_log_L(matrix(n/n_total, nrow = n_total, ncol = n_cat, byrow = T), sample_y)
        R2_ws_ML[i, ]         <- mult_R2(sample_mlogit_fitted, sample_y, L0_sample, n_total, outcomes, n_cat)
        R2_ws_Lasso[i, ]      <- mult_R2(sample_lasso_fitted,  sample_y, L0_sample, n_total, outcomes, n_cat)
        R2_ws_Ridge[i, ]      <- mult_R2(sample_ridge_fitted,  sample_y, L0_sample, n_total, outcomes, n_cat)
        
        # Overall: Out of sample R^2
        L0_pop                <- mult_log_L(matrix(N/N_total, nrow = N_total, ncol = n_cat, byrow = T), pop_y) # Likelihood of the null model, in the development sample
        R2_oos_ML[i, ]        <- mult_R2(pop_mlogit_preds, pop_y, L0_pop, N_total, outcomes, n_cat)
        R2_oos_Lasso[i, ]     <- mult_R2(pop_lasso_preds , pop_y, L0_pop, N_total, outcomes, n_cat)
        R2_oos_Ridge[i, ]     <- mult_R2(pop_ridge_preds , pop_y, L0_pop, N_total, outcomes, n_cat)
        R2_oos_DGM[i, ]       <- mult_R2(DGM_probs       , pop_y, L0_pop, N_total, outcomes, n_cat)
        
    
      ### Discrimination: 
        
          # Pairwise c statistic, With the M-index
        M_ws_ML[ , , i]    <- MyM(sample_mlogit_fitted, sample_y, outcomes)$pairwise_M
        M_ws_Lasso[ , , i] <- MyM(sample_lasso_fitted,  sample_y, outcomes)$pairwise_M
        M_ws_Ridge[ , , i] <- MyM(sample_ridge_fitted,  sample_y, outcomes)$pairwise_M
        
        M_oos_ML[ , , i]    <- MyM(pop_mlogit_preds,     pop_y,    outcomes)$pairwise_M
        M_oos_Lasso[ , , i] <- MyM(pop_lasso_preds,      pop_y,    outcomes)$pairwise_M
        M_oos_Ridge[ , , i] <- MyM(pop_ridge_preds,      pop_y,    outcomes)$pairwise_M
        M_oos_DGM[ , , i]   <- MyM(DGM_probs,            pop_y,    outcomes)$pairwise_M
        
        
          # Pairwise OnevsRest c stat
        OneVsRestC_ws_ML[ , , i]      <- OneVsRestCBi(sample_mlogit_fitted, sample_disease)$i_vs_rest
        OneVsRestC_ws_Lasso[ , , i]   <- OneVsRestCBi(sample_lasso_fitted,  sample_disease)$i_vs_rest
        OneVsRestC_ws_Ridge[ , , i]   <- OneVsRestCBi(sample_ridge_fitted,  sample_disease)$i_vs_rest
        
        OneVsRestC_oos_ML[ , , i]     <- OneVsRestCBi(pop_mlogit_preds, pop_disease)$i_vs_rest
        OneVsRestC_oos_Lasso[ , , i]  <- OneVsRestCBi(pop_lasso_preds,  pop_disease)$i_vs_rest
        OneVsRestC_oos_Ridge[ , , i]  <- OneVsRestCBi(pop_ridge_preds,  pop_disease)$i_vs_rest
        OneVsRestC_oos_DGM[ , , i]    <- OneVsRestCBi(DGM_probs,        pop_disease)$i_vs_rest
        
        
        ### Polytomous Discrimination Index: PDI
        if (computePDI) # To save some time during debugging
        {
        if (mc)
        {
          fitted <- list("ML" = sample_mlogit_fitted, "Lasso" = sample_lasso_fitted, "Ridge" = sample_ridge_fitted)
          preds  <- list("ML" = pop_mlogit_preds,     "Lasso" = pop_lasso_preds,     "Ridge" = pop_ridge_preds    ,  "DGM" = DGM_probs)
          
          ap <- foreach(i = 1:length(fitted), .combine = rbind, .export = c("pdi3i", "PdiDF")) %dopar% (PdiDF(sample_y, fitted[[i]], method_name = names(fitted)[i]))
          te <- foreach(i = 1:length(preds),  .combine = rbind, .export = c("pdi3i", "PdiDF")) %dopar% (PdiDF(pop_y,     preds[[i]], method_name = names(preds)[i]))
          
          PDI_ws_ML[i, ]        <- subset(ap, method_name == "ML"   )[1:5]
          PDI_ws_Lasso[i, ]     <- subset(ap, method_name == "Lasso")[1:5]
          PDI_ws_Ridge[i, ]     <- subset(ap, method_name == "Ridge")[1:5]
          
          PDI_oos_ML[i, ]       <- subset(te, method_name == "ML"   )[1:5]
          PDI_oos_Lasso[i, ]    <- subset(te, method_name == "Lasso")[1:5]
          PDI_oos_Ridge[i, ]    <- subset(te, method_name == "Ridge")[1:5]
          PDI_oos_DGM[i, ]      <- subset(te, method_name == "DGM"  )[1:5]
          
        } else
        {
          PDI_ws_ML[i, ]        <- PdiDF(sample_y, sample_mlogit_fitted)[1:5]
          PDI_ws_Lasso[i, ]     <- PdiDF(sample_y, sample_lasso_fitted)[1:5]
          PDI_ws_Ridge[i, ]     <- PdiDF(sample_y, sample_ridge_fitted)[1:5]
          
          PDI_oos_ML[i, ]       <- PdiDF(pop_y, pop_mlogit_preds)[1:5]
          PDI_oos_Lasso[i, ]    <- PdiDF(pop_y, pop_lasso_preds)[1:5]
          PDI_oos_Ridge[i, ]    <- PdiDF(pop_y, pop_ridge_preds)[1:5]
          PDI_oos_DGM[i, ]      <- PdiDF(pop_y, DGM_probs)[1:5]
        }
        }
        
  
        ### Calibration Slopes using Multinomial calibration.
        if (compute_cal) {# To save some time during debugging
          library(VGAM)
          cal_ML_temp    <- PolCalAll(pop_y, pop_mlogit_preds, i = i, id = id, method = "ML")  
          cal_Lasso_temp <- PolCalAll(pop_y, pop_lasso_preds , i = i, id = id, method = "Lasso")  
          cal_Ridge_temp <- PolCalAll(pop_y, pop_ridge_preds , i = i, id = id, method = "Ridge") 
          detach("package:VGAM", unload = T)
          
          cal_ML_int[ , , i]    <- UnName(cal_ML_temp$intercepts)
          cal_ML_slopes[ , , i] <- UnName(cal_ML_temp$slopes)

          cal_Lasso_int[ , , i]    <- UnName(cal_Lasso_temp$intercepts)
          cal_Lasso_slopes[ , , i] <- UnName(cal_Lasso_temp$slopes)

          cal_Ridge_int[ , , i]    <- UnName(cal_Ridge_temp$intercepts)
          cal_Ridge_slopes[ , , i] <- UnName(cal_Ridge_temp$slopes)
          }
        
        # Assessing the predicted prevalence
        cm_probs_oos_ML[i, ]    <- colMeans(pop_mlogit_preds)
        cm_probs_oos_Lasso[i, ] <- colMeans(pop_lasso_preds)
        cm_probs_oos_Ridge[i, ] <- colMeans(pop_ridge_preds)
        cm_probs_oos_DGM[i, ]   <- colMeans(DGM_probs)
        
        if (i %% 10 == 0 & i != 0) print(paste(i, " iterations complete."))
        if ( (proc.time() - begin)[3] > max.time )  break;
  }
  
  #invisible(close(pb))
  if (mc) stopCluster(cl)

  return(list("C_Pairwise_M" = list("ws_ML"    = M_ws_ML,    "oos_ML"    = M_oos_ML,
                                    "ws_Lasso" = M_ws_Lasso, "oos_Lasso" = M_oos_Lasso,
                                    "ws_Ridge" = M_ws_Ridge, "oos_Ridge" = M_oos_Ridge),
    
              "One_vs_restk_c" = list("ws_ML"    = OneVsRestC_ws_ML,    "oos_ML"    = OneVsRestC_oos_ML,
                                      "ws_Lasso" = OneVsRestC_ws_Lasso, "oos_Lasso" = OneVsRestC_oos_Lasso,
                                      "ws_Ridge" = OneVsRestC_ws_Ridge, "oos_Ridge" = OneVsRestC_oos_Ridge),
              
              "PDI"        = list("ws_ML"      = PDI_ws_ML,        "oos_ML"    = PDI_oos_ML,
                                  "ws_Lasso"   = PDI_ws_Lasso,     "oos_Lasso" = PDI_oos_Lasso,
                                  "ws_Ridge"   = PDI_ws_Ridge,     "oos_Ridge" = PDI_oos_Ridge),
              
              "Brier"      = list("ws_ML"      = Brier_ws_ML    , "oos_ML"    = Brier_oos_ML,
                                  "ws_Lasso"   = Brier_ws_Lasso , "oos_Lasso" = Brier_oos_Lasso,
                                  "ws_Ridge"   = Brier_ws_Ridge , "oos_Ridge" = Brier_oos_Ridge),
              
              "mse_probs"  = list(                               "oos_ML"    = mse_probs_oos_ML    ,
                                                                 "oos_Lasso" = mse_probs_oos_Lasso ,
                                                                 "oos_Ridge" = mse_probs_oos_Ridge ),

              "R2"         = list("ws_ML"      = R2_ws_ML       , "oos_ML"    = R2_oos_ML,
                                  "ws_Lasso"   = R2_ws_Lasso    , "oos_Lasso" = R2_oos_Lasso,
                                  "ws_Ridge"   = R2_ws_Ridge    , "oos_Ridge" = R2_oos_Ridge),
              
              "coefs"      = list("Population" = coefs_pop      , "ML"        = coefs_ML        , 
                                  "Lasso"      = coefs_Lasso    , "Ridge"     = coefs_Ridge     ,
                                  "SE_ML"      = se_ML), 
              
              "cal_int"    = list("oos_ML"     = cal_ML_int     , "oos_Lasso" = cal_Lasso_int,
                                  "oos_Ridge"  = cal_Ridge_int ),
              
              "cal_slopes" = list("oos_ML"     = cal_ML_slopes  , "oos_Lasso" = cal_Lasso_slopes,
                                  "oos_Ridge"  = cal_Ridge_slopes),
              
              "warn_list"  = list("Generic"          = warns_gen       , "ML"               = warns_ML,
                                  "Lasso"            = warns_Lasso     , "Ridge"            = warns_Ridge,
                                  "cal_ML"           = warns_cal_ML    , "cal_Lasso"        = warns_cal_Lasso,
                                  "cal_Ridge"        = warns_cal_Ridge),
              
              "error_list" = list("Generic"          = errors_gen      , "ML"               = errors_ML, 
                                  "Lasso"            = errors_ML       , "Ridge"            = errors_Ridge,
                                  "cal_ML"           = errors_cal_ML   , "cal_Lasso"        = errors_cal_Lasso,
                                  "cal_Ridge"        = errors_cal_Ridge),
         
              "params"     = list("max_i"            = max_i           , "n"                = n, 
                                  "N"                = N               , "n_bi_var"         = n_bi_var,  
                                  "coefs"            = coefs           , "bi_prob"          = bi_prob,
                                  "id"               = id              , "repetition"       = repetition,
                                  "mfor"             = mfor            , "n_cv"             = n_cv,
                                  "min_cor"          = min_cor         , "max_cor"          = max.cor),

              "DGM"        = list("Brier"            = Brier_oos_DGM   , "R2"               = R2_oos_DGM,
                                  "M-index"          = M_oos_DGM       , "OneVsRestC"       = OneVsRestC_oos_DGM,
                                  "PDI"              = PDI_oos_DGM),
              
              "seed"       = list("pop_seed"         = pop_seed        , "sample_seed"      = sample_seed,  
                                  "rep_pop_seed"     = rep_pop_seed    , "rep_sample_seed"  = rep_sample_seed,
                                  "end_seed"         = .Random.seed),
              
              "lambda"     = list("min"              = list(
                                    "Lasso"            = s_lambda_Lasso,    "Ridge"            = s_lambda_Ridge),
                                  "user_sequence"    = list(
                                    "Lasso"            = user_l_lambda,     "Ridge"            = user_r_lambda),
                                  "glmnet_sequence"  = list(
                                    "Lasso"            = s_l_seq_Lasso,     "Ridge"            = s_l_seq_Ridge),
                                  "rank"             = list(
                                    "Lasso"            = lambda_Lasso_rank, "Ridge"            = lambda_Ridge_rank),
                                  "prop_rank"        = list(
                                    "Lasso"            = lambda_Lasso_prop, "Ridge"            = lambda_Ridge_prop),
                                  "max_rank"         = list(
                                    "Lasso"            = lambda_Lasso_max,  "Ridge"            = lambda_Ridge_max),
                                  "min_rank"         = list(
                                    "Lasso"            = lambda_Lasso_min,  "Ridge"            = lambda_Ridge_min))
                                  
              
              #,"DGM_probs" = DGM_probs, "pop_mlogit_preds" = pop_mlogit_preds, "pop_disease" = pop_disease # for debugging only; too big to save each in final simulation
              , "cm_probs" = list("oos_ML"    = cm_probs_oos_ML              ,  "oos_Lasso" = cm_probs_oos_Lasso,
                                  "oos_Ridge" = cm_probs_oos_Ridge           ,  "oos_DGM"   = cm_probs_oos_DGM)
              
              , "cor"      = list("dev"       = dev_mean_cor                 ,  "val"       = val_mean_cor)
              ) 
         )
}