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

# what = the simulation
# n.median.bs = number of bootstraps for median cal slopes /intercepts
LoadData <- function(what = "main", n.median.bs = 1e5)
{
  n_median_bootstraps <- n.median.bs
  
  SE_conv_limit <- log(50)
  bs_test <- 0
  total_reps <- 2000 # Used ONLY for arrays of full data. For all other statistics it is generated from data.
  
  if (what == "main") {
    # Main simulation
    bins    <- c(0)
    cors    <- c(0)
    ids     <- 1:63
    failed  <- c()
    reps    <- 0:9
    effects <- c(1)
    
    sims    <- ids          # sims and n_sims are dummy variables
  } else if (what == "cor") { 
  id      <- 3 # as n_coef, epv and prev are fixed
  bins    <- c(0)
  cors    <- c(0, 2, 3, 5, 7, 9)
  failed  <- c()
  reps    <- 1:10
  effects <- c(1)
  
  sims   <- cors
  } else if (what == "bin") {
    id      <- 3 # as n_coef, epv and prev are fixed
    bins    <- c(0,1)
    cors    <- c(0)
    failed  <- c()
    reps    <- 0:10 # Thanks to my mistake in naming some files, some unexisting files are selected, this throws some warnings which can be ignored
    effects <- c(1)
    
    sims   <- bins
    } else if (what == "effect") {
      id      <- 3 # as n_coef, epv and prev are fixed
      bins    <- c(0)
      cors    <- c(0)
      failed  <- c()
      reps    <- 0:10 # Thanks to my mistake in naming some files, some unexisting files are selected, this throws some warnings which can be ignored
      effects <- c(0, 1)
      
      sims   <- effects
      } else (stop("Nothing to load selected."))
  
  bin    <- bins[1]
  cor_id <- cors[1]
  effect <- effects[1]
  n_sims <- length(sims)
    
    perf_names <- c("ws_ML", "oos_ML", "ws_Lasso", "oos_Lasso", "ws_Ridge", "oos_Ridge")
    par_names  <- c("id", "EPV1", "EPV2", "EPV3", "n_coefs", "prev", "I", "prop NA"
                    , "n_total"
                    , "cor", "bin", "effect")
    
    # Useful functions
    source("Coerce.R")
    source("Merge.R", local = T)
    source("SE.R", local = T) 
    source("Tables.R", local = T)
    
    # For loading.
    out <- list()
    
    # parameter:
    n_pars <- length(par_names)
    
    # out is to be returned:
    out$max_rank   <- out$min_rank <- data.frame(matrix(NA, nrow = n_sims, ncol = 2 + n_pars))
    
    out$Brier <- out$Brier_SE <- out$Brier_pct  <- out$Brier_pct_SE <- data.frame(matrix(NA, nrow = n_sims, ncol = 7 + n_pars))
    
    out$Cox <- out$McFadden <- out$Nagelkerke <- out$Nagelkerke_SE <- out$Nagelkerke_pct <- out$Nagelkerke_pct_SE <- 
      data.frame(matrix(NA, nrow = n_sims, ncol = 7 + n_pars))
    
    out$pdi <- out$pdi_SE <- out$pdi_spec_small <- out$pdi_spec_small_middle <- out$pdiw <- out$pdi_pct <- out$pdi_pct_SE <- 
      out$M <- out$M_SE <- out$OneVsRestC <- data.frame(matrix(NA, nrow = n_sims, ncol = 7 + n_pars))
    
    out$cal_int_median <- out$cal_slo_median <- out$cal_int_mean  <- out$cal_slo_mean  <- data.frame(matrix(NA, nrow = n_sims, ncol = 19 + n_pars))
    
    out$cal_int_median_SE <- out$cal_slo_median_SE <- out$cal_int_mean_SE <- out$cal_slo_mean_SE  <- 
      data.frame(matrix(NA, nrow = n_sims, ncol = 18 + 18 + n_pars))
    
    out$all_cal_slo              <- array(dim = c(n_sims, length(c(par_names, 1:18)), total_reps), dimnames = list("scenario #" = sims, variable = c(par_names, 1:18), i = 1:total_reps))
    out$all_Brier <- out$all_pdi <- array(dim = c(n_sims, length(c(par_names, 1:7)),  total_reps), dimnames = list("scenario #" = sims, variable = c(par_names, 1:7) , i = 1:total_reps))
    
    out$mse_probs <- out$rmse_probs <- data.frame(matrix(NA, nrow = n_sims, ncol = 3 + n_pars))
    coefs_mse_pop <- coefs_mse_DGM <- coefs_bias_list <- coefs_last <- coefs_max_ML_SE <- list()
    out$non_convergence <- rep(NA, n_sims)
    
    out$mean_coefs_rmse_pop <- data.frame(matrix(nrow = n_sims, ncol = 3 + n_pars)) 
    out$mean_coefs_rmse_DGM <- data.frame(matrix(nrow = n_sims, ncol = 4 + n_pars))
    out$mean_coefs_bias <- out$mean_coefs_last <- data.frame(matrix(nrow = n_sims, ncol = 8 + n_pars)) 
    
    out$errors <- out$warns <- data.frame(matrix(NA, nrow = n_sims, ncol = 7))
    out$cal_select <- list()
    out$n_cat <- rep(NA, n_sims)
    
    source("fileParams.R", local = T) # Same file as is used for generating the simulation scenarios.
    
    LoadEH <- function(e){   return(NA) }
    
    n_sim <- 0
    for (s in sims) {    # dummy variables, of which the meanings change in different simulation settings
      n_sim <- n_sim + 1 # dummy variable to keep track of rows of final dataframes
      
      if (what == "main")   id     <- ids[n_sim]     # changes only in main analysis
      if (what == "cor")    cor_id <- cors[n_sim]    # changes only in sensitivity analysis
      if (what == "bin")    bin    <- bins[n_sim]    # ""
      if (what == "effect") effect <- effects[n_sim] # ""
      
      M_list <- M_pop_list <- c_One_list <- c_One_pop_list <- PDI_list <- PDI_pop_list <- Brier_list <- Brier_pop_list <- 
        mse_probs_list <- R2_list <- R2_pop_list <- est_coefs_list <- cal_int_list <- cal_slo_list <- warn_list <- error_list <- 
        max_rank_list <- min_rank_list <- list()
      
      ### Load the files
      for (r in reps) {
        if (any(cors > 0)) {       results_folder <- "cor results"
        } else if (any(bin > 0)) { results_folder <- "bin results"
        } else if (effect != 1)  { results_folder <- "effect results"
        } else                     results_folder <- "results"
        
        results_path <- paste(paste(getwd(), "Output", sep = "/"), results_folder, sep = "/")
        
        file    <- sub("%%%", id, "results_prev@@@_coefs###_epv!!!_id%%%_repetition%#!@_cor@#!_bin#!%" )
        file    <- sub("%#!@", r, file)
        
        n_coefs <- n_coefs_vec[id]
        
        file    <- sub("###", n_coefs, file)
        epv     <- rep(epv1, 9)[id]
        file    <- sub("!!!", epv, file)
        prev    <- 10 * n_mat[id, 3] / n_mat[id, 1]
        file    <- sub("@@@", prev, file)
        
        file    <- sub("@#!", cor_id, file)         
        file    <- sub("#!%", bin, file)        
        
        file <- paste(results_path, file, sep = "/")
        if (effect != 1) file <- paste(file, "_eff0", sep = "")
        file <- paste(file, ".txt", sep = "")
        
        
        spec_pdi_cat <- if (prev == 1.25) 3 else 1 
        
        if (!(any(s == failed))) {
          loaded <- tryCatch(load(file = file), error = LoadEH)
          if (is.na(loaded)) next
          
          ## Put in (temporary) lists
          Brier_list[[r + 1]]      <- results$Brier
          Brier_pop_list[[r + 1]]  <- results$DGM$Brier
          M_list[[r + 1]]          <- results$C_Pairwise_M
          M_pop_list[[r + 1]]      <- results$DGM$`M-index`
          PDI_list[[r + 1]]        <- results$PDI
          PDI_pop_list[[r + 1]]    <- results$DGM$PDI
          c_One_list[[r + 1]]      <- results$One_vs_restk_c
          c_One_pop_list[[r + 1]]  <- results$DGM$OneVsRestC
          R2_list[[r + 1]]         <- results$R2
          R2_pop_list[[r + 1]]     <- results$DGM$R2
          
          cal_int_list[[r + 1]]    <- results$cal_int
          cal_slo_list[[r + 1]]    <- results$cal_slopes
          
          error_list[[r + 1]]      <- results$error_list
          warn_list[[r + 1]]       <- results$warn_list
          
          est_coefs_list[[r + 1]]  <- results$coefs
          
          max_rank_list[[r + 1]]   <- results$lambda$max_rank
          min_rank_list[[r + 1]]   <- results$lambda$min_rank
          
          mse_probs_list[[r + 1]]  <- results$mse_probs
          
        } # end if !failed
      
      } # end reps
      temp_Brier      <- MergeLiLiDf(Brier_list)
      temp_pop_Brier  <- MergeLiDf(Brier_pop_list)
      
      temp_R2         <- MergeLiLiDf(R2_list)
      temp_pop_R2     <- MergeLiDf(R2_pop_list)
      
      temp_M          <- MergeLiLiAr(M_list)
      temp_pop_M      <- MergeLiAr(M_pop_list)
      temp_c_One      <- MergeLiLiAr(c_One_list)
      temp_pop_c_One  <- MergeLiAr(c_One_pop_list)
      
      temp_PDI        <- MergeLiLiDf(PDI_list)
      temp_pop_PDI    <- MergeLiDf(PDI_pop_list)
      
      temp_cal_int    <- MergeLiLiAr(cal_int_list)
      temp_cal_slo    <- MergeLiLiAr(cal_slo_list)
      
      temp_coefs      <- MergeLiLiAr(est_coefs_list)
      
      temp_max_rank   <- MergeLiLiLi(max_rank_list)
      temp_min_rank   <- MergeLiLiLi(min_rank_list)
      
      temp_warns      <- MergeLiLiLi(warn_list)
      temp_errors     <- MergeLiLiLi(error_list)
      
      
      temp_sum_i      <- sum(!is.na(temp_Brier$ws_ML[ ,1]))
      temp_prop_NA    <- mean(is.na(temp_Brier$ws_ML[ ,1]))  
      
      temp_mse_probs  <- MergeLiLiDf(mse_probs_list)
      
      ## Data organization into data.frames.
      type.pred <- if (bin) "Binary" else "Normal"
      parameters <- c(results$params$id[1]                                                                                        # ID of simulation
                      , min(results$params$n) / (length(as.matrix(results$params$coefs)) - nrow(as.matrix(results$params$coefs))) # EPV 1 = Events per parameter in equation (except for the intercept)
                      , min(results$params$n) / (ncol(as.matrix(results$params$coefs)) - 1)                                       # EPV 2 = Events per predictor
                      , mean(results$params$n[-which.max(results$params$n)]) 
                      / (length(as.matrix(results$params$coefs)) - nrow(as.matrix(results$params$coefs)))                         # EPV 3 = Mean events in non-largest cat / # parameters (exc. int.)
                      , n_coefs, prev/10
                      , temp_sum_i, temp_prop_NA
                      , sum(results$params$n)
                      , cor_id
                      , bin
                      , effect
      )   
      outcomes <- seq(from = 0, to = nrow(results$params$coefs))
      pp <- matrix(parameters, ncol = total_reps, nrow = length(parameters))
      out$n_cat[n_sim] <- n_cat    <- nrow(results$params$coefs) + 1
      
      ## Lambdas rank 
      out$max_rank[n_sim, ]  <- c(parameters, unlist(lapply(temp_max_rank, sum)))
      out$min_rank[n_sim, ]  <- c(parameters, unlist(lapply(temp_min_rank, sum)))
      
      ## Brier Score
      temp_Brier$DGM <- data.frame(temp_pop_Brier) 
      
      out$Brier[n_sim, ]    <- c(parameters 
                                 ,lapply(lapply(temp_Brier, colMeans, na.rm = T), sum) )
      out$Brier_SE[n_sim, ] <- c(parameters
                                 , data.frame(lapply(lapply(temp_Brier, rowSums), sd, na.rm = T))/sqrt(temp_sum_i) )
      
      Brier_sums <- lapply(temp_Brier, rowSums)
      out$Brier_pct[n_sim, ]    <- c(parameters, colMeans(Convert(Brier_sums), na.rm = T))
      out$Brier_pct_SE[n_sim, ] <- c(parameters, ColSEMs(Convert(Brier_sums)))
      out$all_Brier[n_sim, , ] <- rbind(pp, t(sapply(temp_Brier, rowSums, na.rm = F)))
      
      ## c_stat
      # Pairwise
      out$M[n_sim, ]   <- c(parameters, data.frame(lapply(temp_M, mean, na.rm = T)), mean(temp_pop_M, na.rm = T)) 
      out$M_SE[n_sim, ]   <- c(parameters, lapply(temp_M, ArrayDimMeansSE, MARGIN = 3, n = temp_sum_i), ArrayDimMeansSE(temp_pop_M, MARGIN = 3, n = temp_sum_i)  ) 
      C_temp       <- data.frame(lapply(temp_c_One, mean, na.rm = T))
      out$OneVsRestC[n_sim, ]   <- c(parameters ,C_temp, mean(temp_pop_c_One, na.rm = T)) 
      
      # PDI
      temp_PDI$DGM            <- data.frame(temp_pop_PDI) 
      PDI                     <- data.frame(lapply(temp_PDI, FUN = colMeans, na.rm = T))
      
      out$pdi[n_sim, ]            <- c(parameters, PDI[4, ]) 
      out$pdi_SE[n_sim, ]         <- c(parameters, sapply(temp_PDI, ColSEMs)[4,])  
      
      out$pdiw[n_sim, ]           <- c(parameters, PDI[5, ]) 
      out$pdi_spec_small[n_sim, ] <- c(parameters, PDI[spec_pdi_cat, ]) 
      out$pdi_spec_small_middle[n_sim, ] <- c(parameters, PDI[which.min(results$params$n),]) 
      
      temp_all_pdi <- sapply(temp_PDI, '[', 4)
      names(temp_all_pdi) <- gsub(".pdi", "", names(temp_all_pdi))
      out$all_pdi[n_sim, , ] <- rbind(pp, t(as.matrix(data.frame(temp_all_pdi))))
      
      out$pdi_pct[n_sim, ]    <- c(parameters, colMeans(Convert(temp_all_pdi), na.rm = T))
      out$pdi_pct_SE[n_sim, ] <- c(parameters, ColSEMs(Convert(temp_all_pdi)))
      
      ## MSE of probabilities
      # Note that contrary to the outcomes, the true values of the probabilities of these outcomes are unknown and therefore estimated.
      # These are therefore estimates as well:
      out$mse_probs[n_sim, ]  <- c(parameters, lapply(lapply(temp_mse_probs, colMeans, na.rm = T), sum))          
      out$rmse_probs[n_sim, ] <- c(parameters, sqrt(sapply(lapply(temp_mse_probs, colMeans, na.rm = T), sum)))    
      
      # R^2
      R2    <- data.frame(lapply(temp_R2, FUN = colMeans, na.rm = T)) 
      R2_SE <- sapply(temp_R2, ColSEMs, na.rm = T)
      
      temp_Nagelkerke <- lapply(temp_R2, `[[`, which(names(temp_R2[[1]]) == "Nagelkerke"))
      temp_Nagelkerke$"DGM" <- temp_pop_R2$Nagelkerke
      out$Nagelkerke_pct[n_sim, ]    <- c(parameters, colMeans(Convert(temp_Nagelkerke), na.rm = T))
      out$Nagelkerke_pct_SE[n_sim, ] <- c(parameters, ColSEMs(Convert(temp_Nagelkerke)))
      
      out$Cox[n_sim, ] <- c(parameters, R2[1, ], mean(temp_pop_R2$Cox, na.rm = T))  
      out$Nagelkerke[n_sim, ] <- c(parameters, R2[2, ], mean(temp_pop_R2$Nagelkerke, na.rm = T)) 
      out$Nagelkerke_SE[n_sim, ] <- c(parameters, R2_SE[2, ], sd(temp_pop_R2$Nagelkerke, na.rm = T)/sqrt(temp_sum_i)) 
      out$McFadden[n_sim, ] <- c(parameters, R2[3, ], mean(temp_pop_R2$McFadden, na.rm = T)) 
      
      ## Coefficients.
      # First organized them into lists, but later decided data frames might be better. So here's both:
      # Lists
      
      coefs_mse_pop[[n_sim]] <- list(
        "ML"    = matrix(apply((temp_coefs$ML[ , , ]    - temp_coefs$Population[ , ,]) ^ 2, MARGIN = c(1,2), mean, na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population)),  # mse of coefs
        "Lasso" = matrix(apply((temp_coefs$Lasso[ , , ] - temp_coefs$Population[ , ,]) ^ 2, MARGIN = c(1,2), mean, na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population)),
        "Ridge" = matrix(apply((temp_coefs$Ridge[ , , ] - temp_coefs$Population[ , ,]) ^ 2, MARGIN = c(1,2), mean, na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population)))
      
      coefs_mse_DGM[[n_sim]] <- list(
        "Val. data ML" = matrix(rowMeans(apply(temp_coefs$Population[ , , ], 3, `-`, results$params$coefs) ^ 2, na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population)),
        "ML"           = matrix(rowMeans(apply(temp_coefs$ML[ , , ],         3, `-`, results$params$coefs) ^ 2, na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population)),
        "Lasso"        = matrix(rowMeans(apply(temp_coefs$Lasso[ , , ],      3, `-`, results$params$coefs) ^ 2, na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population)),
        "Ridge"        = matrix(rowMeans(apply(temp_coefs$Ridge[ , , ],      3, `-`, results$params$coefs) ^ 2, na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population)))
      
      coefs_bias_list[[n_sim]] <- list(
        "Validation" = matrix(rowMeans(apply(temp_coefs$Population[ ,-1, ], 3, `-`, results$params$coefs[ ,-1]), na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population) - 1),
        "ML"         = matrix(rowMeans(apply(temp_coefs$ML[ ,-1, ],         3, `-`, results$params$coefs[ ,-1]), na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population) - 1),
        "Lasso"      = matrix(rowMeans(apply(temp_coefs$Lasso[ ,-1, ],      3, `-`, results$params$coefs[ ,-1]), na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population) - 1),
        "Ridge"      = matrix(rowMeans(apply(temp_coefs$Ridge[ ,-1, ],      3, `-`, results$params$coefs[ ,-1]), na.rm = T), nrow = nrow(temp_coefs$Population), ncol = ncol(temp_coefs$Population) - 1))
      
      coefs_last[[n_sim]] <- list(
        "Validation" = matrix(rowMeans(temp_coefs$Population[ , n_coefs + 1, ], na.rm = T), nrow = nrow(temp_coefs$Population), ncol = 1),
        "ML"         = matrix(rowMeans(temp_coefs$ML        [ , n_coefs + 1, ], na.rm = T), nrow = nrow(temp_coefs$Population), ncol = 1),
        "Lasso"      = matrix(rowMeans(temp_coefs$Lasso[ , n_coefs + 1, ], na.rm = T)     , nrow = nrow(temp_coefs$Population), ncol = 1),
        "Ridge"      = matrix(rowMeans(temp_coefs$Ridge[ , n_coefs + 1, ], na.rm = T)     , nrow = nrow(temp_coefs$Population), ncol = 1))
      
      coefs_max_ML_SE[[n_sim]] <- c(apply(temp_coefs$SE_ML, MARGIN = 3, FUN = max, na.rm = T ))
      out$non_convergence[n_sim] <- mean(coefs_max_ML_SE[[n_sim]] > SE_conv_limit)
      
      # Dataframes 
      out$mean_coefs_rmse_pop[n_sim, ] <- c(parameters, sapply(coefs_mse_pop[[n_sim]], mean) ) # Also include intercept, and sum over outcomes.
      out$mean_coefs_rmse_DGM[n_sim, ] <- c(parameters, sapply(coefs_mse_DGM[[n_sim]], mean) ) 
      out$mean_coefs_bias[n_sim, ]     <- c(parameters, unlist(lapply(coefs_bias_list[[n_sim]], rowMeans, na.rm = T)) ) # = mean bias of coefficients
      out$mean_coefs_last[n_sim, ]     <- c(parameters, unlist(lapply(coefs_bias_list[[n_sim]], rowMeans, na.rm = T)) ) # = bias of last coefficient
      
      # Calibration slopes.
      out$cal_select[[length(out$cal_select) + 1]] <- cal_select
      success_cal_int <- lColApply(temp_cal_int, FUN = CountSuccess, selection = cal_select)
      success_cal_slo <- lColApply(temp_cal_slo, FUN = CountSuccess, selection = cal_select)
      min_success_cal_int <- min(success_cal_int)
      min_success_cal_slo <- min(success_cal_slo)
      
      # All calibration slopes
      out$cal_int_median[n_sim, ] <- c(parameters, lColApply(temp_cal_int, FUN = median, selection = cal_select), min_success_cal_int)
      out$cal_slo_median[n_sim, ] <- c(parameters, lColApply(temp_cal_slo, FUN = median, selection = cal_select), min_success_cal_slo)
      out$cal_int_mean[n_sim, ] <- c(parameters, lColApply(temp_cal_int, FUN = mean, selection = cal_select), min_success_cal_int)
      out$cal_slo_mean[n_sim, ] <- c(parameters, lColApply(temp_cal_slo, FUN = mean, selection = cal_select), min_success_cal_slo)
      
      # SEs of means of calibration slopes only
      out$cal_int_mean_SE[n_sim, ]   <- c(parameters, lColApply(temp_cal_int, FUN = sd, selection = cal_select) / sqrt(success_cal_int), success_cal_int )
      out$cal_slo_mean_SE[n_sim, ]   <- c(parameters, lColApply(temp_cal_slo, FUN = sd, selection = cal_select) / sqrt(success_cal_slo), success_cal_slo )
      
      liar_cs <- lapply(temp_cal_slo, '[', ref = 1:3, slope = c(1,4), iter = 1:2000) # positionally matched to the dimensions, to select slopes only
      
      # Coerce the slopes to a matrix, and add to array
      lima_cs <- lapply(liar_cs, function(X) apply(X, MARGIN = 3, FUN = AsVector, byrow = T))  
      lima_cs <- lapply(lima_cs, AddCalNames, outcomes = outcomes)
      mat_cs  <- MergeLiDf(lima_cs, combine.row.names = T)
      out$all_cal_slo[n_sim, , ] <- rbind(pp, mat_cs)
      
      out$cal_int_median_SE[n_sim, ] <- c(parameters, lColApply(temp_cal_int, FUN = SEMedian, selection = cal_select, n = n_median_bootstraps) , success_cal_int) 
      out$cal_slo_median_SE[n_sim, ] <- c(parameters, lColApply(temp_cal_slo, FUN = SEMedian, selection = cal_select, n = n_median_bootstraps) , success_cal_slo) 
      
      # print(out$cal_slo_median_SE[1:3, ])
      out$errors[n_sim, ] <- lapply(temp_errors, FUN = length)
      out$warns[n_sim, ]  <- lapply(temp_warns,  FUN = length)
      out$last_par <- parameters
      out$n_median_bootstraps <- n_median_bootstraps
      
      print(ToOneString("Loaded id", id, "with correlation #", cor_id, ", binary vars:", ToBool(bin), ", and effect:", effect, sep = " "))
    }
    
    names(out$max_rank) <- names(out$min_rank)   <- c(par_names, "Lasso", "Ridge")
    
    names(out$Brier) <- names(out$Brier_SE) <- names(out$Brier_pct) <- names(out$Brier_pct_SE) <- c(par_names, names(results$Brier), "DGM")
    
    names(out$Cox) <- names(out$Nagelkerke) <- names(out$Nagelkerke_SE) <- names(out$McFadden) <- names(out$Nagelkerke_pct) <- 
      names(out$Nagelkerke_pct_SE) <- c(par_names, names(R2), "DGM")
    
    names(out$pdi) <- names(out$pdi_SE) <- names(out$pdiw) <- names(out$pdi_spec_small) <- names(out$pdi_spec_small_middle) <-
      names(out$pdi_pct) <- c(par_names, names(PDI)) # [-length(names(PDI))]
    names(out$pdi_pct_SE) <- c(par_names, names(PDI))
    
    names(out$M) <- names(out$M_SE) <- names(out$OneVsRestC) <- c(par_names, names(C_temp) , "DGM") 
    
    cal_names <- names(lColApply(temp_cal_int, FUN = mean, selection = cal_select))
    NA_cal_names <- paste("Comp_", cal_names, sep = "")
    names(out$cal_int_median) <- names(out$cal_slo_median) <- names(out$cal_int_mean) <- names(out$cal_slo_mean) <- 
      c(par_names, cal_names, "success")
    
    names(out$cal_int_mean_SE) <- names(out$cal_slo_mean_SE)   <- c(par_names, cal_names, NA_cal_names)
    dimnames(out$all_cal_slo)[[2]] <- unlist(c(par_names, dimnames(mat_cs)[1]))
    dimnames(out$all_cal_slo)[[2]] <- gsub(":", ".", dimnames(out$all_cal_slo)[[2]])
    dimnames(out$all_cal_slo)[[2]] <- gsub("-", ".", dimnames(out$all_cal_slo)[[2]])
    dimnames(out$all_cal_slo)[[2]] <- gsub(" ", "_", dimnames(out$all_cal_slo)[[2]])
    dimnames(out$all_Brier)[[2]] <- unlist(c(par_names, names(temp_Brier)))
    dimnames(out$all_pdi)[[2]]   <- unlist(c(par_names, names(temp_PDI)))
    
    names(out$cal_int_median_SE) <- names(out$cal_slo_median_SE) <- c(par_names, cal_names, NA_cal_names)
    
    names(out$mse_probs) <- names(out$rmse_probs) <- c(par_names, names(results$mse_probs))
    names(out$mean_coefs_rmse_pop) <- c(par_names, names(sapply(coefs_mse_pop[[1]], mean)))
    names(out$mean_coefs_rmse_DGM) <- c(par_names, names(sapply(coefs_mse_DGM[[1]], mean)))
    names(out$mean_coefs_bias) <- names(out$mean_coefs_last) <- c(par_names, names(unlist(lapply(coefs_bias_list[[n_sim]], rowMeans, na.rm = T))) )
    
    names(out$errors)     <- names(results$error_list)
    names(out$warns)      <- names(results$warn_list)

  return(out)
}

n_median_bs <- 1e5

# ## For loading and aggregating the raw data, and bootstrapping the SE of median and saving: :
begin <- proc.time()

# cor_res  <- LoadData(what = "cor", n.median.bs = n_median_bs)
# bin_res  <- LoadData(what = "bin", n.median.bs = n_median_bs) # This trows 2 warnings that can be ignored.
# eff_res  <- LoadData(what = "effect", n.median.bs = n_median_bs) # This trows 2 warnings that can be ignored.
# main_res <- LoadData(n.median.bs = n_median_bs)
# 
# 
# save(main_res, file = ToOneString(getwd(), "Aggregated results", "main_res.RData", sep = "/") )
# save(cor_res,  file = ToOneString(getwd(), "Aggregated results", "cor_res.RData" , sep = "/") )
# save(bin_res,  file = ToOneString(getwd(), "Aggregated results", "bin_res.RData" , sep = "/") )
# save(eff_res,  file = ToOneString(getwd(), "Aggregated results", "eff_res.RData" , sep = "/") )
# 


# # For loading the aggregated data:
source("Coerce.R")
load(file = ToOneString(getwd(), "Aggregated results", "main_res.RData", sep = "/") )
load(file = ToOneString(getwd(), "Aggregated results", "cor_res.RData" , sep = "/") )
load(file = ToOneString(getwd(), "Aggregated results", "bin_res.RData" , sep = "/") )
load(file = ToOneString(getwd(), "Aggregated results", "eff_res.RData" , sep = "/") )

last_par <- main_res$last_par # Makes plotting a lot easier
print(proc.time() - begin)