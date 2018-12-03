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

#### Prediction function using coefs of GLMnet and sample_x

GlmnetPredictList <- function(coefs_list, newdata, add.intercept = T)
{
  if (add.intercept) { newdata = cbind(1, newdata)}
  
  n_cat         <- length(coefs_list)
  y_lin_pred    <- matrix(NA, nrow = nrow(newdata), ncol = n_cat)
  y_probs       <- matrix(0,  nrow = nrow(newdata), ncol = n_cat)
  
  for (cat in 1:n_cat)
  {
    y_lin_pred[ , cat] <-  as.matrix(newdata %*% coefs_list[[cat]]  )
  }
  
  for (cat in 1:n_cat)
  {
    y_probs[ , cat] <- exp(y_lin_pred[,cat]) / (rowSums(exp(y_lin_pred[ , ])) )
  }
  return((y_probs))
}

GlmnetPredictMatCoefAsCol <- function(coefs_mat, newdata, add.intercept = T)
{
  if (add.intercept) { newdata = cbind(1, newdata)}
  n_cat         <- nrow(coefs_mat)
  y_lin_pred    <- matrix(NA, nrow = nrow(newdata), ncol = n_cat)
  y_probs       <- matrix(0,  nrow = nrow(newdata), ncol = n_cat)
  
  for (cat in 1:n_cat)
  {
    y_lin_pred[ , cat] <-  newdata %*% coefs_mat[cat, ]
  }
  
  for (cat in 1:n_cat)
  {
    y_probs[ , cat] <- exp(y_lin_pred[,cat]) / (rowSums(exp(y_lin_pred)) )
  }
  return((y_probs))
}

GlmnetAsRefCat <- function(glmnet_coefs_mat)
{
  return(glmnet_coefs_mat - matrix(glmnet_coefs_mat[1, ], nrow = nrow(glmnet_coefs_mat), ncol = ncol(glmnet_coefs_mat), byrow = T))
}

# 
# # Sample
# n_vec <- c(200, 200, 200)
# coefs <- MakeCoefs(0.330051, 4)
# system.time(retro_data_list3 <- MakeRetroData(n = n_vec, outcome = c(0, 1, 2), 
#                                         coefs = coefs, 
#                                         n_bi_var = 0, bi_prob = c(.5), Brier = T))
# 
# retro_data3       <- data.frame(retro_data_list3$dat)
# sample_x          <- as.matrix(retro_data3[ ,-1])
# sample_y          <- retro_data3[ , 1]
# 
# 
# ############### LASSO
# library(glmnet)
# source("foldid.R")
# source("MakeData.R")
# 
# ncv <- 1
# lambdas_lasso <- rep(NA, ncv)
# for (cv in 1:ncv)
# {
# lasso_cv               <- cv.glmnet(sample_x, sample_y, alpha = 1, nfolds = 10,family = "multinomial")
# lambdas_lasso[cv]      <- lasso_cv$lambda.min
# }
# lambda_lasso <- mean(lambdas_lasso)
# 
# lasso_model          <- glmnet(x = sample_x, y = sample_y, alpha = 1, family = "multinomial")
# 
# sample_lasso_fitted  <- data.frame(predict(lasso_model, newx = sample_x, s = lambda_lasso, type = "response"))
# mean(rowSums((sample_lasso_fitted - retro_data_list3$pop_disease)^2))
# 
# my_lasso_fitted_list      <- GlmnetPredictList(coef(lasso_model, s = lambda_lasso), newdata = sample_x)

# mean(rowSums((sample_lasso_fitted - my_lasso_fitted_list) ^ 2))
# 
# lasso_coef <- t(do.call(cbind, lapply(coef(lasso_model, s = lambda_lasso), matrix)))
# lasso_coef_w_ref_cat <- GlmnetAsRefCat(lasso_coef)
# 
# test_lasso_fitted <- GlmnetPredictMatCoefAsCol(lasso_coef_w_ref_cat, newdata = sample_x)
# lasso_coef_w_ref_cat[2:nrow(lasso_coef_w_ref_cat), ]
# max((test_lasso_fitted - sample_lasso_fitted) ^ 2)
# .Machine$double.eps
# 
# coef_names            <- CreateNames(ncol(as.matrix(coefs)))
# coefs_Lasso_array   <- array(data = NA, dim = c(nrow(coefs), ncol(coefs), 3), dimnames = list(outcome = c(0,1,2)[-1], coef = coef_names, iteration = 1:3)) # Note that glmnet doesnt use a reference category!
# coefs_Lasso_array[ , , 2]  <- GlmnetAsRefCat(t(do.call(cbind, lapply(coef(lasso_model, s = lambda_lasso), matrix))))[2:nrow(lasso_coef), ]
# coefs_Lasso_array[ ,-1, 2] - coefs[, -1]
############# RIDGE
# 
# ridge_cv             <- cv.glmnet(sample_x, sample_y, alpha = 0, nfolds = 10, family = "multinomial", type.measure = "mse")                              # alpha = 0, for ridge
# lambda_ridge         <- ridge_cv$lambda.min
# 
# ridge_model          <- glmnet(x = sample_x, y = sample_y, alpha = 0, family = "multinomial")                             # alpha = 0, for ridge
# 
# sample_ridge_fitted  <- data.frame(predict(ridge_model, newx = sample_x, s = lambda_ridge, type = "response"))
# 
# my_ridge_fitted      <- GlmnetPredictList(coef(ridge_model, s = lambda_ridge), newdata = sample_x)
# 
# mean((sample_ridge_fitted - my_ridge_fitted) ^ 2)

############# Ridge test of foldid
# Validation population:
# retro_data_list4 <- MakeRetroData(n = c(10, 10, 10), outcome = c(0, 1, 2), 
#                                   coefs = rbind(c(0, 0.2, 0.3, 0, 0), c(0, 0.2, 0.3, 0, 0)), 
#                                   n_bi_var = 0, bi_prob = c(.5), Brier = T)
# 
# retro_data4       <- data.frame(retro_data_list4$dat)
# sample_x4          <- as.matrix(retro_data4[ ,-1])
# sample_y4          <- retro_data4[ , 1]
# ######### Default foldid:
# 
# ridge_cv             <- cv.glmnet(sample_x, sample_y, alpha = 0, nfolds = 10, family = "multinomial", type.measure = "mse")                              # alpha = 0, for ridge
# lambda_ridge         <- ridge_cv$lambda.min
# 
# ridge_model          <- glmnet(x = sample_x, y = sample_y, alpha = 0, family = "multinomial")                             # alpha = 0, for ridge
# ridge_coefs          <- GlmnetAsRefCat(t(do.call(cbind, lapply(coef(ridge_model, s = lambda_ridge), matrix))))
# 
# val_ridge_preds      <- GlmnetPredictMatCoefAsCol(ridge_coefs, newdata = sample_x4)
# 
# 
# mean((val_ridge_preds - retro_data_list4$pop_disease) ^ 2) # Brier score
# 
# 
# ######## My foldid function:
# 
# ridge_cv_my_foldid             <- cv.glmnet(sample_x, sample_y, alpha = 0, foldid = FoldidLeaveSetOut(n_vec, 10), family = "multinomial", type.measure = "mse")                              # alpha = 0, for ridge
# lambda_ridge_my_foldid         <- ridge_cv_my_foldid$lambda.min
# 
# ridge_model_my_foldid          <- glmnet(x = sample_x, y = sample_y, alpha = 0, family = "multinomial")                             # alpha = 0, for ridge
# ridge_coefs_my_foldid          <- GlmnetAsRefCat(t(do.call(cbind, lapply(coef(ridge_model_my_foldid, s = lambda_ridge_my_foldid), matrix))))
# 
# val_ridge_preds_my_foldid     <- GlmnetPredictMatCoefAsCol(ridge_coefs_my_foldid, newdata = sample_x4)

# 
# mean((val_ridge_preds_my_foldid - retro_data_list4$pop_disease) ^ 2) # Brier Score
