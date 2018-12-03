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

########################################## My M-index function.  #########################################
# Assumes response is numerical, and predictor columns are sorted in increasing order. 
# e.g. column 1 is outcome 0, column 2 is outcome 1, etc.
# For binomial outcome it breaks down to the AUC or c-statistic.
# Note that probabilities in a column are the (estimated) probabilities that the 
# observation falls into that category. If a binomial outcome is supplied, a n x n_cat matrix of predictions
# must be supplied! 
MyM <- function(preds, outcome, cats = sort(unique(outcome)))
{
  outcome <- as.matrix(outcome)
  preds <- as.matrix(preds)
  
  n_cat <- ncol(preds)
  if (n_cat < 2)  { stop("Prediction matrix should contain at least 2 columns, where each column contains the predicted probabilities for that category.") } 
  if (n_cat != length(cats)) { stop("Number of outcome categories should match number of predictor columns. If predictions are made for an outcome not present in this dataset, the possible outcome categories should be supplied in increasing order to parameter cats.")}

  A_mat <- matrix(NA, nrow = n_cat, ncol = n_cat)
  dimnames(A_mat) <- list("outcome" = cats, "other" = cats)
  for (i in 1:n_cat)
  {
    cat0 <- outcome == cats[i] 
    n0   <- sum(cat0)
    
    for (j in 1:n_cat)
    {
      if (i != j) 
      {
        cat1 <- outcome == cats[j]
        temp_preds <- preds[cat0 | cat1, i] 
        temp_out  <- outcome[cat0 | cat1]
        
        n1   <- sum(cat1)
        
        r <- rank(temp_preds)
        S0 <- sum(as.numeric(r[temp_out == cats[i]]))
        
        A_mat[i, j] <- (S0 - n0 * (n0 + 1)/2)/(as.numeric(n0) * as.numeric(n1))  # = equation 3 from Hand and Till 2001. Except it is not multiplied by 2, as A is not divided by 2 either.
      } 
    }
  }
  return( list("M" = sum(sum(A_mat, na.rm = T)/(n_cat*(n_cat - 1))),  # = left part of equation 7 from Hand and Till 2001
               "pairwise_M" = A_mat ))              # Pairwise c-stats
}

############################### My AUC or C stat implemented as M-index ###################################
# For binomial outcomes. The probabilities are the (estimated) probabilities
# the observation falls in the higher outcome category.

MyC <- function(preds, outcome)
{
  if (is.matrix(preds) || is.data.frame(preds)) { if (ncol(preds) != 1) { stop("Number of outcome columns must be one.")}}
  preds <- as.matrix(preds)
  cats <- sort(unique(outcome))
  n_cat <- length(cats)
  if (n_cat != 2) { stop("Number of outcome categories must be 2.")}
  
    n0   <- sum(outcome == cats[2])
    n1   <- length(outcome) - n0
    
    r <- rank(preds[,1])
    S0 <- sum(as.numeric(r[outcome == cats[2]]))
    
    A <- (S0 - n0 * (n0 + 1)/2)/(as.numeric(n0) * as.numeric(n1))  # = equation 3 from Hand and Till 2001. 
    return(A)  
}


################################### Test it for multinomial outcome #######################################
# library(AUC)
# library(HandTill2001)
# 
# 
# compute_M <- function(response, preds)
#   {
#   
#   M <- multcap(as.factor(response[]), as.matrix(preds[,]));
#   return(auc(M));
#   }

# Sometimes I get errors with the HandTill2001 package, but in the runs it worked, it produced the same results as mine,
# Except mine was 2 to 5 times as quick.


# Sample
# n_vec <- c(2400, 2400, 2400)
# system.time(retro_data_list3 <- MakeRetroData(n = n_vec, outcome = c(0, 1, 2), 
#                                               coefs = rbind(c(0, 0.2, 0.2, 0.5, 0.8), c(0, 0.2, 0.2, 0.5, 0.8)), 
#                                               n_bi_var = 0, bi_prob = c(.5), Brier = T))
# 
# retro_data3       <- data.frame(retro_data_list3$dat)
# sample_x          <- as.matrix(retro_data3[ ,-1])
# sample_y          <- retro_data3[ , 1]
# p   <- retro_data_list3$y_probs
# 
# ############### LASSO
# library(glmnet)
# source("foldid.R")
# 
# 
# lasso_cv          <- cv.glmnet(sample_x, sample_y, alpha = 1, nfolds = 10,family = "multinomial")
# lambda_lasso      <- lasso_cv$lambda.min
# 
# 
# lasso_model          <- glmnet(x = sample_x, y = sample_y, alpha = 1, family = "multinomial")
# 
# sample_lasso_fitted  <- data.frame(predict(lasso_model, newx = sample_x, s = lambda_lasso, type = "response"))

# Ridge
# 
# library(glmnet)
# source("foldid.R")
# 
# 
# ridge_cv          <- cv.glmnet(sample_x, sample_y, alpha = 0, nfolds = 10,family = "multinomial")
# lambda_ridge      <- ridge_cv$lambda.min
# 
# 
# ridge_model          <- glmnet(x = sample_x, y = sample_y, alpha = 0, family = "multinomial")
# ridge_model2         <- ridge_cv$glmnet.fit
# 
# sample_ridge_fitted   <- data.frame(predict(ridge_model,  newx = sample_x, s = lambda_ridge, type = "response"))
# sample_ridge_fitted2  <- data.frame(predict(ridge_model2, newx = sample_x, s = 52, type = "response"))
# 
# MyM(sample_ridge_fitted, sample_y)$M
# MyM(sample_ridge_fitted2, sample_y)$M
### Run the c stats
# max_i <- 100
# my_M   <- rep(NA, max_i)
# M      <- rep(NA, max_i)
# t_my_M <- rep(NA, max_i)
# t_M    <- rep(NA, max_i)
# colnames(sample_lasso_fitted) <-  c("0", "1", "2")
# 
# for (i in 1:max_i)
# {
#   sp <- sample(length(sample_y), size = 100000)
#   t_my_M[i] <- system.time(my_M[i] <- MyM(as.matrix(sample_lasso_fitted[sp,]), sample_y[sp])$M )[1]
#   t_M[i]    <- system.time( M[i]   <- compute_M(sample_y[sp], sample_lasso_fitted[sp,])  )[1]
# }
# head(M)
# head(my_M)
# 
# all.equal(M, my_M)
# mean(t_my_M)
# mean(t_M)/mean(t_my_M)


### One vs rest c-statistic, where outcome is a matrix with columns for each outcome, and 1 denotes presence of that outcome
# and 0 denotes absence of that outcome.
OneVsRestCBi <- function(preds, out_mat)
{
  out_mat      <- as.matrix(out_mat)
  preds        <- as.matrix(preds)
  
  if ( ncol(out_mat) != ncol(preds)) { stop("Number of columns of resonse matrix should be equal to number of columns of prediction matrix")}
  n_cat = ncol(preds)
  
  pairwise_c <- rep(NA, n_cat)
  
  for (i in 1:n_cat)
  {
    pairwise_c[i] <- MyC(preds[,i], out_mat[,i])
  }
  
  return(list("mean_i_vs_rest" = mean(pairwise_c), "i_vs_rest" = pairwise_c))
}

# One vs rest C statistic. Uses exact same formula as OneVsRestCBi, but assumes outcome is a vector of outcome categories.
# Note that it assumes the preds matrix is sorted. e.g. first outcome 0, then 1, then 2... or: first "A", then "B"
# A vector of outcome categories may be supplied to parameter cats, if the sample contains no events of a certain outcome
# category for which predictions are made.
OneVsRestC <- function(preds, outcome, cats = sort(unique(outcome)))
{
  outcome      <- as.matrix(outcome)
  preds        <- as.matrix(preds)
  
  n_cat = ncol(preds)
  if ( length(cats) != n_cat) { warning("Number of columns of prediction matrix should be equal to number of outcome categories.")} # This does not produce an error, because it might be that there are no cases of a certain event, when predictions for a validation set are made.
  
  out_mat <- matrix(0, ncol = n_cat, nrow = nrow(preds)) # convert outcome to matrix
  for (i in 1:n_cat) {out_mat[outcome == cats[i], i] <- 1} 
  
  
  pairwise_c <- rep(NA, n_cat)
  
  for (i in 1:n_cat)
  { if (any(outcome == cats[i])) {
    pairwise_c[i] <- MyC(preds[,i], out_mat[,i])
  }
  }
  
  return(list("mean_i_vs_rest" = sum(pairwise_c, na.rm = T)/n_cat, "i_vs_rest" = pairwise_c))
}
# 
# system.time(one <- OneVsRestC(sample_lasso_fitted, retro_data_list3$pop_disease))
# one


### Test for Binomial Outcome # STill correct?
# n_vec <- c(40000, 40000)
# system.time(retro_data_list3 <- MakeRetroData(n = n_vec, outcome = c(0, 1), 
#                                               coefs = c(0, 0.2, 0.3, 0, 0), 
#                                               n_bi_var = 0, bi_prob = c(.5), Brier = T))
# 
# retro_data3       <- data.frame(retro_data_list3$dat)
# sample_x          <- as.matrix(retro_data3[ ,-1])
# sample_y          <- retro_data3[ , 1]
# 
# 
# # sample and estimate with lasso
# library(glmnet)
# source("foldid.R")
# 
# 
# lasso_cv               <- cv.glmnet(sample_x, sample_y, alpha = 1, nfolds = 10,family = "binomial")
# lambda_lasso           <- lasso_cv$lambda.min
# 
# lasso_model          <- glmnet(x = sample_x, y = sample_y, alpha = 1, family = "binomial")
# 
# sample_lasso_fitted  <- data.frame(predict(lasso_model, newx = sample_x, s = lambda_lasso, type = "response"))
# 
# 
# library(AUC)
# fitted_values <- as.matrix(sample_lasso_fitted, ncol = 1, nrow = length(sample_lasso_fitted))
# 
# max_i  <- 100
# c      <- rep(NA, max_i)
# a      <- rep(NA, max_i)
# t_c    <- rep(NA, max_i)
# t_a    <- rep(NA, max_i)
# for (i in 1:max_i)
# {
#   s <- sample(nrow(fitted_values), size = 100)
#   s_fitted_values <- fitted_values[s,]
#   s_sample_y <- sample_y[s]
#   t_a[i]   <- system.time(a[i] <- AUC::auc(AUC::roc(s_fitted_values, as.factor(s_sample_y))))[1]
#   t_c[i]   <- system.time(c[i] <- MyC(s_fitted_values, s_sample_y))[1] 
# }
#   
# 
# all.equal(a, c)
# mean(t_a)/mean(t_c)


################# Test to see what happens when there are (a lot of) ties #####################
# 
# n_vec <- c(20000, 20000)
# system.time(retro_data_list_ties <- MakeRetroData(n = n_vec, outcome = c(0, 1), 
#                                               coefs = c(0, 0, 0, 0, 0), 
#                                               n_bi_var = 0, bi_prob = c(.5), Brier = T))
# 
# retro_data_ties        <- data.frame(retro_data_list_ties$dat)
# sample_y_ties          <- retro_data_ties[ , 1]
# sample_probs_ties      <- retro_data_list_ties$y_probs
# 
# MyM(sample_probs_ties, sample_y_ties)
# MyC(sample_probs_ties[ ,2], sample_y_ties)
# AUC::auc(AUC::roc(sample_probs_ties[,2], as.factor(sample_y_ties)))
