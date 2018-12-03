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

# Given a vector of numbers of observation of each outcome category, and number of folds for cross-validation,
# this produces a foldid vector for the glmnet package.
FoldidLeaveSetOut <- function(n_vec, n_folds)
{
  if (any(n_folds > n_vec)) { warning("Number of folds exceeds number of events. Leave-Set-Out not possible. Some fold(s) received 0 events.")}
  if (any(n_vec < 1)) {stop("The number of observations of an outcome category cannot be smaller than 1.")}
  foldid <- rep(NA, sum(n_vec))
  start <- 1
  for (i in 1:length(n_vec))
  {
    end <- start + n_vec[i] - 1
    foldid[start:end] <- sample(rep(seq(n_folds), length = n_vec[i]))
    
    start <- end + 1
  }
  
  
  return(foldid)
}

# 
# #### Testing
# outcome <- c(rep(1, 2), rep(2, 3), rep(3, 3)) # outcome vector as in my simulation 
# # sum(outcome == 1) counts how many observations there are with outcome = 1.
# cbind(outcome, foldid = FoldidLeaveSetOut(c(sum(outcome == 1), sum(outcome == 2), sum(outcome == 3)), 2))
# cbind(outcome, foldid = FoldidLeaveSetOut(c(sum(outcome == 1), sum(outcome == 2), sum(outcome == 3)), 3))