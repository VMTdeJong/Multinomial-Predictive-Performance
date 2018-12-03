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

# Function to create variance-covariance matrix, for given correlations
CreateSigma <- function(c = 0, n_pred)
{
  if (length(c) == 1) { c <- rep(c, n_pred ^ 2) } else {
  if (length(c) != factorial(n_pred - 1)) stop("length of c should be 1 or equal factorial(n_pred - 1)")
  }
  
  s <- matrix(nrow = n_pred, ncol = n_pred)
  k <- 1
  
  for (row in 1:n_pred) {
    for (col in 1:row)
    {
      if (row == col) s[row, col] <- 1 else {
        s[row, col] <- c[k]
        s[col, row] <- c[k]
        k <- k + 1
      }
    }
  }
  return(s)
}

# c <- 0.2
# n_pred <- 4
# 
# m <- CreateSigma(c, n_pred)
# m


MakeCorData <- function(n, outcome, coefs, n_bi_var = 0, bi_prob = c(0.5), Brier = FALSE, c = 0)
{  
  if ( !is.matrix(coefs)) { coefs <- matrix(coefs, nrow = 1, ncol = length(coefs))}
  n_var = ncol(as.matrix(coefs))
  n_cont_var = n_var - n_bi_var - 1
  n <- max(n, n_var) # Necessary because mvrnorm() generates an error when n <= number of predictors
  
  ##### Test input.
  n_bi_var = 0
  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5)  abs(x - round(x)) < tol
  if ( n_bi_var   < 0 || !is.wholenumber(n_bi_var)) { stop("n_bi_var must be whole number of 0 or greater")}
  if ( n          < 2 || !is.wholenumber(n))        { stop("n must be whole number of 2 or greater")}
  if ( !isTRUE(all.equal(dim(coefs), c(length(outcome) - 1, n_bi_var + n_cont_var + 1) ) ) )  { stop("Dimensions of coefs should match with predictors.")}
  library(MASS)
  
  
  ##### Sample Predictors.
  xx <- matrix(1, nrow = n, ncol = 1) # Matrix with intercept, to store predictors in.
  
  ### If requested, sample binomial predictors
  if ( n_bi_var > 0 ) 
  { 
    # Commented out, as it does not produce correlated predictors.
    
    #if ( length(bi_prob) != n_bi_var && length(bi_prob) != 1) {
    #  stop("Number of probabilities for binomial variables must match number of binomial variables, or be 1, or not given (in which case .5 will be assumed for all binomial variables).")}  
    #
    #xx <- cbind(xx, matrix(rbinom(n = n * n_bi_var, size = 1, prob = bi_prob), ncol = n_bi_var, nrow = n, byrow = T))
  }
  
  ### If requested, sample correlated normal predictors
  if ( n_cont_var > 0 ) {
    sigma <- CreateSigma(c = c, n_pred = n_cont_var)
    preds <- mvrnorm(n, mu = rep(0, n_cont_var), Sigma = sigma, empirical = TRUE) 
    xx <- cbind(xx, preds)
  }
  
  ##### Sample y.
  ### Start with computing the linear predictors
  n_cat <- length(outcome)
  y_lin_pred <- matrix(0, nrow = n, ncol = n_cat) # exp(0) = 1, but is actually not used (for computational speed)
  
  for (cat in 2:n_cat) {
    y_lin_pred[,cat] <- xx %*% coefs[cat - 1,] 
  }
  
  ### Then compute probability of each outcome
  y_probs <- matrix(0, nrow = n, ncol = cat)
  
  if (n_cat > 2)
  {
    for (cat in 2:n_cat)
    {
      y_probs[ , cat] <- exp(y_lin_pred[,cat]) / (1 + rowSums(exp(y_lin_pred[ , -1])) )
    }
  } else {
    y_probs[ , 2] <- exp(y_lin_pred[,2]) / (1 + exp(y_lin_pred[ , 2]) )
  }
  
  # And finally the reference category must be equal to 1 - P(other categories)
  y_probs[, 1] <- 1 - rowSums(y_probs) 
  
  # Sample y, with computed probabilities
  y <- rep(0, n)
  for (i in 1:n)
  {
    y[i] <- sample(outcome, size = 1, replace = FALSE, prob = y_probs[i, ])
  }
  
  dat <- cbind(y, xx[ , -1]) # return y and x's, without the intercept
  colnames(dat) <- NULL
  colnames(dat) <- colnames(dat, do.NULL = FALSE, prefix = "X")
  
  if (Brier)
  {return(list(dat = dat, y_probs = y_probs)) }
  else {return(list(dat = dat)) }
}


# system.time(data_list2 <- MakeCorData(n = 30000 , outcome = c(0, 1, 2), coefs = rbind(c(0, 1, 1, 1, 1), c(0, 0, 0, 0 , 0) ),   
#                                    Brier = T , c = 0.5))
# 
# data2 <- data_list2$dat
# 
# summary(data2)
# 
# cor(data2[,-1])


MakeRetroCorData <- function(n, outcome, coefs, Brier = T, c = 0)
{
  n_cat   <- length(n)
  more    <- sum(n)
  k       <- 1
  if ( !is.matrix(coefs)) { coefs <- matrix(coefs, nrow = 1, ncol = length(coefs)) }
  ncol    <- ncol(as.matrix(coefs))
  if (nrow(coefs) != length(n) - 1) stop("Wrong number of rows for coefs, for this number of outcome categories: nrow(coefs) == length(n) - 1")
  dat     <- matrix(nrow = 0, ncol = ncol)
  y_probs <- matrix(nrow = 0, ncol = n_cat)
  
  while (more > 0)
  { 
    temp                <- MakeCorData(n = more*k, outcome = outcome, coefs = coefs, Brier = Brier, c = c)
    dat                 <- rbind(dat, temp$dat)
    if (Brier) {y_probs <- rbind(y_probs, temp$y_probs)}
    k                   <- k*2
    
    more <- 0
    for (cat in 1:n_cat)
    {
      this_cat <- max(n[cat] - sum(dat[,1] == outcome[cat]), 0)
      if (this_cat > more) { more <- this_cat }
    }
  }
  
  sel_data    <- matrix(nrow = 0, ncol = ncol)
  sel_y_probs <- matrix(nrow = 0, ncol = n_cat)
  
  for (cat in 1:n_cat)
  {
    # For each category, select up to n[cat] patients.
    sel_data                <- rbind(sel_data,        dat[dat[ , 1] == outcome[cat], ][1:(n[cat]), ])
    if (Brier) {sel_y_probs <- rbind(sel_y_probs, y_probs[dat[ , 1] == outcome[cat], ][1:(n[cat]), ]) }
  }
  colnames(sel_data) <- NULL
  colnames(sel_data) <- colnames(sel_data, do.NULL = FALSE, prefix = "X")
  
  if (Brier)  { 
    sel_disease <- matrix(0, nrow = sum(n), ncol = n_cat)
    for (cat in 1:n_cat) { sel_disease[sel_data[ , 1] == outcome[cat], cat] <- 1} # convert population y
    pop_Brier <- colMeans((sel_disease - sel_y_probs) ^ 2) 
    return( list(dat = sel_data, pop_Brier = pop_Brier, pop_disease = sel_disease, y_probs = sel_y_probs))}
  else {return(list(dat = sel_data) )}
  
}

# coefs <- MakeCoefs(alphas = 0.3308024, n_coefs = 4) 
# system.time(data_list3 <- MakeRetroCorData(n = c(1e4, 1e4, 1e4) , outcome = c(0, 1, 2), coefs = coefs,   
#                                    Brier = T , c = 0.5))
# 
# data3 <- data_list3$dat
# 
# summary(data3)
# 
# cor(data3[,-1])

