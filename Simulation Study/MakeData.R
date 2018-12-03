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

#remove(list = ls())

# Simple function for making coefficients. Assumes 3 outcome categories. 
# First outcomes' coefs are negative. Second's are positive.
MakeCoefs <- function(alphas, n_coefs){
  if (length(alphas) == 1) { alphas <- c(alphas, alphas)}
  if (n_coefs == 0) { coefs <- matrix(alphas, nrow = length(alphas), ncol = 1)
  } else {
    coefs <- rep(c(0.2, 0.2, 0.5, 0.8), n_coefs /4)
    coefs <- rbind(-c(alphas[1], coefs), c(-alphas[2], coefs))
  }
  return(coefs)
}

MakeData <- function(n, outcome, coefs, n_bi_var = 0, bi_prob = c(0.5), Brier = FALSE, rand.order = T)
{  
  if ( !is.matrix(coefs)) { coefs <- matrix(coefs, nrow = 1, ncol = length(coefs))}
  n_var = ncol(as.matrix(coefs))
  n_cont_var = n_var - n_bi_var - 1

  ##### Test input.
  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5)  abs(x - round(x)) < tol
  if ( n_bi_var   < 0 || !is.wholenumber(n_bi_var)) { stop("n_bi_var must be whole number of 0 or greater")}
  if ( n          < 2 || !is.wholenumber(n))        { stop("n must be whole number of 2 or greater")}
  if ( !isTRUE(all.equal(dim(coefs), c(length(outcome) - 1, n_bi_var + n_cont_var + 1) ) ) )  { stop("Dimensions of coefs should match with predictors.")}
  
  ##### Sample Predictors.
  xx <- matrix(1, nrow = n, ncol = 1) # Matrix with intercept, to store predictors in.
  
  ### If requested, sample binomial predictors
  if ( n_bi_var > 0 ) 
    { 
    if ( length(bi_prob) != n_bi_var && length(bi_prob) != 1) {stop("Number of probabilities for binomial variables must match number of binomial variables, or be 1, or not given (in which case .5 will be assumed for all binomial variables).")}  
    binary_vars <- matrix(rbinom(n = n * n_bi_var, size = 1, prob = bi_prob), ncol = n_bi_var, nrow = n, byrow = T) - 1/2
    xx <- cbind(xx, binary_vars)
  }
  
  ### If requested, sample normal predictors
  if ( n_cont_var > 0 ) {
    xx <- cbind(xx, matrix(rnorm( n = n * n_cont_var, mean = 0, sd = 1), ncol = n_cont_var, nrow = n))
  }
  
  ##### Sample y.
  ### Start with computing the linear predictors
  n_cat <- length(outcome)
  y_lin_pred <- matrix(0, nrow = n, ncol = n_cat) # exp(0) = 1, but is actually not used (for computational speed)
  
  for (cat in 2:n_cat) {
    y_lin_pred[,cat] <- xx %*% coefs[cat - 1,] 
  }
  print(y_lin_pred)
  ### Then compute probability of each outcome
  y_probs <- matrix(0, nrow = n, ncol = cat)
  
  if (n_cat > 2)
  {
    for (cat in 2:n_cat) {y_probs[ , cat] <- exp(y_lin_pred[,cat]) / (1 + rowSums(exp(y_lin_pred[ , -1])) ) }
  } else {y_probs[ , 2] <- exp(y_lin_pred[,2]) / (1 + exp(y_lin_pred[ , 2]) ) }
  
  # And finally the reference category must be equal to 1 - P(other categories)
  y_probs[, 1] <- 1 - rowSums(y_probs) 
  
  # Sample y, with computed probabilities
  y <- rep(0, n)
  for (i in 1:n) { y[i] <- sample(outcome, size = 1, replace = FALSE, prob = y_probs[i, ]) }
  
  dat <- cbind(y, xx[ , -1]) # return y and x's, without the intercept
  colnames(dat) <- NULL
  colnames(dat) <- colnames(dat, do.NULL = FALSE, prefix = "X")
  
  if (Brier)
  {return(list(dat = dat, y_probs = y_probs)) }
  else {return(list(dat = dat)) }
}

## Test the function :
# system.time(data_list <- MakeData(n = 30, outcome = c(0, 1, 2), coefs = MakeCoefs(0, n_coefs = 4),
#                             n_bi_var = 0, Brier = T ))

# ## And evaluate:
# require(mlogit)
# test_data       <- data.frame(data_list$dat)
# test_data_long  <- mlogit.data(test_data, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
# mlogit_model    <- mlogit(formula = CreateMlogitFormulaInt(ncol(test_data)), data = test_data_long) 
# mlogit_model
# # # Or :
# mlogit_model2   <- do.call("mlogit", args = list(formula = CreateMlogitFormulaInt(ncol(test_data)), data = test_data_long))
# mlogit_model2$coefficients

#### n = c() of outcome category sizes,
MakeRetroData <- function(n, outcome, coefs, n_bi_var = 0, bi_prob = c(.5), Brier = T)
{
  n_cat   <- length(n)
  more    <- sum(n)
  k       <- 1
  if ( !is.matrix(coefs)) { coefs <- matrix(coefs, nrow = 1, ncol = length(coefs)) }
  ncol    <- ncol(as.matrix(coefs))
  dat     <- matrix(nrow = 0, ncol = ncol)
  y_probs <- matrix(nrow = 0, ncol = n_cat)
  
  while (more > 0)
  { 
    temp                <- MakeData(n = more*k, outcome = outcome, coefs = coefs, n_bi_var = n_bi_var, bi_prob = bi_prob, Brier = Brier)
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
    # For each category, select up to n[cat] observations
    sel_data                <- rbind(sel_data,        dat[dat[ , 1] == outcome[cat], ][1:(n[cat]), ])
    if (Brier) {sel_y_probs <- rbind(sel_y_probs, y_probs[dat[ , 1] == outcome[cat], ][1:(n[cat]), ]) }
  }
  colnames(sel_data) <- NULL
  colnames(sel_data) <- colnames(sel_data, do.NULL = FALSE, prefix = "X")
  
  if (Brier)  { 
    sel_disease <- matrix(0, nrow = sum(n), ncol = n_cat)
    for (cat in 1:n_cat) { sel_disease[sel_data[ , 1] == outcome[cat], cat] <- 1} # convert y
    pop_Brier <- colMeans((sel_disease - sel_y_probs) ^ 2) 
    return( list(dat = sel_data, pop_Brier = pop_Brier, pop_disease = sel_disease, y_probs = sel_y_probs))}
  else {return(list(dat = sel_data) )}

}

# n_vec <- c(100, 100, 100)
# system.time(retro_data_Brier_disease <- MakeRetroData(n = n_vec, outcome = c(0, 1, 2), 
#                                         coefs = rbind(c(0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25, 0.25, 0.8, 0.8, 0.8, 0.8), -c(0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25, 0.25, 0.8, 0.8, 0.8, 0.8)), 
#                                         n_bi_var = 0, bi_prob = c(.5), Brier = T) )
# retro_data_Brier_disease$pop_Brier
# colMeans(retro_data_Brier_disease$y_probs)
# retro_data <- data.frame(retro_data_Brier_disease$dat)

# 
# system.time(retro_data_Brier_disease2 <- MakeRetroData(n = c(400000, 400000, 400000), outcome = c(0, 1, 2), 
#                                                       coefs = rbind(c(0, 0.4, 0.4), c(0, 0.4, 0.4)), 
#                                                       n_bi_var = 0, bi_prob = c(.5), Brier = T) )
# retro_data_Brier_disease2$pop_Brier
# 
# retro_data2 <- retro_data_Brier_disease2$dat
# 
# #retro_data[,1] <- retro_data[, 1] + 2
# require(mlogit)
# table(retro_data$X1)
# retro_data_set_long <- mlogit.data(retro_data, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
# retro_mlogit_model  <- mlogit(formula = CreateMlogitFormulaInt(ncol(retro_data)), data = retro_data_set_long) 
# retro_mlogit_fitted <- retro_mlogit_model$probabilities
# retro_mlogit_model$coefficients
# retro_mlogit_model$logLik[1]
# 
# # retro_data2 <- data.frame(retro_data2)
# table(retro_data2$X1)
# retro_data_set_long2 <- mlogit.data(retro_data2, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
# retro_mlogit_model2  <- mlogit(formula = CreateMlogitFormulaInt(ncol(retro_data2)), data = retro_data_set_long2) 
# retro_mlogit_model2


#### Draw Sample
#dat <- MakeRetroData(n = c(10,20,40), n_bi_var = 2, n_cont_var = 2, outcome = c(0,1,2), coefs = rbind(c(0,0,0,0,0), c(0,0,0,0,0)))

DrawSample <- function(dat, n, N, outcome = sort(unique(dat[,1])))
{
  dat <- as.matrix(dat)
  n_cat <- length(n)
  draw <- matrix(nrow = sum(n), ncol = ncol(dat))
  v <- 0
  for (c in 1:n_cat)
  {
    draw[(1 + v):(n[c] + v), ] <- (dat[dat[,1] == outcome[c], ])[sample(x = N[c], size = n[c], replace = TRUE), ]
    v <- v + n[c]
  }
  return(draw)
}
# 
# system.time(dr <- data.frame(DrawSample(retro_data, n = c(800, 100, 100), N = n_vec  )     ) )
# dr_long <- mlogit.data(dr, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
# dr_mlogit_model <- mlogit(formula = mFormula(X1 ~ 1 | X2 + X3 + X4 + X5), data = dr_long)
# dr_mlogit_model
# 
# s <- matrix(nrow = 5, ncol = length(.Random.seed))
# s[1,] <- .Random.seed
# dr <- data.frame(DrawSample(retro_data, n = c(800, 100, 100), N = n_vec  )     )
# head(dr)
# 
# .Random.seed <- s[1,]
# dr2 <- data.frame(DrawSample(retro_data, n = c(800, 100, 100), N = n_vec  )     )
# head(dr2)
# 
# TestSample <- function(seed)
# {
#   .Random.seed <<- seed
#   return(data.frame(DrawSample(retro_data, n = c(800, 100, 100), N = n_vec  )     ) )
# }
# 
# .Random.seed <- s[1,]
# head(TestSample(s[1,]))
# 
# .Random.seed <- s[1,]
# head(TestSample(s[1,]))
