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

# Computes Log Likelihood for Multinomial Distribution
# Labels has to be a vector of multinomial outcomes.
# Probs has to be matrix or data.frame, where each column contains the probabilies of an outcome.
# Note that the order of the columns has to be in increasing order. E.g. 0, 1, 2,...,10, or "A", "B", "C"
mult_log_L <- function(probs, labels)
{
  L <- 0
  cats <- sort(unique(labels))
  if (length(cats) != ncol(probs)) {stop("Number of categories in labels shoud be equal to number of columns of probs.")}
  if (nrow(probs) != length(labels)) {stop("Number of rows of probs should be equal to length of labels.")}
  for (i in 1:ncol(probs))
  {
    L <- L + sum(log(probs[labels == cats[i], i]))
  }
  
  return(L)
}

# Computes Cox, Nagelkerke and McFaddens R^2.
# Labels has to be a vector of multinomial outcomes.
# Probs has to be matrix or data.frame, where each column contains the probabilities of an outcome.
# Note that the order of the columns has to be in increasing order. E.g. 0, 1, 2,...,10, or "A", "B", "C"
# L0, n_obs, cats and n_cat may be supplied for computational speed, but can also be
# computed using preds and labels alone.
mult_R2 <- function(preds, labels, L0, n_obs, cats, n_cat)
{
  if (missing(n_obs)) {n_obs <- length(labels)}
  if (nrow(as.matrix(preds)) != n_obs) { stop("Number of rows of each preds matrix should be equal to length of labels.")}
  if (missing(cats)) cats <- sort(unique(labels))
  if (missing(n_cat)) n_cat <- length(cats)
  if (missing(L0))
    {
    prevs <- rep(NA, n_cat) # Prevalence
    for (i in 1:n_cat) { prevs[i] <- sum(labels == cats[i])/n_obs } # Prevalence
    L0 <- mult_log_L(matrix(prevs, nrow = n_obs, ncol = n_cat, byrow = T), labels) # Null model Likelihood
  }
  
  # The likelihoods of the model
  L <- mult_log_L(preds, labels)
  
  # And finally the R2:
  Cox <- 1 - exp(-(L - L0) * 2 / n_obs)
  Cox_max <- 1 - exp(2 * n_obs ^ (-1) * L0)
  Nagelkerke <- Cox/Cox_max
  McFadden <- 1 - L / L0
  
  return(data.frame("Cox" = Cox, "Nagelkerke" = Nagelkerke, "McFadden" = McFadden))
}
# 
# set.seed(201601)
# rdbd <- MakeRetroData(n = c(400, 100, 100), outcome = c(0, 1, 2), 
#                                                       coefs = rbind(c(0, 4, 4, 4, 4), c(0, 4, 4, 4, 4)), 
#                                                       n_bi_var = 0, bi_prob = c(.5), Brier = T)
# 
# rd <- data.frame(rdbd$dat)
# 
# require(mlogit)
# rd_long <- mlogit.data(rd, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
# rmm  <- mlogit(formula = CreateMlogitFormulaInt(ncol(rd)), data = rd_long) 
# rmp <- rmm$probabilities
# rmm$logLik[1]
# mult_log_L(rmp, rd[,1])
# system.time(my_R2 <- mult_R2(rmp, rd[,1]))
# my_R2
# write.table(rd, file = "R2test.csv")

