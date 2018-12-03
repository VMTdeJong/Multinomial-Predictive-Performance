# #######################################################################################################################

# Authors: 
# The original code for the PDI was provided to us by Hajime Uno
# Final version provided here is edited by Valentijn M.T. de Jong.
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


#######################################
pdi3i <- function(outcome, prob_vec, c1, c2, c3)
{
  
  p1 <- prob_vec[outcome == c1] ; n1 <- length(p1)
  p2 <- prob_vec[outcome == c2] ; n2 <- length(p2)
  p3 <- prob_vec[outcome == c3] ; n3 <- length(p3)
  
  # pdi1
  wk <- 0
  for (i in 1:n1)
  {
    greater2 <- sum(p2 < p1[i])
    greater3 <- sum(p3 < p1[i])
    ties2    <- sum(p2 == p1[i])
    ties3    <- sum(p3 == p1[i])
    
    wk <- wk + as.numeric(greater2) * as.numeric(greater3) + (ties2*greater3)/2 + (greater2* ties3)/2 + (ties2* ties3)/3
  }
  pdii <- wk /(as.numeric(n1) * as.numeric(n2) * as.numeric(n3))
  
  return(c(pdii, n1, n2, n3))
}

########################################
MyPdi3 <- function(outcome, probs, cats = sort(unique(outcome))){
  
  pdi1 <- pdi3i(outcome, probs[ , 1], cats[1], cats[2], cats[3])[1]
  pdi2 <- pdi3i(outcome, probs[ , 2], cats[2], cats[3], cats[1])[1]
  pdi3_list <- pdi3i(outcome, probs[ , 3], cats[3], cats[1], cats[2])
  pdi3 <- pdi3_list[1]
  
  n <- pdi3_list[2:4]
  Z      <- list()
  Z$pdii <- c(pdi1, pdi2, pdi3)
  Z$pdi  <- mean(Z$pdii)
  Z$pdiw <- t(as.matrix(n)) %*% as.matrix(Z$pdii)/sum(as.matrix(n))
  return(Z)
}

########################################
### PDI that returns dataframe
PdiDF <- function(outcome, probs, cats = sort(unique(outcome)), method_name = "unnamed"){
  
  pdi1 <- pdi3i(outcome, probs[ , 1], cats[1], cats[2], cats[3])[1]
  pdi2 <- pdi3i(outcome, probs[ , 2], cats[2], cats[3], cats[1])[1]
  pdi3_list <- pdi3i(outcome, probs[ , 3], cats[3], cats[1], cats[2])
  pdi3 <- pdi3_list[1]
  
  n <- pdi3_list[2:4]
  
  pdiw <- t(as.matrix(n)) %*% as.matrix(c(pdi1, pdi2, pdi3))/sum(as.matrix(n))
  return(data.frame("pdi1" = pdi1, "pdi2" = pdi2, "pdi3" = pdi3, "pdi" = mean(c(pdi1, pdi2, pdi3)), "pdiw" = pdiw, "method_name" = method_name))
}