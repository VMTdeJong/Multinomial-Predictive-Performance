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

ColSDs <- function(mat, na.rm = T){
  k <- ncol(mat)
  out <- rep(NA, k)
  for (col in 1:k)
  {
    out[col] <- sd(mat[,col], na.rm = na.rm)
  }
  names(out) <- colnames(mat)
  return(out)
}

SEM <- function(vec, na.rm = T) {
  sd(vec, na.rm = na.rm) / sqrt(CountSuccess(vec))
}

#ColSDs(temp_Brier[[1]])

ColSEMs <- function(mat, na.rm = T){ 
  out <- ColSDs(mat, na.rm = na.rm)
  return( out / sqrt(apply(mat, 2, CountSuccess)) )
}

#ColSEMs(temp_Brier[[1]])
# 
# lapply(temp_Brier, FUN = ColSDs, na.rm = T)
# lapply(temp_Brier, FUN = ColSEMs, na.rm = T)


ColSEsMedian <- function(mat, na.rm = T, n = 1e3){
  k <- ncol(mat)
  out <- rep(NA, k)
  for (col in 1:k)
  {
    out[col] <- SEMedian(mat[,col], na.rm = T, n = n)
  }
  return(out)
}

SEMedian <- function(vec, na.rm = T, n = 1e3){
  medians <- rep(NA, n)
  
  for (i in 1:n)
  {
    bs <- sample(vec, length(vec), replace = T)
    medians[i] <- median(bs, na.rm = na.rm)
  }
  if (bs_test == 0) { bs_test <<- n ; print(paste(n, " bootstraps", sep = "") ) }
  return(sd(medians))
} 

SE <- function(vec, na.rm = T, n = 1e3, FUN = mean) 
{
  values <- rep(NA, n)
  
  for (i in 1:n)
  {
    bs <- sample(vec, length(vec), replace = T)
    values[i] <- FUN(bs, na.rm = na.rm)
  }
  return(sd(values))
}

ArrayDimMeansSE <- function(ar, MARGIN, n, na.rm = T) {
  sd(apply(ar, FUN = mean, MARGIN = MARGIN, na.rm = na.rm), na.rm = na.rm)/sqrt(n)
}

# Function for applying a function over columns of an array within a list


CountSuccess <- function(x, ...) sum(!is.na(x))

AsVector <- function(mat, byrow = F, ...) {
  if (!byrow) return(as.vector(mat))
  out <- rep(NA, nrow(mat) * ncol(mat))
  s <- 1
  e <- ncol(mat)
  for (row in 1:nrow(mat)) {
    out[s:e] <- mat[row, ] 
    s <- e + 1
    e <- e + ncol(mat) 
  }
  
  return(out)
}

CalNames <- function(mat, outcomes) {  
  cal_names <- c()
  for (i in 1:length(outcomes)) cal_names <- c(cal_names, paste("ref:", paste(outcomes[i], outcomes[-i], sep = " other:"), sep = "") )
  return(cal_names) 
}

AddCalNames <- function(mat, outcomes, row = T) {
  if (row) {row.names(mat) <- CalNames(mat, outcomes)} else {
    names(mat) <- CalNames(mat, outcomes)
  }
  return(mat)
}

# for applying mean or median to calibration data.
lColApply <- function(l.ar, FUN, selection, ..., na.rm = T) {
  if (missing(selection)) { selection <- 1:dim(l.ar[[1]]) } 
  l_ar_sel <- lapply(l.ar, function(x) x[ , selection, ]) # selects columns
  l_df_sel <- lapply(X = l_ar_sel, FUN = function(X) apply(X, MARGIN = c(1, 2), FUN = FUN, na.rm = na.rm, ...)) # applies function in each array
  outcomes <- as.numeric(row.names(l_df_sel[[1]]))

  slope_names <- lapply(l_df_sel, CalNames, outcomes = outcomes)
  out <- unlist(lapply(l_df_sel, AsVector, byrow = T))
  meth_names <- rep(gsub("oos_", "", names(l_df_sel)), each = length(l_df_sel[[1]]))
  
  names(out) <- paste(meth_names, unlist(slope_names), sep = "-")
  
  return(out)
}

CalApply <- function(cal.list, FUN, selection = c(1,4), which = matrix(c(T, F, T, F, F, T), nrow = 3), na.rm = T)
{
  c2 <- lapply(cal.list, FUN = function(x) x[ , selection, ] )
  c3 <- lapply(c2, function(x) apply(x, FUN = FUN, MARGIN = c(1,2), na.rm = na.rm))
  return(lapply(c3, function(x) x[which]) )
}