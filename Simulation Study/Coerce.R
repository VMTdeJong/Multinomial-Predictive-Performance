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

# Internal function
.OneStr <- function(l, sep) {
  if (length(l) > 0) out <- l[[1]] else return("")
  
  if (length(l) > 1) {
    for (i in 2:length(l)) out <- paste(out, l[[i]], sep = sep)
  }
  return(out)
}

# Joins strings. Returns 1 string
# If argumenets are not string, throws error.
JoinStr <- function(..., sep = "") {
  l <- list(...)
  if (!all(sapply(l, is.character)) ) stop("Arguments must be character / string.") 
  
  return(.OneStr(l, sep))
}

# Coerces arguments to character, and returns strings(s). 
# Due to vectorization in R, returns multiple strings when one argument is a vector.
ToStrings <- function(..., sep = "") {
  l <- list(...)
  l <- lapply(l, as.character)
  
  return(.OneStr(l, sep))
}

# Forces all arguments to become part of the same string.
ToOneString <- function(..., sep = "") {
  l <- list(...)
  l <- as.list(unlist(l))
  l <- lapply(l, as.character)
  
  return(.OneStr(l, sep))
}

# For printing anything as a bool:
ToBool <- function(x) if (x) TRUE else FALSE

AreTrue <- function(vec) {
  out <- rep(F, length(vec))
  for (i in 1:length(vec)) if (isTRUE(vec[i])) out[i] <- T
  return(out)
}

# Combine individual strings of the objects within the arguments.
# Arguments must be matrix or data.frame of same size.
CombIndStrings <- function(..., sep = "", append.by = "", append.to = "")
{
  l <- list(...)
  out <- l[[1]]
  
  for (c in 1:ncol(l[[1]]))
  {
    for (r in 1:nrow(l[[1]]))
    {
      out[r, c] <- ToOneString(append.to, ToOneString(unlist(lapply(l, `[`, r, c)), sep = sep), append.by, sep = "")
    }
  }
  out
}