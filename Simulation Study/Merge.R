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

# Functions for merging (reducing) data structures.
# Name indicates structure. e.g. LiLiAr = List of Lists of Arrays.
MergeLiLiDf <- function(lilidf){
  library(plyr)
  reduced_list <- Reduce(function(...) c(...), lilidf)
  lnames <- names(reduced_list)
  ulnames <- unique(lnames)
  relisted <- plyr::llply(ulnames, function(x) reduced_list[lnames == x])
  
  RedRbind <- function(li){
    data.frame(Reduce(function(...) rbind(...), li))
  }
  lidf <- lapply(relisted, FUN = RedRbind)
  names(lidf) <- ulnames
  lidf
}

MergeLiDf <- function(lidf, combine.row.names = F){
  rn1 <- names(lidf)
  df <- Reduce(function(...) rbind(...), lidf)
  
  if (combine.row.names) {
    rn1 <- rep(gsub("oos_", "", rn1), each = nrow(lidf[[1]]))
    row.names(df) <- paste(rn1, row.names(df), sep = "-")
  }
    
  return(df)
}

MergeLiLiAr <- function(liliar){
  library(plyr)
  library(abind)
  reduced_list <- Reduce(function(...) c(...), liliar)
  lnames <- names(reduced_list)
  ulnames <- unique(lnames)
  relisted <- plyr::llply(ulnames, function(x) reduced_list[lnames == x])
  
  RedAbind <- function(li){
    Reduce(function(...) abind(..., along = 3), li)
  }
  liar <- lapply(relisted, FUN = RedAbind)
  names(liar) <- ulnames
  liar
}

MergeLiAr <- function(liar){
  library(abind)
  Reduce(function(...) abind(..., along = 3), liar)
}

MergeLiLiLi <- function(lilidf){
  library(plyr)
  reduced_list <- Reduce(function(...) c(...), lilidf)
  lnames <- names(reduced_list)
  ulnames <- unique(lnames)
  relisted <- plyr::llply(ulnames, function(x) reduced_list[lnames == x])
  
  RedC <- function(li){
    Reduce(function(...) c(...), li)
  }
  lidf <- lapply(relisted, FUN = RedC)
  names(lidf) <- ulnames
  lidf
}

MergeLiLiLiNum <- function(lilidf){
  library(plyr)
  reduced_list <- Reduce(function(...) c(...), lilidf)
  lnames <- names(reduced_list)
  ulnames <- unique(lnames)
  relisted <- plyr::llply(ulnames, function(x) reduced_list[lnames == x])
  
  RedC <- function(li){
    Reduce(function(...) c(...), li)
  }
  lidf <- lapply(relisted, FUN = RedC)
  lidf <- lapply(lidf, FUN = RedC)
  names(lidf) <- ulnames
  lidf
}

MergeBiasDf <- function(df){
  if (any(names(df) == "ML1")) {
    out <- df[ , 1:length(parameters)]
    out$ML <- rowMeans(abs(cbind(df$ML1, df$ML2)))
    out$Lasso <- rowMeans(abs(cbind(df$Lasso1, df$Lasso2)))
    out$Ridge <- rowMeans(abs(cbind(df$Ridge1, df$Ridge2)))
  } else {
    stop("Wrong dataframe entered.")
  }
  return(out)
}
