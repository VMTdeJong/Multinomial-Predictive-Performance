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

# Marg() Marginalizes dataframes, keeping 1 or 2 margins. Additionally a variable to sum (rather than average) over can be provided.
Marg <- function(df, MARGIN = c("prev", "EPV1"), row.names = T, sort = MARGIN[1], class.name){
  
  by <- unlist(subset(df, select = MARGIN[1]))
  li <- split(df, by)
  
  if (length(MARGIN) > 1) {
    lili <- list()
    for (l in 1:length(li)) {
      by <- unlist(subset(li[[l]], select = MARGIN[2]))  
      lili[[l]] <- split(li[[l]], by)
    }

    li <- Reduce(function(...) c(...), lili)
  } 
  
  if (row.names) {
    row_names <- lapply(lapply(li, subset, select = "id"), as.character)  
  } else row_names <- NULL
  
  m   <- data.frame(t(sapply(li, colMeans)), row.names = row_names)
  m$I <- data.frame(t(sapply(li, colSums)) , row.names = row_names)$I

  if (any(names(df) == "NA_ML"))    { m$NA_ML    <- data.frame(t(sapply(li, colSums)) , row.names = row_names)$NA_ML   }
  if (any(names(df) == "NA_Lasso")) { m$NA_lasso <- data.frame(t(sapply(li, colSums)) , row.names = row_names)$NA_Lasso}
  if (any(names(df) == "NA_Ridge")) { m$NA_ridge <- data.frame(t(sapply(li, colSums)) , row.names = row_names)$NA_Ridge}
  
  # class
  class(m) <- c(class(m), "marg")
  if (!missing(class.name)) { class(m) <- c(class(m), class.name) }
  
  if (!missing(sort)) { return(m[order(subset(m, select = sort)), ])
} else return(m)
}

# 
# Marg(Brier, "EPV1", row.names = F)
# Marg(Brier, "n_coefs", row.names = F)
# Marg(Brier, "prev", row.names = F)
# 
# Marg(Brier, c("n_coefs", "EPV1"))
# Marg(Brier, c("prev"   , "EPV1"))