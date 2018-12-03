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

# # # For on hpc:
c_args  <- commandArgs(trailingOnly = TRUE) # These are the arguments given through LINUX 

if (length(c_args) == 0) {
  id <- 1 
  maxit <- 1000
  reps <- 150
  require(parallel)
} else {
  id      <- as.numeric(c_args[1]) 
  maxit   <- as.numeric(c_args[2]) # Number of iterations each optimizing function is allowed
  reps    <- as.numeric(c_args[3]) # Number of times each optimization is repeated
  lib.loc <- "/path/R/library" # To be replaced
  require(parallel, lib.loc = lib.loc)
}

if (length(c_args) == 0) {
  source("MakeData.R")
  source("FileParams.R")
  
} else {
  proj_loc <- "/path" # To be replaced
  source(paste(proj_loc, "MakeData.R"   , sep = "/"))
  source(paste(proj_loc, "FileParams.R"   , sep = "/"))
  
  proj_loc <- paste(paste(proj_loc, "Optim", sep = "/"), id, sep = "/")
}

# Note: for three outcomes only!!!
Lossf <- function(alphas, n_coefs, n_vec, prevs){
  coefs <- MakeCoefs(alphas = alphas, n_coefs = n_coefs)
  
  dl <- MakeRetroData(n = n_vec, outcome = c(0, 1, 2), 
                      coefs = coefs, 
                      n_bi_var = 0, Brier = T)
  cm <- colMeans(dl$y_probs)
  
  return(sqrt(sum((cm - prevs) ^ 2)))
}

Lossf1 <- function(alpha, n_coefs, n_vec, prevs){
  alphas <- c(alpha, alpha)
  return(Lossf(alphas, n_coefs, n_vec, prevs))
}

MyOptim1 <- function(n_vec, prevs, n_coefs, 
                     par = 0, maxit = 10, upper = 5, lower = -5)
{
  if (!n_coefs %% 4 == 0) {stop("n_coefs must be multiple of 4, or zero.")}
  
  opt <- optim(par, Lossf1, n_coefs = n_coefs, n_vec = n_vec, prevs = prevs, method = "Brent", upper = upper, lower = lower, control = list(maxit = maxit))
  
  return(opt)
}

large_n <- c(3 * 1e5, 3 * 1e5, 3 * 1e5)
small_n <- c(3 * 1e4, 3 * 1e4, 3 * 1e4)
n_coefs_vec <- c(4, 8, 16)

begin <- proc.time()
#### Optimisation using Brent.
## First optimisation. Starting values at 0. 
start_par <- c(0)
optDF <- data.frame(matrix(nrow = 0, ncol = 6))
names(optDF) <- c("rmse", "par1", "par2", "prev_id", "prevs", "n_coefs")


load("OptAlph.txt") # Error is fine.

for (c in 1:3) {
  for (p in 1:3) {
    par   <- start_par
    par   <- MyOptim1(n_vec = round(small_n * n_mat[p * 21, ]/sum(n_mat[p * 21, ])), prevs = prev_mat[p, ], n_coefs = n_coefs_vec[c], par = par, maxit = maxit, upper = 5, lower = -5)$par
    gc()
    opt   <- data.frame("rmse" = testOptim(c(par, par), n_coefs = n_coefs_vec[c], prevs = prev_mat[p, ], n_vec = round(large_n * n_mat[p * 21, ]/sum(n_mat[p * 21, ]))), "par1" = par, "par2" = par, "prev_id" = p, "prevs" = prev_mat[p, 3]/prev_mat[p, 1], "n_coefs" = n_coefs_vec[c])
    optDF <- rbind(optDF, opt)
    gc()
    print(optDF)
  }
}

if (length(c_args) == 0) {
  save(optDF, file = "OptAlph.txt")
} else {
  save(optDF, file = paste(proj_loc, "OptAlph.txt", sep = "/"))
}

t_one_run <- proc.time() - begin
t_one_run

for (i in 2:reps)
{
  for (c in 1:3) {
    for (p in 1:3) {
      par   <- optDF[optDF$prev_id == p & optDF$n_coefs == n_coefs_vec[c], ][which.min(optDF[optDF$prev_id == p & optDF$n_coefs == n_coefs_vec[c], ]$rmse), ]$par1 
      upper <- par + 2
      lower <- par - 2
      par   <- MyOptim1(n_vec = round(small_n * n_mat[p * 21, ]/sum(n_mat[p * 21, ])), prevs = prev_mat[p, ], n_coefs = n_coefs_vec[c], par = par, maxit = maxit, upper = upper, lower = lower)$par
      gc()
      opt   <- data.frame("rmse" = testOptim(c(par, par), n_coefs = n_coefs_vec[c], prevs = prev_mat[p, ], n_vec = round(large_n * n_mat[p * 21, ])/sum(n_mat[p * 21, ])), "par1" = par, "par2" = par, "prev_id" = p, "prevs" = prev_mat[p, 3]/prev_mat[p, 1], "n_coefs" = n_coefs_vec[c])
      gc()
      optDF <- rbind(optDF, opt)
      print(optDF)
    }
  }
  if (length(c_args) == 0) {
    save(optDF, file = "OptAlph.txt")
  } else {
    save(optDF, file = paste(proj_loc, "OptAlph.txt", sep = "/"))
  }
}

end <- proc.time()
optimTime <- end - begin
optimTime

alph_mean <- data.frame(matrix(ncol = 3, nrow = 9))
names(alph_mean) <- c("alph", "prev_id", "n_coefs")
r <- 0
for (c in 1:3) {
  n_coefs <- n_coefs_vec[c]
  for (p in 1:3) {
    r <- r + 1
    prevs <- prev_mat[p, ]
    par_mean <- mean(optDF[optDF$prev_id == p & optDF$n_coefs == n_coefs, ]$par1)
    alph_mean[r, ] <- c(par_mean, p, n_coefs)
  }
} 

if (length(c_args) == 0) {
  save(alph_mean,file = "alph_mean.txt")
} else {
  save(alph_mean,file = paste(proj_loc, "alph_mean.txt",sep = "/"))
}

#load("alph_mean.txt")
