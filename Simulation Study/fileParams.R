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

prev_mat <- rbind(c(1/3, 1/3, 1/3), c(2/20, 9/20, 9/20), c(8/10, 1/10, 1/10))

epv1 <- c(3, 5, 10, 15, 20, 30, 50)
n_mat33 <- epv1 %*% t(c(1, 1, 1)) * 8 
n_mat33 <- rbind(n_mat33, n_mat33 * 2, n_mat33 * 4)

vec45   <- c(1, 4.5, 4.5)
n_mat45 <- epv1  %*% t(vec45) * 8
n_mat45 <- rbind(n_mat45, n_mat45 * 2, n_mat45 * 4)

vec80   <- c(8, 1, 1)
n_mat80 <- epv1 %*% t(vec80) * 8 
n_mat80 <- rbind(n_mat80, n_mat80 * 2, n_mat80 * 4)

n_mat <- rbind(n_mat33, n_mat45, n_mat80)

##### N
##### = Validation sample sizes.
N_mat33 <- matrix(10 ^ 4, nrow = 21, ncol = 3)
N_mat45 <- matrix(vec45 * 3000, nrow = 21, ncol = 3, byrow = T)
N_mat80 <- matrix(vec80 * 3000, nrow = 21, ncol = 3, byrow = T)

N_mat   <- rbind(N_mat33, N_mat45, N_mat80)


##### Coefficients. 
n_coefs_conds <- c(4, 8, 16)
n_coefs_vec <- rep(c(rep(n_coefs_conds[1], length(epv1)), rep(n_coefs_conds[2], length(epv1)), rep(n_coefs_conds[3], length(epv1))), 3)

cal_ec     <- 1 # length(results$params$n) - 2 code wouldn't work cause results is not loaded yet. So I have hard coded it.
cal_select <- seq(from = 1, to = cal_ec * 3 + 1, by = 3) # for selecting the right calibration slope columns in DataPrep
