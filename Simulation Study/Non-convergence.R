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

nrep <- 2000

## Overall non-convergence
AvgNonConv <- function(ll) paste((1 - sum(ll$Brier$I)/(nrow(ll$Brier)*nrep))*100, "%", sep = "") # Average non-convergence in %
MaxNonConv <- function(ll) paste((nrep - min(ll$Brier$I)) / nrep * 100, "%", sep = "")     # Max non-convergence in %


AvgNonConv(main_res) # Average % non-convergence in main
MaxNonConv(main_res) # Maximum % non-convergence in main

AvgNonConv(bin_res)
MaxNonConv(bin_res)

AvgNonConv(cor_res)
MaxNonConv(cor_res)

## Non-convergence of calibration slopes.
AvgNonConvCal <- function(ll) paste((1 - sum(ll$cal_slo_median$success)/(nrow(ll$cal_slo_median)*nrep) ) * 100, "%", sep = "")
MaxNonConvCal <- function(ll) paste((nrep - min(ll$cal_slo_median$success)) / nrep * 100, "%", sep = "")

AvgNonConvCal(main_res) # Average % non-convergence calibration slopes in main
MaxNonConvCal(main_res) # Maximum % non-convergence calibration slopes in main

AvgNonConvCal(bin_res) 
MaxNonConvCal(bin_res)

AvgNonConvCal(cor_res) 
MaxNonConvCal(cor_res)
