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


library(dplyr)
library(tidyr)
library(abind)
source("Coerce.R")
source("Plot.R")

######################### Plotting functions

BoxPlot <- function(data1, data2, 
                    col = c("green", "red", "blue", "green", "red", "blue"), ref.col = "darkgrey",
                    pars = NULL, 
                    xlab = "", ylab = rep("Measure", 2),
                    names = sort(c("Within-sample ML"   , "Out-of-sample ML",  # Sorting appears to be the default in boxplot
                                   "Within-sample lasso", "Out-of-sample lasso",
                                   "Within-sample ridge", "Out-of-sample ridge")),
                    par = F,
                    h = NULL,
                    main = "", main1 = main, main2 = main, # main1 and 2 are obsolete.
                    x.angle = 28,
                    bg = "white",
                    text.col = "black")
  {
  if (par) p <- par(mfrow = c(1,2), bg = bg)

  plotOne <- function(data, i, names) # Note that plotOne takes its other parameters from the BoxPlot environment.
    {
    # Remove the names = rep() part to check the names.
    # Note that the names fall off the picture in R, but not in the eps file :)
    boxplot(formula(Value ~ Method), data = data, pars = pars, xlab = xlab, ylab = ylab[i], col = col, names = rep("", 6), 
       las = 2, main = main[i], col.lab = text.col, col.axis = text.col, col.main = text.col, border = text.col)
    axis(side = 1, col.ticks = text.col, labels = F, col = text.col)
    axis(side = 2, col.ticks = text.col, labels = F, col = text.col)
  abline(h = h, col = ref.col, lwd = 2)
  text(x = 1:6, y = par()$usr[3] - 0.1 * (par()$usr[4] - par()$usr[3]),
  labels = names, srt = x.angle, adj = 1, xpd = T, col = text.col)
  }
  
  plotOne(data1, 1, names)
  plotOne(data2, 2, names)
  
  if (par) par(p)
}

######################### Data Prepping


###### The Data

##### Calibration slopes
# Data
# Change names if different columns are selected above!!!!
selected_slopes <- c(5, 6)
selected_columns <- length(bin_res$last_par) + c(selected_slopes, selected_slopes + 6, selected_slopes + 12)

acs <- bin_res$all_cal_slo[ , c(1:length(bin_res$last_par), selected_columns), ]
acs_norm <- gather(data = as.data.frame(t(acs[1, , ])), key = "Method", value = "Value", ML.ref.2_other.0:Ridge.ref.2_other.1) 
acs_bin  <- gather(data = as.data.frame(t(acs[2, , ])), key = "Method", value = "Value", ML.ref.2_other.0:Ridge.ref.2_other.1) 
acs_bin_wo <- acs_bin[AreTrue(acs_bin$Value < 5 & acs_bin$Value > 0), ]

acs_bin_wo_log <- acs_bin_wo
acs_bin_wo_log$Value <- log(acs_bin_wo$Value)

##### Brier score
# Data
# DGM is removed but is zero anyways
Brier_pct_all <- abind(lapply(lapply(apply(bin_res$all_Brier[, , ], MARGIN = 3, FUN = ShowHow), '[', "df"),'[[', 1), along = 3)
Brier_pct_norm <- gather(data = as.data.frame(t(Brier_pct_all[1, , ])), key = "Method", value = "Value", ws_ML:oos_Ridge, -DGM) 
Brier_pct_bin  <- gather(data = as.data.frame(t(Brier_pct_all[2, , ])), key = "Method", value = "Value", ws_ML:oos_Ridge, -DGM)


##### PDI
# Data
# DGM is removed but is zero anyways
pdi_pct_all <- abind(lapply(lapply(apply(bin_res$all_pdi[, , ], MARGIN = 3, FUN = ShowHow), '[', "df"),'[[', 1), along = 3)
pdi_pct_norm <- gather(data = as.data.frame(t(pdi_pct_all[1, , ])), key = "Method", value = "Value", ws_ML:oos_Ridge, -DGM) 
pdi_pct_bin  <- gather(data = as.data.frame(t(pdi_pct_all[2, , ])), key = "Method", value = "Value", ws_ML:oos_Ridge, -DGM)


