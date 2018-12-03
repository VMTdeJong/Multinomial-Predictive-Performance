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

source("Tables.R")
source("Coerce.R")


write.csv(
  EstSETable(TableNames(main_res$Brier_pct), TableNames(main_res$Brier_pct_SE))
          , file = ToOneString(getwd(), "/Tables/", "MainBrier.csv"), row.names = F)
write.csv(
  EstSETable(TableNames(main_res$pdi_pct), TableNames(main_res$pdi_pct_SE))
          , file = ToOneString(getwd(), "/Tables/", "MainPDI.csv"), row.names = F)
write.csv(EstSETable(TableNames(main_res$Nagelkerke_pct), TableNames(main_res$Nagelkerke_pct_SE))
          , file = ToOneString(getwd(), "/Tables/", "AppNagelkerke.csv"), row.names = F)



sel_col <- c(5, 6)
cal_table <- CalTable(estimates = main_res$cal_slo_median, se = main_res$cal_slo_median_SE, sel.col = sel_col)
write.csv(cal_table
          , file = ToOneString(getwd(), "Tables", "MainCal.csv", sep = "/")
          , row.names = F)
