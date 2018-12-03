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

# Plotting functions.
source("plotcal.R")
source("Plot.R")
source("Marg.R")

##################################################################################################################################
##################################################################################################################################
########################################               Main Analysis                 #############################################
##################################################################################################################################
##################################################################################################################################


#################### Figure titles

coef4title  <- "4 predictors"#. Prevalence marginalized out."
coef8title  <- "8 predictors"#. Prevalence marginalized out."
coef16title <- "16 predictors"#. Prevalence marginalized out."

prev1title   <- "Equal frequencies"  # <- "Frequencies: $33\\%, 33\\%, 33\\%$"#.Number of predictors marginalized out." # 
prev4.5title <- "One small category" # <- "Frequencies: $10\\%, 45\\%, 45\\%$"#.Number of predictors marginalized out." # 
prev8title   <- "One large category" # <- "Frequencies: $80\\%, 10\\%, 10\\%$"#.Number of predictors marginalized out." # 

Nagel_meas <- "Nagelkerke"


path <- paste(getwd(), "Figures/EPS", sep = "/")

sel_col <- c(5, 6)

##################################################################################################################################
########################################          Median calibration slopes          #############################################
##################################################################################################################################

csm <- main_res$cal_slo_median

#, fallback_resolution = 800 # Does not work in this version of R
cairo_ps(paste(path, "Figure_1_VMTdeJong_MainCalSlopes.eps", sep = "/"), family = "serif" 
         , width = 8, height = 11, pointsize = 15)
side_y_label <- "Median calibration slopes"
other_y_label <-  ""
par(mfrow = c(3,3))
cex <- .8

plotcal(csm[csm$n_coefs == 4  & csm$prev == 1 , ], measure = side_y_label , title = coef4title, 
        x.title = "Equal frequencies", sel.col = sel_col, cex = cex, pdf.plot = T) 
plotcal(csm[csm$n_coefs == 8  & csm$prev == 1 , ], measure = other_y_label, title = coef8title, 
        x.title = "Equal frequencies", sel.col = sel_col, cex = cex, pdf.plot = T) 
plotcal(csm[csm$n_coefs == 16 & csm$prev == 1 , ], measure = other_y_label, title = coef16title, 
        x.title = "Equal frequencies", sel.col = sel_col, cex = cex, pdf.plot = T) 

plotcal(csm[csm$n_coefs == 4  & csm$prev == 4.5  , ], measure = side_y_label , title = coef4title, 
        x.title = "One small outcome category", sel.col = sel_col, cex = cex, pdf.plot = T) 
plotcal(csm[csm$n_coefs == 8  & csm$prev == 4.5  , ], measure = other_y_label, title = coef8title, 
        x.title = "One small outcome category", sel.col = sel_col, cex = cex, pdf.plot = T) 
plotcal(csm[csm$n_coefs == 16 & csm$prev == 4.5  , ], measure = other_y_label, title = coef16title, 
        x.title = "One small outcome category", sel.col = sel_col, cex = cex, pdf.plot = T)

plotcal(csm[csm$n_coefs == 4  & csm$prev == 0.125  , ], measure = side_y_label , title = coef4title, 
        x.title = "One large outcome category", sel.col = sel_col, cex = cex, pdf.plot = T) 
plotcal(csm[csm$n_coefs == 8  & csm$prev == 0.125  , ], measure = other_y_label, title = coef8title, 
        x.title = "One large outcome category", sel.col = sel_col, cex = cex, pdf.plot = T) 
plotcal(csm[csm$n_coefs == 16 & csm$prev == 0.125  , ], measure = other_y_label, title = coef16title, 
        x.title = "One large outcome category", sel.col = sel_col, cex = cex, pdf.plot = T)

dev.off()


##################################################################################################################################
########################################                      PDI                    #############################################
##################################################################################################################################

pdimin <- -6
pdimax <- 8

cairo_ps(paste(path, "Figure_2_VMTdeJong_MainPDI.eps", sep = "/"), family = "serif" #, fallback_resolution = 800
         , width = 8, height = 11, pointsize = 15)

par(mfcol = c(3,2))
pdi_mea <- "PDI"
# prevalence marginalized out. Conditional on number of coefficients
P1 <- Marg(main_res$pdi, MARGIN = c("n_coefs", "EPV1"), class.name = "PDI" )

Plot(P1[P1$n_coefs == 4 , ], measure = pdi_mea, title = coef4title , ymin = pdimin, ymax = pdimax, tex = F)
Plot(P1[P1$n_coefs == 8 , ], measure = pdi_mea, title = coef8title , ymin = pdimin, ymax = pdimax, tex = F)
Plot(P1[P1$n_coefs == 16, ], measure = pdi_mea, title = coef16title, ymin = pdimin, ymax = pdimax, tex = F)

# prevalence marginalized out. Conditional on number of coefficients
P2 <- Marg(main_res$pdi, MARGIN = c("prev", "EPV1"), class.name = "PDI" )

Plot(P2[P2$prev == 1.000, ], measure = pdi_mea, title = prev1title  , ymin = pdimin, ymax = pdimax, tex = F)
Plot(P2[P2$prev == 4.500, ], measure = pdi_mea, title = prev4.5title, ymin = pdimin, ymax = pdimax, tex = F)
Plot(P2[P2$prev == 0.125, ], measure = pdi_mea, title = prev8title  , ymin = pdimin, ymax = pdimax, tex = F)

dev.off()


##################################################################################################################################
########################################                 Brier score                 #############################################
##################################################################################################################################

Briermax <-  8
Briermin <- -8
cairo_ps(paste(path, "Figure_3_VMTdeJong_MainBrier.eps", sep = "/"), family = "serif"#, fallback_resolution = 800
         , width = 8, height = 11, pointsize = 15)
par(mfcol = c(3,2))

# prevalence marginalized out. Conditional on number of coefficients
B1 <- Marg(main_res$Brier, MARGIN = c("n_coefs", "EPV1"), class.name = "Brier" )

Plot(B1[B1$n_coefs == 4 , ], measure = "Brier score", title = coef4title,  ymin = Briermin, ymax = Briermax, tex = F)
Plot(B1[B1$n_coefs == 8 , ], measure = "Brier score", title = coef8title,  ymin = Briermin, ymax = Briermax, tex = F)
Plot(B1[B1$n_coefs == 16, ], measure = "Brier score", title = coef16title, ymin = Briermin, ymax = Briermax, tex = F)

# number of coefficients marginalized out. Conditional on prevalence
B2 <- Marg(main_res$Brier, MARGIN = c("prev", "EPV1"), class.name = "Brier" )

Plot(B2[B2$prev == 1.000, ], measure = "Brier score", title = prev1title,   ymin = Briermin, ymax = Briermax, tex = F)
Plot(B2[B2$prev == 4.500, ], measure = "Brier score", title = prev4.5title, ymin = Briermin, ymax = Briermax, tex = F)
Plot(B2[B2$prev == 0.125, ], measure = "Brier score", title = prev8title,   ymin = Briermin, ymax = Briermax, tex = F)

dev.off()

##################################################################################################################################
########################################                Nagelkerke R^2               #############################################
##################################################################################################################################

Nagelmax <-  16
Nagelmin <- -22
cairo_ps(paste(path, "Figure_App_VMTdeJong_MainNagelkerke.eps", sep = "/"), family = "serif" #, fallback_resolution = 800
         , width = 8, height = 11, pointsize = 15)

par(mfcol = c(3,2))

# prevalence marginalized out. Conditional on number of coefficients
Nagel1 <- Marg(main_res$Nagelkerke, MARGIN = c("n_coefs", "EPV1"), class.name = "Nagelkerke" )

Plot(Nagel1[Nagel1$n_coefs == 4 , ], measure = Nagel_meas, title = coef4title, ymin = Nagelmin, ymax = Nagelmax
     , tex = F, legend.m.loc = "bottomright", legend.s.loc = "topright")
Plot(Nagel1[Nagel1$n_coefs == 8 , ], measure = Nagel_meas, title = coef8title, ymin = Nagelmin, ymax = Nagelmax
     , tex = F, legend.m.loc = "bottomright", legend.s.loc = "topright")
Plot(Nagel1[Nagel1$n_coefs == 16, ], measure = Nagel_meas, title = coef16title, ymin = Nagelmin, ymax = Nagelmax
     , tex = F, legend.m.loc = "bottomright", legend.s.loc = "topright")

# number of coefficients marginalized out. Conditional on prevalence
Nagel2 <- Marg(main_res$Nagelkerke, MARGIN = c("prev", "EPV1"), class.name = "Nagelkerke" )

Plot(Nagel2[Nagel2$prev == 1.000, ], measure = Nagel_meas, title = prev1title, ymin = Nagelmin, ymax = Nagelmax
     , tex = F, legend.m.loc = "bottomright", legend.s.loc = "topright")
Plot(Nagel2[Nagel2$prev == 4.500, ], measure = Nagel_meas, title = prev4.5title, ymin = Nagelmin, ymax = Nagelmax
     , tex = F, legend.m.loc = "bottomright", legend.s.loc = "topright")
Plot(Nagel2[Nagel2$prev == 0.125, ], measure = Nagel_meas, title = prev8title, ymin = Nagelmin, ymax = Nagelmax
     , tex = F, legend.m.loc = "bottomright", legend.s.loc = "topright")

dev.off()

##################################################################################################################################
##################################################################################################################################
########################################          Sensitivity Analysis: Correlations #############################################
##################################################################################################################################
##################################################################################################################################

cairo_ps(paste(path, "Figure_4_VMTdeJong_SensCor.eps", sep = "/"), family = "serif" #, fallback_resolution = 800
         , width = 8, height = 11, pointsize = 15)
par(mfrow = c(3,1))

########################################          Median calibration slopes          #############################################

sens.csm <- cor_res$cal_slo_median
side_y_label <- "Median calibration slopes"
other_y_label <-  ""


plotcal(sens.csm, measure = side_y_label , title = "", sel.col = sel_col, x.var = "cor", ymax = 1.2
        , legend.type.pos = "topleft", pdf.plot = T, small.margins = F, meas.margin = 0, cex = 1.2, cex.text = 1.2, 
        cex.legend = 1/1.2) 

########################################                      PDI                    #############################################

Plot(cor_res$pdi, measure = "PDI", title = "", x.var = "cor"
     , legend.m.x = 0, legend.m.y = 2.5
     , legend.s.x = 0, legend.s.y = -.7
     , ymin = -3, ymax = 3.5, tex = F)

########################################                 Brier score                 #############################################

Plot(cor_res$Brier, measure = "Brier score", title = "", x.var = "cor"
     , legend.m.x = 0, legend.m.y = 1.8
     , legend.s.x = 0, legend.s.y = -.7
     , ymin = -2.5, ymax = 2.5, tex = F)


dev.off()

##################################################################################################################################
##################################################################################################################################
########################################          Sensitivity Analysis: Type: Binary #############################################
##################################################################################################################################
##################################################################################################################################
source("BoxPlot.R")



# unique(acs_norm$Method)
cal_plot_names <- sort(paste(rep(c("ML,", "Lasso,", "Ridge,"), each = 2) , rep(c("cat. 3 vs 1", "cat. 3 vs 2"), 3), sep = " "))

cairo_ps(paste(path, "Figure_5_VMTdeJong_SensBin.eps", sep = "/"), family = "serif" #, fallback_resolution = 800
         , width = 8, height = 11, pointsize = 15)

p <- par(mfrow = c(3, 2))
BoxPlot(acs_norm, acs_bin, pars = list(ylim = c(0, 5)), ylab = rep("Calibration slope", 2), names = cal_plot_names, h = 1, 
        main = c("Normal predictors", "Binary predictors"), col = c("green", "green", "red", "red", "blue", "blue"))

# p <- par(mfrow = c(1, 2))
BoxPlot(pdi_pct_norm, pdi_pct_bin, pars = list(ylim = c(-30, 25), cex = 1.2), xlab = "", ylab = "PDI, % difference of reference", 
        h = 0)
BoxPlot(Brier_pct_norm, Brier_pct_bin, pars = list(ylim = c(-14,10), cex = 1.2), xlab = "", 
        ylab = "Brier score, % difference of reference", h = 0)

par(p)

dev.off()



