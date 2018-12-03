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

plotcal <- function(cal, lwd = 1, col = c("red", "green", "blue", "darkgrey"), pch = c(0, 1, 2, 8), 
                    measure = "Calibration Slopes", cex = 1.2, cex.text = cex, title = "", cex.legend = 1, 
                    sub = "", legend = T, ymin = 0.6, ymax = 1.2, yby = .1, 
                    x.var = "EPV1", pdf.plot = F, 
                    y.axis = T, include.ref = T, meas.margin = -1.5,
                    legend.pos = "bottomright", x.title = "", sel.col = NA, legend.type.pos = "topright",
                    default_y = 'ML-ref:0 other:1',
                    bg = "white", text.col = "black", small.margins = "T")
{
  if (is.na(ymin)) { ymin <- round(min(cal[ , -c(1:length(parameters))], na.rm = T), digits = 3)}
  if (is.na(ymax)) { ymax <- round(max(cal[ , -c(1:length(parameters))], na.rm = T), digits = 3)
  yby <- round((ymax - ymin) / 6, digits = 3) } 
  op <- par(mar = c(5,6 + meas.margin,4,2) + 0.1, bg = bg)
  if (small.margins) par(mai = c(.6, 0.6, 0.6, 0.4))
  lty <- c(2, 2, 2, 1, 2, 1) # Comparisons of 2 vs 3: solid line, others: dotted line
  lty_leg <- c("1 vs 2", "1 vs 3", "2 vs 1", "2 vs 3", "3 vs 1", "3 vs 2")
  dummy_y <- unlist(subset(cal, select = default_y))
  
  # Select the correct x-variable and labels.
  x_var <- unlist(subset(cal, select = x.var))
  
  if (x.var == "EPV1" || x.var == "EPV2") {
    xlab    <- "Events per variable"
    xlabels <- cal$EPV1
    add.n.total <- T
    bty <- "c"
  } else {
    add.n.total <- F
    bty <- "L"
    if (x.var == "cor") {
      xlab    <- "Correlations between predictors"
      xlabels <- cal$cor / 10
    } else if (x.var == "bin") {
      xlab    <- "Type of predictors"
      xlabels <- rep("Normal", nrow(cal))
      xlabels[cal$bin == 1] <- "Binary"
    } else if (x.var == "effect") {
      xlab    <- "Predictive value"
      xlabels <- rep("No", nrow(cal))
      xlabels[cal$effect == 1] <- "Yes"
    }
  }
  
  # Plot the frame
  plot(x_var, dummy_y, type = "n", col = col[1], lwd = lwd, 
       ylim = c(ymin, ymax), ylab = "", 
       col.lab = text.col,
       xlab = "",
       xaxt = "n", yaxt = "n", bty = bty, cex.lab = cex, pch = pch[1], cex = cex, cex.main = cex,
       sub = sub)
  
  # As reference: Data Generating Mechansim (By definition equal to 1)
  points(x_var, rep(1, length(x_var)), type = "b", col = col[4] , lwd = lwd, pch = pch[4])
  
  # Function for selecting columns of the data and plotting them
  .Points <- function(x_var, df, sel.meth, sel.col, lty,...) {
    plot_data <- df[ , grepl(sel.meth, names(df))] # Selects correct method.
    if (!any(is.na(sel.col))) 
    {
      plot_data <- plot_data[ , sel.col] # And the right columns
      lty <- lty[sel.col] # and the right line type
    }
    
    for (col in 1:ncol(plot_data)) {
      points(x_var, plot_data[ , col], lty = lty[col], ...)
    }
  } 
  
  .Points(x_var, cal, "ML"   , sel.col, lty = lty, type = "b", col = col[1], lwd = lwd, pch = pch[1])
  .Points(x_var, cal, "Lasso", sel.col, lty = lty, type = "b", col = col[2], lwd = lwd, pch = pch[2])
  .Points(x_var, cal, "Ridge", sel.col, lty = lty, type = "b", col = col[3], lwd = lwd, pch = pch[3])
  
  if (legend)
  { 
    if (include.ref) {
      legend_names <- c("ML", "Lasso", "Ridge", "Reference")
    } else {legend_names <- c("ML", "Lasso", "Ridge")}
    legend(legend.pos, legend = legend_names, col = col, bty = "n", title = "Method", lwd = lwd, lty = 1, pch = pch, text.col = text.col, cex = cex * cex.legend) 
    
    legend(legend.type.pos, legend = lty_leg[sel.col], col = text.col, bty = "n", title = "Comparison", lwd = lwd, lty = lty[sel.col], text.col = text.col, cex = cex * cex.legend)
  }
  
  # Axes:
  # x - axis
  title(xlab = xlab, line = 2, cex.lab = cex.text)
  axis(1, at = c(x_var), cex.axis = cex, labels = xlabels, col = text.col, col.axis = text.col) 

  
  # y - axis
  if (y.axis) { axis(2, at = seq(ymin, ymax, by = yby ),las = 2, cex.axis = cex, col = text.col, col.axis = text.col) }
  title(ylab = measure, line = 3.8 + meas.margin, cex.lab = cex, col.lab = text.col)
  
  # Top x - axis
  title(main = title, adj = 0, line = 3.2, col.main = text.col, cex.main = cex)
  if (add.n.total) {
    axis(side = 3, at = x_var, labels = round(cal$n_total), cex.axis = cex, col = text.col, col.axis = text.col)
    if (pdf.plot) {
      mtext("Total sample size", side = 3, line = 2.0, cex = cex * cex.text, adj = 1, col = text.col)
      mtext(x.title, side = 4, line = .3, cex = cex * cex.text, adj = 1, font = 2, col = text.col)
    } else {
      mtext("Total sample size", side = 3, line = 2.5, cex = 0.8, adj = 1, col = text.col) # For some reason pdf handles cex for mtext different than tikz
      mtext(x.title, side = 4, line = .5, cex = 0.8, adj = 1, font = 2, col = text.col)
    }
  } 
  par(op)
}


PlotAllCal <- function(cal,...) {
  for (i in 1:9) {
    plotcal(cal[(1:7) + (i - 1) * 7, ],...)
  }
}

# PlotAllCal(cal_slo_median)

# plotcal2 <- function(cal, lwd = 1, col = c("red", "green", "blue", "darkgrey"), pch = c(0, 1, 2, 8), 
#                      measure = "Calibration Slopes", cex = 1.2, title = "", 
#                      sub = "", legend = T, ymin = 0.6, ymax = 1.1, yby = .1, 
#                      xlab = "Events per variable",  x.meas = "EPV1", pdf.plot = F, 
#                      add.n.total = T, bty = "c", y.axis = T, include.ref = T, meas.margin = 0,
#                      legend.pos = "bottomright", x.title = "", sel.col = NA)
# {
#   if (is.na(ymin)) { ymin <- round(min(cal[ , -c(1:length(parameters))], na.rm = T), digits = 3)}
#   if (is.na(ymax)) { ymax <- round(max(cal[ , -c(1:length(parameters))], na.rm = T), digits = 3)
#   yby <- round((ymax - ymin) / 6, digits = 3) } 
#   x_var <- unlist(subset(cal, select = x.meas))
#   
#   if (any(names(cal) == "bin")) { 
#     xlabels <- rep("Normal", nrow(cal)) # Normal
#     xlabels[cal$bin == 1] <- "Binary"
#   }
#   
#   if (xlab == "bin") xlab <- "Type of predictors"
#   if (substring(xlab, 1,1) == "C" || substring(xlab, 1,1) == "c") xlabels <- cal$cor / 10 else if (xlab != "bin" && xlab != "Type of predictors") {xlabels <- cal$EPV1}
#   
#   op <- par(mar = c(5,6 + meas.margin,4,2) + 0.1)
#   
#   plot(x_var, cal$oos_ML1, type = "n", col = col[1], lwd = lwd, 
#        ylim = c(ymin, ymax), ylab = "", 
#        xlab = xlab,
#        xaxt = "n", yaxt = "n", bty = bty, cex.lab = cex, pch = pch[1], cex = cex, 
#        sub = sub)
#   
#   
#   # as Reference: Data Generating Mechansim (By definition equal to 1)
#   points(x_var, rep(1, length(x_var)), type = "b", col = col[4] , lwd = lwd, pch = pch[4])
#   
#   # Function for selecting columns of the data and plotting them
#   .Points <- function(x_var, df, select,...) {
#     plot_data <- df[ , grepl(select, names(df))] # Selects the column of the correct method.
#     for (col in 1:ncol(plot_data)) {
#       if (!any(is.na(sel.col))) if (!any(sel.col == col)) next # if NA, ignore selection method.
#       points(x_var, plot_data[ , col], ...)
#       }
#   } 
#   
#   .Points(x_var, cal, "ML"   , type = "b", col = col[1], lwd = lwd, pch = pch[1])
#   .Points(x_var, cal, "Lasso", type = "b", col = col[2], lwd = lwd, pch = pch[2])
#   .Points(x_var, cal, "Ridge", type = "b", col = col[3], lwd = lwd, pch = pch[3])
#   
#   
#   if (legend)
#   { 
#     if (include.ref) {
#       legend_names <- c("ML", "lasso", "ridge", "reference")
#     } else {legend_names <- c("ML", "lasso", "ridge")}
#     legend(legend.pos, legend = legend_names, col = col, bty = "n", title = "Method", lwd = lwd, lty = 1, pch = pch) 
#   }
#   
#   # Axes:
#   # x - axis
#   axis(1, at = c(x_var), cex.axis = cex, labels = xlabels) 
#   
#   # y - axis
#   if (y.axis) { axis(2, at = seq(ymin, ymax, by = yby ),las = 2, cex.axis = cex) }
#   title(ylab = measure, line = 3.8 + meas.margin, cex.lab = cex)
#   
#   
#   # Top x - axis
#   title(main = title, adj = 0, line = 2.8)
#   if (add.n.total) {
#     axis(side = 3, at = x_var, labels = round(cal$n_total), cex.axis = cex)
#     if (pdf.plot) {
#       mtext("Total sample size", side = 3, line = 2.5, cex = cex , adj = 1)
#       mtext(x.title, side = 4, line = .5, cex = cex, adj = 1)
#     } else {
#       mtext("Total sample size", side = 3, line = 2.5, cex = 0.8, adj = 1) # For some reason pdf handles cex for mtext different than tikz
#       mtext(x.title, side = 4, line = .5, cex = 0.8, adj = 1)
#     }
#   } 
#   par(op)
# }