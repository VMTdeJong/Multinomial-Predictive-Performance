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

Plot <- function(df, lwd = 1, col = c("red", "green", "blue", "darkgrey"), 
                  pch = c(0, 1, 2, 8), ylab = "", cex = 1.2, sub = "", legend = T,
                  title = "",
                  y_lab_digits = 1, ylab_length = 7, ymin = NA, ymax = NA, measure,
                  how = c("pct.inc", "raw", "subtract.ref", "pct"),
                  pdf.plot = F, y.axis = T,
                  meas.margin = 0, x.var = "EPV1", DGM.mean = F, 
                  legend.m.loc = "topright", legend.s.loc = "bottomright",
                  legend.m.x = NA, legend.m.y = NA,
                  legend.s.x = NA, legend.s.y = NA, legend.outside = F,
                  bg = "white", text.col = "black", tex = T)
{
  if (missing(measure)) {
    if (!is.na(class(df)[3])) { measure <- class(df)[3]
    } else {measure <- ""}
  }
  
  if (any(names(df) == "ws_ML")) ws_sample <- T else ws_sample <- F
  if (x.var != "EPV1" && x.var != "EPV2" && x.var != "EPV3") add.n.total <- F
  op <- par(mar = c(5,6 + meas.margin,4,2) + 0.1, bg = bg)
  if (DGM.mean) df$DGM <- rep(mean(df$DGM), length(df$DGM))

  # Select the correct x-variable and labels.
  x_var <- df[[x.var]] 
  if (x.var == "EPV1" || x.var == "EPV2" || x.var == "EPV3") {
    xlab    <- expression("Events per variable")
    xlabels <- df[[x.var]]
    add.n.total <- T
    bty <- "c"
  } else {
    add.n.total <- F
    bty <- "L"
    if (x.var == "cor") {
      xlab    <- expression("Correlations between predictors")
      xlabels <- df$cor / 10
    } else if (x.var == "bin") {
      xlab    <- expression("Type of predictors")
      xlabels <- rep("Normal", nrow(df))
      xlabels[df$bin == 1] <- "Binary"
    } else if (x.var == "effect") {
      xlab    <- expression("Predictive value")
      xlabels <- rep("No", nrow(df))
      xlabels[df$effect == 1] <- "Yes"
    }
  }
  
temp_data <- ShowHow(df, how, measure, ws_sample, tex = tex)
df <- temp_data$df
if (how[1] != "raw") measure <- temp_data$measure
  
  if (is.na(ymin)) { ymin <- min(df[ , -c(1:length(last_par))], na.rm = T)} 
  if (is.na(ymax)) { ymax <- max(df[ , -c(1:length(last_par))], na.rm = T)} # last_par is assigned in LoadData.R
  
  plot(x_var, df$oos_ML, type = "b", col = col[1], lwd = lwd, 
       ylim = c(ymin, ymax), ylab = "", xlab = xlab, col.lab = text.col, main = "", 
       xaxt = "n", yaxt = "n", bty = bty, cex.lab = cex, pch = pch[1],  sub = comment(df))
  
  points(x_var, df$oos_Lasso, type = "b", col = col[2], pch = pch[2], lwd = lwd)
  points(x_var, df$oos_Ridge, type = "b", col = col[3], pch = pch[3], lwd = lwd)
  
  if (ws_sample) {
    points(x_var, df$ws_ML, type = "b", col = col[1], pch = pch[1], lty = 3, lwd = lwd)
    points(x_var, df$ws_Lasso, type = "b", col = col[2], pch = pch[2], lty = 3, lwd = lwd)
    points(x_var, df$ws_Ridge, type = "b", col = col[3], pch = pch[3], lty = 3, lwd = lwd)
  }
  
  par(xpd = legend.outside)
  if (legend) {
    if (!is.na(legend.m.x) && !is.na(legend.m.y)) {
      legend(x = legend.m.x, y = legend.m.y, legend = c("ML", "Lasso", "Ridge", "Reference"), col = col, text.col = text.col,
             bty = "n", title = "Method", lwd = lwd, lty = 1, pch = pch)
    } else {
      legend(legend.m.loc, legend = c("ML", "Lasso", "Ridge", "Reference"), col = col, text.col = text.col,
             bty = "n", title = "Method", lwd = lwd, lty = 1, pch = pch)
    }
    
    if (!is.na(legend.s.x) && !is.na(legend.s.y)) {
      legend(x = legend.s.x, y = legend.s.y, legend = c("Within-sample", "Out-of-sample"), bty = "n", lty = c(3, 1), lwd = lwd, text.col = text.col, col = text.col)
    } else {
      legend(legend.s.loc, legend = c("Within-sample", "Out-of-sample"), bty = "n", lty = c(3, 1), lwd = lwd, text.col = text.col, col = text.col) 
    }
  }
  
  # Axes:
  # x - axis
  axis(1, at = c(x_var), cex.axis = cex, labels = xlabels, col = text.col, col.axis = text.col)
  
  # y - axis
  y_seq <- round(seq(ymin, ymax, length.out = ylab_length), digits = y_lab_digits)
    
  if (y.axis) { axis(2, at = y_seq, las = 2, cex.axis = cex, col = text.col, col.axis = text.col)
    title(ylab = measure, line = 3.8 + meas.margin, cex.lab = cex, col.lab = text.col)
  }
  
  # Top x-axis
  title(main = title, adj = 0, line = 2.8, col.main = text.col)
  if (add.n.total) {
    axis(side  = 3, at = x_var, labels = round(df$n_total), cex.axis = cex, col = text.col, col.axis = text.col) 
    if (pdf.plot) {
      mtext("Total sample size", side = 3, line = 2.5,cex = cex, adj = 1, col = text.col)    
    } else {
      mtext("Total sample size", side = 3, line = 2.5,cex = 0.8, adj = 1, col = text.col) # For some reason pdf handles cex for mtext different than tikz
    }
  }
  
  points(x_var, df$DGM, type = "b", col = col[4], pch = pch[4], lwd = lwd) # DGM
  par(op)
}

ShowHow <- function(df, how = c("pct.inc", "raw", "subtract.ref", "pct"), measure = "", ws_sample = T, tex = T) {
  df <- as.data.frame(df)
  how <- as.character(how[1])
  if (how == "subtract.ref") {
    if (ws_sample) {
      df$ws_ML     <- df$ws_ML     - df$DGM      
      df$ws_Lasso  <- df$ws_Lasso  - df$DGM
      df$ws_Ridge  <- df$ws_Ridge  - df$DGM
    }
    df$oos_ML    <- df$oos_ML    - df$DGM
    df$oos_Lasso <- df$oos_Lasso - df$DGM
    df$oos_Ridge <- df$oos_Ridge - df$DGM
    
    df$DGM <- 0
  } else if (how == "pct") {
    df$DGM <- df$DGM / 100
    if (ws_sample) {
      df$ws_ML     <- df$ws_ML     / df$DGM
      df$ws_Lasso  <- df$ws_Lasso  / df$DGM
      df$ws_Ridge  <- df$ws_Ridge  / df$DGM
    }
    df$oos_ML    <- df$oos_ML / df$DGM
    df$oos_Lasso <- df$oos_Lasso / df$DGM
    df$oos_Ridge <- df$oos_Ridge / df$DGM
    df$DGM <- 100
    if (!measure == "") {if (tex) measure <- paste(measure, ", $\\%$ of reference", sep = "") else measure <- paste(measure, ", % of reference", sep = "")}
    
  } else if (how == "pct.inc") {
    Pct <- function(n, d) (n - d) / d * 100
    if (ws_sample) {
      df$ws_ML     <- Pct(df$ws_ML, df$DGM)
      df$ws_Lasso  <- Pct(df$ws_Lasso, df$DGM)
      df$ws_Ridge  <- Pct(df$ws_Ridge, df$DGM)
    }
    df$oos_ML    <- Pct(df$oos_ML, df$DGM)
    df$oos_Lasso <- Pct(df$oos_Lasso, df$DGM)
    df$oos_Ridge <- Pct(df$oos_Ridge, df$DGM)
    df$DGM       <- Pct(df$DGM, df$DGM)
    
    if (!measure == "") {if (tex) measure <- paste(measure, ", $\\%$ difference of reference", sep = "") else 
    {if (measure == "Nagelkerke")
        measure <- expression(Nagelkerke ~ R^{2} ~ ", % difference of reference") else measure <- paste(measure, ", % difference of reference", sep = "")}
    }
  }
  
  return(list(df = df, measure = measure))
}


