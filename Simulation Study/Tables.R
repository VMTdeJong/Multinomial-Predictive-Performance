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

source("Coerce.R")

EstSETable <- function(est, CE, panel.by = NA, comb.at = "w ML", ref = "Ref.", digits = 2,
                       med.val = .05, med.icon = "'", large.val = .10, large.icon = "*"){
  v    <- which(names(est) == comb.at):ncol(est) 
  pars <- est[ , -v]
  est  <- est[ , v]
  CE   <- CE[ , v]
  app  <- F
  
  if (!is.na(ref)) 
    if (any(names(est) == ref)) 
    {
      ref_col <- est[ref]
      est <- est[names(est) != ref]
      app <- T
    }
  
  print(ToOneString("Minimum SE: ", format(round(min(as.numeric(as.matrix(CE))), digits = digits), digits = digits), 
                    ". Maximum SE: ", format(round(max(as.numeric(as.matrix(CE))), digits = digits), digits = digits),"."))
  
  CE[CE < med.val]    <- " "
  CE[CE >= large.val] <- large.icon
  CE[CE >= med.val & CE < large.val]   <- med.icon

  comb_values <- CombIndStrings(format(est, digits = digits), CE) # , sep = " (", append.by = ")")
  comb_values <- cbind(pars, comb_values)
  if (app) comb_values <- cbind(comb_values, ref_col)
  return(comb_values)
}

SelectCal <- function(data, sel.meth = c("ML", "Lasso", "Ridge"), sel.col = NA) {
  out <- list()
  
  for (m in 1:length(sel.meth))
  {
    slopes <- data[ , grepl(sel.meth[m], names(data))]     # Selects correct method.
    if (!any(is.na(sel.col))) slopes <- slopes[ , sel.col] # And the right columns
    
    
    for (col in 1:ncol(slopes)) {
      selection <- names(slopes)[col]
      out[selection] <- slopes[selection]
    }
  }
  
  return(data.frame(out))
}

CalTable <- function(estimates, se, pars = c("prev", "n_coefs", "EPV1", "n_total"), sel.meth = c("ML", "Lasso", "Ridge"), 
                     sel.col = NA, digits = 3, 
                     large.val = .005, med.val = .0025, ...)
{
  par_data <- data.frame(matrix(nrow = nrow(estimates), ncol = 0))
  par_data[pars] <- estimates[pars]
  
  est_data <- round(SelectCal(data = estimates, sel.meth = sel.meth, sel.col = sel.col), digits = digits) 
  se_data  <- round(SelectCal(data = se       , sel.meth = sel.meth, sel.col = sel.col), digits = digits) 
  
  comb_data <- EstSETable(est_data, se_data, large.val = large.val, med.val = med.val, comb.at = "ML.ref.2.other.0", digits = digits)
  
  cbind(par_data, comb_data)
}

Convert <- function(datalist, n, FUN = sd, ref = "DGM", meth = "pct.bias", na.rm = T)
{
  ref_values      <- datalist[[ref]]
  mean_ref_value  <- mean(ref_values, na.rm = na.rm)
  data            <- datalist[names(datalist)[names(datalist) != ref]]
  
  if (meth == "pct.bias" || meth == "prop.bias" || meth == "bias") data <- lapply(data, `-`, mean_ref_value)
  if (meth == "pct.bias" || meth == "prop.bias")                   data <- lapply(data, `/`, mean_ref_value)
  if (meth == "pct.bias")                                          data <- lapply(data, `*`, 100)
  data[[ref]] <- ref_values
  data.frame(data)
}


TableNames <- function(df, digits = 2, panel.by = NA, scientific = F){
  tempprev <- df$prev
  df <- format(round(df, digits = digits), scientific = scientific)

  df_table <- data.frame("Freq." = rep(0, length(df$EPV1)), "Predictors" = df$n_coefs, "EPV" = df$EPV1, 
                         "Total sample size" = df$n_total)
  df_table$"Freq."[tempprev == 1]     <- c("1/3, 1/3, 1/3.")    # c("$33\\%$") # c("$33\\%, 33\\%, 33\\%$")
  df_table$"Freq."[tempprev == 4.500] <- c("2/20, 9/20, 9/20.") # c("$45\\%$") # c("$10\\%, 45\\%, 45\\%$")
  df_table$"Freq."[tempprev == 0.125] <- c("8/10, 1/10, 1/10.") # c("$80\\%$") # c("$80\\%, 10\\%, 10\\%$")

  
  if (any(names(df) == "oos_ML")) {   # For all tables, except calibration
    df_table$"w ML"     <- df$ws_ML     
    df_table$"w lasso"  <- df$ws_Lasso  
    df_table$"w ridge"  <- df$ws_Ridge  
    
    df_table$"o ML"     <- df$oos_ML    
    df_table$"o lasso"  <- df$oos_Lasso 
    df_table$"o ridge"  <- df$oos_Ridge 
    df_table$"Ref."     <- df$DGM
    
  } else {# For calibration tables
    df_table$"ML"       <- df$ML
    df_table$"lasso"    <- df$Lasso
    df_table$"ridge"    <- df$Ridge
    
    df_table$"NA ML"       <- df$NA_ML
    df_table$"NA lasso"    <- df$NA_Lasso
    df_table$"NA ridge"    <- df$NA_Ridge
  }
  
return(df_table)
}

PanelBy <- function(df, panel.by = NA)
{
  if (!is.na(panel.by)) {
    labs <- unlist(unique(subset(df, select = panel.by)))
    panel_list <- list()
    
    for (lab in labs)
    {
      panel_list[[length(panel_list) + 1]] <- df[subset(df, select = panel.by) == lab, ] 
    }
    return(Reduce(cbind, panel_list))
  }
  return(df)
}

FinalTable <- function(df1, df2) PanelBy(EstSETable(TableNames(df1), TableNames(df2)), panel.by = "Freq.")
