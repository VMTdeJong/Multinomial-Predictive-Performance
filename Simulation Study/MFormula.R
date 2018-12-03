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

## Create Mlogit formula, using logical vector
CreateMlogitFormulaLogical <- function(logic_vec) {
  if (logic_vec[1]) {mfor_str <- "X1 ~ 1"} else {mfor_str <- "X1 ~ 0"}
  
  if (any(logic_vec)) {
    
    mfor_str <- paste(mfor_str, " |")
    this_term <- " X@@"
    
    for (v in 2:length(logic_vec))
    {
      if (logic_vec[v])
      {
        this_term <- gsub("@@", v , this_term)
        mfor_str <- paste(mfor_str, this_term)
        this_term <- " + X@@"
      }
    }
  }  
  return(mFormula(as.formula(mfor_str)))
}

## Create mlogit formula of specified length, using integer (calls function above)
CreateMlogitFormulaInt <- function(int) {
  logic_vec <- rep(TRUE, int)
  
  return(CreateMlogitFormulaLogical(logic_vec))
}

# Create names
CreateNamesLogi <- function(logic_vec) {
  if (logic_vec[1]) {name <- c("intercept")} else {name <- character()}
  
  if (any(logic_vec)) {
    
    this_term <- "X@@"
    
    for (v in 2:length(logic_vec))
    {
      if (logic_vec[v])
      {
        this_term <- gsub("@@", v , this_term)
        name <- c(name, this_term)
        this_term <- "X@@"
      }
    }
  }  
  return(name)
}

CreateNames <- function(int) {
  logic_vec <- rep(TRUE, int)
  
  return(CreateNamesLogi(logic_vec))
}
### Testing:
# 
# 
# test_data_list       <- MakeRetroData(n = c(10, 10, 10), outcome = c(0, 1, 2), 
#                               coefs = rbind(c(-0, -1, 1, 1, -1), c(0, 1, -1, -1, 1)), 
#                               n_bi_var = 0, bi_prob = c(.5), Brier = T)
# 
# test_data            <- test_data_list$dat
# 
# 
# 
# require(mlogit)
# test_data            <- data.frame(test_data)
# test_data_long       <- mlogit.data(test_data, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
# 
# 
# mfor                 <- CreateMlogitFormulaInt(5, envi = environment()) #### <---- gebruikt envi = environment() zodat hij altijd gevonden wordt
# mlogit_model         <- mlogit(formula = mfor, data = test_data_long)
# mlogit_model
# 
# ## Predict
# sample_mlogit_fitted <- predict(mlogit_model, newdata = test_data_long, type = "probs")
# 
# test_data_pop_list   <- MakeRetroData(n = c(100, 100, 100), outcome = c(0, 1, 2), 
#                                  coefs = rbind(c(-0, -1, 1, 1, -1), c(0, 1, -1, -1, 1)), 
#                                  n_bi_var = 0, bi_prob = c(.5), Brier = T)
# test_data_pop        <- data.frame(test_data_pop_list$dat)
# 
# test_data_pop_long   <- mlogit.data(test_data_pop, shape = "wide", choice = "X1", chid.var = "individual", alt.var = "outcome")
# 
# pop_mlogit_preds     <- predict(mlogit_model, newdata = test_data_pop_long,    type = "probs")
# 
# 
# 
# ## Notes for subset selection
# selection_crit <- cbind(c("backwards", 0.5), c("forward", 0.05))
# 
# if (length(selection_crit) > 0) { if (any(selection_crit == "backwards")) {print("backwards is true")}}
# 
# ##
# mFormula(X1 ~ 1 | X2 + X3 + X4 + X5)
