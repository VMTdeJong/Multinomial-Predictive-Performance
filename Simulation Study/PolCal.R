# #######################################################################################################################

# This code is based on code published in an appendix of 
#   Van Hoorde K, Vergouwe Y, Timmerman D, Van Huffel S, Steyerberg EW, Van Calster B. Assessing calibration of
#   multinomial risk prediction models. Statistics in medicine. 2014;33(15):2585–2596.


# #######################################################################################################################
# Code edited by: Valentijn M.T. de Jong.
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




Polcal <- function(outcome, k, p, LP, r = 1, estimates = FALSE, dfr = 2, doplot = F, plotoverall = F, datapoints = F,
                   smoothing = F, smoothpar = 1, intercept = T, slope = T, test = FALSE){
  # NOTE: This function is written for a multinomial outcome with three categories.
  #       	  If there are more than three categories, the code can be easily adapted.
  #       	  Comments are added to state where and how the code should be changed if needed.
  
  ##################################################################################
  # outcome 	  column vector containing the outcome for every case, with values 1 to k (i.c. k=3)
  # k 		      number of outcome categories (i.c. 3)
  # p 		      matrix with the probabilities of the prediction model, ordered from prob(cat. 1) to 
  # 		        prob(cat. k)
  # LP		      matrix with all the linear predictors with respect to the chosen reference category, 
  #		          ordered (e.g. LP2vs1 and LP3vs1)
  # r 		      reference category (default: category 1)
  # estimates   indicates whether the coefficients of the parametric recalibration framework are desired 
  #		          (default=FALSE)
  # dfr   		  degrees of freedom for the non-parametric calibration (default=2)
  # plotoverall indicates whether overall (non-)parametric calibration plots are constructed 
  #		          (default=TRUE)
  # datapoints	indicates whether the individual datapoints are shown in the overall (non-)parametric 
  #		          calibration plots (default = TRUE)
  # smoothing 	indicates whether a smoothed line (using cubic splines) is added to the calibration plots 
  #		          (default=TRUE)
  # smoothpar 	smoothing parameter for the smoothed line (default=1)
  # intercept 	indicates whether calibration intercepts are desired (default=FALSE)
  # slope 		  indicates whether calibration slopes are desired (default=FALSE)
  # test 		    indicates whether statistical tests for calibration are desired (default=FALSE)
  ##################################################################################
  
  # for this function you need to use library (VGAM)
  library(VGAM)
  
  # checks
  if (k != length(table(outcome))) {stop('The number of categories in outcome does not equal the specified number of categories.')}
  if (dim(p)[2] != k) {stop('The number of columns in p does not equal the specified number of categories.')}
  if (dim(LP)[2] != k - 1) {stop('The number of columns in LP does not equal the specified number of categories minus 1.')}
  if (!r %in% 1:k) {stop('The given reference category (r) is not possible.')}     
  if (!is.matrix(p)) {stop('p is not a matrix.')}
  if (!is.matrix(LP)) {stop('LP is not a matrix.')}
  if (isTRUE(plotoverall) && !isTRUE(datapoints) && !isTRUE(smoothing)) {stop('For the overall (non-)parametric calibration plots either 
                                                                            datapoints or smoothed lines should be requested.')}
  
  # if tests for perfect calibration are requested, automatically calibration intercepts and calibration slopes 
  # are given
  if (isTRUE(test)) {intercept <- slope <- TRUE}
  
  # probabilities
  probs <- split(p,col(p))    
  
  # linear predictors necessary for non-parametric calibration plot - give a name to each linear predictor 
  # seperately
  lps <- split(LP,col(LP))
  for (i in 1:(k - 1)) {assign(paste("lp", i, sep = ""), unlist(lps[[i]]))}
  
  ########################################
  # estimation of calibration intercepts 
  # cf. section 2.2.3. and 2.2.4.        
  ########################################
  
  if (isTRUE(intercept)) {int <- vgam(outcome ~ 1, offset = LP, family = multinomial(refLevel = r))
  coeffint <- coefficients(int)
  se <-  sqrt(diag(vcov(int)))
  ci1i <- cbind(LL1 = coeffint[1] - qnorm(0.975) * se[1], UL1 = coeffint[1] + qnorm(0.975) * se[1])
  ci2i <- cbind(LL2 = coeffint[2] - qnorm(0.975) * se[2], UL2 = coeffint[2] + qnorm(0.975) * se[2])
  estint <- c(coeffint[1],ci1i,coeffint[2],ci2i)
  names(estint) <- paste('CALINT',c('int1','LLint1','ULint1','int2','LLint2','ULint2'), sep = '.')}
  
  
  ####################################
  # estimation of calibration slopes 
  # cf. section 2.2.3. and 2.2.4.    
  ####################################
  
  if (isTRUE(slope)) {
    i <- diag(k - 1)
    i2 <- rbind(1, 0)
    i3 <- rbind(0, 1)
    clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3)
    slopes <- vgam(outcome ~ lp1 + lp2, family = multinomial(refLevel = r), constraints = clist)
    
    coeffslopes <- coefficients(slopes)[k:length(coefficients(slopes))]
    se <- sqrt(diag(vcov(slopes)))
    ci1s <- cbind(LL1 = coeffslopes[1] - qnorm(0.975) * se[3], UL1 = coeffslopes[1] + qnorm(0.975) * se[3])
    ci2s <- cbind(LL2 = coeffslopes[2] - qnorm(0.975) * se[4], UL2 = coeffslopes[2] + qnorm(0.975) * se[4])
    estslopes <- c(coeffslopes[1],ci1s,coeffslopes[2],ci2s)
    names(estslopes) <- paste('CALSLOPES', c('lp1','LLlp1','ULlp1','lp2','LLlp2','ULlp2'), sep = '.')
  }
  
  # Printing of results
  # The probabilities of calibration intercepts and slopes are only shown when the hypothesis of perfect 
  # calibration is rejected.
  
  results <- list(if (isTRUE(estimates)) {est} else{'Not requested'}, if (isTRUE(intercept)) {estint} else{'Not requested'}, if(isTRUE(slope)) {estslopes} else{'Not requested'}, if (isTRUE(test)) {c(deviancewithout, devint, devslopes)} else{'Not requested'}, if (isTRUE(test)) {c(poverall, if (poverall < 0.05) { c(pint, pslopes)})} else{'Not requested'})
  names(results) <- c("Coefficients of parametric recalibration framework","Calibration Intercepts with 95% CI","Calibration Slopes with 95% CI","Deviances","P-values")
  n <- 1:5
  selection <- c(isTRUE(estimates),isTRUE(intercept),isTRUE(slope),isTRUE(test),isTRUE(test))
  results[n[selection]]
  
  }


# Wrapper function that applies Polcal for each category as reference once
# Note: ref is the index of the column, it is not the actual outcome.
PolCalAll <- function(outcome, p.mat, ref = 1:ncol(p.mat), i = 0, id = 0, method = "NA") { # The i, id and method are not necessary. THey are only there to register for which sample in the simulation the code failed to compute a slope (This appeared to be the case mostly when the lasso removed all predictors from the model. The cal slope would then be inf or -inf.) 
  p.mat <- as.matrix(p.mat)
  if (!is.numeric(ref)) stop("ref should be the indexes (i.e. integers) of the column(s) chosen as reference.")
  ref <- sort(ref)
  K <- ncol(p.mat)
  
  if (K != 3) stop("Predicted probability matrix p.mat must have 3 columns, see Polcal() .")
  
  cal_ints   <- data.frame(matrix(nrow = length(ref), ncol = 6))
  cal_slopes <- data.frame(matrix(nrow = length(ref), ncol = 6))
  
  row.names(cal_ints)   <- ref
  row.names(cal_slopes) <- ref
  ref_num <- 1
  for (r in ref) {
    
  # Compute LP
  lp <- matrix(nrow = nrow(p.mat), ncol = K - 1)
  cat <- 1
  for (k in 1:K) {
    if (k != r) {
      lp[ , cat] <- log(p.mat[ , k] /  p.mat[ , r])
      cat <- cat + 1
    }
  }
  
#  # Some code for if you'd like to return a square matrix, with missing diagonal entries.
#   lp <- matrix(nrow = nrow(p.mat), ncol = ncol(p.mat))
#   
#   for (k in 1:K) {
#     if (k != r) {
#       lp[ , k] <- log(p.mat[ , k] / p.mat[ , r])
#     }
#   }
  
  # Error handler
  PolcalEH <- function(e) {
    print("Error caught in Polcal.")
    l <- list(e = e, outcome = outcome, K = K, p.mat = p.mat, lp = lp, r = r, id = id, i = i, method = method)
    save(l, file = 
           paste(paste(paste(paste(paste(getwd(), paste(method, "error_Polcal_id", sep = "_"), sep = "/"), id, sep = ""), "it", sep = ""), i, sep = ""), ".RData", sep = "")
         )
    return(NA) 
  }
  
  # Compute calibraton with each category as reference once.
  if (!is.na( tryCatch(cal <- 
             Polcal(outcome = outcome, k = K, p = p.mat, LP = lp, r = r, doplot = F, plotoverall = F)
             , error = PolcalEH)
  )[1] ) {
    cal_ints[ref_num, ]   <- cal$`Calibration Intercepts with 95% CI`
    cal_slopes[ref_num, ] <- cal$`Calibration Slopes with 95% CI`
  }
  ref_num <- ref_num + 1
  }
  
  colnames(cal_ints)   <- PolCalIntNames(K)
  colnames(cal_slopes) <- PolCalSlopesNames(K)
  
  return(list("intercepts" = cal_ints, "slopes" = cal_slopes))
}

PolCalIntNames <- function(k) {
  c <- k - 1
  out <- list()
  for (i in 1:c) {
    out <- c(out, paste(paste('CalInt',c('alpha','LL','UL') , sep = ""), i, sep = "")  )    
  }
 return(unlist(out)) 
}
  
PolCalSlopesNames <- function(k) {
  c <- k - 1
  out <- list()
  for (i in 1:c) {
    out <- c(out, paste(paste('CalSlope',c('Beta','LL','UL') , sep = ""), i, sep = "")  )    
  }
  return(unlist(out)) 
}  

UnName <- function(mat) {
  return(matrix(unlist(mat), nrow = nrow(mat), ncol = ncol(mat)))
}