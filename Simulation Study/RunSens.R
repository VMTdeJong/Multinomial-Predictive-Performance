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

c_args        <- commandArgs(trailingOnly = T) # These are the arguments given through LINUX 

if (length(c_args) == 0) {
  id <- 3                                # id of the simulation. Defines EPV, relative prevalences and # predictors. Each combination of paramters has its own id.
  repetition <- 1                        # To increase the number of iterations, the simulation can be repeated. This is the number of the repetition (Starts at 0)
  max_i <- 2                             # Number of iterations.
  cor_id <- 0                            # id for correlation. cor = cor_id / 10
  bin    <- 1                            # Are all predictors binary?  0 = false, all else is true
  effect <- 1                            # Multiplies size of predictors.
  source("MakeData.R")
  source("MakeCorData.R")
  source("MFormula.R")
  source("foldid.R")
  source("net predict.R")
  source("M-index.R")
  source("pseudo R2.R")
  source("PolCal.R")
  source("PDI.R")
  source("Simulation.R")
  proj_loc <- getwd()
  
} else {
  id            <- as.numeric(c_args[1]) # id of the simulation. Defines all further parameters for the main simulation. Each combination of paramters has its own id.
  repetition    <- as.numeric(c_args[2]) # To increase the number of iterations, the simulation can be preated. This is the number of the repetition (Starts at 0)
  max_i         <- as.numeric(c_args[3]) # Number of iterations.
  cor_id        <- as.numeric(c_args[4]) # Only for sensitivity analysis non-zero. id for correlation. 0 = 0, 1 = 0.2, 2 = 0.5, ...
  bin           <- as.numeric(c_args[5]) # Only for sensitivity analysis non-zero. Are all predictors binary?  0 = false, all else is true
  effect        <- if (is.na(as.numeric(commandArgs(trailingOnly = T)[6]))) 1 else as.numeric(c_args[6]) 
  
  proj_loc <- "/hpc/shared/julius_bs/VMTdeJong/MultEPVSens" # Path to my directory, and this project.
  source(paste(proj_loc, "MakeData.R"   , sep = "/"))
  source(paste(proj_loc, "MakeCorData.R", sep = "/"))
  source(paste(proj_loc, "MFormula.R"   , sep = "/"))
  source(paste(proj_loc, "foldid.R"     , sep = "/"))
  source(paste(proj_loc, "net predict.R", sep = "/"))
  source(paste(proj_loc, "M-index.R"    , sep = "/"))
  source(paste(proj_loc, "pseudo R2.R"  , sep = "/"))
  source(paste(proj_loc, "PolCal.R"     , sep = "/"))
  source(paste(proj_loc, "PDI.R"        , sep = "/"))
  source(paste(proj_loc, "Simulation.R" , sep = "/"))
}

print(paste("command args: ", c_args))
print(paste("id = ", id))
print(paste("rep = ", repetition))
print(paste("max_i = ", max_i))
print(paste("cor = ", cor_id))
print(paste("bin = ", bin))



################################ Load all necessary packages immediately ####################
# Require returns TRUE if successful, or if already loaded.
packages_loaded <- require(parallel) && require(glmnet) && require(mlogit) && require(MASS) 
# packages_loaded <- packages_loaded && require(doParallel) && require(foreach)

################################# Draw seed. ################################################
if (cor_id + bin > 0) {
  seed <- 1260 + repetition + 20 * cor_id + 200 * bin + (1 - effect) * 2000
} else {
  seed <- id * repetition   
}

RNGkind("L'Ecuyer-CMRG")
set.seed(8387691)

RNGstreams <- matrix(nrow = max(seed, 2)  , ncol = length(.Random.seed))
RNGstreams[1, ] <- .Random.seed

for (i in 2:nrow(RNGstreams)) {
  RNGstreams[i, ] <- nextRNGStream(RNGstreams[i - 1,])
}

seed <- RNGstreams[seed, ]

prev_id <- if (id <= 21) 1 else if (id <= 42) 2 else 3

##### n 
##### = development sample size, n_vec = vector of development sample category sizes.
if (length(c_args) == 0) { source("fileParams.R") } else source(paste(proj_loc, "fileParams.R" , sep = "/"))

n_vec   <- n_mat[id, ]
n_vec

##### N
##### = Validation sample sizes.
N_vec   <- N_mat[id, ]
N_vec

##### Coefficients.
n_coefs <- n_coefs_vec[id]
n_coefs

## Intercept
if (bin > 0) {# When predictors are binary
  alpha_file <- paste(paste("OptimBin/alpha for bin", bin, sep = ""), ".RDATA", sep = "")
  load(file = paste(proj_loc, alpha_file, sep = "/"))
  alpha <- mean(alpha)
}
##Correlations among predictors
if (cor_id > 0) {
  alpha_file <- paste(paste("OptimCor/alpha for cor0", cor_id, sep = ""), ".RDATA", sep = "")
  load(file = paste(proj_loc, alpha_file, sep = "/" ))
  alpha <- mean(alpha)
  
} else {if (bin == 0) {
  load(file = paste(proj_loc, "alph_mean.txt"   , sep = "/"))  
  alpha <- alph_mean[alph_mean$prev_id == prev_id & alph_mean$n_coefs == n_coefs, 1]
}
}


if (cor_id < 10 && cor_id >= 0) {
  max_cor <- min_cor <- cor_id / 10 
} else {
  stop("Wrong cor_id")
}

coefs <- effect * MakeCoefs(alphas = alpha, n_coefs = n_coefs) # The second alpha is added inside the function. # default for effect is 1
coefs

##### Proportion(s) for binary variables
bi_prob <- 0.5 # Ignored when n_bi_var = 0

if (bin) { n_bi_var <- n_coefs} else {n_bi_var <- 0}
################################# Save all parameters in a log file ####################################
out_loc      <- paste(proj_loc, "Output"   , sep = "/")  # Path for saving
results_loc  <- paste(out_loc , "results"  , sep = "/")  # path for results
scenario_loc <- paste(out_loc , "scenarios", sep = "/")  # path for scenarios

### Parameters of scenario
scenario <- list(max_i = max_i, n = n_vec, N = N_vec, coefs = coefs, n_bi_var = n_bi_var, bi_prob = bi_prob, id = id, seed = seed)

### Name
# Full name, to be filled in.

scenario_name <- "prev@@@_coefs###_epv!!!_id%%%_repetition&&&_cor@#!_bin#!%" 
if (!effect) scenario_name <- paste(scenario_name, "_eff0", sep = "" )
scenario_name <- paste(scenario_name, ".txt", sep = "")
# Prevalence ratio
prev_ratio <- as.character(10*n_vec[3] / n_vec[1]) 
scenario_name <- sub("@@@", prev_ratio, scenario_name)

# Number of coefs
#n_coefs <- ncol(coefs) - 1 # old
scenario_name <- sub("###", as.character(n_coefs), scenario_name)

# EPV
epv <- min(n_vec) / (n_coefs * nrow(coefs))
scenario_name <- sub("!!!", as.character(epv), scenario_name)

# id and repetition
scenario_name <- sub("%%%", as.character(id), scenario_name)
scenario_name <- sub("&&&", as.character(repetition), scenario_name)

# correlation
scenario_name <- sub("@#!", as.character(cor_id), scenario_name)

# Predictors binary?
scenario_name <- sub("#!%", as.character(bin), scenario_name)


scenario_name_params     <- paste("param", scenario_name, sep = "_")              # File name
scenario_path_name_params <- paste(scenario_loc, scenario_name_params, sep = "/") #} # Path and file 
print(scenario_path_name_params)



### And finally save it. (This is done before running so that the parameters are saved even when an unexpected error is encountered.)
save(scenario, file = scenario_path_name_params)

#################################         Run the scenario!        ####################################

if (packages_loaded) {
  results <- Simulation(max_i = max_i 
                        ,n = n_vec
                        ,N = N_vec 
                        ,coefs = coefs
                        ,n_bi_var = n_bi_var 
                        ,bi_prob = bi_prob
                        ,id = id                 # only a reference number
                        ,repetition = repetition # Again only reference.
                        ,seed = seed
                        
                        ,nfolds = 10
                        ,n_cv = 1    
                        ,cores = 1               # ignored, as < 3 ( = 1, as only 1 core per simulation should be used on hpc)
                        ,computePDI = T          # default = T
                        ,compute_cal = T         # default = T
                        ,n_lambda_Lasso = 100    # Seems like the default option, but actually forces glmnet to run all 100 (default) lambdas.
                        ,n_lambda_Ridge = 200    # Adds more lambdas for Ridge, so that the optimal can be chosen
                        ,proj_loc = proj_loc
                        ,min.cor = min_cor
                        ,max.cor = max_cor)
} else {
  results <- "packages could not be loaded"
}

################################# Save all results in a file!   ####################################

scenario_name_results      <- paste("results", scenario_name, sep = "_")           # File name
#if (length(c_args) == 0) { 
#  scenario_path_name_results <- scenario_name_results 
#} else {
scenario_path_name_results <- paste(results_loc, scenario_name_results, sep = "/") # } # Path and file combined

print(scenario_path_name_results)

save(results, file = scenario_path_name_results)
