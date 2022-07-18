## Predictive performance of multinomial logistic prediction models
* Author of this file and of the final version of the code: Valentijn M.T. de Jong
* File last modified: July 2022
* Code last modified: < April 2018
* This is code for a simulation study presented in: de Jong VMT, Eijkemans MJC, van Calster B, et al. Sample size considerations and predictive performance of multinomial logistic prediction models.Statistics in Medicine. 2019;1–19.https://doi.org/10.1002/sim.8063

* For questions, email: 
    * V.M.T.deJong-2@umcutrecht.nl
    * Valentijn.M.T.de.Jong@gmail.com

## How to...
1. [Get started](#start)
2. [Run the simulation in R](#inR)
3. [Run the simulation from command line](#cmd)
4. [Reproduce the data](#rep-data)
5. [Load the data.](#load)
6. [Reproduce the plots in the manuscript](#plots)
7. [Rereproduce the tables in the manuscript](#tables)
8. [Inspect non-convergence](#non-conv)

#### 1. Get started <a name="start"></a>
    a. All simulations were carried out in R 3.2.2
    b. The .R files contain elaborate comments on how to use each function. Therefore, this file is only there to direct users to the appropriate .R files in the correct order, and which parameters to use.
    c. Note that many of the .R files require other files to be present in the same directory, or in the directory they were put.
    d. The used R packages are: install.packages(c("glmnet", "mlogit", "MASS", "abind", "tidyr", "VGAM", "DBI", "lazyeval")) # The last two are not necessary for some versions of R. Note that this simulation requires a previous version of mlogit" (otherwise an error with dfidx::dfidx will show), which can be downloaded as follows:
	      install.packages("devtools")
	      require(devtools)
	      install_version("mlogit", version = "0.4.2", repos = "http://cran.us.r-project.org")
    e. The full data files are very large, and therefore not included in this data archive. Instead, .RDATA files are incuded, which include all summary data in the manuscript, and more. The estimated alpha parameters are located in .txt files, which are automatically loaded by the scripts.

#### 2. Run the simualation in R <a name="inR"></a>
    a. Make sure all of the project's .R files and alph_mean.txt are in the correct directory
    b. Open RunSens.R
    c. Change the variables to whatever you desire:
	      * id         # This sets all parameters for main sim (possible values: 1:63)
	      * repetition # This adjusts the seed                 (possible values: all >= 0)
 	      * max_i      # This sets the number of iterations.   (possible values >= 2)
	      * cor_id	   # correlation = cor_id/10
	      * bin        # 1 = binary, 0 = normal predictors
	      * effect     # multiplies coefficients. 1 = as specified in manuscript. 
	      * upper.folder # name of the folder the results are stored in. Change only for sensitivity analyses.
    d. Run it.

#### 3. Run the simulation from command line in linux <a name="cmd"></a>
    a. Go through 1, 2a and 2b
    b. Call RunSens.R . Supply values for the parameters id, repetition and max_i, cor_id, bin and effect e.g.:	/your/path/RunSens.R 3 $SGE_TASK_ID 200 2 0
    c. Alternative for b.: Use the sh files in the sh files folder.

#### 4. Reproduce the data. <a name="rep-data"></a>
    a. See 1, as well as 2.
    b. Make sure the file alph_mean.txt is in the same directory as the R files, and OptimCor and OptimBin contain their files with intercepts to skip this step. Or, if you want to reproduce the estimated intecepts. This can be done by running: optimFinal.R, OptimizeAlphaForCor and OptimizeBin, with the following parameters:
  	    * id = 1 		# Increment with 1 when more iterations are to be added
  		  * maxit = 1000	# Sets the number of iterations the algorithm is allowed to run
  		  * reps = 150 	# Number of repetitions of each scenario (increase to increase accuracy)
    c. Run all scenarios 2000 times. A scenario is selected by setting id to an appropriate value (1:63). See RunSens.R for how exactly this selects a scenario. One might want to run all 2000 in one go. Then set max_i to 2000, and leaverepetition at 0. If one wants to split up each scenario in batches, increment repetition with 1 for every time a scenario is restarted. This ensures a new seed is selected, i.e. independent of all previous seeds (thanks to the parallel package).
    d. For exact replication, run every main scenario (id = 1:63) with max_i = 200, and repetition 0 till 9 and every sensitivity analysis repetion 1 till 10. This can be achieved by using the sh files in linux (see step 3), by running the files in the folder "sh files", after editing them so that they conform to your own setup.

#### 5. Load the data. <a name="load"></a>
    a. See 1.
    b. If you want to load the raw data, you'll have to reproduce them first, see 2 and 4. The Aggregate data is provided with this archive and can easily be loaded in c. and d.
    c. Place all files in the correct directory (if this archive is unchanged, this should be automatic).
    d. Run "LoadData.R". Note that: 
	      * Data from each scenario must be present. 
	      * Loading of raw (instead of aggregated) data may take several minutes. Loading of aggregate data should take a few seconds max.
	      * 4 warnings of missing files may appear. Sorry. They are not really missing, but misnamed.

#### 6. Reproduce the plots in the manuscript <a name="plots"></a>
    a. See 1, 2 and 4. for generating data.
    b. See 1 and 5 for loading of data.
    c. Run "Plots 2nd sub.R"

#### 7. Reproduce the tables in the manuscript. <a name="tables"></a>
    a. See 1, 2, 4. for generating data.
    b. See 1 and 5 for loading of data.
    c. Run "Manuscript Tables.R"

#### 8. Inspect non-convergence <a name="non-conv"></a>
    a. See 1, 2, 4. for generating data.
    b. See 1 and 5 for loading of data.
    c. The proportion of non-convergence can be found in "Non-convergence.R"

#### The end
