#$ -S /bin/bash
#$ -cwd
#$ -M V.M.T.deJong-2@umcutrecht.nl
#$ -m beas
module load R
Rscript  /hpc/shared/julius_bs/VMTdeJong/MultEPVSens/RunSens.R $SGE_TASK_ID 0 200 0 0