#!/bin/bash
#SBATCH -J IndBMSFA     # Job name
#SBATCH --array=1
#SBATCH -N 1          
#SBATCH -n 1         # Tasks per node that run concurrently
#SBATCH -c 1          
#SBATCH -t 10:00:00   
#SBATCH --mem=200G       
#SBATCH -e IndBMSFA-%a.err           
#SBATCH -o IndBMSFA-%a.out           

module load r/4.4.0-yycctsj          # Load the R module

echo "Starting job $SLURM_ARRAY_sbaTASK_ID on $HOSTNAME with seed $SEED"

# Call your R script and pass the seed to the function
# Rscript -e "source('./main_code/sim_scenario3.R'); \
#              result <- sim_scenario3(($SLURM_ARRAY_TASK_ID) * 3); \
#              saveRDS(result, paste0('./RDS/sc3/sc3_', $SLURM_ARRAY_TASK_ID, '.rds'))"

Rscript -e "source('./main_code/run_curatedData.R')"
