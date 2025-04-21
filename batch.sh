#!/bin/bash
#SBATCH -J nutrition_PFA     # Job name
#SBATCH --array=1
#SBATCH -N 1          
#SBATCH -n 1         # Tasks per node that run concurrently
#SBATCH -c 1          
#SBATCH -t 5:00:00   
#SBATCH --mem=25G       
#SBATCH -e nutrition_PFA-%a.err           
#SBATCH -o nutrition_PFA-%a.out           

module load r/4.4.0-yycctsj          # Load the R module

echo "Starting job $SLURM_ARRAY_sbaTASK_ID on $HOSTNAME with seed $SEED"

# Call your R script and pass the seed to the function
Rscript -e "source('./main_code/sim_scenario4.R'); \
             result <- sim_scenario4_mis(($SLURM_ARRAY_TASK_ID) * 3); \
             saveRDS(result, paste0('./RDS/sc4_mis/sc4_mis_', $SLURM_ARRAY_TASK_ID, '.rds'))"

#Rscript -e "source('./main_code/run_nutrition.R')"
