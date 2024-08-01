#!/bin/bash

#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 30:00:00

run_R.sh

module load r/4.4.0-yycctsj
Rscript ./main_code/multicore.R 
