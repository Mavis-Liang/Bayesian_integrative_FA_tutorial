#!/bin/bash

#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH -t 30:00:00

module load r/4.2.2-z6qdiis
Rscript ./main_code/multicore.R
