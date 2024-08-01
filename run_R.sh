#!/bin/bash

#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH -t 30:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=liangxw27@gmail.com

module load r/4.4.0-yycctsj
Rscript ./random_code/sim_Tetris.R
