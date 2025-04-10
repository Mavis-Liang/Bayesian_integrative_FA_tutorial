#!/bin/bash

#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liangxw27@gmail.com

module load r/4.4.0-yycctsj
Rscript -e "source('./functions/gen_scenario2.R'); \
            source('./functions/run_Tetris.R'); \
            set.seed(1)
            true_data <- gen_scenario2(S=4, N_s=c(35, 35, 35, 35), P=15, J_s=c(1, 1, 1, 1), K=4)
            new_control <- function(nrun = 100, burn = 80, thin = 1,
                        nu = 3, 
                        a1 = 2.1, b1 = 1,
                        a2 = 3.1, b2 = 1,
                        apsi = 1, bpsi = 0.3,
                        alpha = 10, beta = 0.5)
                {
                    return(list(nrun = nrun, burn = burn, thin = thin,
                            nu = nu, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
                            apsi = apsi, bpsi = bpsi, alpha = alpha, beta = beta))
                }
            peakRAM::peakRAM({
                fit_Tetris_fixT <- run_Tetris(true_data, big_T=true_data$big_T)})"
            
