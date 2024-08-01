source("Tetris.R")
source("functions/gen_senerioBMSFA.R")
packages <- c('MSFA', 'peakRAM', 'BFR.BE', 'tidyverse', 'matlab', 'MatrixCorrelation', "SUFA")
lapply(packages, library, character.only = TRUE)

#run_fixed <- tetris(sim_data_test$Y_list,alpha=5,beta=1,fixed=TRUE,A_fixed=A)
#Lambda <- getLambda(run_fixed,A)

sim_data_small <- gen_senerioBMSFA(S=3, N=100, P=15, K=2, c(1,1,1,1), "dense")
run_small <- tetris(sim_data_small$Y_list,alpha=5,beta=1)
A_small <- choose.A(run_small,alpha_IBP=5,S=3)
save(A_small, file = "./RDS/A_sim_small.rds")

sim_data_test1 <- readRDS("./RDS/sim_data_test2B_8_1.rds")
run <- tetris(sim_data_test1$Y_list,alpha=5,beta=1)
A <- choose.A(run,alpha_IBP=5,S=4)
save(A, file = "./RDS/A_sim_2.rds")