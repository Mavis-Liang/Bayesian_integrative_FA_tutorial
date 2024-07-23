library(SUFA)
library(mvtnorm)
library(tidyverse)
library(MASS)
source("./functions/gen_senerioSS.R")
burnin=500
sim_data_test <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
res_ss<-fit_SUFA(sim_data_test$Y_list,qmax=5,nthreads = 5,nrun = 7.5e3,
              nleapfrog = 4, leapmax = 9, del_range = c(0,.01))
Phi_est <- lam.est(res_ss$Lambda,burn = burnin)
Phi_true <- sim_data_test$Phi
MatrixCorrelation::RV(Phi_est, Phi_true)

