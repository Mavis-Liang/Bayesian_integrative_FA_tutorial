source("./FBPFA-PFA.R")
source("./FBPFA-PFA with fixed latent dim.R")
source("./Tetris.R")
source("./functions/run_Tetris.R")
source("./functions/post_PFA.R")
source("./functions/post_BMSFA.R")
source("./functions/post_stackFA.R")
source("./functions/post_SUFA.R")
source("./functions/post_MOMSS.R")
source("./functions/post_Tetris.R")
source("./functions/post_IndFA.R")
source("./functions/measurements.R")
# --------------------------- Scenario BMSFA ---------------------------
source("./functions/gen_senerioBMSFA.R")
set.seed(6)
sim_B_mid <- gen_senerioBMSFA(S=4, N=500, P=20, K=4, c(2,2,1,1), "dense")
saveRDS(sim_B_mid, "./RDS/sim_B_mid.rds")

# Stack FA
fit_stackFA_senB_mid <- MSFA::sp_fa(sim_B_mid$Y_mat, k = 6, trace = FALSE,
                              control = list(nrun = 5000, burn = 4000))
print("Stack FA done")
# Ind FA
fit_indFA_senB_mid <-
  lapply(1:length(sim_B_mid$n_s), function(s){
    j_s = c(6, 6, 5, 5)
    MSFA::sp_fa(sim_B_mid$Y_list[[s]], k = j_s[s], trace = FALSE,
                control = list(nrun = 5000, burn = 4000))
  })
print("Ind FA done")

# MOMSS
library(BFR.BE)
fit_MOMSS_senB_mid <- BFR.BE::BFR.BE.EM.CV(x = sim_B_mid$Y_mat, v = NULL, 
                                             b = sim_B_mid$A, q = 4, scaling = FALSE)
print("MOMSS done")

# BMSFA
fit_BMSFA_senB_mid <- MSFA::sp_msfa(sim_B_mid$Y_list, k = 4, j_s = c(2, 2,1,1),
                                outputlevel = 1, scaling = FALSE, trace = FALSE,
                                control = list(nrun = 5000, burn = 4000))
print("BMSFA done")
# PFA
fit_PFA_senB_mid <- PFA(Y=t(sim_B_mid$Y_mat), 
                        latentdim = 4,
                        grpind = rep(1:length(sim_B_mid$n_s), 
                                     times = sim_B_mid$n_s))
print("PFA done")
# SUFA
fit_SUFA_senB_mid <- SUFA::fit_SUFA(sim_B_mid$Y_list, qmax=4,nthreads = 5,
                                nrun = 7.5e3,
                                nleapfrog = 4, leapmax = 9, 
                                del_range = c(0,.01))
print("SUFA done")
# Tetris
# Don't run this code block if you don't have the time to wait
fit_Tetris_senB_mid <- run_Tetris(sim_B_mid)
saveRDS(fit_Tetris_senB_mid, "./RDS/fit_Tetris_senB_mid.rds")
fit_Tetris_senB_mid <- readRDS("./RDS/fit_Tetris_senB_mid.rds")
print("Tetris done")

# Results
postfit_scB_mid_list <- list(stackFA = post_stackFA(fit_stackFA_senB_mid, S=4), 
                               IndFA = post_IndFA(fit_indFA_senB_mid),
                               MOMSS = post_MOMSS(fit_MOMSS_senB_mid),
                               BMSFA = post_BMSFA(fit_BMSFA_senB_mid),
                               PFA = post_PFA(fit_PFA_senB_mid),
                               SUFA = post_SUFA(fit_SUFA_senB_mid),
                               Tetris = post_Tetris(fit_Tetris_senB_mid))
saveRDS(postfit_scB_mid_list, "./RDS/postfit_scB_mid_list.rds")
#metrics_table(postfit_scB_small_list, sim_B_small)

# Plot
#plot_loading(postfit_scB_small_list, sim_B_small)
#plot_cov(postfit_scB_small_list, sim_B_small)