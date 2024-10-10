source("./FBPFA-PFA.R")
source("./FBPFA-PFA with fixed latent dim.R")
source("./Tetris.R")
library(BFR.BE)
source("./functions/run_Tetris.R")
source("./functions/post_PFA.R")
source("./functions/post_BMSFA.R")
source("./functions/post_stackFA.R")
source("./functions/post_SUFA.R")
source("./functions/post_MOMSS.R")
source("./functions/post_Tetris.R")
source("./functions/post_IndFA.R")
source("./functions/measurements.R")
# --------------------------- Scenario PFA (scenario 3)---------------------------
# This function produce the result of a single simulation of the scenario 3
# includes data generation, fitting the 7 models, post processing for the point estimates, and FN and RV metrics
# This function is designed for use in parallel computing. 
source("./functions/gen_scenario3.R")
sim_scenario3 <- function(seed=1){
  set.seed(seed)
  # # Data generation 
  true_data <- gen_scenario3(4, N_s=c(100, 100, 100, 100), P=40, K=4)
  # Fitting
  # Stack FA
  profile_stackFA <- peakRAM::peakRAM({
    fit_stackFA <- MSFA::sp_fa(true_data$Y_mat, k = 4,
                               control = list(nrun = 10000, burn = 8000))
  })
  
  # Ind FA
  profile_indFA <- peakRAM::peakRAM({
    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        MSFA::sp_fa(true_data$Y_list[[s]], k = 4,
                    control = list(nrun = 10000, burn = 8000))
      })
  })
  
  # MOMSS
  profile_MOMSS <- peakRAM::peakRAM({
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = NULL, 
                                      b = true_data$A, q = 4, scaling = FALSE)
  })
  
  # BMSFA
  profile_BMSFA <- peakRAM::peakRAM({
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 4, j_s = c(1,1,1,1), # j_s can't be 0
                               outputlevel = 1, scaling = FALSE,
                               control = list(nrun = 10000, burn = 8000))
  })
  
  # PFA
  profile_PFA <- peakRAM::peakRAM({
    fit_PFA <- PFA(Y=t(true_data$Y_mat), latentdim = 4, grpind = true_data$grpind,
                   Thin = 5, Total_itr = 10000, burn = 8000)
  })
  
  # SUFA
  profile_SUFA <- peakRAM::peakRAM({
    fit_SUFA <- SUFA::fit_SUFA(true_data$Y_list, qmax=4,nrun = 10000)
  })
  # SUFA now also support specifying the number of study-specific factors. J_s = c(0, 1, 1, 1). J_s can be 0.
  
  # Tetris is runned on separate process, since it is slow and not parallelizable for large cores
  # Tetris
  #  profile_Tetris <- peakRAM::peakRAM({
  #   fit_tetris <- run_Tetris(true_data)
  # })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(stackFA = post_stackFA(fit_stackFA, S=4), 
                       IndFA = post_IndFA(fit_indFA),
                       MOMSS = post_MOMSS(fit_MOMSS),
                       BMSFA = post_BMSFA(fit_BMSFA),
                       PFA = post_PFA(fit_PFA),
                       SUFA = post_SUFA(fit_SUFA))
  # Time and peak RAM profile
  profile_df <- rbind(profile_stackFA[,c(2, 4)], profile_indFA[,c(2, 4)], 
                      profile_MOMSS[,c(2, 4)],
                      profile_BMSFA[,c(2, 4)], profile_PFA[,c(2, 4)], 
                      profile_SUFA[,c(2, 4)]
  ) %>% 
    t() %>% `colnames<-`(c("stackFA", "IndFA", "MOMSS", "BMSFA",
                           "PFA", "SUFA"))
  
  # Results
  final_list <- list(seed = seed,
                     true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

