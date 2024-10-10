source("./FBPFA-PFA.R")
source("./FBPFA-PFA with fixed latent dim.R")
#source("./Tetris.R")
library(BFR.BE)
#source("./functions/run_Tetris.R")
source("./functions/post_PFA.R")
source("./functions/post_BMSFA.R")
source("./functions/post_stackFA.R")
source("./functions/post_SUFA.R")
source("./functions/post_MOMSS.R")
#source("./functions/post_Tetris.R")
source("./functions/post_IndFA.R")
source("./functions/measurements.R")

# --------------------------- Scenario MOMSS (scenario 1)---------------------------
# This function produce the result of a single simulation of the scenario 1
# includes data generation, fitting the 7 models, post processing for the point estimates, and FN and RV metrics
# This function is designed for use in parallel computing. 
source("./functions/gen_scenario1.R")
sim_scenario1 <- function(i){
  set.seed(i)
  # Data generation
  true_data <- gen_scenario1(S=4, N_s=c(100, 100, 100, 100), P=40, K=4, Q=2)
  
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
  # Stack FA
  profile_stackFA <- peakRAM::peakRAM({
    fit_stackFA <- MSFA::sp_fa(true_data$Y_mat, k = 4,
                               control = list(nrun = 10000, burn = 8000))
  })
  # Ind FA
  profile_indFA <- peakRAM::peakRAM({
    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        j_s = c(4, 4, 4, 4)
        MSFA::sp_fa(true_data$Y_list[[s]], k = j_s[s], trace = FALSE,
                    control = list(nrun = 10000, burn = 8000))
      })
  })
  # MOMSS
  profile_MOMSS <- peakRAM::peakRAM({
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = true_data$X, 
                                      b = true_data$A, q = 4, scaling = FALSE)
  })
  # BMSFA
  profile_BMSFA <- peakRAM::peakRAM({
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 4, j_s = c(1,1,1,1),
                               outputlevel = 1, scaling = FALSE, trace = FALSE,
                               control = list(nrun = 10000, burn = 8000))
  })
  # PFA
  profile_PFA <- peakRAM::peakRAM({
    fit_PFA <- PFA(Y=t(true_data$Y_mat), 
                   latentdim = 4,
                   grpind = rep(1:length(true_data$N_s), 
                                times = true_data$N_s))
  })
  # SUFA
  profile_SUFA <- peakRAM::peakRAM({
    fit_SUFA <- SUFA::fit_SUFA(true_data$Y_list, qmax=4,
                               nrun = 10000)
  })
  # Tetris
  # Don't run this code block if you don't have the time to wait
  # profile_Tetris <- peakRAM::peakRAM({
  #   fit_Tetris <- run_Tetris(true_data)
  # })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(stackFA = post_stackFA(fit_stackFA, S=4), 
                       IndFA = post_IndFA(fit_indFA),
                       MOMSS = post_MOMSS(fit_MOMSS),
                       BMSFA = post_BMSFA(fit_BMSFA),
                       PFA = post_PFA(fit_PFA),
                       SUFA = post_SUFA(fit_SUFA)#,
                       #Tetris = post_Tetris(fit_tetris)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(profile_stackFA[,c(2, 4)], profile_indFA[,c(2, 4)], 
                      profile_MOMSS[,c(2, 4)],
                      profile_BMSFA[,c(2, 4)], profile_PFA[,c(2, 4)], 
                      profile_SUFA[,c(2, 4)]#,
                      #profile_Tetris[,c(2, 4)]
                      ) %>% 
    t() %>% `colnames<-`(c("stackFA", "IndFA", "MOMSS", "BMSFA",
                           "PFA", "SUFA"#, "Tetris"
                           ))
  
  # Results
  final_list <- list(true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

# Usage
# result <- sim_scenario2(1)
# print(result$metrics)
test <- sim_scenario1(1)
