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
# --------------------------- Scenario SUFA (scenario 6 [scenario 3 in the manuscript])---------------------------
# This function produce the result of a single simulation of the scenario 6
# includes data generation, fitting the 7 models, post processing for the point estimates, and FN and RV metrics
# This function is designed for use in parallel computing. 
# There is a smaller version of this scenario for a quick start. See sim_scBMSFA.R
source("./functions/gen_scenario6.R")
sim_scenario6 <- function(seed=1){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario6(S=4, N_s=c(100, 100, 100, 100), P=40, K=4, c(1,1,1,1))
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()
  
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
  #Stack FA
  profile_stackFA <- peakRAM::peakRAM({
    fit_stackFA <- MSFA::sp_fa(true_data$Y_mat, k = 4, scaling = FALSE, centering = TRUE, 
                                      control = list(nrun = 10000, burn = 8000))
  })
  # # Ind FA
  profile_indFA <- peakRAM::peakRAM({
    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        j_s = c(5, 5, 5, 5)
        MSFA::sp_fa(true_data$Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE, 
                    control = list(nrun = 10000, burn = 8000))
      })
  })

  # PFA
  profile_PFA <- peakRAM::peakRAM({
  fit_PFA <- PFA(Y=t(Y_mat_scaled), 
                         latentdim = 4,
                         grpind = rep(1:length(true_data$N_s), 
                                      times = true_data$N_s),
                Thin = 5,
                Total_itr = 10000, burn = 8000)
  })

  # MOMSS
  profile_MOMSS <- peakRAM::peakRAM({
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = NULL, 
                                             b = true_data$M, q = 4, scaling = FALSE)
  })
  
  
  # # SUFA (Now SUFA also support specifying the number of study-specific factors.)
  profile_SUFA_fixJs <- peakRAM::peakRAM({
    fit_SUFA_fixJs <- SUFA::fit_SUFA(Y_list_scaled, qmax=4, nrun = 10000, qs = c(1, 1, 1, 1))
  })
  # # SUFA 
  profile_SUFA <- peakRAM::peakRAM({
  fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=4, nrun = 10000)
  })

  # # BMSFA
  profile_BMSFA <- peakRAM::peakRAM({
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 4, j_s = c(1, 1, 1, 1),
                                      outputlevel = 1, scaling = FALSE, centering = TRUE, trace = FALSE,
                                      control = list(nrun = 10000, burn = 8000))
  })
  
  # # Tetris
   profile_Tetris_fixT <- peakRAM::peakRAM({
     fit_Tetris_fixT <- run_Tetris(true_data, big_T=true_data$big_T)
  })
  # Don't run this code block if you don't have the time to wait
  # profile_Tetris <- peakRAM::peakRAM({
  #   fit_Tetris <- run_Tetris(true_data)
  # })

  ######################################################################
  # Post processing for each methods
  postfit_list <- list(stackFA = post_stackFA(fit_stackFA, S=length(true_data$N_s)), 
                       IndFA = post_IndFA(fit_indFA),
                       MOMSS = post_MOMSS(fit_MOMSS),
                       BMSFA = post_BMSFA(fit_BMSFA),
                        PFA = post_PFA(fit_PFA),
                       SUFA = post_SUFA(fit_SUFA),
                       SUFA_fixJs = post_SUFA(fit_SUFA_fixJs),
                       Tetris = post_Tetris(profile_Tetris),
                       Tetris_fixT = post_Tetris(fit_Tetris_fixT)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(profile_stackFA[,c(2, 4)], profile_indFA[,c(2, 4)], 
                      profile_MOMSS[,c(2, 4)],
                      profile_BMSFA[,c(2, 4)], 
                       profile_PFA[,c(2, 4)], 
                      profile_SUFA[,c(2, 4)], profile_SUFA_fixJs[,c(2, 4)],
                      profile_Tetris[,c(2, 4)],
                      profile_Tetris_fixT[, c(2, 4)]
                      ) %>% 
    t() %>% `colnames<-`(c("stackFA", "IndFA", "MOMSS", "BMSFA",
                            "PFA",
                           "SUFA", "SUFA_fixJs", "Tetris_fixT"
                           )
                           )
  
  # Results
  final_list <- list(seed = seed,
                     true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

# Usage
# result <- sim_scenario6(1)
# print(result$metrics)


sim_scenario6_mis <- function(seed=1){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario6(S=4, N_s=c(100, 100, 100, 100), P=40, K=4, c(1,1,1,1))
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()
  
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
  # Stack FA
  profile_stackFA <- peakRAM::peakRAM({
    fit_stackFA <- MSFA::sp_fa(true_data$Y_mat, k = 6, scaling = FALSE, centering = TRUE, 
                                      control = list(nrun = 10000, burn = 8000))
  })
  # Ind FA
  profile_indFA <- peakRAM::peakRAM({
    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        j_s = c(6, 6, 6, 6)
        MSFA::sp_fa(true_data$Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE, 
                    control = list(nrun = 10000, burn = 8000))
      })
  })

  # PFA
  profile_PFA <- peakRAM::peakRAM({
  fit_PFA <- PFA(Y=t(Y_mat_scaled), 
                         latentdim = 6,
                         Cutoff = 0.001,
                         grpind = rep(1:length(true_data$N_s), 
                                      times = true_data$N_s),
                Thin = 5,
                Total_itr = 10000, burn = 8000)
  })

  # MOMSS
  profile_MOMSS <- peakRAM::peakRAM({
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = NULL, 
                                             b = true_data$M, q = 6, scaling = FALSE)
  })
  
  
  # SUFA 
  profile_SUFA <- peakRAM::peakRAM({
  fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=6, nrun = 10000)
  })

  # BMSFA
  profile_BMSFA <- peakRAM::peakRAM({
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 6, j_s = c(2, 2, 2, 2),
                                      outputlevel = 1, scaling = FALSE, centering = TRUE, trace = FALSE,
                                      control = list(nrun = 10000, burn = 8000))
  })

  ######################################################################
  # Post processing for each methods
  postfit_list <- list(stackFA = post_stackFA(fit_stackFA, S=length(true_data$N_s)), 
                       IndFA = post_IndFA(fit_indFA),
                       PFA = post_PFA(fit_PFA),
                       MOMSS = post_MOMSS(fit_MOMSS, version = 2),
                       BMSFA = post_BMSFA(fit_BMSFA),
                       SUFA = post_SUFA(fit_SUFA))
  
  # Results
  final_list <- list(seed = seed,
                     true = true_data,
                     point_est = postfit_list)
  return(final_list)
}