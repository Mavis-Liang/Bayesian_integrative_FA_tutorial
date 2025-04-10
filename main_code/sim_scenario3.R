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
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()
  # Fitting
  # Stack FA
  profile_stackFA <- peakRAM::peakRAM({
    fit_stackFA <- MSFA::sp_fa(true_data$Y_mat, k = 4,scaling = FALSE, centering = TRUE, trace = FALSE,
                               control = list(nrun = 10000, burn = 8000))
  })
  
  # Ind FA
  profile_indFA <- peakRAM::peakRAM({
    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        MSFA::sp_fa(true_data$Y_list[[s]], k = 4,scaling = FALSE, centering = TRUE, trace = FALSE,
                    control = list(nrun = 10000, burn = 8000))
      })
  })
  
  # PFA
  profile_PFA <- peakRAM::peakRAM({
    fit_PFA <- PFA(Y=t(Y_mat_scaled), latentdim = 4, grpind = true_data$grpind,
                   Thin = 5, Total_itr = 10000, burn = 8000)
  })
  
  # MOMSS
  profile_MOMSS <- peakRAM::peakRAM({
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = NULL, 
                                      b = true_data$M, q = 4, scaling = FALSE)
  })
  
  # SUFA (Now SUFA also support specifying the number of study-specific factors.)
  profile_SUFA_fixJs <- peakRAM::peakRAM({
    fit_SUFA_fixJs <- SUFA::fit_SUFA(true_data$Y_list, qmax=4, nrun = 10000, qs = c(1, 1, 1, 1))
  })
  
  # SUFA
  profile_SUFA <- peakRAM::peakRAM({
    fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=4,nrun = 10000)
  })
  # SUFA now also support specifying the number of study-specific factors. J_s = c(0, 1, 1, 1). J_s can be 0.
  
  # BMSFA
  profile_BMSFA <- peakRAM::peakRAM({
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 4, j_s = c(1,1,1,1), # j_s can't be 0
                               outputlevel = 1, scaling = FALSE, centering = TRUE, trace = FALSE,
                               control = list(nrun = 10000, burn = 8000))
  })

  # Tetris with fixed T
  T_mat <- matrix(0, nrow = 4, ncol = 4 + 1 + 1 + 1 + 1)
  T_mat[,1:4] <- T_mat[1, 5] <- T_mat[2, 6] <- T_mat[3, 7] <- T_mat[4, 8] <- 1
  profile_Tetris_fixT <- peakRAM::peakRAM({
    fit_Tetris_fixT <- run_Tetris(true_data, big_T=T_mat)
  })
  
  # Tetris is runned on separate process, since it is slow and not parallelizable for large cores
  # Tetris
  #  profile_Tetris <- peakRAM::peakRAM({
  #   fit_tetris <- run_Tetris(true_data)
  # })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(stackFA = post_stackFA(fit_stackFA, S=4), 
                       IndFA = post_IndFA(fit_indFA),
                       PFA = post_PFA(fit_PFA),
                       MOMSS = post_MOMSS(fit_MOMSS),
                       SUFA_fixJs = post_SUFA(fit_SUFA_fixJs),
                       SUFA = post_SUFA(fit_SUFA),
                       BMSFA = post_BMSFA(fit_BMSFA),
                       Tetris_fixT = post_Tetris(fit_Tetris_fixT)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(profile_stackFA[,c(2, 4)], profile_indFA[,c(2, 4)], 
                      profile_PFA[,c(2, 4)], profile_MOMSS[,c(2, 4)],
                      profile_SUFA_fixJs[,c(2, 4)],profile_SUFA[,c(2, 4)],
                      profile_BMSFA[,c(2, 4)], profile_Tetris_fixT[,c(2, 4)]
                      
  ) %>% 
    t() %>% `colnames<-`(c("stackFA", "IndFA","PFA",  "MOMSS", "SUFA_fixJs",
                           "SUFA", "BMSFA", "Tetris_fixT"))
  
  # Results
  final_list <- list(seed = seed,
                     true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}



sim_scenario3_mis <- function(seed=1){
  set.seed(seed)
  # # Data generation 
  true_data <- gen_scenario3(4, N_s=c(100, 100, 100, 100), P=40, K=4)
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()
  # Fitting
  # Stack FA
    fit_stackFA <- MSFA::sp_fa(true_data$Y_mat, k = 6,scaling = FALSE, centering = TRUE, trace = FALSE,
                               control = list(nrun = 10000, burn = 8000))
  
  # Ind FA

    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        MSFA::sp_fa(true_data$Y_list[[s]], k = 6,scaling = FALSE, centering = TRUE, trace = FALSE,
                    control = list(nrun = 10000, burn = 8000))
      })

  # PFA
    fit_PFA <- PFA(Y=t(Y_mat_scaled), latentdim = 6, grpind = true_data$grpind,
                   Thin = 5, Total_itr = 10000, burn = 8000)

  
  # MOMSS
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = NULL, 
                                      b = true_data$M, q = 6, scaling = FALSE)

  
  # SUFA

    fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=6,nrun = 10000)
 
  
  # BMSFA
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 6, j_s = c(2,2,2,2), # j_s can't be 0
                               outputlevel = 1, scaling = FALSE, centering = TRUE, trace = FALSE,
                               control = list(nrun = 10000, burn = 8000))

  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(stackFA = post_stackFA(fit_stackFA, S=4), 
                       IndFA = post_IndFA(fit_indFA),
                       PFA = post_PFA(fit_PFA),
                       MOMSS = post_MOMSS(fit_MOMSS, version = 2),
                       SUFA = post_SUFA(fit_SUFA),
                       BMSFA = post_BMSFA(fit_BMSFA)
                       )
  
  
  # Results
  final_list <- list(seed = seed,
                     point_est = postfit_list)
  return(final_list)
}
