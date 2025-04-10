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
source("./functions/post_MOMSS.R")
source("./functions/measurements.R")

# --------------------------- Scenario mimic genomic data (scenario 5 # [Scenario 6 in the manuscript])---------------------------
# This function produce the result of a single simulation of the scenario 4
# includes data generation, fitting the 7 models, post processing for the point estimates, and FN and RV metrics
# This function is designed for use in parallel computing. 
# Tetris is run on a separate file.
source("./functions/gen_scenario5.R")
sim_scenario5_IndFA <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()
  
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
  
  # Ind FA
  profile_indFA <- peakRAM::peakRAM({
    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        j_s = c(17, 17, 17, 17)
        MSFA::sp_fa(true_data$Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE,
                    control = list(nrun = 10000, burn = 8000))
      })
  })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list( IndFA = post_IndFA(fit_indFA)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(#profile_stackFA[,c(2, 4)]#, 
                      profile_indFA[,c(2, 4)]) %>% 
    t() %>% as.data.frame()  %>% `colnames<-`(c(#"stackFA"#, 
                            "IndFA"
                           ))
  
  # Results
  final_list <- list(true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}
sim_scenario5_SUFAfixJ <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()
  
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
  
  
  # SUFA (Now SUFA also support specifying the number of study-specific factors.)
  profile_SUFA_fixJs <- peakRAM::peakRAM({
    fit_SUFA_fixJs <- SUFA::fit_SUFA(Y_list_scaled, qmax=15, nrun = 10000, 
                                    qs = c(2, 2, 2, 2)) #"Identifiability condition for `infomration switching` not statisfied!"
  })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list( SUFA_fixJs = post_SUFA(fit_SUFA_fixJs)#,
                       )
  # Time and peak RAM profile
  profile_df <- rbind(#profile_stackFA[,c(2, 4)]#, 
                      #profile_indFA[,c(2, 4)]#, 
                      #profile_PFA[,c(2, 4)], 
                       #profile_MOMSS[,c(2, 4)],
                       profile_SUFA_fixJs[,c(2, 4)]#,
                      # profile_SUFA[,c(2, 4)],
                      #profile_BMSFA[,c(2, 4)]#, 
                      #profile_Tetris_fixT[,c(2, 4)]
                      #,profile_Tetris[,c(2, 4)]
                      ) %>% 
    t() %>% as.data.frame()  %>% `colnames<-`(c(#"stackFA"#, 
                            #"IndFA"#,
                            #"PFA", 
                            #"MOMSS",
                             "SUFA_fixJs"#, 
                             #"SUFA"#, 
                           #"BMSFA"#, 
                           #"Tetris_fixT"
                           ))
  
  # Results
  final_list <- list(#true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

sim_scenario5_2 <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()

   
  # MOMSS
  profile_MOMSS <- peakRAM::peakRAM({
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = NULL, 
                                      b = true_data$M, q = 15, scaling = FALSE)
  }) 

  # # BMSFA
  profile_BMSFA <- peakRAM::peakRAM({
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 15, j_s = c(2, 2, 2, 2),
                               outputlevel = 1, scaling = FALSE, centering = TRUE, 
                               control = list(nrun = 10000, burn = 8000))
  })

  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(#stackFA = post_stackFA(fit_stackFA, S=length(true_data$N_s)), 
                       #IndFA = post_IndFA(fit_indFA),
                       #PFA = post_PFA(fit_PFA),
                       MOMSS = post_MOMSS(fit_MOMSS),
                        #SUFA_fixJs = post_SUFA(fit_SUFA_fixJs)#,
                      #  SUFA = post_SUFA(fit_SUFA),
                       BMSFA = post_BMSFA(fit_BMSFA)#,
                      #  Tetris_fixT = post_Tetris(fit_Tetris_fixT)
                       #,Tetris = post_Tetris(fit_tetris)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(#profile_stackFA[,c(2, 4)], profile_indFA[,c(2, 4)], 
                      #profile_PFA[,c(2, 4)], 
                       profile_MOMSS[,c(2, 4)],
                       #profile_SUFA_fixJs[,c(2, 4)]#,
                      # profile_SUFA[,c(2, 4)],
                      profile_BMSFA[,c(2, 4)]#, 
                      #profile_Tetris_fixT[,c(2, 4)]
                      #,profile_Tetris[,c(2, 4)]
                      ) %>% 
    t() %>% `colnames<-`(c(#"stackFA", "IndFA",#"PFA", 
                            "MOMSS",
                             #"SUFA_fixJs"#, 
                             #"SUFA", 
                           "BMSFA"#, 
                           #"Tetris_fixT"
                           ))
  
  # Results
  final_list <- list(#true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

sim_scenario5_3 <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()
  
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
 
  # # SUFA
  profile_SUFA <- peakRAM::peakRAM({
    fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=15,nrun = 10000)
  })
  profile_Tetris_fixT <- peakRAM::peakRAM({
    fit_Tetris_fixT <- run_Tetris(true_data, big_T=true_data$big_T)# A function made to wrap Tetris
  })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(SUFA = post_SUFA(fit_SUFA),
                        Tetris_fixT = post_Tetris(fit_Tetris_fixT)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(profile_SUFA[,c(2, 4)],
                      # profile_BMSFA[,c(2, 4)]#, 
                      profile_Tetris_fixT[,c(2, 4)]
                      #,profile_Tetris[,c(2, 4)]
                      ) %>% 
    t() %>% `colnames<-`(c( "SUFA", 
                          #  "BMSFA"#, 
                           "Tetris_fixT"
                           ))
  
  # Results
  final_list <- list(#true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

sim_scenario5_TetrisfixT <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
  profile_Tetris_fixT <- peakRAM::peakRAM({
    fit_Tetris_fixT <- run_Tetris(true_data, big_T=true_data$big_T)# A function made to wrap Tetris
  })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(Tetris_fixT = post_Tetris(fit_Tetris_fixT)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(profile_Tetris_fixT[,c(2, 4)]) %>% 
    t() %>% as.data.frame()  %>% `colnames<-`(c( "Tetris_fixT"))
  
  # Results
  final_list <- list(#true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

sim_scenario5_Tetris <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
  profile_Tetris <- peakRAM::peakRAM({
    fit_Tetris <- run_Tetris(true_data)# A function made to wrap Tetris
  })
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(Tetris_fixT = post_Tetris(fit_Tetris)
                       )
  # Time and peak RAM profile
  profile_df <- rbind(profile_Tetris[,c(2, 4)]) %>% 
    t() %>% as.data.frame()  %>% `colnames<-`(c( "Tetris"))
  
  # Results
  final_list <- list(#true = true_data,
                     point_est = postfit_list, 
                     metrics = metrics_tables(postfit_list, true_data), 
                     efficiency = list(profile_df))
  return(final_list)
}

sim_scenario5_StackFA_mis <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  
  # Ind FA
  fit_stackFA <- MSFA::sp_fa(true_data$Y_mat, k = 20, scaling = FALSE, centering = TRUE, trace = TRUE,
                               control = list(nrun = 10000, burn = 8000))
  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list( stackFA = post_stackFA(fit_stackFA, S=length(true_data$N_s))
                       )
  
  # Results
  final_list <- list(point_est = postfit_list)
  return(final_list)
}


sim_scenario5_IndFA_mis <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  
    fit_indFA <-
      lapply(1:length(true_data$N_s), function(s){
        j_s = c(20, 20, 20, 20)
        MSFA::sp_fa(true_data$Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE,
                    control = list(nrun = 10000, burn = 8000))
      })

  
  ######################################################################
  # Post processing for each methods
  postfit_list <- list(IndFA = post_IndFA(fit_indFA))
  
  # Results
  final_list <- list(point_est = postfit_list)
  return(final_list)
}


sim_scenario5_MOMSS_BMSFA_mis <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()

   
  # MOMSS
  profile_MOMSS <- peakRAM::peakRAM({
    fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = true_data$Y_mat, v = NULL, 
                                      b = true_data$M, q = 20, scaling = FALSE)
  }) 

  # # BMSFA
  profile_BMSFA <- peakRAM::peakRAM({
    fit_BMSFA <- MSFA::sp_msfa(true_data$Y_list, k = 20, j_s = c(4, 4, 4, 4),
                               outputlevel = 1, scaling = FALSE, centering = TRUE, 
                               control = list(nrun = 10000, burn = 8000))
  })

  # Post processing for each methods
  postfit_list <- list(MOMSS = post_MOMSS(fit_MOMSS, version = 2),
                       BMSFA = post_BMSFA(fit_BMSFA)
                       )
  
  # Results
  final_list <- list(point_est = postfit_list)
  return(final_list)
}

sim_scenario5_SUFA_mis <- function(seed=i){
  set.seed(seed)
  # Data generation
  true_data <- gen_scenario5()
  # Centering the data
  Y_list_scaled <- lapply(
    true_data$Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
  )
  
  # Fitting (peakRAM is used to measure the peak RAM usage and time consumed of each model)
 
  # # SUFA
  fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=20,nrun = 10000)

  
    # Post processing for each methods
  postfit_list <- list(SUFA = post_SUFA(fit_SUFA))
  
  # Results
  final_list <- list(point_est = postfit_list)
  return(final_list)
}

