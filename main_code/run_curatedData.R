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

# Whole data or training set
   Y_list <- readRDS("./RDS/CuratedOvarian_processed.rds")
  #  Y_list <- readRDS("./RDS/curated_train.rds")
  # Centering the data
  Y_list_scaled <- lapply(
    Y_list, function(x) scale(x, center = TRUE, scale = TRUE)
  )
  Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()

  # ------------------------------- Stack FA --------------------------------------
  # Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
  # fit_stackFA <- MSFA::sp_fa(Y_mat, k = 20, scaling = TRUE, centering = TRUE, 
  #                              control = list(nrun = 10000, burn = 8000))
  # res_stackFA = post_stackFA(fit_stackFA, S=4)
  # saveRDS(res_stackFA, "./RDS/results_curated/RCuratedOvarian_stackFA_scaled.rds")

  # ------------------------------- Stack FA 2 --------------------------------------
  # Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
  # fit_stackFA <- MSFA::sp_fa(Y_mat, k = 2, scaling = TRUE, centering = TRUE, 
  #                              control = list(nrun = 10000, burn = 8000))
  # res_stackFA = post_stackFA(fit_stackFA, S=4)
  # saveRDS(res_stackFA, "./RDS/results_curated/RCuratedOvarian_stackFA2_scaled.rds")

# ------------------------------- Ind FA ----------------------------------------
# fit_indFA <-
#       lapply(1:4, function(s){
#         j_s = c(20, 20, 20, 20)
#         MSFA::sp_fa(Y_list[[s]], k = j_s[s], scaling = TRUE, centering = TRUE,
#                     control = list(nrun = 10000, burn = 8000))
#       }) 
# res_IndFA = post_IndFA(fit_indFA)
# saveRDS(res_IndFA, "./RDS/results_curated/RCuratedOvarian_IndFA_scaled.rds")

# ------------------------------- Ind FA 2 ----------------------------------------
# fit_indFA <-
#       lapply(1:4, function(s){
#         j_s = c(1, 1, 4, 1)
#         MSFA::sp_fa(Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE,
#                     control = list(nrun = 10000, burn = 8000))
#       }) 
# res_IndFA = post_IndFA(fit_indFA)
# saveRDS(res_IndFA, "./RDS/results_curated/RCuratedOvarian_IndFA2.rds")

# ------------------------------- MOM-SS ----------------------------------------
# Convert Y_list into Y_mat
# Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
# # Construct the membership matrix
# N_s <- sapply(Y_list, nrow)
# M_list <- list()
#   for(s in 1:4){
#     M_list[[s]] <- matrix(1, nrow = N_s[s], ncol = 1)
#   }
# M <- as.matrix(bdiag(M_list))

# fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = Y_mat, v = NULL, b = M, q = 20, scaling = TRUE)
# res_MOMSS = post_MOMSS(fit_MOMSS, version = 2)
# saveRDS(res_MOMSS, "./RDS/results_curated/RCuratedOvarian_MOMSS_scaled.rds")

# ------------------------------- BMSFA ------------------------------------------
# fit_BMSFA <- MSFA::sp_msfa(Y_list, k = 20, j_s = c(4, 4, 4, 4),
#                                outputlevel = 1, scaling = TRUE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_BMSFA <- post_BMSFA(fit_BMSFA)
# saveRDS(res_BMSFA, "./RDS/results_curated/RCuratedOvarian_BMSFA.rds")

# ------------------------------- BMSFA 2------------------------------------------
# fit_BMSFA <- MSFA::sp_msfa(Y_list, k = 6, j_s = c(4, 4, 4, 4),
#                                outputlevel = 1, scaling = TRUE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_BMSFA <- post_BMSFA(fit_BMSFA)
# saveRDS(res_BMSFA, "./RDS/results_curated/RCuratedOvarian_BMSFA2_scaled.rds")
  
  # ------------------------------- SUFA --------------------------------------
  # fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=20,nrun = 10000)
  # res_SUFA = post_SUFA(fit_SUFA)
  # saveRDS(res_SUFA, "./RDS/results_curated/RCuratedOvarian_SUFA_scaled.rds")
  
  # ------------------------------- PFA ----------------------------------------
  N_s <- sapply(Y_list, nrow)
  fit_PFA <- PFA(Y=t(Y_mat_scaled), 
                         latentdim = 20,
                         grpind = rep(1:4, 
                                      times = N_s),
                Thin = 5,
                Cutoff = 0.001,
                Total_itr = 10000, burn = 8000)
  
  res_PFA = post_PFA(fit_PFA)
  saveRDS(res_PFA, "./RDS/results_curated/RCuratedOvarian_PFA_scaled.rds")

  # ------------------------------- Tetris ----------------------------------------
  # set_alpha <- ceiling(1.25*4)
  # fit_Tetris <- tetris(Y_list_scaled, alpha=set_alpha, beta=1, nprint = 200, nrun=10000, burn=8000)
  # print("Start to choose big_T. It might take a long time.")
  # big_T <- choose.A(fit_Tetris, alpha_IBP=set_alpha, S=4)
  # run_fixed <- tetris(Y_list_scaled, alpha=set_alpha, beta=1, 
  #                     fixed=TRUE, A_fixed=big_T, nprint = 200, nrun=10000, burn=8000)
  # res_Tetris = post_Tetris(run_fixed)
  # saveRDS(res_Tetris, "./RDS/results_curated/RCuratedOvarian_Tetris_scaled.rds")


# Fitting on the training data
# # ------------------------------- PFA ----------------------------------------
  # N_s <- sapply(Y_list, nrow)
  # fit_PFA <- PFA(Y=t(Y_mat_scaled), 
  #                        latentdim = 20,
  #                        grpind = rep(1:4, 
  #                                     times = N_s),
  #               Cutoff = 0.001,
  #               Thin = 5,
  #               Total_itr = 10000, burn = 8000)
  
  # res_PFA = post_PFA(fit_PFA)
  # saveRDS(res_PFA, "./RDS/results_curated/RCuratedOvarian_PFA_train.rds")

  #   # ------------------------------- Stack FA --------------------------------------
  # Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
  # fit_stackFA <- MSFA::sp_fa(Y_mat, k = 2, scaling = FALSE, centering = TRUE, 
  #                              control = list(nrun = 10000, burn = 8000))
  # res_stackFA = post_stackFA(fit_stackFA, S=4)
  # saveRDS(res_stackFA, "./RDS/results_curated/RCuratedOvarian_stackFA_train.rds")

  # ------------------------------- Ind FA 2 ----------------------------------------
# fit_indFA <-
#       lapply(1:4, function(s){
#         j_s = c(4, 7, 6, 6)
#         MSFA::sp_fa(Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE,
#                     control = list(nrun = 10000, burn = 8000))
#       }) 
# res_IndFA = post_IndFA(fit_indFA)
# saveRDS(res_IndFA, "./RDS/results_curated/RCuratedOvarian_IndFA2_train.rds")

# #------------------------------- MOM-SS ----------------------------------------
# # # Convert Y_list into Y_mat
# Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
# # Construct the membership matrix
# N_s <- sapply(Y_list, nrow)
# M_list <- list()
#   for(s in 1:4){
#     M_list[[s]] <- matrix(1, nrow = N_s[s], ncol = 1)
#   }
# M <- as.matrix(bdiag(M_list))

# fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = Y_mat, v = NULL, b = M, q = 20, scaling = FALSE)
# res_MOMSS = post_MOMSS(fit_MOMSS, version = 2)
# saveRDS(res_MOMSS, "./RDS/results_curated/RCuratedOvarian_MOMSS_train.rds")

# # # ------------------------------- BMSFA ------------------------------------------
# fit_BMSFA <- MSFA::sp_msfa(Y_list, k = 6, j_s = c(4, 4, 4, 4),
#                                outputlevel = 1, scaling = FALSE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_BMSFA <- post_BMSFA(fit_BMSFA)
# saveRDS(res_BMSFA, "./RDS/results_curated/RCuratedOvarian_BMSFA_train.rds")

# #   # ------------------------------- SUFA --------------------------------------
#   fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=20,nrun = 10000)
#   res_SUFA = post_SUFA(fit_SUFA)
#   saveRDS(res_SUFA, "./RDS/results_curated/RCuratedOvarian_SUFA_train.rds")
# #   # ------------------------------- Tetris ----------------------------------------
#   set_alpha <- ceiling(1.25*4)
#   fit_Tetris <- tetris(Y_list, alpha=set_alpha, beta=1, nprint = 200, nrun=10000, burn=8000)
#   print("Start to choose big_T. It might take a long time.")
#   big_T <- choose.A(fit_Tetris, alpha_IBP=set_alpha, S=4)
#   run_fixed <- tetris(Y_list, alpha=set_alpha, beta=1, 
#                       fixed=TRUE, A_fixed=big_T, nprint = 200, nrun=10000, burn=8000)
#   # res_Tetris = post_Tetris(run_fixed)
#   saveRDS(run_fixed, "./RDS/results_curated/Tetris_train_run_fixed.rds")


