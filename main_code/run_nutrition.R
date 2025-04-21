source("./FBPFA-PFA.R")
source("./FBPFA-PFA with fixed latent dim.R")
source("./functions/post_PFA.R")
source("./functions/measurements.R")
source("./Tetris.R")
source("./functions/post_BMSFA.R")
library(BFR.BE)
source("./functions/post_MOMSS.R")
source("./functions/post_stackFA.R")
source("./functions/post_SUFA.R")
source("./functions/post_IndFA.R")
source("./functions/post_Tetris.R")

library(Matrix) # for bdiag function


Y_list <- readRDS("./RDS/nutrition_processed.rds")
#Y_list <- readRDS("./RDS/nutrition_train.rds") # training set

# Centering the data
Y_list_scaled <- lapply(
Y_list, function(x) scale(x, center = TRUE, scale = FALSE)
)
Y_mat_scaled <- Y_list_scaled %>% do.call(rbind, .) %>% as.matrix()

# ------------------------------- Stack FA --------------------------------------
# Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
# fit_stackFA <- MSFA::sp_fa(Y_mat, k = 6, scaling = FALSE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_stackFA = post_stackFA(fit_stackFA, S=6)
# saveRDS(res_stackFA, "./RDS/results_nutrition/Rnutrition_stackFA.rds")


# ------------------------------- Ind FA ----------------------------------------
# fit_indFA <-
#       lapply(1:6, function(s){
#         j_s = c(8, 8, 8, 8, 8, 8)
#         MSFA::sp_fa(Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE,
#                     control = list(nrun = 10000, burn = 8000))
#       }) 
# res_IndFA = post_IndFA(fit_indFA)
# saveRDS(res_IndFA, "./RDS/results_nutrition/Rnutrition_IndFA.rds")

# ------------------------------- MOM-SS ----------------------------------------
# Convert Y_list into Y_mat
# Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
# # Construct the membership matrix
# N_s <- sapply(Y_list, nrow)
# M_list <- list()
#   for(s in 1:6){
#     M_list[[s]] <- matrix(1, nrow = N_s[s], ncol = 1)
#   }
# M <- as.matrix(bdiag(M_list))

# fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = Y_mat, v = NULL, b = M, q = 6, scaling = FALSE)
# res_MOMSS = post_MOMSS(fit_MOMSS, version = 2)
# saveRDS(res_MOMSS, "./RDS/results_nutrition/Rnutrition_MOMSS.rds")


# # ------------------------------- SUFA ------------------------------------------
# fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=6, nrun = 10000)
# res_SUFA = post_SUFA(fit_SUFA)
# saveRDS(res_SUFA, "./RDS/results_nutrition/Rnutrition_SUFA.rds")

# ------------------------------- BMSFA -----------------------------------------
# fit_BMSFA <- MSFA::sp_msfa(Y_list, k = 6, j_s = c(2, 2, 2, 2, 2, 2),
#                                outputlevel = 1, scaling = FALSE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_BMSFA <- post_BMSFA(fit_BMSFA)
# saveRDS(res_BMSFA, "./RDS/results_nutrition/Rnutrition_BMSFA.rds")

# ------------------------------- PFA -------------------------------------------
N_s <- sapply(Y_list, nrow)
fit_PFA <- PFA(Y=t(Y_mat_scaled), 
                        latentdim = 6,
                        grpind = rep(1:6, 
                                    times = N_s),
                Thin = 5,
                Cutoff = 0.001,
                Total_itr = 10000, burn = 8000)
  
res_PFA = post_PFA(fit_PFA)
saveRDS(res_PFA, "./RDS/results_nutrition/Rnutrition_PFA.rds")

# ----------------------------- Tetris -----------------------------------------
# set_alpha <- ceiling(1.25*6)
# fit_Tetris <- tetris(Y_list, alpha=set_alpha, beta=1, nprint = 200, nrun=10000, burn=8000)
# print("Start to choose big_T. It might take a long time.")
# big_T <- choose.A(fit_Tetris, alpha_IBP=set_alpha, S=6)
# run_fixed <- tetris(Y_list, alpha=set_alpha, beta=1, 
#                     fixed=TRUE, A_fixed=big_T, nprint = 200, nrun=10000, burn=8000)
# res_Tetris = post_Tetris(run_fixed)
# saveRDS(res_Tetris, "./RDS/results_nutrition/Rnutrition_Tetris.rds")

# ----------------------------- StackFA 2 --------------------------------------
# Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
# fit_stackFA <- MSFA::sp_fa(Y_mat, k = 4, scaling = FALSE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_stackFA = post_stackFA(fit_stackFA, S=6)
# saveRDS(res_stackFA, "./RDS/results_nutrition/Rnutrition_stackFA_2.rds")

# ------------------------------- Ind FA 2 ----------------------------------------
# fit_indFA <-
#       lapply(1:6, function(s){
#         j_s = c(4, 5, 5, 5, 4, 4)
#         MSFA::sp_fa(Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE,
#                     control = list(nrun = 10000, burn = 8000))
#       }) 
# res_IndFA = post_IndFA(fit_indFA)
# saveRDS(res_IndFA, "./RDS/results_nutrition/Rnutrition_IndFA_2.rds")

# ------------------------------- BMSFA 2 -----------------------------------------
# fit_BMSFA <- MSFA::sp_msfa(Y_list, k = 4, j_s = c(2, 2, 2, 2, 2, 2),
#                                outputlevel = 1, scaling = FALSE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_BMSFA <- post_BMSFA(fit_BMSFA)
# saveRDS(res_BMSFA, "./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")


# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
# ----------------------------Run training set--------------------------------------
# Stack FA
# Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
# fit_stackFA <- MSFA::sp_fa(Y_mat, k = 4, scaling = FALSE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_stackFA = post_stackFA(fit_stackFA, S=6)
# saveRDS(res_stackFA, "./RDS/results_nutrition/Rnutrition_stackFA_train.rds")

# # Ind FA
# fit_indFA <-
#       lapply(1:6, function(s){
#         j_s = c(4, 5, 5, 5, 4, 4)
#         MSFA::sp_fa(Y_list[[s]], k = j_s[s], scaling = FALSE, centering = TRUE,
#                     control = list(nrun = 10000, burn = 8000))
#       })
# res_IndFA = post_IndFA(fit_indFA)
# saveRDS(res_IndFA, "./RDS/results_nutrition/Rnutrition_IndFA_train.rds")

# # # BMSFA
# fit_BMSFA <- MSFA::sp_msfa(Y_list, k = 4, j_s = c(2, 2, 2, 2, 2, 2),
#                                outputlevel = 1, scaling = FALSE, centering = TRUE, 
#                                control = list(nrun = 10000, burn = 8000))
# res_BMSFA <- post_BMSFA(fit_BMSFA)
# saveRDS(res_BMSFA, "./RDS/results_nutrition/Rnutrition_BMSFA_train.rds")

# # PFA
# N_s <- sapply(Y_list, nrow)
# fit_PFA <- PFA(Y=t(Y_mat_scaled), 
#                         latentdim = 6,
#                         grpind = rep(1:6, 
#                                     times = N_s),
#                 Cutoff = 0.001,
#                 Thin = 5,
#                 Total_itr = 10000, burn = 8000)
  
# res_PFA = post_PFA(fit_PFA)
# saveRDS(res_PFA, "./RDS/results_nutrition/Rnutrition_PFA_train.rds")

# # # MOM-SS
# # # Convert Y_list into Y_mat
# Y_mat =  Y_list %>% do.call(rbind, .) %>% as.matrix()
# # Construct the membership matrix
# N_s <- sapply(Y_list, nrow)
# M_list <- list()
#   for(s in 1:6){
#     M_list[[s]] <- matrix(1, nrow = N_s[s], ncol = 1)
#   }
# M <- as.matrix(bdiag(M_list))

# fit_MOMSS <- BFR.BE::BFR.BE.EM.CV(x = Y_mat, v = NULL, b = M, q = 6, scaling = FALSE)
# res_MOMSS = post_MOMSS(fit_MOMSS, version = 2)
# saveRDS(res_MOMSS, "./RDS/results_nutrition/Rnutrition_MOMSS_train.rds")

# # # SUFA 
# fit_SUFA <- SUFA::fit_SUFA(Y_list_scaled, qmax=6, nrun = 10000)
# res_SUFA = post_SUFA(fit_SUFA)
# saveRDS(res_SUFA, "./RDS/results_nutrition/Rnutrition_SUFA_train.rds")

# Tetris 
# set_alpha <- ceiling(1.25*6)
# fit_Tetris <- tetris(Y_list, alpha=set_alpha, beta=1, nprint = 200, nrun=10000, burn=8000)
# print("Start to choose big_T. It might take a long time.")
# big_T <- choose.A(fit_Tetris, alpha_IBP=set_alpha, S=6)
# run_fixed <- tetris(Y_list, alpha=set_alpha, beta=1, 
#                     fixed=TRUE, A_fixed=big_T, nprint = 200, nrun=10000, burn=8000)
# res_Tetris = post_Tetris(run_fixed)
# saveRDS(res_Tetris, "./RDS/results_nutrition/Rnutrition_Tetris_train.rds")