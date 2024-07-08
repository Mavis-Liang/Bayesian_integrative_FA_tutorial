packages <- c('MSFA', 'peakRAM', 'BFR.BE', 'tidyverse', 'matlab', 'MatrixCorrelation')
lapply(packages, library, character.only = TRUE)
library(devtools)
source("./functions/gen_senerioSS.R")
source("./functions/gen_senerioBMSFA.R")
source("./functions/calculateRV.R")

################################# Functions #################################
# fit the models
fitting <- function(sim_data){
  sen1_Y_mat <- sim_data$Y_mat
  sen1_A <- sim_data$A
  sen1_X <- sim_data$X
  sen1_Y_list <- sim_data$Y_list
  
  profile_MOMSS <- peakRAM({
    #Fit the MOM-SS model
    result_MOMSS_sen1 <- BFR.BE.EM.CV(x = sen1_Y_mat, v = sen1_X, 
                                      b = sen1_A, q = 5, scaling = FALSE)
  })
  
  profile_BMSFA <- peakRAM({
    # #Fit the BMSFA model
    result_BMSFA_sen1 <- sp_msfa(sen1_Y_list, k = 5, j_s = c(1,1,1,1),
                                 outputlevel = 1, scaling = FALSE,
                                 control = list(nrun = 5000, burn = 4000))
  })
  return(list(result_MOMSS = result_MOMSS_sen1, result_BMSFA = result_BMSFA_sen1,
              profile_MOMSS = profile_MOMSS, profile_BMSFA = profile_BMSFA))
}

# measurements
gen_MOMSS_MSE <- function(sim_data, MOMSS_result){
  est_F <- MOMSS_result$Ez
  est_alpha <- MOMSS_result$Theta[, c(3, 4, 5, 6)]
  est_beta <- MOMSS_result$Theta[, c(1, 2)]
  est_Phi <- MOMSS_result$M
  est_Y <- 
    sim_data$A %*% t(est_alpha) + 
    sim_data$X %*% t(est_beta) +
    est_F %*% t(est_Phi)
  MSE <- mean((sim_data$Y_mat - est_Y)^2)
  return(MSE)
}
gen_BMSFA_MSE <- function(sim_data, BMSFA_result){
  ## Assuming the number of common factors and specific factors are known
  est_Phi <- sp_OP(BMSFA_result$Phi, trace=FALSE)$Phi
  est_Lambda <- lapply(BMSFA_result$Lambda, function(x) sp_OP(x, trace=FALSE)$Phi)
  est_Psi <- lapply(BMSFA_result$psi, function(x) sp_OP(x, trace=FALSE)$Phi)# Without diagonalizing
  
  est_Psi_inv <- lapply(est_Psi, function(x) diag(as.vector(1/x))) # NOTE: Inverse are done simplified
  
  est_Y <- NULL
  for (j in 1:4) {
    proj <- solve(t(cbind(est_Phi, est_Lambda[[j]])) %*% est_Psi_inv[[j]] %*% cbind(est_Phi, est_Lambda[[j]]))%*%
      t(cbind(est_Phi, est_Lambda[[j]])) %*%
      est_Psi_inv[[j]]
    est_factors <- sim_data$Y_list[[j]] %*% t(proj)
    
    est_F <- est_factors[, 1:5]
    est_L <- est_factors[, 6]
    est_Y <- rbind(est_Y, est_F %*% t(est_Phi) + est_L %*% t(est_Lambda[[j]]))
  }
  return(mean((sim_data$Y_mat - est_Y)^2))
}

est_perform <- function(true_data, result_BMSFA, result_MOMSS){
  # Scale MOMSS's M
  MOMSS_Phi_scaled <- scale(result_MOMSS$M)
  # RV and FN for SigmaPhi
  RV_MOMSS_Phi <- RV(MOMSS_Phi_scaled, true_data$Phi)
  FN_MOMSS_Phi <- frobenius.norm(tcrossprod(MOMSS_Phi_scaled) - tcrossprod(true_data$Phi))
  
  est_Phi <- sp_OP(result_BMSFA$Phi, trace=FALSE)$Phi

  RV_BMSFA_Phi <- RV(est_Phi, true_data$Phi)
  FN_BMSFA_Phi <- frobenius.norm(tcrossprod(est_Phi) - tcrossprod(true_data$Phi))
  MSE_MOMSS <- gen_MOMSS_MSE(true_data, result_MOMSS)
  MSE_BMSFA <- gen_BMSFA_MSE(true_data, result_BMSFA)
  
  # Output
  return(list(RV_MOMSS_Phi=RV_MOMSS_Phi,
              FN_MOMSS_Phi = FN_MOMSS_Phi,
              RV_BMSFA_Phi = RV_BMSFA_Phi,
              FN_BMSFA_Phi = FN_BMSFA_Phi,
              MSE_MOMSS = MSE_MOMSS,
              MSE_BMSFA = MSE_BMSFA))
}

computing_perform <- function(profiles_list){
  output <- list()
  
  for (i in seq_along(profiles_list)) {
    profile <- profiles_list[[i]]
    name <- names(profiles_list)[i]
    
    output[[paste0("run_time_", name)]] <- profile[,'Elapsed_Time_sec']
    output[[paste0("peak_RAM_", name)]] <- profile[,'Peak_RAM_Used_MiB']
  }
  
  return(output)
}

################################# Example ###########################
# sim_data_test <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
# result_test <- fitting(sim_data_test)
# 
# gen_MOMSS_MSE(sim_data_test, result_test$result_MOMSS)
# 
# gen_BMSFA_MSE(sim_data_test, result_test$result_BMSFA)
# 
# est_perform(sim_data_test, result_test$result_BMSFA, result_test$result_MOMSS)
# computing_perform(result_test[c("profile_MOMSS", "profile_BMSFA")])





  
