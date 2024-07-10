source("./functions/calculateRV.R")
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
  est_Phi <- BMSFA_result$est_Phi
  est_Lambda <- BMSFA_result$est_Lambda
  est_Psi <- BMSFA_result$est_Psi
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

est_perform <- function(true_data, result_BMSFA, result_MOMSS, result_PFA){
  # Scale MOMSS's M
  MOMSS_Phi <- result_MOMSS$M
  # MOMSS
  RV_MOMSS_Phi <- RV(MOMSS_Phi, true_data$Phi)
  FN_MOMSS_Phi <- frobenius.norm(tcrossprod(MOMSS_Phi) - tcrossprod(true_data$Phi))
  RV_MOMSS_SigmaPhi <- calculateRV(tcrossprod(MOMSS_Phi), tcrossprod(true_data$Phi))
  MSE_MOMSS <- gen_MOMSS_MSE(true_data, result_MOMSS)
  
  # BMSFA
  RV_BMSFA_Phi <- RV(result_BMSFA$est_Phi, true_data$Phi)
  FN_BMSFA_Phi <- frobenius.norm(tcrossprod(result_BMSFA$est_Phi) - tcrossprod(true_data$Phi))
  RV_BMSFA_SigmaPhi <- calculateRV(tcrossprod(result_BMSFA$est_Phi), tcrossprod(true_data$Phi))
  MSE_BMSFA <- gen_BMSFA_MSE(true_data, result_BMSFA)
  
  # PFA
  RV_PFA_SigmaPhi <- calculateRV(result_PFA$est_SigmaPhi, tcrossprod(true_data$Phi))
  # Output
  return(list(RV_MOMSS_Phi=RV_MOMSS_Phi,
              FN_MOMSS_Phi = FN_MOMSS_Phi,
              RV_BMSFA_Phi = RV_BMSFA_Phi,
              FN_BMSFA_Phi = FN_BMSFA_Phi,
              RV_MOMSS_SigmaPhi = RV_MOMSS_SigmaPhi,
              RV_BMSFA_SigmaPhi = RV_BMSFA_SigmaPhi,
              RV_PFA_SigmaPhi = RV_PFA_SigmaPhi,
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