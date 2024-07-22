source("./functions/calculateRV.R")
gen_MOMSS_MSE <- function(sim_data, MOMSS_result){
  est_F <- MOMSS_result$Ez
  est_alpha <- MOMSS_result$Theta[, (ncol(sim_data$Beta) + 1) : ncol(MOMSS_result$Theta)]
  est_beta <- MOMSS_result$Theta[, 1:ncol(sim_data$Beta)]
  est_Phi <- MOMSS_result$M
  est_Y <- 
    sim_data$A %*% t(est_alpha) + 
    sim_data$X %*% t(est_beta) +
    est_F %*% t(est_Phi)
  MSE <- norm(sim_data$Y_mat - est_Y, "F")^2/nrow(sim_data$Y_mat)
  return(MSE)
}

gen_BMSFA_MSE <- function(sim_data, BMSFA_result){
  est_Phi <- BMSFA_result$est_Phi
  est_Lambda <- BMSFA_result$est_Lambda
  est_Psi <- BMSFA_result$est_Psi
  est_Psi_inv <- lapply(est_Psi, function(x) diag(as.vector(1/x))) # NOTE: Inverse are done simplified
  
  SSE <- NULL
  for (s in 1:length(sim_data$Y_list)) {
    combined_loadings <- cbind(est_Phi, est_Lambda[[s]])
    proj <- solve(t(combined_loadings) %*% est_Psi_inv[[s]] %*% combined_loadings)%*%
      t(combined_loadings) %*%
      est_Psi_inv[[s]]
    est_factors <- proj %*% t(sim_data$Y_list[[s]])
    
    est_Y <- t(combined_loadings %*% est_factors)
    
    SSE_S <- norm(sim_data$Y_list[[s]] - est_Y, "F")^2
    SSE <- c(SSE, SSE_S)
  }
  return(sum(SSE)/nrow(sim_data$Y_mat))
}

est_perform <- function(true_data, result_all){
  result_MOMSS <- result_all$result_MOMSS
  result_BMSFA <- result_all$result_BMSFA
  result_PFA <- result_all$result_PFA
  result_SUFA <- result_all$result_SUFA
  # Scale MOMSS's M
  MOMSS_Phi <- result_MOMSS$M
  #------------ MOMSS --------------
  RV_MOMSS_Phi <- MatrixCorrelation::RV(MOMSS_Phi, true_data$Phi)
  #FN_MOMSS_Phi <- frobenius.norm(tcrossprod(MOMSS_Phi) - tcrossprod(true_data$Phi))
  RV_MOMSS_SigmaPhi <- MatrixCorrelation::RV(tcrossprod(MOMSS_Phi), 
                                             tcrossprod(true_data$Phi))
  MSE_MOMSS <- gen_MOMSS_MSE(true_data, result_MOMSS)
  
  #------------ BMSFA --------------
  RV_BMSFA_Phi <- MatrixCorrelation::RV(result_BMSFA$est_Phi, true_data$Phi)
  #FN_BMSFA_Phi <- frobenius.norm(tcrossprod(result_BMSFA$est_Phi) - tcrossprod(true_data$Phi))
  RV_BMSFA_SigmaPhi <- MatrixCorrelation::RV(tcrossprod(result_BMSFA$est_Phi), 
                                             tcrossprod(true_data$Phi))
  MSE_BMSFA <- gen_BMSFA_MSE(true_data, result_BMSFA)
  
  #------------- PFA ---------------
  RV_PFA_Phi <- MatrixCorrelation::RV(result_PFA$est_Phi, true_data$Phi)
  RV_PFA_SigmaPhi <- MatrixCorrelation::RV(result_PFA$est_SigmaPhi, 
                                           tcrossprod(true_data$Phi))
  
  
  #------------- SUFA ---------------
  RV_SUFA_Phi <- MatrixCorrelation::RV(result_SUFA$est_Phi, true_data$Phi)
  #FN_SUFA_Phi <- frobenius.norm(tcrossprod(result_SUFA$est_Phi) - tcrossprod(true_data$Phi))
  RV_SUFA_SigmaPhi <- MatrixCorrelation::RV(result_SUFA$est_SigmaPhi, 
                                            tcrossprod(true_data$Phi))
  # Output
  return(list(RV_MOMSS_Phi=RV_MOMSS_Phi,
              RV_MOMSS_SigmaPhi = RV_MOMSS_SigmaPhi,
              #FN_MOMSS_Phi = FN_MOMSS_Phi,
              RV_BMSFA_Phi = RV_BMSFA_Phi,
              #FN_BMSFA_Phi = FN_BMSFA_Phi,
              RV_BMSFA_SigmaPhi = RV_BMSFA_SigmaPhi,
              RV_PFA_Phi = RV_PFA_Phi,
              RV_PFA_SigmaPhi = RV_PFA_SigmaPhi,
              RV_SUFA_Phi = RV_SUFA_Phi,
              #FN_SUFA_Phi = FN_SUFA_Phi,
              RV_SUFA_SigmaPhi = RV_SUFA_SigmaPhi,
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