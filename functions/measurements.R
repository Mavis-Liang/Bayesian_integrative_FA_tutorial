source("./functions/RV_pd.R")
library(tidyverse)
measure_common <- function(sim_data, post_output){
  # Obtain True common covariance matrix and loading
  true_Phi <- sim_data$Phi
  true_SigmaPhi <- sim_data$SigmaPhi
  
  # Obtain estimated common covariance matrix and loading
  est_Phi <- post_output$Phi
  est_SigmaPhi <- post_output$SigmaPhi
  
  # Metrics
  if(is.null(true_Phi) | is.null(est_Phi)){ # e.g. Ind FA does not have Phi
    FN_Phi <- NA
    FN_Phi_filled <- NA
    RV_Phi <- NA
    RV_Phi_filled <- NA
    FN_PhiPhiT <- NA
    FN_SigmaPhi <- NA
    RV2_SigmaPhi <- NA
    RV_SigmaPhi <- NA
  } else {
  # Frobenius Norm for Phi (case 1: assigning NA) e.g. in Tetris
  if(dim(true_Phi)[2] == dim(est_Phi)[2]){ # Same K dimension
    FN_Phi <- frobenius.norm(est_Phi - true_Phi)
  } else { # Different K dimensions
    FN_Phi <- NA
  }
  
  # Frobenius Norm and RV for Phi (case 2: filling smaller Phi with 0s) e.g. in Tetris
  if(dim(true_Phi)[2] > dim(est_Phi)[2]){ # true_Phi has larger K
    est_Phi_filled <- cbind(est_Phi, matrix(0, nrow = nrow(est_Phi), ncol = dim(true_Phi)[2] - dim(est_Phi)[2]))
    FN_Phi_filled <- frobenius.norm(est_Phi_filled - true_Phi)
    RV_Phi_filled <- MatrixCorrelation::RV(est_Phi_filled, true_Phi)
  } else if(dim(true_Phi)[2] < dim(est_Phi)[2]){ # est_Phi has larger K
    true_Phi_filled <- cbind(true_Phi, matrix(0, nrow = nrow(true_Phi), ncol = dim(est_Phi)[2] - dim(true_Phi)[2]))
    FN_Phi_filled <- frobenius.norm(est_Phi - true_Phi_filled)
    RV_Phi_filled <- MatrixCorrelation::RV(est_Phi, true_Phi_filled)
  } else { # Same K dimension
    FN_Phi_filled <- frobenius.norm(est_Phi - true_Phi)
    RV_Phi_filled <- MatrixCorrelation::RV(est_Phi, true_Phi)
  }
  
  # Frobenius Norm for Phi*Phi^T and SigmaPhi, and RV for SigmaPhi
  RV_Phi <- MatrixCorrelation::RV(est_Phi, true_Phi)
  FN_PhiPhiT <- frobenius.norm(tcrossprod(est_Phi) - tcrossprod(true_Phi))
  FN_SigmaPhi <- frobenius.norm(est_SigmaPhi - true_SigmaPhi)
  RV2_SigmaPhi <- RV_pd(est_SigmaPhi, true_SigmaPhi)
  RV_SigmaPhi <- MatrixCorrelation::RV(est_SigmaPhi, true_SigmaPhi)
}
  # Return results
  return(list(FN_Phi = FN_Phi, FN_Phi_filled = FN_Phi_filled, 
              FN_PhiPhiT = FN_PhiPhiT, FN_SigmaPhi = FN_SigmaPhi, 
              RV_Phi = RV_Phi, RV_Phi_filled = RV_Phi_filled, 
              RV2_SigmaPhi = RV2_SigmaPhi, RV_SigmaPhi = RV_SigmaPhi))
}

# This is a helper function for measure_specific()
measure_specific_unmatch_Lambda <- function(est_LambdaList, true_LambdaList, S){
  # Particularly in Tetris, the number of columns in Lambda can be different
  # If dimensions are different, then NA is assigned
  FN_Lambda <- sapply(1:S, function(i) {
    if(dim(true_LambdaList[[i]])[2] == dim(est_LambdaList[[i]])[2]){
      frobenius.norm(est_LambdaList[[i]] - true_LambdaList[[i]])
    } else {
      NA
    }
  })
  RV_Lambda <- sapply(1:S, function(i) {
      MatrixCorrelation::RV(est_LambdaList[[i]], true_LambdaList[[i]], center = FALSE)
  })
  # Filling zeros for either true or estimated matrix with lower dimensions
  FN_Lambda_filled <- sapply(1:S, function(i) {
    if(dim(true_LambdaList[[i]])[2] > dim(est_LambdaList[[i]])[2]){
      est_Lambda_filled <- cbind(est_LambdaList[[i]], matrix(0, nrow = nrow(est_LambdaList[[i]]), ncol = dim(true_LambdaList[[i]])[2] - dim(est_LambdaList[[i]])[2]))
      frobenius.norm(est_Lambda_filled - true_LambdaList[[i]])
    } else if(dim(true_LambdaList[[i]])[2] < dim(est_LambdaList[[i]])[2]){
      true_Lambda_filled <- cbind(true_LambdaList[[i]], matrix(0, nrow = nrow(true_LambdaList[[i]]), ncol = dim(est_LambdaList[[i]])[2] - dim(true_LambdaList[[i]])[2]))
      frobenius.norm(est_LambdaList[[i]] - true_Lambda_filled)
    } else {
      frobenius.norm(est_LambdaList[[i]] - true_LambdaList[[i]])
    }
  })
  RV_Lambda_filled <- sapply(1:S, function(i) {
    if(dim(true_LambdaList[[i]])[2] > dim(est_LambdaList[[i]])[2]){
      est_Lambda_filled <- cbind(est_LambdaList[[i]], matrix(0, nrow = nrow(est_LambdaList[[i]]), ncol = dim(true_LambdaList[[i]])[2] - dim(est_LambdaList[[i]])[2]))
      MatrixCorrelation::RV(est_Lambda_filled, true_LambdaList[[i]], center = FALSE)
    } else if(dim(true_LambdaList[[i]])[2] < dim(est_LambdaList[[i]])[2]){
      true_Lambda_filled <- cbind(true_LambdaList[[i]], matrix(0, nrow = nrow(true_LambdaList[[i]]), ncol = dim(est_LambdaList[[i]])[2] - dim(true_LambdaList[[i]])[2]))
      MatrixCorrelation::RV(est_LambdaList[[i]], true_Lambda_filled, center = FALSE)
    } else {
      MatrixCorrelation::RV(est_LambdaList[[i]], true_LambdaList[[i]], center = FALSE)
    }
  })
  FN_Lambda_mean = mean(FN_Lambda, na.rm = TRUE) 
  RV_Lambda_mean = mean(RV_Lambda, na.rm = TRUE)
  FN_Lambda_mean_filled = mean(FN_Lambda_filled, na.rm = TRUE)
  RV_Lambda_mean_filled = mean(RV_Lambda_filled, na.rm = TRUE)
  # Mean of the metrics
  return(list(FN_Lambda_mean = FN_Lambda_mean, RV_Lambda_mean = RV_Lambda_mean,
              FN_Lambda_mean_filled = FN_Lambda_mean_filled, RV_Lambda_mean_filled = RV_Lambda_mean_filled))
}




measure_specific <- function(sim_data, post_output){
  S <- length(sim_data$SigmaLambdaList)
  # Obtain True study-specific covariance and loadings
  true_LambdaList <- sim_data$LambdaList
  true_SigmaLambdaList <- sim_data$SigmaLambdaList
  
  # Obtain estimated study-specific covariance and loadings
  est_LambdaList <- post_output$LambdaList
  est_SigmaLambdaList <- post_output$SigmaLambdaList
  
  # Metrics
  # ----------- Lambda ------------
  if(is.null(true_LambdaList) | is.null(est_LambdaList)){ # e.g. PFA does not have Lambda
    FN_Lambda_mean <- NA
    RV_Lambda_mean <- NA
    FN_Lambda_mean_filled <- NA
    RV_Lambda_mean_filled <- NA
  } else {
    measure_lambda <- measure_specific_unmatch_Lambda(est_LambdaList, true_LambdaList, S)
    FN_Lambda_mean <- measure_lambda$FN_Lambda_mean
    RV_Lambda_mean <- measure_lambda$RV_Lambda_mean
    FN_Lambda_mean_filled <- measure_lambda$FN_Lambda_mean_filled
    RV_Lambda_mean_filled <- measure_lambda$RV_Lambda_mean_filled
  }
  
  # ----------- SigmaLambda ------------
  if(!is.null(true_SigmaLambdaList) & !is.null(est_SigmaLambdaList)){ # e.g. MOMSS does not have SigmaLambda
    FN_SigmaLambda_mean <- sapply(1:S, function(s) {
      frobenius.norm(est_SigmaLambdaList[[s]] - true_SigmaLambdaList[[s]])
    }) %>% mean(na.rm = TRUE)
    RV_SigmaLambda_mean <- sapply(1:S, function(s) {
      MatrixCorrelation::RV(est_SigmaLambdaList[[s]], true_SigmaLambdaList[[s]],
                            center = FALSE)
    }) %>% mean(na.rm = TRUE)
    RV2_SigmaLambda_mean <- sapply(1:S, function(s) {
      RV_pd(est_SigmaLambdaList[[s]], true_SigmaLambdaList[[s]])
    }) %>% mean(na.rm = TRUE)
  } else {
    FN_SigmaLambda_mean <- NA
    RV_SigmaLambda_mean <- NA
    RV2_SigmaLambda_mean <- NA
  }
  
  return(list(FN_Lambda_mean = FN_Lambda_mean, RV_Lambda_mean = RV_Lambda_mean,
              FN_Lambda_mean_filled = FN_Lambda_mean_filled, RV_Lambda_mean_filled = RV_Lambda_mean_filled,
              FN_SigmaLambda_mean = FN_SigmaLambda_mean, RV_SigmaLambda_mean = RV_SigmaLambda_mean,
              RV2_SigmaLambda_mean = RV2_SigmaLambda_mean))
}



measure_marginal <- function(sim_data, post_output){
  # Obtain True marginal covariance matrices
  true_SigmaMarginal <- sim_data$SigmaMarginal
  # Obtain estimated marginal covariance matrices
  est_SigmaMarginal <- post_output$SigmaMarginal
  
  # Metrics
  S <- length(true_SigmaMarginal)
  FN_SigmaMarginal_mean <- sapply(1:S, function(s) {
    frobenius.norm(est_SigmaMarginal[[s]] - true_SigmaMarginal[[s]])
  }) %>% mean(na.rm = TRUE)
  RV_SigmaMarginal_mean <- sapply(1:S, function(s) {
    MatrixCorrelation::RV(est_SigmaMarginal[[s]], true_SigmaMarginal[[s]],
                          center = FALSE)
  }) %>% mean(na.rm = TRUE)
  RV2_SigmaMarginal_mean <- sapply(1:S, function(s) {
    RV_pd(est_SigmaMarginal[[s]], true_SigmaMarginal[[s]])
  }) %>% mean(na.rm = TRUE)
  
  return(list(FN_SigmaMarginal_mean = FN_SigmaMarginal_mean,
              RV_SigmaMarginal_mean = RV_SigmaMarginal_mean,
              RV2_SigmaMarginal_mean = RV2_SigmaMarginal_mean))
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