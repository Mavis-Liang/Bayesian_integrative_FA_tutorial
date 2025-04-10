if (!require("caret", quietly = TRUE)) {
  install.packages("caret")
}
if (!require("MASS", quietly = TRUE)) {
  install.packages("MASS")
}
if (!require("matlab", quietly = TRUE)) {
  install.packages("matlab")
}
if (!require("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
# [Scenario 5 in the manuscript][Mimicing the ovarian cancer gene expression data, model based on Tetris]
gen_scenario5 <- function(S = 4, N_s = c(157, 195, 285, 117), P = 1060) {
  N <- sum(N_s)
  # Define "big-T" matrix (K = 15, J_s = c(2, 2, 2, 2), partial factors = 2 + 1)
  K <- 15 + 4 * 2 + 2 + 1
  big_T <- matrix(0, nrow = S, ncol = K)
  # define common factors
  big_T[, 1:15] <- 1
  # define partially-shared factors
  late <- c(1,2,4)
  Affymetrix <- c(2, 3, 4)
  big_T[late,16] <- big_T[late, 17] <- big_T[Affymetrix, 18]  <- 1
  # define study-specific factors
  big_T[1, 19:20] <- big_T[2, 21:22] <- big_T[3, 23:24] <- big_T[4, 25:26] <- 1
  
  # Define T_s
  T_s <- lapply(1:S, function(s) {
    T_s <- diag(big_T[s, ])
    T_s
  })
  
  
  # Loading matrix
  sparsity <- 0.8
  Phi_long <- as.vector(matlab::zeros(P, K))
  noZERO_count <- P * K * (1 - sparsity)
  noZEROs <- runif(noZERO_count, 0.2, 0.6)
  sign <- sample(x = length(noZEROs), 
                 size = (length(noZEROs) / 2))# Randomly assign negative sign
  noZEROs[sign] <- noZEROs[sign] * (-1)
  position_noZERO <- sample(x = K * P, size = length(noZEROs))
  Phi_long[position_noZERO] <- noZEROs
  Phi_star <- matrix(Phi_long, P, K)
  
  # Generate Sigma_s
  Psi_list <- list()
  Psi_list <- lapply(1:S, function(s) {
    Psi_s <- diag(runif(P, 0, 0.01), P)
    Psi_s
  })
  Sigma_list <- list()
  Sigma_list <- lapply(1:S, function(s) {
    Sigma_s <- Phi_star %*% T_s[[s]] %*% t(Phi_star) + Psi_list[[s]]
    Sigma_s
  })
  
  # Generate Y_s
  Y_list <- lapply(1:S, function(s) {
    Sigma_s <- Sigma_list[[s]]
    Y_s <- MASS::mvrnorm(n = N_s[[s]], mu = rep(0, P), Sigma = Sigma_s)
  })
  
  Y_mat <- do.call(rbind, Y_list)
  M_list <- lapply(1:S, function(s) {
    matrix(1, nrow = N_s[s], ncol = 1)
  })
  M <- Matrix::bdiag(M_list) %>% as.matrix()
  Phi <- as.matrix(Phi_star[,colSums(big_T)==S])
  SigmaPhi = tcrossprod(Phi)
  P = diag((colSums(big_T) == S)*1)
  LambdaList <- lapply(1:S, function(s){
    T_s <- diag(big_T[s,])
    Lambda_s <- Phi_star %*% (T_s - P)
    Lambda_s <- Lambda_s[,-which(colSums(Lambda_s == 0) == nrow(Lambda_s))]
    Lambda_s <- matrix(Lambda_s, nrow=nrow(Phi_star))}
  )
  SigmaLambdaList <- lapply(1:S, function(s){
    tcrossprod(LambdaList[[s]])})
  Psi_mat <- lapply(Psi_list,diag) %>% do.call(cbind, .)
  
  return(list(Y_mat=Y_mat, Y_list=Y_list, N_s=N_s, M=M, 
              Phi=Phi, SigmaPhi = tcrossprod(Phi), 
              LambdaList=LambdaList,
              SigmaLambdaList = SigmaLambdaList,
              SigmaMarginal = Sigma_list,
              Psi_list=Psi_list, 
              Psi_mat=Psi_mat, big_T=big_T))
}
