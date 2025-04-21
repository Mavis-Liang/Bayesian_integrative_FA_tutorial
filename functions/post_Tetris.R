source("./functions/run_Tetris.R")# Contains some user defined functions
post_Tetris <- function(fit){
  # Estimated common covariance
  A <- fit$A[[1]]
  Lambda <- getLambda(fit,A)
  S <- dim(A)[1]
  est_Phi <- as.matrix(Lambda[,colSums(A)==S])
  est_SigmaPhi <- tcrossprod(est_Phi)
  # Estimated study-specific covariance
  P = diag((colSums(A) == S)*1)
  T_s <- list()
  est_LambdaList <- list()
  for(s in 1:S){
    T_s[[s]] <- diag(A[s,])
    Lambda_s <- Lambda %*% (T_s[[s]] - P)
    Lambda_s <- Lambda_s[,-which(colSums(Lambda_s == 0) == nrow(Lambda_s))]
    Lambda_s <- matrix(Lambda_s, nrow=nrow(Lambda))
    est_LambdaList[[s]] <- Lambda_s}
  est_SigmaLambdaList <- lapply(1:S, function(s){
    tcrossprod(est_LambdaList[[s]])})
  
  # Estimated marginal covariance
  Psi <- list()
  est_SigmaMarginal <- lapply(1:S, function(s){
    Psi[[s]] <- diag(Reduce("+", fit$Psi[[s]])/length(fit$Psi[[s]])) # Psi is not returned. Modify if needed.
    Sigma_s <- Lambda %*% T_s[[s]] %*% t(Lambda) + Psi[[s]]
    })
  
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi,
              LambdaList = est_LambdaList, SigmaLambdaList = est_SigmaLambdaList,
              Psi = Psi, T_s = T_s,
              SigmaMarginal = est_SigmaMarginal))
}

