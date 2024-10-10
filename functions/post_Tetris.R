source("Tetris.R")
post_Tetris <- function(fit){
  # Estimated common covariance
  A <- fit$A[[1]]
  Lambda <- getLambda(fit,A)
  S <- dim(A)[1]
  est_Phi <- Lambda[,colSums(A)==S]
  est_SigmaPhi <- tcrossprod(est_Phi)
  # Estimated study-specific covariance
  A_common = diag((colSums(A) == S)*1)
  est_LambdaList <- lapply(1:S, function(s){
    T_s <- diag(A[s,])
    Lambda_s <- Lambda %*% (T_s - A_common)
    Lambda_s <- Lambda_s[,-which(colSums(Lambda_s == 0) == nrow(Lambda_s))]
    Lambda_s <- matrix(Lambda_s, nrow=nrow(Lambda))}
    )
  est_SigmaLambdaList <- lapply(1:S, function(s){
    tcrossprod(est_LambdaList[[s]])})
  
  # Estimated marginal covariance
  est_SigmaMarginal <- lapply(1:S, function(s){
    T_s <- diag(A[s,])
    Psi_s <- Reduce("+", fit$Psi[[s]])/length(fit$Psi[[s]])
    Sigma_s <- Lambda %*% T_s %*% t(Lambda) + diag(Psi_s)
    })
  
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi,
              LambdaList = est_LambdaList, SigmaLambdaList = est_SigmaLambdaList,
              SigmaMarginal = est_SigmaMarginal))
}


