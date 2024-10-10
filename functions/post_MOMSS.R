post_MOMSS <- function(fit){
  est_Phi <- fit$M
  est_SigmaPhi <- tcrossprod(est_Phi)
  
  # Marginal covariance
  S <- dim(fit$sigma)[2]
  est_PsiList <- est_SigmaMarginal <-  list()
  for(s in 1:S){
    est_PsiList[[s]] <- fit$sigma[,s]
    est_SigmaMarginal[[s]] <- est_SigmaPhi + diag(fit$sigma[,s])
  }
  
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi, 
              SigmaMarginal = est_SigmaMarginal))
}

