post_MOMSS <- function(fit, version = 1){
  est_Phi <- fit$M
  if (version==2){est_Phi <- fit$Mpost}
  est_SigmaPhi <- tcrossprod(est_Phi)
  
  # Marginal covariance
  S <- dim(fit$sigma)[2]
  est_PsiList <- est_SigmaMarginal <-  list()
  for(s in 1:S){
    est_PsiList[[s]] <- fit$sigma[,s]
    est_SigmaMarginal[[s]] <- est_SigmaPhi + diag(fit$sigma[,s])
  }
  # last S columns of fit$Theta are the study-specific intercepts
  est_alphas <- fit$Theta[, (dim(fit$Theta)[2]-S+1):dim(fit$Theta)[2]]
  # The rest are coeficients for the known covariates
  est_B <- fit$Theta[, 1:(dim(fit$Theta)[2]-S)]
  
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi, Psi = est_PsiList, alpha = est_alphas, B = est_B,
              SigmaMarginal = est_SigmaMarginal))
}

