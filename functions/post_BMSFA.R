library(MSFA)
library(tidyverse)
post_BMSFA <- function(fit){
  # Common covariance matrix and loading
  est_Phi <- sp_OP(fit$Phi, trace=FALSE)$Phi
  est_SigmaPhi <- tcrossprod(est_Phi)
  
  # Study-specific covariance matrices and loadings
  est_LambdaList <- lapply(fit$Lambda, function(x) sp_OP(x, trace=FALSE)$Phi)
  est_SigmaLambdaList <- lapply(est_LambdaList, function(x) tcrossprod(x))
  
  # Marginal covariance matrices
  S <- length(est_SigmaLambdaList)
  # Get point estimate of each Psi_s
  est_PsiList <- lapply(1:S, function(s) {
    apply(fit$psi[[s]], c(1, 2), mean)
  })
  est_margin_cov <- lapply(1:S, function(s) {
    est_SigmaPhi + est_SigmaLambdaList[[s]] + diag(est_PsiList[[s]] %>% as.vector())
  })
  
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi,
         LambdaList = est_LambdaList, SigmaLambdaList = est_SigmaLambdaList,
         PsiList = est_PsiList,
         SigmaMarginal = est_margin_cov))
}



