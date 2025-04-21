library(SUFA)
#-----------------------Post processing-----------------------
post_SUFA <- function(fit){
  all <- dim(fit$Lambda)[3]
  burnin <- floor(all * 0.8) # We will use the last 20% samples
  # shared and study-specific loading matrices
  loadings <- lam.est.all(fit, burn = burnin)
  # Obtain common covariance matrix and loading from fitting
  est_Phi <- loadings$Shared
  est_SigmaPhi <- SUFA_shared_covmat(fit, burn = burnin)
  est_Psi <- diag(colMeans(fit$residuals))
  # Study-specific loadings
  est_LambdaList <- loadings$Study_specific
  
  # Obtain study-specific covariance matrices
  S <- length(fit$A)
  marginal_cov <- sufa_marginal_covs(fit, burn = burnin)
  est_SigmaLambdaList <- list()
  for (s in 1:S) {
    est_SigmaLambdaList[[s]] <- marginal_cov[,,s] - est_SigmaPhi
  }
  
  return(list(SigmaPhi = est_SigmaPhi, Phi = est_Phi, 
              SigmaLambdaList = est_SigmaLambdaList,
              LambdaList = est_LambdaList, 
              Psi = est_Psi,
              SigmaMarginal = lapply(1:S, function(s) marginal_cov[,,s])
              ))
}
