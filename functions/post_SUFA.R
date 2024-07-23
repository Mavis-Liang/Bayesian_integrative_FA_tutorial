library(SUFA)
post_SUFA <- function(fit){
  burnin <- 500
  est_Phi <- lam.est(fit$Lambda,burn = burnin)
  est_SigmaPhi <- SUFA_shared_covmat(fit,burn = burnin)
  return(list(est_Phi = est_Phi, est_SigmaPhi = est_SigmaPhi))
}