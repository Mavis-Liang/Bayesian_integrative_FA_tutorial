library(SUFA)
post_SUFA <- function(fit){
  est_Phi <- lam.est(fit$Lambda,burn = burnin)
  est_SigmaPhi <- SUFA_shared_covmat(fit,burn = burnin)
  return(list(est_Phi = est_Phi, est_SigmaPhi = est_SigmaPhi))
}