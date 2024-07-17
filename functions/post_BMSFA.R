post_BMSFA <- function(fit){
  ## Assuming the number of common factors and specific factors are known
  est_Phi <- sp_OP(fit$Phi, trace=FALSE)$Phi
  est_Lambda <- lapply(fit$Lambda, function(x) sp_OP(x, trace=FALSE)$Phi)
  est_Psi <- lapply(fit$psi, function(x) rowMeans(x))
  return(list(est_Phi = est_Phi, est_Lambda = est_Lambda, est_Psi = est_Psi))
}
