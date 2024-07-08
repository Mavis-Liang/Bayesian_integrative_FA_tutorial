post_BMSFA <- function(result_BMSFA){
  ## Assuming the number of common factors and specific factors are known
  est_Phi <- sp_OP(result_BMSFA$Phi, trace=FALSE)$Phi
  est_Lambda <- lapply(result_BMSFA$Lambda, function(x) sp_OP(x, trace=FALSE)$Phi)
  est_Psi <- lapply(result_BMSFA$psi, function(x) sp_OP(x, trace=FALSE)$Phi)# Without diagonalizing
  return(list(est_Phi = est_Phi, est_Lambda = est_Lambda, est_Psi = est_Psi))
}
