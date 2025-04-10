post_stackFA <- function(fit, S){
  est_Phi <- MSFA::sp_OP(fit$Lambda, trace=FALSE)$Phi
  est_SigmaPhi <- tcrossprod(est_Phi)
  est_SigmaMarginal <-  lapply(1:S, function(s)
    apply(fit$Sigma, c(1, 2), mean)
  )
  Psi_chain <- list()
  for(i in 1:dim(fit$Sigma)[3]){
    Psi_chain[[i]] <- fit$Sigma[, , i] - tcrossprod(fit$Lambda[, , i])
  }
  est_Psi <- Reduce('+', Psi_chain)/length(Psi_chain)
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi, Psi = est_Psi,
              SigmaMarginal = est_SigmaMarginal))
}