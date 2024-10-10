post_stackFA <- function(fit, S){
  est_Phi <- MSFA::sp_OP(fit$Lambda, trace=FALSE)$Phi
  est_SigmaPhi <- tcrossprod(est_Phi)
  est_SigmaMarginal <-  lapply(1:S, function(s)
    apply(fit$Sigma, c(1, 2), mean)
  )
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi,
              SigmaMarginal = est_SigmaMarginal))
}