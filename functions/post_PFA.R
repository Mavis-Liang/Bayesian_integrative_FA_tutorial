post_PFA <- function(fit){
  sigmap <- Reduce('+', fit$Latentsigma[200:499])/length(fit$Latentsigma[200:499])
  lambdap <- (Reduce('+', fit$Loading[200:499])/length(fit$Latentsigma[200:499])) %*% diag(sigmap)
  Sigma_e <- Reduce('+', fit$Errorsigma[200:499])/length(fit$Errorsigma[200:499])
  SigmaPhi <- lambdap %*% diag(sigmap^2) %*% t(lambdap) + diag(Sigma_e^2)
  return(list(est_Phi=lambdap, est_SigmaPhi=SigmaPhi))
}
