post_PFA <- function(fit){
  sigma2p <- Reduce('+', fit$Latentsigma[200:499])/length(fit$Latentsigma[200:499])
  lambdap <- (Reduce('+', fit$Loading[200:499])/length(fit$Latentsigma[200:499])) %*% diag(sigma2p)
  return(list(est_SigmaPhi=lambdap))
}
