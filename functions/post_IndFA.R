post_IndFA <- function(fit_list){
  # Estimated study-specific covariance and loading
  S <- length(fit_list)
  est_LambdaList <- lapply(1:S, function(s){
    MSFA::sp_OP(fit_list[[s]]$Lambda, trace=FALSE)$Phi
  })
  est_SigmaLambdaList <- lapply(est_LambdaList, function(x) tcrossprod(x))
  
  # Marginal covariance matrices
  est_SigmaMarginal <- lapply(1:S, function(s) {
    fit <- fit_list[[s]]
    apply(fit$Sigma, c(1, 2), mean)
  })
  
  return(list(LambdaList = est_LambdaList, SigmaLambdaList = est_SigmaLambdaList,
              SigmaMarginal = est_SigmaMarginal))
}
