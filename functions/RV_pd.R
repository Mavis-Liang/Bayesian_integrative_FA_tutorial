# RV when X, y are semi-positive definite matrices
RV_pd <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  # Calculate the trace of the cross-products
  trace_XY <- sum(diag(t(X) %*% Y))
  trace_XX <- sum(diag(t(X) %*% X))
  trace_YY <- sum(diag(t(Y) %*% Y))
  # Calculate the RV coefficient
  RV_coefficient <- trace_XY / sqrt(trace_XX * trace_YY)
  return(RV_coefficient)
}

# # # Example
# true_Phi <- sim_data_test$Phi
# est_Phi <- (result_test$result_MOMSS)$M
# 
# # RV of the loadings
# RV(true_Phi, est_Phi)
# calculateCovRV(tcrossprod(true_Phi), tcrossprod(est_Phi))
# 
# # RV of covariances
# RV(tcrossprod(true_Phi), tcrossprod(est_Phi))

