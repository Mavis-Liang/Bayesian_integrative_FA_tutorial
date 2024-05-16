calculateRV <- function(X, Y) {
  X <- X + matrix(rnorm(n = nrow(X) * ncol(X), mean = 0, sd = 0.0001), nrow = nrow(X), ncol = ncol(X))
  Y <- Y + matrix(rnorm(n = nrow(Y) * ncol(Y), mean = 0, sd = 0.0001), nrow = nrow(Y), ncol = ncol(Y))
  X <- scale(X)
  Y <- scale(Y)
  # Calculate the trace of the crossproducts
  trace_XY <- sum(diag(X %*% Y))
  trace_XX <- sum(diag(X %*% X))
  trace_YY <- sum(diag(Y %*% Y))
  # Calculate the RV coefficient
  RV_coefficient <- trace_XY / sqrt(trace_XX * trace_YY)
  return(RV_coefficient)
}