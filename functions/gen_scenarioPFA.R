gen_scenarioPFA <- function(S, N, P, K){
  # Number of observations in each study
  # data generation
  n <- N
  p  <- P
  r  <- K
  grp <- S
  # Number of observations in each study
  n_s <- rmultinom(1, N, prob = rep(1/S, S))
  grpind = rep(1:length(n_s), times = n_s)
  
  
  eta0 <- matrix(rnorm(r*n), r, n)# standard normal
  
  # generate common loading BMSFA's code (a smaller loading than PFA)
  grid <- seq(-1, 1,length.out = p)
  lambda0 <- matrix(grid, nrow = p, ncol = r)
  rate<-trunc(p/(r*2))
  for(k in 2:r){
    lambda0[,k]<-grid[c((k*rate):p, 1:(k*rate-1))]
  }

  Y <- matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)#matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n) #t(rmvnorm(n, sigma = solve(pdmat)))#
  #Plambda0 <- diag(p) - lambda0 %*% solve(crossprod(lambda0)) %*% t(lambda0)
  Qlist <- matrix(array(diag(p)), p^2, grp)
  po <- 1
  Qmean <- array(1*diag(p))
  
  #####################Generating the perturbation matrices################
  for(i in 2:grp){
    Qlist[, i] <- array(matrix(rnorm(p^2, Qmean, sd = sqrt(0.01)), p, p))
  }
  
  QYpr <- function(i, mat = Y){
    temp <- matrix(Qlist[, grpind[i]], p, p)
    return(temp%*%mat[, i])
  }
  
  Y <- parallel::mcmapply(1:n, FUN = QYpr, MoreArgs = list(mat=Y)) #matrix(Q %*% array(Y), p, n)
  
  # Summarize the data for the use of other methods
  A_list <- list()
  for(s in 1:S){
    A_list[[s]] <- matrix(1, nrow = n_s[s], ncol = 1)
  }
  A <- as.matrix(bdiag(A_list))
  
  # Common covariance, study-specific covariance, Y in list format
  common_cov = lambda0 %*% t(lambda0) + diag(p)
  Y_list <- spec_cov <- marginal_cov <-  list()
  for(s in 1:S){
    Y_list[[s]] <- t(Y)[which(A[,s]==1),]
    #Q_s_inv <- solve(matrix(Qlist[, s], p, p))
    Q_s <- matrix(Qlist[, s], p, p)
    marginal_cov[[s]] <- Q_s%*%common_cov%*%t(Q_s)
    spec_cov[[s]] <- marginal_cov[[s]] - common_cov
    }
  
  return(list(Y_mat=t(Y), Y_list = Y_list, grpind=grpind, A = A,n_s = n_s,
              Phi = lambda0, SigmaPhi = common_cov, 
              SigmaMarginal = marginal_cov,
              SigmaLambdaList = spec_cov))
}




