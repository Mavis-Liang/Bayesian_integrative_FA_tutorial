# This code is adapted from Roy's code https://github.com/royarkaprava/Perturbed-factor-model/blob/master/LRFpertgrp1.R.
library(Matrix) # for bdiag function
library(MASS) # for mvrnorm function
library(matlab)# for zeros function
library(tidyverse) # for %>% operator
gen_scenario3 <- function(S, N_s, P, K, sparsity = 0.4){
  # Number of observations in each study
  # data generation
  N <- sum(N_s)
  S <- S
  # Number of observations in each study
  grpind = rep(1:length(N_s), times = N_s)
  
  
  eta0 <- matrix(rnorm(K*N), K, N)# standard normal
  
  # generate common loading
  Phi_long <- as.vector(zeros(P, K))
  noZERO_count <- P * K * (1 - sparsity)
  noZEROs <- runif(noZERO_count, 0.6, 1)
  sign <- sample(x = length(noZEROs), 
                 size = (length(noZEROs) / 2))# Randomly assign negative sign
  noZEROs[sign] <- noZEROs[sign] * (-1)
  position_noZERO <- sample(x = K * P, size = length(noZEROs))
  Phi_long[position_noZERO] <- noZEROs
  Phi <- matrix(Phi_long, P, K)

  Y <- matrix(rnorm(P*N, mean = Phi %*% eta0, sd = 0.5), P, N)
  Qlist <- matrix(array(diag(P)), P^2, S)
  po <- 1
  Qmean <- array(1*diag(P))
  
  #####################Generating the perturbation matrices################
  for(i in 2:S){
    Qlist[, i] <- array(matrix(rnorm(P^2, Qmean, sd = sqrt(0.01)), P, P))
  }
  
  QYpr <- function(i, mat = Y){
    temp <- matrix(Qlist[, grpind[i]], P, P)
    return(temp%*%mat[, i])
  }
  
  Y <- parallel::mcmapply(1:N, FUN = QYpr, MoreArgs = list(mat=Y)) 
  
  # Summarize the data for the use of other methods
  M_list <- list()
  for(s in 1:S){
    M_list[[s]] <- matrix(1, nrow = N_s[s], ncol = 1)
  }
  M <- as.matrix(bdiag(M_list))
  
  # Common covariance, study-specific covariance, Y in list format
  common_cov = Phi %*% t(Phi) + diag(P)
  Y_list <- spec_cov <- marginal_cov <-  spec_load <- list()
  for(s in 1:S){
    Y_list[[s]] <- t(Y)[which(M[,s]==1),]
    #Q_s_inv <- solve(matrix(Qlist[, s], P, P))
    Q_s <- matrix(Qlist[, s], P, P)
    Q_s_inv = solve(Q_s)
    spec_load[[s]] = (Q_s_inv - diag(P)) %*% common_cov
    marginal_cov[[s]] <- Q_s%*%common_cov%*%t(Q_s)
    spec_cov[[s]] <- marginal_cov[[s]] - common_cov
    }
  
  return(list(Y_mat=t(Y), Y_list = Y_list, grpind=grpind, M = M,N_s = N_s, Q_s=Q_s, LambdaList = spec_load,
              Phi = Phi, SigmaPhi = common_cov, 
              SigmaMarginal = marginal_cov,
              SigmaLambdaList = spec_cov))
}



