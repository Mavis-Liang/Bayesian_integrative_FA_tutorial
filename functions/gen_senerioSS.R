library(Matrix)
library(MASS)
library(matlab)

# gen data with known covariates, and without study-specific factors
gen_senerioSS <- function(S, N, P, Q, K, genPhi = "dense"){
  
  # Number of observations in each study
  n_s <- rmultinom(1, N, prob = rep(1/S, S))
  
  # common variables and parameters
  alpha <- rmvnorm(P, -(1:S), 0.2 * diag(1:S))
  X <- matrix(runif(N*Q, 0, 1),
              nrow = N, ncol = Q)
  Beta = rmvnorm(P, -(1:Q), diag(1:Q))
  
  # Generate Phi
  if(genPhi == "dense"){
    grid <- seq(-1, 1,length.out = P)
    Phi <- matrix(grid, nrow = P, ncol = K)
    rate<-trunc(P/(K*2))
    for(k in 2:K){
      Phi[,k]<-grid[c((k*rate):P, 1:(k*rate-1))]
    }
  } else if(genPhi == "sparse") {
    Ph <- as.vector(zeros(P, K))
    noZERO_count <- P * K * 0.3
    noZEROs <- runif(noZERO_count, 0.6, 1)
    sign <- sample(x = length(noZEROs), 
                   size = (length(noZEROs) / 2))# Randomly assign negative sign
    noZEROs[sign] <- noZEROs[sign] * (-1)
    position_noZERO <- sample(x = K * P, size = length(noZEROs))
    Ph[position_noZERO] <- noZEROs
    Phi <- matrix(Ph, P, K)
  }

  F_matrix <- rmvnorm(N, numeric(K), diag(K))  
  Psi <- matrix(0.2 * (1:S), nrow = P, ncol = S, byrow = TRUE)
  
  # study-specific variables and parameters
  A_list <- E_list <-  list()
  for(s in 1:S){
    A_list[[s]] <- matrix(1, nrow = n_s[s], ncol = 1)
    E_list[[s]] <- mvrnorm(n_s[s], rep(0, P), 
                           diag(Psi[, s], nrow = P))
  }
  E <- do.call(rbind, E_list)
  A <- bdiag(A_list) %>% as.matrix()
  
  # Generate matrix Y
  Y_mat <- A%*%t(alpha)+ X%*%t(Beta) + F_matrix%*%t(Phi)+ E
  
  # Reorgainze the data for BMSFA
  Y_list <- list()
  for(s in 1:S){
    Y_list[[s]] <- Y_mat[which(A[,s]==1),]
  }
  
  return(list(Y_mat=Y_mat, Y_list=Y_list, A=A, X=X, Beta=Beta, Phi=Phi, Psi=Psi))
}

