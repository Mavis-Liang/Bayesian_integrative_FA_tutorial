library(Matrix)
library(MASS)
library(matlab)
library(mvtnorm)

# gen data with known covariates, and without study-specific factors
gen_scenario1 <- function(S, N_s, P, Q, K, genPhi = "sparse"){
  # Number of total observations
  N <- sum(N_s)
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
    sparsity <- 0.4
    Phi_long <- as.vector(zeros(P, K))
    noZERO_count <- P * K * (1 - sparsity)
    noZEROs <- runif(noZERO_count, 0.6, 1)
    sign <- sample(x = length(noZEROs), 
                   size = (length(noZEROs) / 2))# Randomly assign negative sign
    noZEROs[sign] <- noZEROs[sign] * (-1)
    position_noZERO <- sample(x = K * P, size = length(noZEROs))
    Phi_long[position_noZERO] <- noZEROs
    Phi <- matrix(Phi_long, P, K)
  }

  F_matrix <- rmvnorm(N, numeric(K), diag(K))  
  Psi <- matrix(0.2 * (1:S), nrow = P, ncol = S, byrow = TRUE)
  Psi_list <- lapply(1:S, function(s){
    diag(Psi[, s], nrow = P)
  })
  
  # study-specific variables and parameters
  A_list <- E_list <-  list()
  for(s in 1:S){
    A_list[[s]] <- matrix(1, nrow = N_s[s], ncol = 1)
    E_list[[s]] <- mvrnorm(N_s[s], rep(0, P), 
                           diag(Psi[, s], nrow = P))
  }
  E <- do.call(rbind, E_list)
  A <- bdiag(A_list) %>% as.matrix()
  
  # Generate matrix Y
  Y_mat <- A%*%t(alpha)+ X%*%t(Beta) + F_matrix%*%t(Phi)+ E
  
  # Reorgainze the data for other methods
  Y_list <- list()
  for(s in 1:S){
    Y_list[[s]] <- Y_mat[which(A[,s]==1),]
  }
  
  SigmaPhi <- tcrossprod(Phi)
  SigmaMarginal <- lapply(1:S, function(s){
    SigmaPhi + Psi_list[[s]]
  })
  
  return(list(Y_mat=Y_mat, Y_list=Y_list, N_s=N_s, A=A, 
              X=X, Beta=Beta, Phi=Phi, SigmaPhi = SigmaPhi, 
              SigmaMarginal=SigmaMarginal,
              Psi_list=Psi_list, Psi_mat=Psi))
}

