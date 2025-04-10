library(Matrix) # for bdiag function
library(MASS) # for mvrnorm function
library(matlab)# for zeros function
library(tidyverse) # for %>% operator

# [Scenario 3 in the manuscript][Based on SUFA]
gen_scenario6 <- function(S, N_s, P, K, J_s, 
                             genPhi = "sparse", sparsity = 0.4){
  # Total number of observations
  N <- sum(N_s)

  # common variables and parameters
  # Generate Phi
  if(genPhi == "dense"){# Lack of randomness
    grid <- seq(-1, 1,length.out = P)
    Phi <- matrix(grid, nrow = P, ncol = K)
    rate<-trunc(P/(K*2))
    for(k in 2:K){
      Phi[,k]<-grid[c((k*rate):P, 1:(k*rate-1))]
    }
  } else if(genPhi == "sparse") {
    Phi_long <- as.vector(zeros(P, K))
    noZERO_count <- P * K * (1 - sparsity)
    noZEROs <- runif(noZERO_count, 0.6, 1)
    sign <- sample(x = length(noZEROs), 
                   size = (length(noZEROs) / 2))# Randomly assign negative sign
    noZEROs[sign] <- noZEROs[sign] * (-1)
    position_noZERO <- sample(x = K * P, size = length(noZEROs))
    Phi_long[position_noZERO] <- noZEROs
    Phi <- matrix(Phi_long, P, K)
  } else if(genPhi == "sparse_block") {
    # Create Phi that nonzero values appears in blocks
    per_blocks <- 5 # must be an integer and can be divided by P
    block_count <- P * K / per_blocks # Total blocks in Phi
    nonzero_block_counts <- P * K * (1-sparsity) / per_blocks # Number of nonzero blocks
    nonzero_blocks <- sample(1:block_count, nonzero_block_counts, replace = FALSE) # Randomly choose nonzero blocks
    negative_block <- sample(c(-1, 1), nonzero_block_counts, replace = TRUE) # Randomly assign negative sign to nonzero blocks
    nonzeros <- runif(nonzero_block_counts * per_blocks, 0.6, 1) # Generate nonzero values
    Phi_long <- as.vector(zeros(P, K)) # Initialize Phi: columns of Phi are stacked
    for(i in 1:nonzero_block_counts){
      begin <- (nonzero_blocks[i]-1)*per_blocks+1 # Begin index of nonzero blocks in Phi
      end <- nonzero_blocks[i]*per_blocks # End index of nonzero blocks in Phi
      begin_nonzero <- (i-1)*per_blocks+1 # Begin index of nonzero values
      end_nonzero <- i*per_blocks # End index of nonzero values
      Phi_long[begin:end] <- nonzeros[begin_nonzero:end_nonzero] * negative_block[i]
    }
    Phi <- matrix(Phi_long, P, K)
  }

  Psi <- diag(runif(P, 0, 1), P)

  # study-specific variables and parameters
  Lambda_list <- M_list  <- A_s_list <- Psi_list <- Y_list <- Sigma_list <- SigmaLambda_list<-  list()
  for(s in 1:S){
    A_s_list[[s]] <- matrix (rnorm(K*J_s[s],sd=0.4),  nrow=K,ncol=J_s[s])
    Lambda_list[[s]] <- Phi %*% A_s_list[[s]]

    Psi_list[[s]] <- Psi
    # Covariance for the marginal distribution of Y
    SigmaLambda_list[[s]] <- tcrossprod(Lambda_list[[s]])
    Sigma_list[[s]] <- tcrossprod(Phi)  + SigmaLambda_list[[s]]  + Psi_list[[s]]
    Y_list[[s]] <- mvrnorm(N_s[s], rep(0, length = P), Sigma_list[[s]])
    M_list[[s]] <- matrix(1, nrow = N_s[s], ncol = 1)
  }
  
  M <- bdiag(M_list) %>% as.matrix()

  # Store Y_list in Y_mat
  Y_mat <- do.call(rbind, Y_list)
  # Store Psi_list in column-wise Psi matrix
  Psi_mat <- lapply(Psi_list,diag) %>% do.call(cbind, .)
  # Create T in Tetris
  T_mat <- matrix(0, nrow = S, ncol = K + sum(J_s))
  T_mat[,1:K] <- 1
  T_mat[1, (K + 1) :  (K + J_s[1])] <- 1
  for(s in 2:S){
      T_mat[s, (K + sum(J_s[1:s-1]) + 1) :  (K + sum(J_s[1:s-1]) + J_s[s])] <- 1
  }
  return(list(Y_mat=Y_mat, Y_list=Y_list, N_s=N_s, M=M, 
              Phi=Phi, SigmaPhi = tcrossprod(Phi) + Psi, ## Use the definition in SUFA!
              LambdaList=Lambda_list,
              SigmaLambdaList = SigmaLambda_list,
              SigmaMarginal = Sigma_list,
              Psi_list=Psi_list, 
              Psi_mat=Psi_mat, big_T=T_mat))
}