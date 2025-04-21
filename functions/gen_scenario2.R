library(Matrix) # for bdiag function
library(MASS) # for mvrnorm function
library(matlab)# for zeros function
library(tidyverse) # for %>% operator
# gen data with known covariates, and without study-specific factors
gen_scenario2 <- function(S, N_s, P, K, J_s, 
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

  # study-specific variables and parameters
  Lambda_list <- M_list <- Psi_list <- Y_list <- Sigma_list <- SigmaLambda_list<-  list()
  for(s in 1:S){
      Lam <- as.vector(zeros(P, J_s[s]))
      notZERO_count <- P * J_s[s] * (1 - sparsity)
      notZERO_value <- runif(notZERO_count, 0.6, 1)
      sign <- sample(x = length(notZERO_value), 
                     size = (length(notZERO_value) / 2))# Randomly assign negative sign
      notZERO_value[sign] <- notZERO_value[sign] * (-1)
      position_notZERO <- sample(x = J_s[s] * P, 
                                     size = length(notZERO_value))
      Lam[position_notZERO] <- notZERO_value
      Lambda_list[[s]] <- matrix(Lam, P, J_s[s])

    Psi_list[[s]] <- diag(runif(P, 0, 1), P)
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
              Phi=Phi, SigmaPhi = tcrossprod(Phi), 
              LambdaList=Lambda_list,
              SigmaLambdaList = SigmaLambda_list,
              SigmaMarginal = Sigma_list,
              Psi_list=Psi_list, 
              Psi_mat=Psi_mat, big_T=T_mat))
}

# d <- gen_scenario2(4, N_s=c(35, 35, 35, 35), P=100, K=6, J_s=c(2, 1, 1, 1))
# plot_single(d$Phi) + ggtitle(TeX("$\\Phi$")) 
# grid.arrange(plot_single(d$LambdaList[[1]])+ ggtitle(TeX("$\\Lambda_1$"))+ theme(legend.position = "none"),
#              plot_single(d$LambdaList[[2]])+ ggtitle(TeX("$\\Lambda_2$"))+ theme(legend.position = "none"),
#              plot_single(d$LambdaList[[3]])+ ggtitle(TeX("$\\Lambda_3$"))+ theme(legend.position = "none"),
#              plot_single(d$LambdaList[[4]])+ ggtitle(TeX("$\\Lambda_4$"))+ theme(legend.position = "none"),
#              nrow = 1)
# grid.arrange(plot_single(d$Psi_list[[1]]) + ggtitle(TeX("$\\Psi_1$"))+ theme(legend.position = "none", axis.text.x  = element_blank(), axis.text.y  = element_blank()),
#              plot_single(d$Psi_list[[2]]) + ggtitle(TeX("$\\Psi_2$"))+ theme(legend.position = "none", axis.text.x  = element_blank(), axis.text.y  = element_blank()),
#              plot_single(d$Psi_list[[3]]) + ggtitle(TeX("$\\Psi_3$"))+ theme(legend.position = "none", axis.text.x  = element_blank(), axis.text.y  = element_blank()),
#              plot_single(d$Psi_list[[4]]) + ggtitle(TeX("$\\Psi_4$"))+ theme(legend.position = "none", axis.text.x  = element_blank(), axis.text.y  = element_blank()),
#              nrow = 2)
# plot_single(d$SigmaPhi) + ggtitle(TeX("$\\Sigma_{\\Phi}$")) + theme(axis.text.x  = element_blank())

