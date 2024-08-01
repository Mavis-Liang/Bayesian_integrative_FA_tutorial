library(Matrix)
library(MASS)
library(matlab)

# gen data with known covariates, and without study-specific factors
gen_senerioBMSFA <- function(S, N, P, K, j_s=c(1,1,1,1), 
                             genPhi = "dense"){
  # Number of observations in each study
  n_s <- rmultinom(1, N, prob = rep(1/S, S))
  
  # common variables and parameters
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

  # study-specific variables and parameters
  Lambda_list <- A_list <- Psi_list <- Y_list <- Sigma_list <- SigmaLambda_list<-  list()
  for(s in 1:S){
    # Study-specific loading matrix
      Lam <- as.vector(zeros(P, j_s[s]))
      notZERO_count <- P * j_s[s] * 0.3
      notZERO_value <- runif(notZERO_count, 0.6, 1)
      sign <- sample(x = length(notZERO_value), 
                     size = (length(notZERO_value) / 2))# Randomly assign negative sign
      notZERO_value[sign] <- notZERO_value[sign] * (-1)
      position_notZERO <- sample(x = j_s[s] * P, 
                                     size = length(notZERO_value))
      Lam[position_notZERO] <- notZERO_value
      Lambda_list[[s]] <- matrix(Lam, P, j_s[s])

    Psi_list[[s]] <- diag(runif(P, 0, 1), P)
    # Covariance for the marginal distribution of Y
    SigmaLambda_list[[s]] <- tcrossprod(Lambda_list[[s]])
    Sigma_list[[s]] <- tcrossprod(Phi)  + SigmaLambda_list[[s]]  + Psi_list[[s]]
    Y_list[[s]] <- mvrnorm(n_s[s], rep(0, length = P), Sigma_list[[s]])
    A_list[[s]] <- matrix(1, nrow = n_s[s], ncol = 1)
  }
  
  A <- bdiag(A_list) %>% as.matrix()

  # Store Y_list in Y_mat
  Y_mat <- do.call(rbind, Y_list)
  # Store Psi_list in column-wise Psi matrix
  Psi_mat <- lapply(Psi_list,diag) %>% do.call(cbind, .)
  return(list(Y_mat=Y_mat, Y_list=Y_list, n_s=n_s, A=A, 
              Phi=Phi, SigmaPhi = tcrossprod(Phi), 
              SigmaLambdaStack = do.call(rbind, SigmaLambda_list),
              SigmaLambda_list = SigmaLambda_list,
              Psi_list=Psi_list, 
              Psi_mat=Psi_mat))
}



