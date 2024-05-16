library(Matrix)
library(MASS)
library(matlab)

# gen data with known covariates, and without study-specific factors
gen_senerioBMSFA <- function(S, N, P, K, j_s=c(1,1,1,1), 
                             genPhi = "dense", genLambda = "dense"){
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
  
  F_matrix <- rmvnorm(N, numeric(K), diag(K)) 
  Psi <- matrix(0.2 * (1:S), nrow = P, ncol = S, byrow = TRUE)
  
  # study-specific variables and parameters
  L_list <- Lambda_list <- A_list <- E_list <- list()
  for(s in 1:S){
    L_list[[s]] <- matrix(runif(n_s[s] * j_s[s], -1, 1), n_s[s], j_s[s])
    
    if (genLambda == "dense") {
      grid <- seq(-1, 1,length.out = P)
      Lambda_list[[s]] <- matrix(grid, nrow = P, ncol = j_s[s])
      rate<-trunc(P/(j_s[s]*2))
      if(j_s[s] > 1){
        for(js in 2:j_s[s]){
          (Lambda_list[[s]])[,js]<-grid[c((js*rate):P, 1:(js*rate-1))]
        }
      }
    } else if (genLambda == "sparse") {
      Lam <- as.vector(zeros(P, j_s[s]))
      noZERO_count <- P * j_s[s] * 0.3
      noZERO_value <- runif(noZERO_count, 0.6, 1)
      sign <- sample(x = length(noZERO_value), 
                     size = (length(noZERO_value) / 2))# Randomly assign negative sign
      noZERO_value[sign] <- noZERO_value[sign] * (-1)
      position_noZERO <- sample(x = j_s[s] * P, 
                                     size = length(noZERO_value))
      Lam[position_noZERO] <- noZERO_value[[s]]
      Lambda_list[[s]] <- matrix(Lam, P, j_s[s])
    }
    A_list[[s]] <- matrix(1, nrow = n_s[s], ncol = 1)
    E_list[[s]] <- mvrnorm(n_s[s], rep(0, P), diag(Psi[, s], nrow = P))
  }
  
  A <- bdiag(A_list) %>% as.matrix()
  E <- do.call(rbind, E_list)
  L <- bdiag(L_list) %>% as.matrix()
  Lambda <- do.call(cbind, Lambda_list)
  Y_mat <- F_matrix %*% t(Phi) + L %*% t(Lambda) + E
  
  # Reorgainze the data for BMSFA
  Y_list <- list()
  for(s in 1:S){
    Y_list[[s]] <- Y_mat[which(A[,s]==1),]
  }
  return(list(Y_mat=Y_mat, Y_list=Y_list, A=A, Phi=Phi, Psi=Psi))
}



