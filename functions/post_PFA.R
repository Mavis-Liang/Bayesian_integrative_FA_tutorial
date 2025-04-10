source("./FBPFA-PFA.R")## It's neccessary because it seems the "FBPFA-PFA with fixed latent dim.R" depends on this.
source("./FBPFA-PFA with fixed latent dim.R")
library(tidyverse)
post_PFA <- function(fit){
  p <- nrow(fit$Loading[[1]])
  k <- ncol(fit$Loading[[1]])
  npost <- length(fit$Loading)
  Q_list <- fit$Pertmat
  S <- dim(Q_list[[1]])[2]
  posteriorPhis <- array(0, dim = c(p, k, npost))
  posteriorLams <- list()
  
  #---estimated common loadings and shared variance---
  # Element-wise multiplication
  for(s in 1:S){
  posteriorLams[[s]] <- array(0, dim = c(p, k, 1000))
  for(i in 1:npost){
    posteriorPhis[,,i] <- fit$Loading[[i]] %*% diag(fit$Latentsigma[[i]])
    posteriorLams[[s]][,,i] <- (solve(matrix(Q_list[[i]][, s], p, p)) - diag(p)) %*% posteriorPhis[,,i]
  }
}

  # Varimax for loadings
  est_Phi <- MSFA::sp_OP(posteriorPhis, itermax = 10, trace = FALSE)$Phi
  est_speLoad <- lapply(posteriorLams, function(x) MSFA::sp_OP(x, itermax = 10, trace = FALSE)$Phi)

  
  # Estimated shared covariance, study-specific covariance matrix, and SigmaMarginal
  sharevar <- list()
  est_SigmaLambdaList <- list()
  est_SigmaMarginal <- list()
  est_Psi <- list()
  for(s in 1:S){ # Loop over each study
  post_SigmaLambda_s <- list()
  post_SigmaMarginal_s <- list()
  Psi <- list()
  for(i in 1:npost){ # Loop over each posterior sample
    sharevar[[i]] <- fit$Loading[[i]]%*%diag(fit$Latentsigma[[i]]^2)%*%t(fit$Loading[[i]]) + 
      diag(fit$Errorsigma[[i]]^2) # Get the shared variance
    Q_temp_inv <- solve(
      matrix(Q_list[[i]][, s], p, p)
      )
    post_SigmaMarginal_s[[i]] <- Q_temp_inv%*%sharevar[[i]]%*%t(Q_temp_inv)
    post_SigmaLambda_s[[i]] <- post_SigmaMarginal_s[[i]] - sharevar[[i]]
    Psi[[i]] <- diag(fit$Errorsigma[[i]]^2)
  }
  est_SigmaMarginal[[s]] <- Reduce('+', post_SigmaMarginal_s)/length(post_SigmaMarginal_s)
  est_SigmaLambdaList[[s]] <- Reduce('+', post_SigmaLambda_s)/length(post_SigmaLambda_s)
  }
  est_Psi <- Reduce('+', Psi)/length(Psi)
  est_SigmaPhi <- Reduce('+', sharevar)/length(sharevar)
  est_Q <- Reduce('+', Q_list)/length(Q_list)
  est_Q_list <- lapply(1:S, function(s) matrix(est_Q[, s], p, p))
  
  # Return the results
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi, Psi = est_Psi, Q = est_Q_list, LambdaList=est_speLoad,
              SigmaLambdaList = est_SigmaLambdaList,
              SigmaMarginal = est_SigmaMarginal))
}

