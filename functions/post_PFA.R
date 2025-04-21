source("./FBPFA-PFA.R")## It's neccessary because it seems the "FBPFA-PFA with fixed latent dim.R" depends on this.
source("./FBPFA-PFA with fixed latent dim.R")
library(tidyverse)
post_PFA <- function(fit) {
  library(tidyverse)

  # Determine posterior dimension (number of factors per sample)
  k_vec <- sapply(fit$Loading, ncol)
  mode_k <- as.numeric(names(sort(table(k_vec), decreasing = TRUE)[1]))
  
  # Filter posterior samples to those with mode_k
  keep_idx <- which(k_vec == mode_k)
  fit$Loading <- fit$Loading[keep_idx]
  fit$Latentsigma <- fit$Latentsigma[keep_idx]
  fit$Errorsigma <- fit$Errorsigma[keep_idx]
  fit$Pertmat <- fit$Pertmat[keep_idx]
  
  npost <- length(fit$Loading)
  p <- nrow(fit$Loading[[1]])
  k <- mode_k
  S <- dim(fit$Pertmat[[1]])[2]
  
  posteriorPhis <- array(0, dim = c(p, k, npost))
  posteriorLams <- vector("list", S)

  for(s in 1:S){
    posteriorLams[[s]] <- array(0, dim = c(p, k, npost))
    for(i in 1:npost){
      posteriorPhis[,,i] <- fit$Loading[[i]] %*% diag(fit$Latentsigma[[i]])
      posteriorLams[[s]][,,i] <- (solve(matrix(fit$Pertmat[[i]][, s], p, p)) - diag(p)) %*% posteriorPhis[,,i]
    }
  }

  # Varimax rotation
  est_Phi <- MSFA::sp_OP(posteriorPhis, itermax = 10, trace = FALSE)$Phi
  est_speLoad <- lapply(posteriorLams, function(x) MSFA::sp_OP(x, itermax = 10, trace = FALSE)$Phi)

  # Estimated covariance components
  sharevar <- list()
  est_SigmaLambdaList <- vector("list", S)
  est_SigmaMarginal <- vector("list", S)
  est_Psi_list <- list()

  for(s in 1:S){
    post_SigmaLambda_s <- vector("list", npost)
    post_SigmaMarginal_s <- vector("list", npost)
    Psi <- vector("list", npost)

    for(i in 1:npost){
      sharevar[[i]] <- fit$Loading[[i]] %*% diag(fit$Latentsigma[[i]]^2) %*% t(fit$Loading[[i]]) + 
        diag(fit$Errorsigma[[i]]^2)
      Q_temp_inv <- solve(matrix(fit$Pertmat[[i]][, s], p, p))
      post_SigmaMarginal_s[[i]] <- Q_temp_inv %*% sharevar[[i]] %*% t(Q_temp_inv)
      post_SigmaLambda_s[[i]] <- post_SigmaMarginal_s[[i]] - sharevar[[i]]
      Psi[[i]] <- diag(fit$Errorsigma[[i]]^2)
    }

    est_SigmaMarginal[[s]] <- Reduce('+', post_SigmaMarginal_s) / npost
    est_SigmaLambdaList[[s]] <- Reduce('+', post_SigmaLambda_s) / npost
    est_Psi_list[[s]] <- Reduce('+', Psi) / npost
  }

  est_Psi <- Reduce('+', est_Psi_list) / S
  est_SigmaPhi <- Reduce('+', sharevar) / npost
  est_Q <- Reduce('+', fit$Pertmat) / npost
  est_Q_list <- lapply(1:S, function(s) matrix(est_Q[, s], p, p))

  return(list(
    Phi = est_Phi,
    SigmaPhi = est_SigmaPhi,
    Psi = est_Psi,
    Q = est_Q_list,
    LambdaList = est_speLoad,
    SigmaLambdaList = est_SigmaLambdaList,
    SigmaMarginal = est_SigmaMarginal,
    mode_k = mode_k,
    kept_samples = length(keep_idx)
  ))
}


