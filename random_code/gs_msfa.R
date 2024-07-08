# adapted code from github.com/rdevito/MSFA

gs_msfa <- function(X_s,  k,  j_s, trace = TRUE, nprint = 1000,
                    outputlevel = 1, control = list(...), ...){
  #### read data ####
  S <- length(X_s)
  p <- dim(X_s[[1]])[2]                    #variables
  Y_s <- list()                         #centered data
  n_s <- c()                            #sample size
  
  
  time_start <- Sys.time()
  
  ################ setting up priors and initialize
  control <- do.call("sp_msfa_control", control)
  nrun <- control$nrun
  thin <- control$thin
  burn <- control$burn
  sp  <- (nrun - burn) / thin
  a_psi_s <- b_psi_s <-  c()    ####gamma hyperparameters for Psi_s
  nu_s <- c()                # gamma hyperpameters for omega_s
  a1_s <- b1_s <- c()      # gamma hyperparameters for delta_1
  a2_s <- b2_s <- c()      #gamma hyperparameters for delta_l^1
  
  #### priors
  apsi <- control$apsi
  bpsi<- control$bpsi
  nu <- control$nu
  nus <- control$nus
  a1 <- control$a1
  b1 <- control$b1
  a2 <- control$a2
  b2 <-  control$b2
  a1s <- control$a1s
  b1s <- control$b1s
  a2s <- control$a2s
  b2s <-  control$b2s
  
  ####initial values
  psi_s <- Psi_s <- Lambda_s <- l_s <- list()
  Meta_s <-  Veta_s <- f_s <- MetaP_s <- VetaP_s <- list()
  
  ####prior setup
  omegajh_s <- delta_s <- tauh_s <- Plam_s <- list()
  
  #### output
  Lambdaout <- psiout <- l_sout <- f_sout  <- SigmaLambda <- list()
  
  for (s in 1:S){
    n_s[s] <- nrow(X_s[[s]])
    Y_s[[s]] <- scale(X_s[[s]], center = TRUE, scale = FALSE)
    a_psi_s[s] <- apsi
    b_psi_s[s] <- bpsi
    nu_s[s] <- nus
    a1_s[s] <- a1s
    b1_s[s] <- b1s
    a2_s[s] <- a2s
    b2_s[s] <- b2s
    
    #### initial values ####
    ###for covariance error terms
    psi_s[[s]] <- rgamma(p, shape = a_psi_s[s], scale = 1 / b_psi_s[s])
    Psi_s[[s]] <- diag(1 / psi_s[[s]])
    ###for study-specific f.l.
    Lambda_s[[s]] <- zeros(p, j_s[s])
    l_s[[s]] <- matrix(rnorm(n_s[s] * j_s[s]), n_s[s], j_s[s])
    Meta_s[[s]] <- zeros(n_s[s], j_s[s])
    Veta_s[[s]] <- eye(j_s[s])
    ###for common f.l.
    f_s[[s]] <- matrix(rnorm(n_s[s]*k), n_s[s], k)
    MetaP_s[[s]] <- zeros(n_s[s], k)
    VetaP_s[[s]] <- eye(k)
    ###prior for study-specific f.l
    omegajh_s[[s]] <- matrix(rgamma(p * j_s[s], shape = nu_s[s] / 2,
                                    scale = 2 / nu_s[s]), p, j_s[s])
    delta_s[[s]] <- rgamma(j_s[s], shape=c(a1_s[s], rep(a2_s[s], j_s[s] - 1)),
                           scale = c(b1_s[s], rep(b2_s[s], j_s[s] - 1)))
    tauh_s[[s]] <- cumprod(delta_s[[s]])
    Plam_s[[s]] <- matvec(omegajh_s[[s]], tauh_s[[s]])
    #### save lambda and Psi####
    if(outputlevel == 1) {
      Lambdaout[[s]] <- array(0, dim=c(p, j_s[s], sp))
      psiout[[s]] <- array(0, dim=c(p, 1, sp))
      #l_sout[[s]] <- zeros(n_s[s], j_s[s])
      #f_sout[[s]] <- zeros(n_s[s], k)
    }
    if(outputlevel == 2) {
      #Lambdaout[[s]] <- zeros(p, j_s[s])
      #psiout[[s]] <- zeros(p, 1)
      l_sout[[s]] <- array(0, dim=c(n_s[s], j_s[s], sp))
      f_sout[[s]] <- array(0, dim=c(n_s[s], k, sp))
      #SigmaLambda[[s]] <- zeros(p, p)
    }
    if(outputlevel == 3) {
      SigmaPhi <- zeros(p, p)
      SigmaLambda[[s]] <- zeros(p, p)
      Lambdaout[[s]] <- zeros(p, j_s[s])
      #psiout[[s]] <- zeros(p, 1)
      #l_sout[[s]] <- zeros(n_s[s], j_s[s])
      #f_sout[[s]] <- zeros(n_s[s], k)
    }
    if(outputlevel == 4) {
      SigmaPhi <- zeros(p, p)
      #psiout[[s]] <- zeros(p, 1)
      #l_sout[[s]] <- zeros(n_s[s], j_s[s])
      #f_sout[[s]] <- zeros(n_s[s], k)
    }
  }
  #Phi
  Phi <- zeros(p, k)
  omegajh <- matrix(rgamma(p * k, shape = nu / 2, scale = 2 / nu), p, k)
  delta <- rgamma(k, shape = c(a1, rep(a2, k-1)), scale=c(b1, rep(b2, k - 1)))
  tauh <- cumprod(delta)
  Pphi <- matvec(omegajh, tauh)
  if(outputlevel == 1) Phiout <- array(0, dim=c(p, k, sp))
  if(outputlevel == 3) Phiout <- zeros(p, k)
  if(outputlevel == 4) {
    Phiout <- zeros(p, k)
    PhioutOP <- array(0, dim=c(p, k, sp))
  }
  #### start posterior sampling
  for(r in 1:nrun)
  {
    # if((Sys.time() - time_start) > 1 * 24 * 60 * 60){break} #see if algorithm has been runing > 1 day
    f_s2 <- l_s2 <- list()
    for (s in 1:S){
      ###Step 1: common latent factors
      LmsgP1 <- vecmat(psi_s[[s]], Phi)
      VetaP1 <- diag(k) + t(LmsgP1) %*% Phi
      TP1 <- chol(VetaP1)
      qrTP1 <- qr(TP1)
      QP1 <- qr.Q(qrTP1)
      RP1 <- qr.R(qrTP1)
      SP1 <- solve(RP1)
      VetaP11 <- tcrossprod(SP1)
      MetaP1 <- (Y_s[[s]] %*% LmsgP1 - l_s[[s]] %*% t(Lambda_s[[s]]) %*% LmsgP1) %*% VetaP11
      xP1 <- matrix(rnorm(n_s[s] * k), nrow = n_s[s], ncol = k)
      f_s[[s]] <- MetaP1 + xP1 %*% t(SP1)
      f_s2[[s]] <- crossprod(f_s[[s]])
      
      #Step 2: study-specific latent factors
      Lmsg1 <- vecmat(psi_s[[s]], Lambda_s[[s]])
      Veta1 <- diag(j_s[s]) + t(Lmsg1) %*% Lambda_s[[s]]
      T1 <- chol(Veta1)
      qrT1 <- qr(T1)
      Q1 <- qr.Q(qrT1)
      R1 <- qr.R(qrT1)
      S1 <- solve(R1)
      Veta11 <- tcrossprod(S1)
      Meta1 <- (Y_s[[s]] %*% Lmsg1 -  f_s[[s]] %*% t(Phi) %*% Lmsg1 ) %*% Veta11
      x1 <- matrix(rnorm(n_s[s] * j_s[s]), nrow = n_s[s], ncol = j_s[s])
      l_s[[s]] <- Meta1 + x1 %*% t(S1)
      l_s2[[s]] <- crossprod(l_s[[s]])
    }
    
    ### Step 3: common factor loadings
    for(i in 1:p){
      q1 <- b1list <- list()
      for(s in 1:S){
        q1[[s]] <- psi_s[[s]][i] * f_s2[[s]]
        b1list[[s]] <- psi_s[[s]][i]*(t(f_s[[s]]) %*% Y_s[[s]][, i]) - psi_s[[s]][i] * (t(f_s[[s]]) %*% l_s[[s]] %*% Lambda_s[[s]][i,])
      }
      q1_sum <- Reduce('+', q1)
      bphi <- Reduce('+', b1list)
      Qphi <- diag(Pphi[i,], k) + q1_sum
      Lphi <-  t(chol(Qphi))
      zphi <- rnorm(k)
      vphi <- forwardsolve(Lphi, bphi)
      mphi <-  backsolve(t(Lphi), vphi)
      yphi <- backsolve(t(Lphi), zphi)
      Phi[i,] <- t(yphi + mphi)
    }
    
    ### Step 4: specific factor loadings, with constraints
    for(s in 1:S){
      for(i in 1:p){
        Qlam <- diag(Plam_s[[s]][i,], j_s[s]) + psi_s[[s]][i] * l_s2[[s]]
        blam <- psi_s[[s]][i] * (t(l_s[[s]]) %*% Y_s[[s]][, i]) - psi_s[[s]][i] * (t(l_s[[s]]) %*% f_s[[s]] %*% Phi[i,])
        Llam <- t(chol(Qlam))
        zlam <- rnorm(j_s[s])
        vlam <- forwardsolve(Llam, blam)
        mlam <- backsolve(t(Llam), vlam)
        ylam <- backsolve(t(Llam), zlam)
        zlam4 <- zlam
        Lambda_s[[s]][i,] <- t(ylam + mlam)
      }
    }
    
    ### Step 5: omegajh for phi
    for(h in 1:k){
      omegajh[, h] <- rgamma(p, shape = (nu + 1) / 2, rate = (nu + tauh[h] * Phi[,h]^2) / 2)
    }
    mat <- omegajh * Phi^2
    ad <- a1 + 0.5 * p * k
    bd <- b1 + 0.5 * (1 / delta[1]) * sum(tauh * colSums(mat))
    delta[1] <- rgamma(1, shape = ad, scale = 1 / bd)
    tauh <- cumprod(delta)
    if (k>1){
      for(h in 2:k)
      {
        ad <- a2 + 0.5 * p * (k-h+1)
        bd <- b2 + 0.5 * (1 / delta[h]) * sum(tauh[h:k] * colSums(as.matrix(mat[, h:k])))
        delta[h] <- rgamma(1, shape = ad, scale = 1 / bd)
        tauh <- cumprod(delta)
      }
    }
    Pphi <- matvec(omegajh, tauh)
    ###################### omegajh for Lambda_s
    for(s in 1:S){
      for(h in 1:j_s[s]){
        omegajh_s[[s]][, h] <- rgamma(p, shape= (nu_s[[s]] + 1) / 2,
                                      rate = (nu_s[[s]] + tauh_s[[s]][h] * Lambda_s[[s]][,h]^2) / 2)
      }
      mat_s <- omegajh_s[[s]] * Lambda_s[[s]]^2
      ad <- a1_s[s] + 0.5 * p*j_s[s]
      bd <- b1_s[s] + 0.5 * (1 / delta_s[[s]][1]) * sum(tauh_s[[s]] * colSums(mat_s))
      delta_s[[s]][1] <- rgamma(1, shape = ad, scale = 1 / bd)
      tauh_s[[s]] <- cumprod(delta_s[[s]])
      if (j_s[s]>1){
        for(h in 2:j_s[s]){
          ad <- a2_s[s] + 0.5 * p * (j_s[s] - h + 1)
          bd <- b2_s[s] + 0.5 * (1 / delta_s[[s]][h]) * sum(tauh_s[[s]][h:j_s[s]] * colSums(as.matrix(mat_s[, h:j_s[s]])))
          delta_s[[s]][h] <- rgamma(1, shape = ad, scale = 1 / bd)
          tauh_s[[s]] <- cumprod(delta_s[[s]])
        }
      }
      Plam_s[[s]] <- matvec(omegajh_s[[s]], tauh_s[[s]])
    }
    
    ### Step 6: Psi_s
    for (s in 1:S){
      Ytil_s <- Y_s[[s]] - (l_s[[s]] %*% t(Lambda_s[[s]]) + f_s[[s]] %*% t(Phi))
      psi_s[[s]] <- rgamma(p, shape = a_psi_s[s] + 0.5 * n_s[s], rate = b_psi_s[s] + 0.5 * colSums(Ytil_s^2))
      Psi_s[[s]] <- diag(1 / psi_s[[s]])
    }
    
    if(r > burn){
      neff <- (r - burn) / thin
      teff <- (nrun - burn) / thin
      if(outputlevel == 1) Phiout[, , neff] <- Phi
      #if(outputlevel == 2) Phiout <- Phiout + Phi / teff
      if(outputlevel == 3) {
        Phiout <- Phiout + Phi / teff
        SigmaPhi <- SigmaPhi + tcrossprod(Phi) / teff
      }
      if(outputlevel == 4) {
        PhioutOP[, , neff] <- Phi
        #Phiout[, , neff] <- PhioutOP + Phi / teff
        SigmaPhi <- SigmaPhi + tcrossprod(Phi) / teff
      }
      
      if(outputlevel==1)
      {
        for(s in 1:S){
          Lambdaout[[s]][, , neff] <- Lambda_s[[s]]
          psiout[[s]][, , neff] <- 1 / psi_s[[s]]
          #l_sout[[s]] <- l_sout[[s]] + l_s[[s]] / teff
          #f_sout[[s]] <- f_sout[[s]] + f_s[[s]] / teff
        }
      }
      if(outputlevel==2){
        for(s in 1:S){
          #Lambdaout[[s]] <- Lambdaout[[s]] + Lambda_s[[s]]  / teff
          #psiout[[s]] <- psiout[[s]] + (1 / psi_s[[s]]) / teff
          l_sout[[s]][, , neff] <- l_s[[s]]
          f_sout[[s]][, , neff] <- f_s[[s]]
          #SigmaLambda[[s]] <- SigmaLambda[[s]] + tcrossprod(Lambda_s[[s]]) / teff
        }
      }
      if(outputlevel==3){
        for(s in 1:S){
          Lambdaout[[s]] <- Lambdaout[[s]] + Lambda_s[[s]]  / teff
          #psiout[[s]] <- psiout[[s]] + (1 / psi_s[[s]]) / teff
          #l_sout[[s]] <- l_sout[[s]] + l_s[[s]] / teff
          #f_sout[[s]] <- ff_sout[[s]] + f_s[[s]] / teff
          SigmaLambda[[s]] <- SigmaLambda[[s]] + tcrossprod(Lambda_s[[s]]) / teff
        }
      }
    }
    if (trace & r %% nprint == 0) cat("r=",r,"\n")
  }
  ##### Save and exit
  if(outputlevel < 3)  {
    SigmaPhi <- NULL
    SigmaLambda <- NULL
  }
  if(outputlevel == 1)  {
    l_sout <- NULL
    f_sout <- NULL
    PhioutOP <- NULL
  }
  if(outputlevel == 2)  {
    Phiout <- NULL
    Lambdaout <- NULL
    psiout <- NULL
    PhioutOP <- NULL
  }
  if(outputlevel == 3)  {
    Phiout <- NULL
    Lambdaout <- NULL
    psiout <- NULL
    l_sout <- NULL
    f_sout <- NULL
    PhioutOP <- NULL
  }
  if(outputlevel == 4)  {
    Lambdaout <- NULL
    psiout <- NULL
    l_sout <- NULL
    f_sout <- NULL
    SigmaLambda <- NULL
  }
  out <- list(Phi = Phiout, PhiOP = PhioutOP, Lambda = Lambdaout, psi = psiout, l_s = l_sout, f_s = f_sout,
              SigmaPhi = SigmaPhi, SigmaLambda = SigmaLambda)
  return(out)
}