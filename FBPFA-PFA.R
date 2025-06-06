library(mvtnorm)
library(MASS)
library(pracma)
library(expm)
library(Rcpp)
library(RcppArmadillo)


################################################Function#################################################### 
#' @title The function is to fit PFA, FBPFA and Common Factor model with heterscedastic latent factors
#' @description Takes the variable (Y), group indices, some additional information to generate posterior draws
#' of the parameters in the model Q_jY_{ij} = Lambda*Eta_{ij}+E_{ij}
#' @references Roy et. al. (2019)
#'     "Perturbed factor analysis: Improving generalizability across studies" arXiv 
#'
#' @param Y is the data matrix with each column, representing one replication
#' @param d is the parameter of the prior of the latent factor's variance such that eta_{ijl}~N(0,sigma_{l}) and sigma_l~IG(d, 0.1)
#' @param grpind is the vector of group indices. Default is NULL, meaning no multi-group study. 
#' @param measureerror is the indicator, mentioning if it is the measurement error model in Section 2.2 of the paper
#' @param FB is the indicator if perturbation parameter alpha is provided or if it is a fully Bayesian approach
#' @param ini.PCA is the indicator to initialize using PCA or not 
#' @param Cutoff is the truncation value to remove columns from the loading matrix
#' @param alph is the preset value of perturbation parameter alpha if it is not FB
#' @param no.core is to specify number of cores for mcmapply() in linux, windows it should be 1
#' @param Thin is the an integer to specify thinig of postburn sample
#' @param burn is the burn length of posterior samples, multiple of thin
#' @param Total_itr is the total number of iterations of MCMC

#' @return List of posterior samples.  = lambda_p,  = eta_p, Errorsigma = sigma1_p,  = sigma2_p, = alph_p \cr
#' \item{Loading}{The posterior samples of Lambda}
#' \item{Latent}{The posterior samples of latent variables}
#' \item{Errorsigma}{The posterior samples of standard deviations of error}
#' \item{Latentsigma}{The posterior samples of standard deviations of latent variable}
#' \item{Pertmat}{The posterior samples of perturbation matrices each entry of the list is a matrix of dimension p^2 X (number of groups) i.e. ith column of the matrix represents vectorized perturation matrix of ith group}
#' \item{Alpha}{The posterior samples of perturbation parameter alpha}
#'
#' @export
#' @examples
Rcpp::sourceCpp('PFA.cpp')
PFA <- function(Y=Y, d = 100, grpind = NULL, measureerror = F, FB=T, ini.PCA=T, Cutoff = 1e-2, alph= 0.0001, no.core = 1, Thin= 10, burn = 10, Total_itr = 5000){
  QYpr <- function(g, mat = Y, vec = grpind, Ql = Qlist){
    temp <- matrix(Ql[, g], p, p)
    return(temp%*%mat[, which(grpind==g)])
  }
  
  n <- ncol(Y)
  p <- nrow(Y)
  grp <- 1
  if(measureerror){grpind <- 1:n}
  
  if(length(grpind)){grp <- length(unique(grpind))}
  if(length(grpind)==0){grpind <- rep(1, n)}
  
  Qlist <- matrix(array(diag(p)), p^2, grp)
  Qlistmat <- Qlist
  Qmean <- array(diag(p))
  
  # parameters initialization
  
  r = p
  nu <- 1
  a1 <- 2.1
  a2 <- 3.1
  eta.var2 <- rep(0, r)
  var <- eta.var2
  delta <- rnorm(r)
  eta.var2[r:1] <- rep(1, r)
  phi <- matrix(rgamma(p * r, nu / 2, nu / 2), p, r)
  
  psi <- c(rgamma(1,1,1), rgamma(r-1, 2, 1))
  tau <- exp(cumsum(log(psi)))
  
  eta    <- matrix(rnorm(r*n), r, n)
  lambda <- Y %*% ginv(crossprod(Y)) %*% t(Y)
  eta <- ginv(crossprod(lambda)) %*% (crossprod(lambda, Y))
  lambdaginv <- ginv(crossprod(lambda)) %*% t(lambda)
  
  gamma <- lambdaginv 
  
  if(ini.PCA){
    svY    <- svd(Y)
    lambda <- svY$u %*% diag(svY$d)
    lambda <- lambda[, 1:r]
    eta    <- t(svY$v)
    eta    <- eta[1:r, ]
    
    tau <- (1/svY$d[1:r])^2
    psi <- log(tau)
    psi <- c(psi[1], psi[2:r]-psi[1:(r-1)])
    psi <- exp(psi)
    
    gamma <- matrix(0, r, p)
  }
  
  sigma1 <- apply(Y - lambda %*% eta, 1, sd) 
  
  sigma2 <- apply(eta - gamma %*% Y, 1, sd)
  
  Yhat  <- lambda %*% eta
  for(i in 1:p){
    al       <- 0.1 + n /2
    be       <- 0.1 + sum((Y[i, ]-Yhat[i, ])^2)/2
    sigma1[i] <- sqrt(1/rgamma(1, al, be))
  }
  
  etahat  <- gamma %*% Y
  
  for(i in 1:r){
    al       <- 0.1 + n /2
    be       <- 0.1 + sum((eta[i, ]-etahat[i, ])^2)/2
    sigma2[i] <- sqrt(1/rgamma(1, al, be))
  }
  
  eta_p   <- list()
  lambda_p <- list()
  sigma1_p <- list()
  sigma2_p <- list()
  itr <- 0
  sds1 <- 1e-3
  sds2 <- 1e-3
  sdlam <- 1e-1
  sdgam <- 1e-3
  R <- 100
  incre <- 4
  secparam <- 10*(1:r)
  ard1 <- 0
  ard2 <- 0
  arl <- 0
  arg <- 0
  L2 <- 10
  epsilon2 <- t(rmvnorm(n, sigma = diag(sigma2)))
  po <- 1
  Qmean <- array(1*diag(p))
  Q_p <- list()
  alph <- alph
  alph_p <- rep(0, Total_itr)
  Qiup <- function(i, le){
    Qtemp    <- matrix(Qlistmat[, i], p, p)
    index <- which(grpind==i)
    YG   <- Y[, index]
    YG2  <- YG^2
    sigmaSq <- sigma1^2
    alphaf  <- alph
    leG     <- le[, index]
    
    QiupC(Qtemp, YG, YG2, leG, sigmaSq, alphaf)
    
    return(Qtemp)
  }
  
  QiupLR <- function(i, le){
    Qtemp    <- matrix(Qlistmat[, i], p, p)
    index <- which(grpind==i)
    YG   <- Y[, index]
    YG2  <- YG^2
    sigmaSq <- sigma1^2
    alphaf  <- alph
    leG     <- le[, index]
    
    QiupC(Qtemp, YG, YG2, leG, sigmaSq, alphaf)
    
    return(Qtemp)
  }
  
  
  Qmeanmat <- matrix(array(1*diag(p)), p^2, grp)
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  QY <- Y
  while (itr < Total_itr) {
    itr <- itr + 1
    temp <- parallel::mcmapply(1:grp, FUN = QYpr, SIMPLIFY=F, MoreArgs = list(mat=Y, vec = grpind, Ql = Qlist))
    
    for(g in 1:grp){
      QY[, which(grpind==g)] <- temp[[g]]
    }
    Yhatred <- QY - lambda %*% epsilon2
    
    for(i in 1:p){
      al       <- 0.1 + n /2
      be       <- 0.1 + sum((Yhatred[i, ])^2)/2
      sigma1[i] <- sqrt(1/rgamma(1, al, be))
    }
    
    for(i in 1:r){
      al       <- d + n /2
      be       <- 0.1 + sum((epsilon2[i, ])^2)/2
      sigma2[i] <- sqrt(1/rgamma(1, al, be))
    }
    
    var.pm <- ginv(crossprod(lambda / sigma1) + diag(1 / sigma2^2))
    var.pm <- (var.pm + t(var.pm)) / 2
    
    mean.etai <- t(lambda/sigma1^2) %*% (QY)
    mean.etai <- var.pm %*% mean.etai
    
    temp <- apply(mean.etai, 2, FUN=function(x){rmvnorm(1, x, var.pm)})
    epsilon2  <- temp
    
    eta <- epsilon2
    
    for(i in 1:r){
      mean.lami <- rowSums((QY - lambda[, -i] %*% eta[ - i, ])*matrix(eta[i, ], p, n, byrow=T)/sigma1^2)
      var.lami  <- 1/(sum(eta[i, ]^2)/sigma1^2 + phi[, i]*tau[i])
      mean.lami <- var.lami * mean.lami
      
      lambda[, i] <- rnorm(p, mean.lami, sqrt(var.lami))
      
      phi[, i]    <- rgamma(p, (nu + 1) / 2, ((lambda[, i]^2) * tau[i] + nu) / 2)
    }
    
    temp1 <- colSums((lambda ^ 2) * phi)
    
    tauprime1 <- tau / (psi[1])
    
    psi[1] <- rgamma(1, a1 + p*r/2, 1 + sum(temp1 * tauprime1) / 2)
    tau    <- tauprime1 * psi[1]
    
    for(i in 2:r){
      tauprime1 <- tau / (psi[i])
      psi[i] <- rgamma(1, a2 + p*(r - i + 1) / 2, 1 + sum(temp1[i:r] * tauprime1[i:r]) / 2)
      tau    <- tauprime1 * psi[i]
    }
    
    if(grp>1){
      Qlist[, 2:grp] <- parallel::mcmapply(2:grp, FUN = Qiup, MoreArgs = list(le=lambda%*%eta), mc.cores = no.core)
      
      Qlistmat <- (Qlist)
    }
    
    if(measureerror){
      Qlist[, 1:grp] <- parallel::mcmapply(2:grp, FUN = Qiup, MoreArgs = list(le=lambda%*%eta), mc.cores = no.core)
      
      Qlistmat <- (Qlist)
    }
    
    if(FB){
      #update alph
      al       <- 0.1 + (grp-1)*p*p /2
      be       <- 0.1 + sum((Qlistmat-Qmeanmat)^2)/2
      alph       <- (1/rgamma(1, al, be))
    }
    
    eta.var2[r:1] <- sigma2[r:1]^2
    
    if(itr > burn){
      BI <- burn / Thin
      if(itr %% Thin == 0){
        sigma1_p[[itr/Thin - BI]] <- sigma1
        sigma2_p[[itr/Thin - BI]] <- sigma2
        alph_p[itr/Thin - BI]     <- alph
        eta_p[[itr/Thin - BI]]    <- eta
        Q_p[[itr/Thin - BI]]      <- Qlist
        lambda_p[[itr/Thin - BI]] <- lambda 
      } 
    }
    
    if(itr %% R==0){
      u <- runif(1)
      if(u < exp(-1 - itr * 5 * 10^(-4) )){
        temp    <- colMeans(abs(lambda))
        c <- (which(temp < Cutoff))
        if(r-length(c)<3){c <- NULL}
        if(r < 3) {c <- NULL}
        if(length(c)> 0){
          r <- r - length(c)
          lambda <- lambda[, -c]
          phi <- phi[, -c]
          tau   <- tau[-c]
          psi   <- psi[-c]
          eta  <- eta[-c, ]
          epsilon2  <- epsilon2[-c, ]
          eta.var2 <- eta.var2[ -c]
          sigma2 <- eta.var2
          gamma    <- gamma[-c, ]
          delta    <- delta[ - c]
        }
      }
      R     <- R + incre
    }
    #Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
  }
  close(pb)
  
  return(list(Loading = lambda_p, Latent = eta_p, Errorsigma = sigma1_p, Latentsigma = sigma2_p, Pertmat = Q_p, Alpha = alph_p))
}
