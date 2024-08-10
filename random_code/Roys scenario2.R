library(mvtnorm)
library(MASS)
library(pracma)
library(expm)
source("./FBPFA-PFA.R")## It's neccessary because it seems the "FBPFA-PFA with fixed latent dim.R" depends on this.
source("./FBPFA-PFA with fixed latent dim.R")

QYpr <- function(i, mat = Y){
  temp <- matrix(Qlist[, ceiling(i/50)], p, p)
  return(temp%*%mat[, i])
}

set.seed(1)
# data generation
n <- 500
p  <- 128
r  <- 5
grp <- 10

eta0 <- matrix(rnorm(r*n), r, n)# standard normal

# generate common loading
lambda0 <- matrix(0, p, r)
lambda0[1:64, 1] <- rnorm(64, mean = 5, sd = 1)
lambda0[c(1+0:29, 65+0:29), 2] <- rnorm(60, mean = 5, sd = 1)
temp <- c(0:12 + 1, 0:12+30, 0:12+65, 0:12+96)
lambda0[temp, 3] <- rnorm(52, mean = 5, sd = 1)
temp <- c(0:4 + 1, 0:4 + 14, 0:4+30,0:4 + 43)
temp <- c(temp, 0:4 + 65, 0:4 + 78, 0:4+96,0:4 + 109)
lambda0[temp, 4] <- rnorm(40, mean = 5, sd = 1)
temp <- c(0:1 + 1, 0:1 + 6, 0:1+14, 0:1 + 19, 0:1 + 30, 0:1+35, 0:1+43, 0:1+48)
temp <- c(temp, 0:1 + 65, 0:1 + 70, 0:1+78, 0:1 + 83, 0:1 + 96, 0:1+101, 0:1+109, 0:1+114)
lambda0[temp, 5] <- rnorm(32, mean = 5, sd = 1)
image(t(lambda0))


Y <- matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)#matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n) #t(rmvnorm(n, sigma = solve(pdmat)))#
#Plambda0 <- diag(p) - lambda0 %*% solve(crossprod(lambda0)) %*% t(lambda0)
Qlist <- matrix(array(diag(p)), p^2, grp)
po <- 1
Qmean <- array(1*diag(p))

#####################Generating the perturbation matrices################
for(i in 2:grp){
  Qlist[, i] <- array(matrix(rnorm(p^2, Qmean, sd = sqrt(0.01)), p, p))
}


Y <- parallel::mcmapply(1:n, FUN = QYpr, MoreArgs = list(mat=Y)) #matrix(Q %*% array(Y), p, n)

# Fitting: PFA
fit_PFA_senPFA <- PFA(Y=Y, 
                      latentdim = 50,
                      grpind = rep(1:10, times = 50))

############## Results ################
# True common covariance
true_comCov <- lambda0 %*% t(lambda0) + diag(p)
# Estimated common covariance
sharevar <- list()# store posterio samples
for(i in 1:1000){
  sharevar[[i]] <- fit_PFA_senPFA$Loading[[i]]%*%
    diag(fit_PFA_senPFA$Latentsigma[[i]]^2)%*%
    t(fit_PFA_senPFA$Loading[[i]]) + 
    diag(fit_PFA_senPFA$Errorsigma[[i]]^2)
}
est_comCov <- Reduce('+', sharevar)/length(sharevar) # Point estimate
# RV of common covariance
MatrixCorrelation::RV(true_comCov, est_comCov)# PFA's

