packages <- c('MSFA', 'peakRAM', 'BFR.BE', 'tidyverse', 'matlab', 'MatrixCorrelation')
lapply(packages, library, character.only = TRUE)
library(devtools)
source("./FBPFA-PFA.R")## It's neccessary because it seems the "FBPFA-PFA with fixed latent dim.R" depends on this.
source("./FBPFA-PFA with fixed latent dim.R")
source("./functions/gen_senerioSS.R")
source("./functions/gen_senerioBMSFA.R")
source("./functions/calculateRV.R")
source("./functions/post_BMSFA.R")
source("./functions/post_PFA.R")
source("./functions/measurements.R")


set.seed(6)
sim_data_test <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
fit <- PFA(Y=t(sim_data_test$Y_mat), 
           latentdim = 5,
           grpind = rep(1:length(sim_data_test$n_s), times = sim_data_test$n_s))
Phi <- post_PFA(fit)
fields::image.plot(t(lambdap))
save(fit, file = "./random_code/fit_PFA.rds")
RV(Phi$est_Phi, sim_data_test$Phi)

################### Roy's example ###################

QYprG <- function(i, mat = Y, vec = grpind, Ql = Qlist){
  p <- nrow(Y)
  temp <- matrix(Ql[, vec[i]], p, p)
  return(temp%*%mat[, i])
}
n <- 500
p  <- 21
r  <- 5
grp <- 10


eta0 <- matrix(rnorm(r*n), r, n) 

lambda0           <- matrix(0, p, r)
lambda0[1:5, 1]   <- rnorm(5, mean = 5, sd = 1)
lambda0[5:9, 2]   <- rnorm(5, mean = 5, sd = 1)
lambda0[9:13, 3]  <- rnorm(5, mean = 5, sd = 1)
lambda0[13:17, 4] <- rnorm(5, mean = 5, sd = 1)
lambda0[17:21, 5] <- rnorm(5, mean = 5, sd = 1)

Y <- matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)#matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)
Qlist <- matrix(array(diag(p)), p^2, grp)
Qmean <- array(1*diag(p))

set.seed(100)
Qlist <- matrix(array(diag(p)), p^2, grp)
Qmean <- array(1*diag(p))
set.seed(1)

#Generating the perturbation matrices
for(i in 2:grp){
  Qlist[, i] <- array(matrix(rnorm(p^2, Qmean, sd = sqrt(0.0001)), p, p))
}
grpind = rep(1:10, each = 50)
Y <- parallel::mcmapply(1:n, FUN = QYprG, MoreArgs = list(mat=Y)) #matrix(Q %*% array(Y), p, n)
fit <- PFA(Y,latentdim = 5, grpind = rep(1:10, each = 50), ini.PCA = T)
sigma2p <- Reduce('+', fit$Latentsigma[200:499])/length(fit$Latentsigma[200:499])
lambdap <- (Reduce('+', fit$Loading[200:499])/length(fit$Latentsigma[200:499])) %*% diag(sigma2p)
RV(lambdap, lambda0)
source("./FBPFA-PFA.R")
PFA(Y, grpind = rep(1:10, each = 50), ini.PCA = T)
# OP
sigmas  <- fit$Latentsigma
lambdals <- fit$Loading

posteriorlam <- array(0, dim = c(21, 5, length(200:499)))
 
for(i in 1:length(200:499)){
  posteriorlam[,,i] <- lambdals[[i]] %*% diag(sigmas[[i]])
}
posteriorlam <- MSFA::sp_OP(posteriorlam, itermax = 10)

RV(posteriorlam$Phi, lambda0)


