source("random_code/simulate_data_fxns.R")
S=5 ## no. of studies
d= 200 ## observed dimension
q= 20 ## latent dimension of the shared subspace
n=d;ns=sapply(rpois(S,n/S),max,n/S) ## sample-sizes in each study
qs=sapply(rpois(S,floor(q/S))+1,min,floor(q/S))
val=2
lam=simulate_lambda_sparse (k = q,p=d,pr = .25,val = val)
blank_rows=which(rowSums(lam)==0)
for(bl in blank_rows){
  lam[bl,sample.int(q,5)]= runif(5,-val,val)
}
library(ggplot2)
heatplot(lam)
diag_sd=sqrt(.5)
null_lam=MASS::Null(lam) 
Y=replicate(S,list(1))
lambda=replicate(S,list(1))
covmat=replicate(S,list(1))
for(s in 1:S){
  eta=matrix(rnorm(ns[s]*q), nrow=ns[s],ncol=q)
  
  lambda_s=matrix (rnorm(q*qs[s],sd=.25),  nrow=q,ncol=qs[s])
  lambda[[s]]=lam %*% lambda_s + matrix(rnorm(d*qs[s],sd=.1),nrow = d,ncol=qs[s])
  
  ls=matrix(rnorm(ns[s]*qs[s]), nrow=ns[s],ncol=qs[s])
  
  Y[[s]] =tcrossprod(eta,lam)+tcrossprod(ls,lambda[[s]]) + matrix(rnorm(d*ns[s],sd=diag_sd),nrow=ns[s],ncol=d)
}

covmat_shared=tcrossprod(lam) + diag_sd*diag(d) ##True shared covariance matrix
time.taken=system.time(res<-fit_SUFA(Y,qmax=25,nthreads = 5,nrun = 7.5e3,
                                     nleapfrog = 4, leapmax = 9, del_range = c(0,.01)))
covmat_shared_est=SUFA_shared_covmat(res,burn = burnin)

## ----plot_sigmas, fig.dim = c(8, 5)-------------------------------------------
cormat_true=cov2cor(covmat_shared)
cormat_est=cov2cor(covmat_shared_est)
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----plot_sigmas_sparse, fig.dim = c(8, 5)------------------------------------
cormat_true[abs(cormat_true)<.1] =NA
cormat_est[abs(cormat_est)<.1] =NA
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----remove_ggplots, echo=FALSE-----------------------------------------------
rm(sig_true,sig_est,cormat_true,cormat_est)

## ----posterior_lambda---------------------------------------------------------
est_lam=lam.est(res$Lambda,burn = burnin)

## ----plot_lambdass, fig.dim = c(6, 5)-----------------------------------------
lam.sp=lam; lam.sp[lam==0]=NA
est_lam.sp=est_lam; est_lam.sp[est_lam==0]=NA

lam_true.gg=heatplot(lam.sp)+labs(title="True")
lam_est.gg=heatplot(est_lam.sp)+labs(title="Estimated")
ggpubr::ggarrange(lam_true.gg,lam_est.gg,ncol = 2)

####-------------ScenarioSS---------------------------------------------------
library(tidyverse)
library(MASS)
sim_data_test <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
res_ss<-fit_SUFA(sim_data_test$Y_list,qmax=5,nthreads = 5,nrun = 7.5e3,
              nleapfrog = 4, leapmax = 9, del_range = c(0,.01))
Phi_est <- lam.est(res_ss$Lambda,burn = burnin)
Phi_true <- sim_data_test$Phi
MatrixCorrelation::RV(Phi_est, Phi_true)

SigmaPhi_est <- SUFA_shared_covmat(res_ss,burn = burnin)
SigmaPhi_true <- tcrossprod(sim_data_test$Phi)
calculateRV(SigmaPhi_est, SigmaPhi_true)










