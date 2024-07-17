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
source("./functions/post_SUFA.R")


# fit the models
fitting <- function(sim_data){
  Y_mat <- sim_data$Y_mat
  A <- sim_data$A
  X <- sim_data$X
  Y_list <- sim_data$Y_list
  n_s <- sim_data$n_s
  set_K <- 5
  
  profile_PFA <- peakRAM({
    # #   #Fit the PFA model
    fit_PFA <- PFA(Y=t(Y_mat), 
                   latentdim = set_K,
                   grpind = rep(1:length(n_s), times = n_s))
  })
  
  profile_BMSFA <- peakRAM({
    # #Fit the BMSFA model
    fit_BMSFA <- sp_msfa(Y_list, k = set_K, j_s = c(1,1,1,1),
                                 outputlevel = 1, scaling = FALSE,
                                 control = list(nrun = 5000, burn = 4000))
  }) 
  profile_MOMSS <- peakRAM({
    #Fit the MOM-SS model
    fit_MOMSS <- BFR.BE.EM.CV(x = Y_mat, v = X, 
                                      b = A, q = set_K, scaling = FALSE)
  })
  
  profile_SUFA <- peakRAM({
    #Fit the SUFA model
    fit_SUFA <- fit_SUFA(Y_list, qmax = set_K, nthreads = 5, nrun = 5000,
                            nleapfrog = 4, leapmax = 9, del_range = c(0, .01))
  })
  return(list(result_MOMSS = fit_MOMSS, result_BMSFA = post_BMSFA(fit_BMSFA),
              result_PFA = post_PFA(fit_PFA),
              result_SUFA = post_SUFA(fit_SUFA),
              profile_MOMSS = profile_MOMSS, profile_BMSFA = profile_BMSFA,
              profile_PFA = profile_PFA,
              profile_SUFA = profile_SUFA
              ))
}

################################# Example ###########################
# set.seed(6)
# sim_data_test <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
# result_test <- fitting(sim_data_test)
# # # 
# gen_MOMSS_MSE(sim_data_test, result_test$result_MOMSS)
# # # 
# gen_BMSFA_MSE(sim_data_test, result_test$result_BMSFA)
# # # 
# est_perform(sim_data_test, result_test)
# computing_perform(result_test[c("profile_MOMSS", "profile_BMSFA", "profile_PFA", "profile_SUFA")])
# fit <- PFA(t(sim_data_test$Y_mat), grpind = rep(1:length(sim_data_test$n_s), times = sim_data_test$n_s), ini.PCA = T)
# SigmaPhi <- post_PFA(fit)
# fields::image.plot(t(lambdap))




  
