packages <- c('MSFA', 'peakRAM', 'BFR.BE', 'tidyverse', 'matlab', 'MatrixCorrelation')
lapply(packages, library, character.only = TRUE)
library(devtools)
source("./functions/gen_senerioSS.R")
source("./functions/gen_senerioBMSFA.R")
source("./functions/calculateRV.R")
source("./functions/post_BMSFA.R")
source("./functions/post_PFA.R")
source("./functions/measurements.R")


# fit the models
fitting <- function(sim_data){
  Y_mat <- sim_data$Y_mat
  A <- sim_data$A
  X <- sim_data$X
  Y_list <- sim_data$Y_list
  n_s <- sim_data$n_s
  
  profile_BMSFA <- peakRAM({
    # #Fit the BMSFA model
    result_BMSFA <- sp_msfa(Y_list, k = 5, j_s = c(1,1,1,1),
                                 outputlevel = 1, scaling = FALSE,
                                 control = list(nrun = 5000, burn = 4000))
  }) 
  profile_MOMSS <- peakRAM({
    #Fit the MOM-SS model
    result_MOMSS <- BFR.BE.EM.CV(x = Y_mat, v = X, 
                                      b = A, q = 5, scaling = FALSE)
  })
   profile_PFA <- peakRAM({
  #   #Fit the PFA model
     fit_PFA <- PFA(t(Y_mat), grpind = rep(1:length(n_s), times = n_s), ini.PCA = T)
   })
  return(list(result_MOMSS = result_MOMSS, result_BMSFA = post_BMSFA(result_BMSFA),
              result_PFA = post_PFA(fit_PFA),
              profile_MOMSS = profile_MOMSS, profile_BMSFA = profile_BMSFA,
              profile_PFA = profile_PFA
              ))
}

################################# Example ###########################
# set.seed(6)
# sim_data_test <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
#result_test <- fitting(sim_data_test)
# # 
# # gen_MOMSS_MSE(sim_data_test, result_test$result_MOMSS)
# # 
# # gen_BMSFA_MSE(sim_data_test, result_test$result_BMSFA)
# # 
#est_perform(sim_data_test, result_test$result_BMSFA, result_test$result_MOMSS, result_test$result_PFA)
#computing_perform(result_test[c("profile_MOMSS", "profile_BMSFA", "profile_PFA")])
# fit <- PFA(t(sim_data_test$Y_mat), grpind = rep(1:length(sim_data_test$n_s), times = sim_data_test$n_s), ini.PCA = T)
# SigmaPhi <- post_PFA(fit)
# fields::image.plot(t(lambdap))




  
