library(foreach)
library(doParallel)
library(iterators)
library(peakRAM)
library(devtools)
packages <- c('MSFA', 'peakRAM', 'BFR.BE', 'tidyverse', 'matlab', 'MatrixCorrelation',"Rcpp","SUFA")
lapply(packages, library, character.only = TRUE)


# Detect the number of cores
numCores <- parallel::detectCores()

# Register the parallel backend
cl <- makeCluster(numCores - 3)
registerDoParallel(cl)

###############################Senerio 1 Dense Phi##############################
sen1_dense <- foreach(i = 1:50, .combine = 'rbind',
                      # "QiupC" is cpp function sourced in PFA.R, and should be import separately.
                   .packages = packages, .noexport = c("QiupC")) %dopar% {
                     # Source C++ functions
                     Rcpp::sourceCpp('PFA.cpp')
                     
                     # Source necessary R scripts
                     source("./FBPFA-PFA.R") ## It's neccessary because it seems the "FBPFA-PFA with fixed latent dim.R" depends on this.
                     source("./FBPFA-PFA with fixed latent dim.R")
                     source("./functions/gen_senerioSS.R")
                     source("./functions/gen_senerioBMSFA.R")
                     source("./functions/calculateRV.R")
                     source("./functions/post_BMSFA.R")
                     source("./functions/post_PFA.R")
                     source("./functions/measurements.R")
                     source("./functions/post_SUFA.R")
                     source("./functions/sim_once.R")
                                   data_sen1 <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
                                   results <- fitting(data_sen1)

                                   # RV and FN for Phi
                                    metrices <- est_perform(data_sen1, results)
                                   # Computing performance
                                    comp_metrices <- computing_perform(results[c("profile_MOMSS", "profile_BMSFA", 
                                                                                 "profile_PFA", "profile_SUFA")])

                                   # Output
                                   c(metrices, comp_metrices)
                                 }
saveRDS(sen1_dense, "./sen1_dense.rds")
# Stop the parallel backend after the job is done
stopImplicitCluster()
###############################################################################

#########################Senerio 1 Sparse Phi#################################
# registerDoParallel(cl)
# sen1_sparse <- foreach(i = 1:50, .combine = 'rbind', 
#                       .packages = packages) %dopar% {
#                                       data_sen1 <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5, genPhi = "sparse")
#                                       results <- fitting(data_sen1)
#                                       
#                                       # RV and FN for Phi
#                                       metrices <- est_perform(data_sen1, results$result_BMSFA, results$result_MOMSS, results$result_PFA)
#                                       # Computing performance
#                                       comp_metrices <- computing_perform(results[c("profile_MOMSS", "profile_BMSFA", "profile_PFA")])
#                                       
#                                       # Output
#                                       c(metrices, comp_metrices)
#                                     }
# saveRDS(sen1_sparse, "./sen1_sparse.rds")
# # Stop the parallel backend after the job is done
# stopImplicitCluster()
###############################################################################



