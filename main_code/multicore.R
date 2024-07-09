library(foreach)
library(doParallel)
library(iterators)
library(peakRAM)
library(devtools)
source("./functions/gen_senerioSS.R")
source("./functions/gen_senerioBMSFA.R")
source("./main_code/sim_once.R")
source("./functions/measurements.R")
packages <- c('MSFA', 'peakRAM', 'BFR.BE', 'tidyverse', 'matlab', 'MatrixCorrelation') ## MSFA used is the Mavis version


# Detect the number of cores
numCores <- parallel::detectCores()

# Register the parallel backend
cl <- makeCluster(numCores - 3)
registerDoParallel(cl)

###############################Senerio 1 Dense Phi##############################
#sen1_dense <- foreach(i = 1:50, .combine = 'rbind', 
 #                  .packages = packages) %dopar% {
  #                                 
   #                                data_sen1 <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
    #                               results <- fitting(data_sen1)
                                   
                                   # RV and FN for Phi
     #                               metrices <- est_perform(data_sen1, results$result_BMSFA, results$result_MOMSS)
                                   # Computing performance
      #                              comp_metrices <- computing_perform(results[c("profile_MOMSS", "profile_BMSFA")])
                                   
                                   # Output
       #                            c(metrices, comp_metrices)
        #                         }
#saveRDS(sen1_dense, "./sen1_dense.rds")
###############################################################################

##########################Senerio 1 Sparse Phi#################################
sen1_sparse <- foreach(i = 1:50, .combine = 'rbind', 
                      .packages = packages) %dopar% {
                                      data_sen1 <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5, genPhi = "sparse")
                                      results <- fitting(data_sen1)
                                      
                                      # RV and FN for Phi
                                      metrices <- est_perform(data_sen1, results$result_BMSFA, results$result_MOMSS)
                                      # Computing performance
                                      comp_metrices <- computing_perform(results[c("profile_MOMSS", "profile_BMSFA")])
                                      
                                      # Output
                                      c(metrices, comp_metrices)
                                    }
saveRDS(sen1_sparse, "./sen1_sparse.rds")
###############################################################################

# Stop the parallel backend after the job is done
stopImplicitCluster()
