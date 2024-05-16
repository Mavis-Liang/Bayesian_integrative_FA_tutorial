library(foreach)
library(doParallel)
library(iterators)
library(peakRAM)
library(devtools)
source("gen_senerioSS.R")
source("gen_senerioBMSFA.R")
source("calculateRV.R")

# Detect the number of cores
numCores <- parallel::detectCores()

# Register the parallel backend
cl <- makeCluster(numCores - 3)
registerDoParallel(cl)

###############################Senerio 1 Dense Phi##############################
time_start <- Sys.time()
sen1_dense <- foreach(i = 1:50, .combine = 'rbind', 
                   .packages = c('MSFA', 'peakRAM','BFR.BE',
                                 'tidyverse','matlab','MatrixCorrelation')) %dopar% {
                                   data_sen1 <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5)
                                   sen1_Y_mat <- data_sen1$Y_mat
                                   sen1_A <- data_sen1$A
                                   sen1_X <- data_sen1$X
                                   sen1_Y_list <- data_sen1$Y_list
                                   
                                   profile_MOMSS <- peakRAM({
                                     #Fit the MOM-SS model
                                     result_MOMSS_sen1 <- BFR.BE.EM.CV(x = sen1_Y_mat, v = sen1_X, 
                                                                       b = sen1_A, q = 5, scaling = FALSE)
                                   })
                                   
                                    profile_BMSFA <- peakRAM({
                                     # #Fit the BMSFA model
                                     result_BMSFA_sen1 <- sp_msfa(sen1_Y_list, k = 5, j_s = c(1,1,1,1),
                                                                  outputlevel = 3, scaling = FALSE,
                                                                  control = list(nrun = 5000, burn = 4000))
                                   })
                                   
                                    # RV and FN for Phi
                                    RV_MOMSS_Phi <- RV(result_MOMSS_sen1$M, 
                                                            data_sen1$Phi)
                                    FN_MOMSS_Phi <- frobenius.norm(tcrossprod(result_MOMSS_sen1$M)
                                                   -tcrossprod(data_sen1$Phi))
                                    
                                    RV_BMSFA_Phi <- calculateRV(result_BMSFA_sen1$SigmaPhi, 
                                                                tcrossprod(data_sen1$Phi))
                                    FN_BMSFA_Phi <- frobenius.norm(result_BMSFA_sen1$SigmaPhi 
                                                   - tcrossprod(data_sen1$Phi))
                                   
                                   # Output
                                   list(run_time_MOMSS = profile_MOMSS[,'Elapsed_Time_sec'],
                                        peak_RAM_MOMSS = profile_MOMSS[,'Peak_RAM_Used_MiB'],
                                        run_time_BMSFA = profile_BMSFA[,'Elapsed_Time_sec'],
                                        peak_RAM_BMSFA = profile_BMSFA[,'Peak_RAM_Used_MiB'],
                                        RV_MOMSS_Phi=RV_MOMSS_Phi,
                                        FN_MOMSS_Phi = FN_MOMSS_Phi,
                                        RV_BMSFA_Phi = RV_BMSFA_Phi,
                                        FN_BMSFA_Phi = FN_BMSFA_Phi)
                                 }

saveRDS(sen1_dense, "sen1_dense.rds")
time_end1 <- Sys.time() - time_start
###############################################################################

##########################Senerio 1 Sparse Phi#################################
time_start <- Sys.time()
sen1_sparse <- foreach(i = 1:50, .combine = 'rbind', 
                      .packages = c('MSFA', 'peakRAM','BFR.BE',
                                    'tidyverse','matlab','MatrixCorrelation')) %dopar% {
                                      data_sen1 <- gen_senerioSS(S=4, N=500, P=50, Q=2, K=5, genPhi = "sparse")
                                      sen1_Y_mat <- data_sen1$Y_mat
                                      sen1_A <- data_sen1$A
                                      sen1_X <- data_sen1$X
                                      sen1_Y_list <- data_sen1$Y_list
                                      
                                      profile_MOMSS <- peakRAM({
                                        #Fit the MOM-SS model
                                        result_MOMSS_sen1 <- BFR.BE.EM.CV(x = sen1_Y_mat, v = sen1_X, 
                                                                          b = sen1_A, q = 5, scaling = FALSE)
                                      })
                                      
                                      profile_BMSFA <- peakRAM({
                                        # #Fit the BMSFA model
                                        result_BMSFA_sen1 <- sp_msfa(sen1_Y_list, k = 5, j_s = c(1,1,1,1),
                                                                     outputlevel = 3, scaling = FALSE,
                                                                     control = list(nrun = 5000, burn = 4000))
                                      })
                                      
                                      # RV and FN for Phi
                                      RV_MOMSS_Phi <- RV(result_MOMSS_sen1$M, 
                                                         data_sen1$Phi)
                                      FN_MOMSS_Phi <- frobenius.norm(tcrossprod(result_MOMSS_sen1$M)
                                                                     -tcrossprod(data_sen1$Phi))
                                      
                                      RV_BMSFA_Phi <- calculateRV(result_BMSFA_sen1$SigmaPhi, 
                                                                  tcrossprod(data_sen1$Phi))
                                      FN_BMSFA_Phi <- frobenius.norm(result_BMSFA_sen1$SigmaPhi 
                                                                     - tcrossprod(data_sen1$Phi))
                                      
                                      # Output
                                      list(run_time_MOMSS = profile_MOMSS[,'Elapsed_Time_sec'],
                                           peak_RAM_MOMSS = profile_MOMSS[,'Peak_RAM_Used_MiB'],
                                           run_time_BMSFA = profile_BMSFA[,'Elapsed_Time_sec'],
                                           peak_RAM_BMSFA = profile_BMSFA[,'Peak_RAM_Used_MiB'],
                                           RV_MOMSS_Phi=RV_MOMSS_Phi,
                                           FN_MOMSS_Phi = FN_MOMSS_Phi,
                                           RV_BMSFA_Phi = RV_BMSFA_Phi,
                                           FN_BMSFA_Phi = FN_BMSFA_Phi)
                                    }
saveRDS(sen1_sparse, "sen1_sparse.rds")
Sys.time() - time_start

###############################################################################

########################Senerio 2 Sparse Phi Sparse Lambda#####################
time_start <- Sys.time()
sen1_sparse <- foreach(i = 1:50, .combine = 'rbind', 
                       .packages = c('MSFA', 'peakRAM','BFR.BE',
                                     'tidyverse','matlab','MatrixCorrelation')) %dopar% {
                                       data_sen2 <- gen_senerioBMSFA(S=4, N=500, P=50, K=5, j_s=c(1,1,1,1), 
                                                                     genPhi = "sparse", genLambda = "sparse")
                                       sen2_Y_mat <- data_sen2$Y_mat
                                       sen2_A <- data_sen2$A
                                       sen2_Y_list <- data_sen2$Y_list
                                       
                                       # profile_MOMSS <- peakRAM({
                                       #   #Fit the MOM-SS model
                                       #   result_MOMSS_sen2 <- BFR.BE.EM.CV(x = sen2_Y_mat,
                                       #                                     b = sen2_A, q = 5, scaling = FALSE)
                                       # })
                                       
                                       profile_BMSFA <- peakRAM({
                                         # #Fit the BMSFA model
                                         result_BMSFA_sen2 <- sp_msfa(sen2_Y_list, k = 5, j_s = c(1,1,1,1),
                                                                      outputlevel = 3, scaling = FALSE,
                                                                      control = list(nrun = 5000, burn = 4000))
                                       })
                                       
                                       # RV and FN for Phi
                                       # RV_MOMSS_Phi <- RV(result_MOMSS_sen2$M, 
                                       #                    data_sen2$Phi)
                                       # FN_MOMSS_Phi <- frobenius.norm(tcrossprod(result_MOMSS_sen2$M)
                                       #                                -tcrossprod(data_sen2$Phi))
                                       
                                       RV_BMSFA_Phi <- calculateRV(result_BMSFA_sen2$SigmaPhi, 
                                                                   tcrossprod(data_sen2$Phi))
                                       FN_BMSFA_Phi <- frobenius.norm(result_BMSFA_sen2$SigmaPhi 
                                                                      - tcrossprod(data_sen2$Phi))
                                       
                                       # Output
                                       list(#run_time_MOMSS = profile_MOMSS[,'Elapsed_Time_sec'],
                                           # peak_RAM_MOMSS = profile_MOMSS[,'Peak_RAM_Used_MiB'],
                                            run_time_BMSFA = profile_BMSFA[,'Elapsed_Time_sec'],
                                            peak_RAM_BMSFA = profile_BMSFA[,'Peak_RAM_Used_MiB'],
                                            #RV_MOMSS_Phi=RV_MOMSS_Phi,
                                            #FN_MOMSS_Phi = FN_MOMSS_Phi,
                                            RV_BMSFA_Phi = RV_BMSFA_Phi,
                                            FN_BMSFA_Phi = FN_BMSFA_Phi)
                                     }
saveRDS(sen2_sparse, "sen2_sparse.rds")
time_end3 <- Sys.time() - time_start
###############################################################################

########################Senerio 2 Dense Phi Dense Lambda#####################
time_start <- Sys.time()
sen1_dense <- foreach(i = 1:50, .combine = 'rbind', 
                       .packages = c('MSFA', 'peakRAM','BFR.BE',
                                     'tidyverse','matlab','MatrixCorrelation')) %dopar% {
                                       data_sen2 <- gen_senerioBMSFA(S=4, N=500, P=50, K=5, j_s=c(1,1,1,1), 
                                                                     genPhi = "dense", genLambda = "dense")
                                       sen2_Y_mat <- data_sen2$Y_mat
                                       sen2_A <- data_sen2$A
                                       sen2_Y_list <- data_sen2$Y_list
                                       
                                       profile_MOMSS <- peakRAM({
                                         #Fit the MOM-SS model
                                         result_MOMSS_sen2 <- BFR.BE.EM.CV(x = sen2_Y_mat,
                                                                           b = sen2_A, q = 5, scaling = FALSE)
                                       })
                                       
                                        profile_BMSFA <- peakRAM({
                                          # #Fit the BMSFA model
                                          result_BMSFA_sen2 <- sp_msfa(sen2_Y_list, k = 5, j_s = c(1,1,1,1),
                                                                       outputlevel = 3, scaling = FALSE,
                                                                       control = list(nrun = 5000, burn = 4000))
                                        })

                                       # RV and FN for Phi
                                       RV_MOMSS_Phi <- RV(result_MOMSS_sen2$M,
                                                          data_sen2$Phi)
                                       FN_MOMSS_Phi <- frobenius.norm(tcrossprod(result_MOMSS_sen2$M)
                                                                      -tcrossprod(data_sen2$Phi))
                                       
                                       RV_BMSFA_Phi <- calculateRV(result_BMSFA_sen2$SigmaPhi, 
                                                                   tcrossprod(data_sen2$Phi))
                                       FN_BMSFA_Phi <- frobenius.norm(result_BMSFA_sen2$SigmaPhi 
                                                                      - tcrossprod(data_sen2$Phi))
                                       
                                       # Output
                                       list(#run_time_MOMSS = profile_MOMSS[,'Elapsed_Time_sec'],
                                            #peak_RAM_MOMSS = profile_MOMSS[,'Peak_RAM_Used_MiB'],
                                            run_time_BMSFA = profile_BMSFA[,'Elapsed_Time_sec'],
                                            peak_RAM_BMSFA = profile_BMSFA[,'Peak_RAM_Used_MiB'],
                                            #RV_MOMSS_Phi=RV_MOMSS_Phi,
                                            #FN_MOMSS_Phi = FN_MOMSS_Phi,
                                            RV_BMSFA_Phi = RV_BMSFA_Phi,
                                            FN_BMSFA_Phi = FN_BMSFA_Phi)
                                     }
saveRDS(sen2_dense, "sen2_dense.rds")
time_end4 <- Sys.time() - time_start



# Stop the parallel backend after the job is done
stopImplicitCluster()
