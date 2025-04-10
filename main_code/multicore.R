packages <- c('BFR.BE', 'MSFA',"SUFA", 'peakRAM', 'tidyverse', 'matlab', 
              'MatrixCorrelation',"Rcpp")
lapply(packages, library, character.only = TRUE)
library(foreach)
library(doParallel)
library(iterators)
library(peakRAM)
library(devtools)

# Running scenario 1 ~ 3 50 times in parallel. Tetris excluded (runned separatly in sim_Tetris.R file).
scenario <- 3
source(paste0("./main_code/sim_scenario", scenario, ".R"))

registerDoParallel(50)

multicore_result <- foreach(i = 1:50, .errorhandling = "pass",
               .packages = packages, .noexport = c("QiupC"# "QiupC" is cpp function sourced in PFA.R, and should be import separately.
               )) %dopar% {
                 
                 # Source necessary R scripts
                 source("./FBPFA-PFA.R") ## It's neccessary because it seems the "FBPFA-PFA with fixed latent dim.R" depends on this.
                 source("./FBPFA-PFA with fixed latent dim.R")
                 source("./Tetris.R")
                 curr <- case_when(
                   scenario == 1 ~ sim_scenario1(seed = i*3),
                   scenario == 2 ~ sim_scenario2(seed = i*3),
                   scenario == 3 ~ sim_scenario3(seed = i*3),
                 )
                 saveRDS(curr, paste0("./RDS/sc", scenario, "/sc", scenario, "_", i, ".rds"))
                 return(curr)
               }
stopImplicitCluster()
saveRDS(multicore_result, "./RDS/sc3_multicore_noTetris.rds")
############################################################################################