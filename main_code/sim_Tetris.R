# This is the separate code to run Tetris on the replicated simulated data, since it runs slow and
# is not parallelizable on large cores.
source("./Tetris.R")
source("./functions/post_Tetris.R")
source("./functions/run_Tetris.R") # don't source Tetris.R again. Control parameters are already defined in run_Tetris.R. Also, source this after sourcing post_Tetris.R because it contains Tetris.R.
source("./functions/measurements.R")
sim_Tetris <- function(i, scenario){
  file_name <- paste0("./RDS/sc", scenario, "/sc", scenario, "_", i, ".rds")
  curr_data <- readRDS(file_name)
  new_control <- function(nrun = 2000, burn = 1500, thin = 1,
                        nu = 3, 
                        a1 = 2.1, b1 = 1,
                        a2 = 3.1, b2 = 1,
                        apsi = 1, bpsi = 0.3,
                        alpha = 10, beta = 0.5)
    {
      return(list(nrun = nrun, burn = burn, thin = thin,
                  nu = nu, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
                  apsi = apsi, bpsi = bpsi, alpha = alpha, beta = beta))
    }
  profile_Tetris <- peakRAM::peakRAM({
    fit_Tetris <- run_Tetris(curr_data$true)
  })
  postfit_tetris <- post_Tetris(fit_Tetris)
  metrics_tables(list(Tetris=postfit_tetris), curr_data$true)
  final_list <- list(seed = curr_data$seed,
                     true = curr_data$true,
                     point_est = list(Tetris = postfit_tetris),
                     metrics = metrics_tables(list(Tetris=postfit_tetris), curr_data$true), 
                     efficiency = list(profile_Tetris[,c(2, 4)]))
  return(final_list)
}

