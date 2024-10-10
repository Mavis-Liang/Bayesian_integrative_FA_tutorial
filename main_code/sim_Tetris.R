# This is the separate code to run Tetris on the replicated simulated data, since it runs slow and
# is not parallelizable on large cores.
source("./Tetris.R")
source("./functions/run_Tetris.R")
source("./functions/post_Tetris.R")
source("./functions/measurements.R")
sim_Tetris <- function(i){
  file_name <- paste0("./RDS/sc2/sc2_", i, ".rds")
  curr_data <- readRDS(file_name)
  profile_Tetris <- peakRAM::peakRAM({
    fit_Tetris <- run_Tetris(curr_data$true, big_T=(curr_data$true)$big_T)
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

