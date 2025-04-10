# sc1
plot_covs(fit_sc1$point_est, fit_sc1$true)
ggplot2::ggsave("./Figs/sc1_covs.png",  width = 6, height = 7, units = "in", dpi = 300)
# sc2
plot_covs(fit_sc2$point_est, fit_sc2$true)
ggplot2::ggsave("./Figs/sc2_covs.png",  width = 8, height = 7, units = "in", dpi = 300)
# sc3
plot_covs(fit_sc3$point_est, fit_sc3$true)
ggplot2::ggsave("./Figs/sc3_covs.png",  width = 6, height = 7, units = "in", dpi = 300)



# ---------------------- Scenario 3 --------------------------------------
# Create a vector of file paths
file_paths <- paste0("./RDS/sc3/sc3_", 1:5, ".rds")
Tetris_file_paths <- paste0("./RDS/sc3/Tetris/sc3_Tetris_", 1:5, ".rds")
# Load all RDS files into a list
sc3_list <- lapply(file_paths, readRDS)
Tetris_sc3_list <- lapply(Tetris_file_paths, readRDS)

# A single run and heatmap
postfit_list <- c(sc3_list[[1]]$point_est, Tetris_sc3_list[[1]]$point_est)
source("./functions/plot_covs.R")
plot_covs(postfit_list, sc3_list[[1]]$true)
ggplot2::ggsave("./Figs/sc3_covs.png",  width = 6, height = 7, units = "in", dpi = 300)
source("./functions/plot_loadings.R")
plot_loadings(postfit_list, sc3_list[[1]]$true)
ggplot2::ggsave("./Figs/sc3_loadings.png", width = 6, height = 7, units = "in", dpi = 300)
plot_single(sc3_list[[1]]$true$Phi)
ggplot2::ggsave("./Figs/sc3_Phi.png", width = 6, height = 7, units = "in", dpi = 300)

# RV_Phi
# Extract the RV_Phi metrics from all loaded files
RV_Phi_list <- lapply(sc3_list, function(x) x[["metrics"]][["common"]]["RV_Phi",])
RV_Phi_combined <- do.call(rbind, RV_Phi_list) %>% as.data.frame()
RV_Phi_combined$seed <- (1:50)*3
RV_Phi_long <- pivot_longer(as.data.frame(RV_Phi_combined), cols=-seed, 
                            names_to = "method", values_to = "RV_Phi")
#unlist
RV_Phi_long$RV_Phi <- sapply(RV_Phi_long$RV_Phi, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})