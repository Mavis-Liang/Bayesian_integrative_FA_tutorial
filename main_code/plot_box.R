# ---------------------- Scenario 1 --------------------------------------
file_paths <- paste0("./RDS/sc1/sc1_", 1:50, ".rds")
sc1_list <- lapply(file_paths, readRDS)
# RV_Phi
RV_Phi_list_1 <- lapply(sc1_list, function(x) x[["metrics"]][["common"]]["RV_Phi",])
RV_Phi_combined_1 <- do.call(rbind, RV_Phi_list_1) %>% as.data.frame()
RV_Phi_combined_1$seed <- (1:50)*3
RV_Phi_long_1 <- pivot_longer(as.data.frame(RV_Phi_combined_1), cols=-seed, 
                            names_to = "method", values_to = "RV_Phi")
#unlist
RV_Phi_long_1$RV_Phi <- sapply(RV_Phi_long_1$RV_Phi, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# RV of Lambda
RV_Lambda_combined_1 <- data.frame(matrix(NA, nrow = 50, ncol = 7))
colnames(RV_Lambda_combined_1) <- names(RV_Phi_combined_1)
RV_Lambda_combined_1$seed <- (1:50)*3
RV_Lambda_long_1 <- pivot_longer(as.data.frame(RV_Lambda_combined_1), cols=-seed, 
                               names_to = "method", values_to = "RV_Lambda")

# RV of SigmaMarginal
RV_SigmaMarginal_list_1 <- lapply(sc1_list, function(x) x[["metrics"]][["marginal"]]["RV_SigmaMarginal_mean",])
RV_SigmaMarginal_combined_1 <- do.call(rbind, RV_SigmaMarginal_list_1) %>% as.data.frame()
RV_SigmaMarginal_combined_1$seed <- (1:50)*3
RV_SigmaMarginal_long_1 <- pivot_longer(as.data.frame(RV_SigmaMarginal_combined_1), cols=-seed, 
                                      names_to = "method", values_to = "RV_SigmaMarginal")
# Unlist
RV_SigmaMarginal_long_1$RV_SigmaMarginal <- sapply(RV_SigmaMarginal_long_1$RV_SigmaMarginal, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# Time
Time_list_1 <- lapply(sc1_list, function(x) x[["efficiency"]][[1]]["Elapsed_Time_sec",])
Time_combined_1 <- do.call(rbind, Time_list_1) %>% as.data.frame()/60
Time_combined_1$seed <- (1:50)*3
Time_long_1 <- pivot_longer(as.data.frame(Time_combined_1), cols=-seed, 
                          names_to = "method", values_to = "Time_min")
# Unlist
Time_long_1$Time_min <- sapply(Time_long_1$Time_min, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# Memory
RAM_list_1 <- lapply(sc1_list, function(x) x[["efficiency"]][[1]]["Peak_RAM_Used_MiB",])
RAM_combined_1 <- do.call(rbind, RAM_list_1) %>% as.data.frame()
RAM_combined_1$seed <- (1:50)*3
  
RAM_long_1 <- pivot_longer(as.data.frame(RAM_combined_1), cols=-seed, 
                         names_to = "method", values_to = "RAM_MiB")
# Unlist
RAM_long_1$RAM_MiB <- sapply(RAM_long_1$RAM_MiB, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# Merge
df_1 <- merge(RV_Phi_long_1, RV_Lambda_long_1, by = c("seed", "method"))
df_1 <- merge(df_1, RV_SigmaMarginal_long_1, by = c("seed", "method"))
df_1 <- merge(df_1, Time_long_1, by
              = c("seed", "method"))
df_1 <- merge(df_1, RAM_long_1, by = c("seed", "method"))
df_long_1 <- pivot_longer(df_1, cols = -c(seed, method), 
                        names_to = "metric", values_to = "value")
df_long_1$scenario <- "Scenario 1"

# ---------------------- Scenario 2 --------------------------------------
file_paths <- paste0("./RDS/sc2/sc2_", 1:50, ".rds")
sc2_list <- lapply(file_paths, readRDS)
# RV_Phi
RV_Phi_list_2 <- lapply(sc2_list, function(x) x[["metrics"]][["common"]]["RV_Phi",])
RV_Phi_combined_2 <- do.call(rbind, RV_Phi_list_2) %>% as.data.frame()
RV_Phi_combined_2$seed <- (1:50)*3
RV_Phi_long_2 <- pivot_longer(as.data.frame(RV_Phi_combined_2), cols=-seed, 
                            names_to = "method", values_to = "RV_Phi")
#unlist
RV_Phi_long_2$RV_Phi <- sapply(RV_Phi_long_2$RV_Phi, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# RV_Lambda
# Extract the RV_SigmaLambda metrics from all loaded files
RV_Lambda_list_2 <- lapply(sc2_list, function(x) x[["metrics"]][["specific"]]["RV_Lambda_mean",])
RV_Lambda_combined_2 <- do.call(rbind, RV_Lambda_list_2) %>% as.data.frame()
RV_Lambda_combined_2$seed <- (1:50)*3
RV_Lambda_long_2 <- pivot_longer(as.data.frame(RV_Lambda_combined_2), cols=-seed, 
                                    names_to = "method", values_to = "RV_Lambda")
# Unlist
RV_Lambda_long_2$RV_Lambda <- sapply(RV_Lambda_long_2$RV_Lambda, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# RV of SigmaMarginal
RV_SigmaMarginal_list_2 <- lapply(sc2_list, function(x) x[["metrics"]][["marginal"]]["RV_SigmaMarginal_mean",])
RV_SigmaMarginal_combined_2 <- do.call(rbind, RV_SigmaMarginal_list_2) %>% as.data.frame()
RV_SigmaMarginal_combined_2$seed <- (1:50)*3
RV_SigmaMarginal_long_2 <- pivot_longer(as.data.frame(RV_SigmaMarginal_combined_2), cols=-seed, 
                                      names_to = "method", values_to = "RV_SigmaMarginal")
# Unlist
RV_SigmaMarginal_long_2$RV_SigmaMarginal <- sapply(RV_SigmaMarginal_long_2$RV_SigmaMarginal, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# Time
Time_list_2 <- lapply(sc2_list, function(x) x[["efficiency"]][[1]]["Elapsed_Time_sec",])
Time_combined_2 <- do.call(rbind, Time_list_2) %>% as.data.frame()/60
Time_combined_2$seed <- (1:50)*3
Time_long_2 <- pivot_longer(as.data.frame(Time_combined_2), cols=-seed, 
                          names_to = "method", values_to = "Time_min")
# Memory
RAM_list_2 <- lapply(sc2_list, function(x) x[["efficiency"]][[1]]["Peak_RAM_Used_MiB",])
RAM_combined_2 <- do.call(rbind, RAM_list_2) %>% as.data.frame()
RAM_combined_2$seed <- (1:50)*3
RAM_long_2 <- pivot_longer(as.data.frame(RAM_combined_2), cols=-seed, 
                         names_to = "method", values_to = "RAM_MiB")

# Merge
df_2 <- merge(RV_Phi_long_2, RV_Lambda_long_2, by = c("seed", "method"))
df_2 <- merge(df_2, RV_SigmaMarginal_long_2, by = c("seed", "method"))
df_2 <- merge(df_2, Time_long_2, by = c("seed", "method"))
df_2 <- merge(df_2, RAM_long_2, by = c("seed", "method"))
df_long_2 <- pivot_longer(df_2, cols = -c(seed, method), 
                        names_to = "metric", values_to = "value")
df_long_2$scenario <- "Scenario 2"

# ---------------------- Scenario 3 --------------------------------------
# Create a vector of file paths
file_paths <- paste0("./RDS/sc3/sc3_", 1:50, ".rds")

# Load all RDS files into a list
sc3_list <- lapply(file_paths, readRDS)

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
ggplot(RV_Phi_long, aes(x = method, y = RV_Phi)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of RV_Phi by method",
       x = "Method",
       y = "RV_Phi")


# RV of SigmaPhi
# Extract the RV_SigmaPhi metrics from all loaded files
RV_SigmaPhi_list <- lapply(sc3_list, function(x) x[["metrics"]][["common"]]["RV_SigmaPhi",])
RV_SigmaPhi_list_Tetris <- lapply(Tetris_sc3_list, function(x) x[["metrics"]][["common"]]["RV_SigmaPhi",])
RV_SigmaPhi_list <- mapply(c, RV_SigmaPhi_list, RV_SigmaPhi_list_Tetris, SIMPLIFY = FALSE)
RV_SigmaPhi_combined <- do.call(rbind, RV_SigmaPhi_list) %>% as.data.frame()
RV_SigmaPhi_combined$seed <- (1:5)*3
RV_SigmaPhi_long <- pivot_longer(as.data.frame(RV_SigmaPhi_combined), cols=-seed, 
                                 names_to = "method", values_to = "RV_SigmaPhi")
# Unilist
RV_SigmaPhi_long$RV_SigmaPhi <- sapply(RV_SigmaPhi_long$RV_SigmaPhi, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# RV of Lambda
# Fill with NA
RV_Lambda_combined <- data.frame(matrix(NA, nrow = 50, ncol = 8))
colnames(RV_Lambda_combined) <- names(RV_SigmaPhi_combined)
RV_Lambda_combined$seed <- (1:5)*3
RV_Lambda_long <- pivot_longer(as.data.frame(RV_Lambda_combined), cols=-seed, 
                               names_to = "method", values_to = "RV_Lambda")

# RV of SigmaLambda
# Extract the RV_SigmaLambda metrics from all loaded files
RV_SigmaLambda_list <- lapply(sc3_list, function(x) x[["metrics"]][["specific"]]["RV_SigmaLambda_mean",])
RV_SigmaLambda_combined <- do.call(rbind, RV_SigmaLambda_list) %>% as.data.frame()
RV_SigmaLambda_combined$seed <- (1:50)*3
RV_SigmaLambda_long <- pivot_longer(as.data.frame(RV_SigmaLambda_combined), cols=-seed, 
                                    names_to = "method", values_to = "RV_SigmaLambda")
# Unlist
RV_SigmaLambda_long$RV_SigmaLambda <- sapply(RV_SigmaLambda_long$RV_SigmaLambda, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# RV of SigmaMarginal
RV_SigmaMarginal_list <- lapply(sc3_list, function(x) x[["metrics"]][["marginal"]]["RV_SigmaMarginal_mean",])
RV_SigmaMarginal_combined <- do.call(rbind, RV_SigmaMarginal_list) %>% as.data.frame()
RV_SigmaMarginal_combined$seed <- (1:50)*3
RV_SigmaMarginal_long <- pivot_longer(as.data.frame(RV_SigmaMarginal_combined), cols=-seed, 
                                      names_to = "method", values_to = "RV_SigmaMarginal")
# Unlist
RV_SigmaMarginal_long$RV_SigmaMarginal <- sapply(RV_SigmaMarginal_long$RV_SigmaMarginal, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# Time
Time_list <- lapply(sc3_list, function(x) x[["efficiency"]][[1]]["Elapsed_Time_sec",])
Time_combined <- do.call(rbind, Time_list) %>% as.data.frame()/60
Time_combined$seed <- (1:50)*3
Time_long <- pivot_longer(as.data.frame(Time_combined), cols=-seed, 
                          names_to = "method", values_to = "Time_min")
# Unlist
Time_long$Time_min <- sapply(Time_long$Time_min, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# Peak RAM
RAM_list <- lapply(sc3_list, function(x) x[["efficiency"]][[1]]["Peak_RAM_Used_MiB",])
RAM_combined <- do.call(rbind, RAM_list) %>% as.data.frame()
RAM_combined$seed <- (1:50)*3
RAM_long <- pivot_longer(as.data.frame(RAM_combined), cols=-seed, 
                         names_to = "method", values_to = "RAM_MiB")
# Unlist
RAM_long$RAM_MiB <- sapply(RAM_long$RAM_MiB, function(x) {
  # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
  if (is.logical(x)) return(NA)
  # Return the numeric value if it's a list with one element
  return(as.numeric(x))
})

# Combine RV_Phi, RV_Lambda, RV_SigmaMarginal, Time, and RAM, matching seed and method
df <- merge(RV_Phi_long, RV_Lambda_long, by = c("seed", "method"))
df <- merge(df, RV_SigmaMarginal_long, by = c("seed", "method"))
df <- merge(df, Time_long, by = c("seed", "method"))
df <- merge(df, RAM_long, by = c("seed", "method"))
df_long <- pivot_longer(df, cols = -c(seed, method), 
                        names_to = "metric", values_to = "value")
df_long$scenario <- "Scenario 3"


# Box plot
method_colors <- c(
  "stackFA" = "#66C2A5",   # Bright turquoise
  "IndFA" = "#FC8D62",     # Bright orange for better visibility
  "MOMSS" = "#8DA0CB",     # Light blue
  "BMSFA" = "#E78AC3",     # Bright pink
  "PFA" = "#A6D854",       # Bright green
  "SUFA" = "#FFD92F",      # Bright yellow
  "SUFA_fixJs" = "#E5C494",  # Light brown
  "Tetris_fixT" = "#B3B3B3"  # Light grey
)
# Ordering
df_long_1$method <- factor(df_long_1$method, levels = c("stackFA", "IndFA", "MOMSS", "BMSFA", "PFA", "SUFA"))
df_long_1$metric <- factor(df_long_1$metric, levels = c("RV_Phi", "RV_Lambda", "RV_SigmaMarginal", "Time_min", "RAM_MiB"))

#library(ggh4x)
plot_1 <- 
ggplot(df_long_1, aes(x = method, y = value, colour = method)) +
  geom_boxplot() +
  # ggh4x::facet_grid2(scenario ~ metric, scales = "free", 
  #                    axes = "x", independent = "all") +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey95"),  # Keep major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines (optional)
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "grey95", color = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Rotate x-axis labels
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),  # Smaller legend keys
        legend.text = element_text(size = 8),  # Smaller legend text
        legend.title = element_text(size = 9)) +  # Place the legend on the right-hand side
  labs(#title = "Boxplot of Metrics by Method",
       x = "",
       y = "",
       fill = "method") +
  scale_color_manual(values = method_colors) 



df_long_2$method <- factor(df_long_2$method, levels = c("stackFA", "IndFA", "MOMSS", "BMSFA", "PFA", "SUFA",                                                         "SUFA_fixJs", "Tetris_fixT"))
df_long_2$metric <- factor(df_long_2$metric, levels = c("RV_Phi", "RV_Lambda", "RV_SigmaMarginal", "Time_min", "RAM_MiB"))
plot_2 <- 
  ggplot(df_long_2, aes(x = method, y = value, colour = method)) +
  geom_boxplot() +
  # ggh4x::facet_grid2(scenario ~ metric, scales = "free", 
  #                    axes = "x", independent = "all") +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey95"),  # Keep major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines (optional)
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "grey95", color = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Rotate x-axis labels
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),  # Smaller legend keys
        legend.text = element_text(size = 8),  # Smaller legend text
        legend.title = element_text(size = 9)) +  # Place the legend on the right-hand side
  labs(#title = "Boxplot of Metrics by Method",
    x = "",
    y = "",
    fill = "method") +
  scale_color_manual(values = method_colors) +  # Use consistent color mapping
  guides(colour = guide_legend(override.aes = list(size = 0.5)))



df_long$method <- factor(df_long$method, levels = c("stackFA", "IndFA", "MOMSS", "BMSFA", "PFA", "SUFA"))
df_long$metric <- factor(df_long$metric, levels = c("RV_Phi", "RV_Lambda", "RV_SigmaMarginal", "Time_min", "RAM_MiB"))
plot_3 <- 
  ggplot(df_long, aes(x = method, y = value, colour = method)) +
  geom_boxplot() +
  # ggh4x::facet_grid2(scenario ~ metric, scales = "free", 
  #                    axes = "x", independent = "all") +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey95"),  # Keep major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines (optional)
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "grey95", color = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Rotate x-axis labels
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),  # Smaller legend keys
        legend.text = element_text(size = 8),  # Smaller legend text
        legend.title = element_text(size = 9)  # Smaller legend title
        ) +  # Place the legend on the right-hand side
  labs(#title = "Boxplot of Metrics by Method",
    x = "Method",
    y = "",
    fill = "method") +
  scale_color_manual(values = method_colors) +  # Use consistent color mapping
  guides(colour = guide_legend(override.aes = list(size = 0.5)))



# Scenario setting
sc1_settings <- data.frame(
  "x" = rep(0, 8),
  "y" = 7:0/7,
  "text" = c(
    TeX("Scenario 1", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 100", output="character"),
    TeX("$N_s = (25, 25, 25, 5)$", output = "character"),
    TeX("K = 6", output="character"),
    TeX("and", output="character"),
    TeX("2 knwon covariates", output="character"),
    TeX("$\\alpha_s \\sim MVN$", output="character")
  )
) %>%
  ggplot(aes(x = x, y = y)) +
  geom_text(aes(label = text), parse = TRUE, size = 10/.pt) +
  ylim(c(-0.1,1.1))  +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

sc2_settings <- data.frame(
  "x" = rep(0, 8),
  "y" = 7:0/7,
  "text" = c(
    TeX("Scenario 2", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 100", output="character"),
    TeX("$N_s = (35, 35, 35, 35)$", output = "character"),
    TeX("K = 6", output="character"),
    TeX(" and ", output="character"),
    TeX("Study-specific factors", output="character"),
    TeX("$J_s = (2, 1, 1, 1)$", output="character")
  )
) %>%
  ggplot(aes(x = x, y = y)) +
  geom_text(aes(label = text), parse = TRUE, size = 10/.pt) +
  ylim(c(-0.1,1.1))  +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

sc3_settings <- data.frame(
  "x" = rep(0, 7),
  "y" = 6:0/6,
  "text" = c(
    TeX("Scenario 3", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 100", output="character"),
    TeX("$N_s = (35, 35, 35, 35)$", output = "character"),
    TeX("K = 6", output="character"),
    TeX(" and ", output="character"),
    TeX("Perturbation", output="character")
  )
) %>%
  ggplot(aes(x = x, y = y)) +
  geom_text(aes(label = text), parse = TRUE, size = 10/.pt) +
  ylim(c(-0.1,1.1))  +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

library(patchwork)
layout <- "
ABBBBB
CDDDDD
EFFFFF
"

sc1_settings + plot_1 + sc2_settings + plot_2 + sc3_settings + plot_3 +
  plot_layout(design=layout)

