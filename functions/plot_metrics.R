library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(latex2exp)

plot_metrics <- function(multicore_output) {
  # Convert to data frame
  df <- as.data.frame(multicore_output)
  
  # Boxplot - runtime
  p1 <- df %>% 
    select(run_time_profile_MOMSS, run_time_profile_BMSFA) %>% 
    rename_with(~ gsub("run_time_profile_", "", .x), everything()) %>% 
    pivot_longer(everything(), names_to = "method", values_to = "value") %>%
    unnest(value) %>%
    ggplot(aes(x = method, y = value)) +
    geom_boxplot() +
    labs(title = "Runtime comparison") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    labs(x = "Method", y = "Runtime (seconds)")
  
  # Boxplot - peak_RAM
  p2 <- df %>% 
    select(peak_RAM_profile_MOMSS, peak_RAM_profile_BMSFA) %>% 
    rename_with(~ gsub("peak_RAM_profile_", "", .x), everything()) %>% 
    pivot_longer(everything(), names_to = "method", values_to = "value") %>%
    unnest(value) %>%
    ggplot(aes(x = method, y = value)) +
    geom_boxplot() +
    labs(title = "Peak RAM comparison") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    labs(x = "Method", y = "Peak RAM (MB)")
  
  # Boxplot - RV of Phi
  p3 <- df %>% 
    select(RV_MOMSS_Phi, RV_BMSFA_Phi) %>% 
    rename_with(~ gsub("RV_(.*)_Phi", "\\1", .x), everything()) %>% 
    pivot_longer(everything(), names_to = "method", values_to = "value") %>%
    unnest(value) %>%
    ggplot(aes(x = method, y = value)) +
    geom_boxplot() +
    ggtitle(TeX("RV of $\\Phi$ comparison")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    labs(x = "Method", y = "RV")
  
  # Boxplot - FN of Phi
  p4 <- df %>% 
    select(FN_MOMSS_Phi, FN_BMSFA_Phi) %>% 
    rename_with(~ gsub("FN_(.*)_Phi", "\\1", .x), everything()) %>% 
    pivot_longer(everything(), names_to = "method", values_to = "value") %>%
    unnest(value) %>%
    ggplot(aes(x = method, y = value)) +
    geom_boxplot() +
    ggtitle(TeX("FN of $\\Phi$ comparison")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    labs(x = "Method", y = "FN")
  
  # Boxplot - RV of Sigma Phi
  p5 <- df %>% 
    select(RV_MOMSS_SigmaPhi, RV_BMSFA_SigmaPhi) %>% 
    rename_with(~ gsub("RV_(.*)_SigmaPhi", "\\1", .x), everything()) %>% 
    pivot_longer(everything(), names_to = "method", values_to = "value") %>%
    unnest(value) %>%
    ggplot(aes(x = method, y = value)) +
    geom_boxplot() +
    ggtitle(TeX("RV of $\\Sigma_\\Phi$ comparison")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    labs(x = "Method", y = "RV")
  
  # Boxplot - MSE
  p6 <- df %>% 
    select(MSE_MOMSS, MSE_BMSFA) %>% 
    rename_with(~ gsub("MSE_", "", .x), everything()) %>% 
    pivot_longer(everything(), names_to = "method", values_to = "value") %>%
    unnest(value) %>%
    ggplot(aes(x = method, y = value)) +
    geom_boxplot() +
    ggtitle("MSE comparison") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    labs(x = "Method", y = "MSE")
  
  # Arrange all plots in a grid
  return(grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3))
}

# Example usage
# plot_metrics(readRDS("RDS/sen1_dense_7.9.rds"))

