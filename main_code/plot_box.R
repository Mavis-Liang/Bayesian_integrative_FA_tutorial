library(patchwork)
library(latex2exp)
library(tidyverse)
# This function returns a long format dataframe of a specific metric from RDS files
extract_metric <- function(scenario, metric){
  # Load all RDS files into a list
  file_paths <- paste0("./RDS/sc", scenario, "/sc", scenario, "_", 1:50, ".rds")
  sc_list <- lapply(file_paths, readRDS)
  # Move efficiency table under the metrics sublist
  sc_list <- lapply(sc_list, function(x) {
    x$metrics$efficiency <- x$efficiency[[1]]
    x$efficiency <- NULL
    return(x)
  })
  # Identify which type of the metric belongs to
  if(metric %in% rownames(sc_list[[1]]$metrics$common)){
    type = "common"
  } else if(metric %in% rownames(sc_list[[1]]$metrics$specific)){
    type = "specific"
  } else if(metric %in% rownames(sc_list[[1]]$metrics$marginal)){
    type = "marginal"
  } else if(metric %in% rownames(sc_list[[1]]$metrics$efficiency)){
    type = "efficiency"
  } else{
    stop("The metric is not found in the metrics list")
  }
  # Extract the metrics from all loaded files
  metric_list <- lapply(sc_list, function(x) x[["metrics"]][[type]][metric,])
  metric_long <- bind_rows(metric_list) %>% as.data.frame() %>% 
    mutate(seed = (1:50)*3) %>%
    pivot_longer(cols = -seed, names_to = "method", values_to = "metric_value") %>% 
    mutate(metric_value = sapply(metric_value, function(x) {
      # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
      if (is.logical(x)) return(NA)
      # Return the numeric value if it's a list with one element
      return(as.numeric(x))
    })) %>% 
    dplyr::rename(!!sym(metric) := metric_value)

  # Extract Tetris data if the files exist
  if(file.exists(paste0("./RDS/sc", scenario, "/Tetris"))){
    Tetris_file_paths <- paste0("./RDS/sc", scenario, "/Tetris/sc", scenario, "_Tetris_", 1:50, ".rds")
    Tetris_sc_list <- lapply(Tetris_file_paths, readRDS)
    Tetris_sc_list <- lapply(Tetris_sc_list, function(x) {
      x$metrics$efficiency <- x$efficiency[[1]] %>% t()
      x$efficiency <- NULL
      return(x)
    })
    metric_list_Tetris <- lapply(Tetris_sc_list, function(x) x[["metrics"]][[type]][metric,])
    metric_long_Tetris <- do.call(rbind, metric_list_Tetris) %>% as.data.frame() %>%
      mutate(seed = (1:50)*3,
             method = "Tetris") %>%
      dplyr::rename(metric_value = V1) %>%
      dplyr::select(seed, method, metric_value) %>%
      mutate(metric_value = sapply(metric_value, function(x) {
        # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
        if (is.logical(x)) return(NA)
        # Return the numeric value if it's a list with one element
        return(as.numeric(x))
      }))%>% 
      dplyr::rename(!!sym(metric) := metric_value) ## !! can achieve the value stored in metric
    return(bind_rows(metric_long, metric_long_Tetris))
  } else{
    metrc_long_Tetris <- data.frame(matrix(NA, nrow = 50, ncol = 3)) %>% 
      `colnames<-`(c("seed", "method", metric)) %>%
      mutate(seed = (1:50)*3,
             method = "Tetris") %>%
      mutate(!!sym(metric) := NA)
    return(rbind(metric_long, metrc_long_Tetris))
  }
}

# This function create an NULL dataframe
method_all <- c("stackFA", "IndFA", "PFA", "MOMSS", "SUFA_fixJs", "SUFA", "BMSFA", "Tetris_fixT", "Tetris")
create_null_df <- function(name){
  df <- data.frame(matrix(NA, nrow = 50, ncol = length(method_all))) %>% 
    `colnames<-`(method_all) %>% 
    mutate(seed = (1:50)*3) %>% 
    pivot_longer(cols = -seed, names_to = "method", values_to = name)
  return(df)
}

method_colors <- c(
  "Stack FA" = "#67001f",  
  "Ind FA" = "#b2182b",     
  "PFA" = "#d6604d",       
  "MOM-SS" = "#f4a582",
  "SUFA_fixJs" = "#fddbc7",  
  "SUFA" = "#92c5de",
  "BMSFA" = "#4393c3",     
  "Tetris_fixT" = "#2166ac",
  "Tetris" = "#053061"
)

method_all <- c("stackFA", "IndFA", "PFA", "MOMSS", "SUFA_fixJs", "SUFA", 
                "BMSFA", "Tetris_fixT", "Tetris")
method_lab <- c("Stack FA", "Ind FA", "PFA", "MOM-SS", "SUFA_fixJs", "SUFA", 
                "BMSFA", "Tetris_fixT", "Tetris")
main_metrics <- c("RV_Phi", "RV_Lambda_mean", "RV_SigmaMarginal_mean",
                  "Elapsed_Time_sec", "Peak_RAM_Used_MiB")
supp_metrics <- c("FN_Phi", "FN_SigmaPhi", "RV_SigmaPhi", "FN_Lambda_mean",
                  "FN_SigmaLambda_mean", "RV_SigmaLambda_mean", "FN_SigmaMarginal_mean")
labels_main <- c(TeX("RV of $(\\Phi,\\widehat{\\Phi})$"), 
                          TeX("RV of $(\\Lambda_s, \\widehat{\\Lambda_s})$"), 
                          TeX("RV of $(\\Sigma_{s}, \\widehat{\\Sigma_{s}})$"), 
                          "Time (sec)",
                          "RAM (MiB)")
labels_supp <- c(TeX(paste0("FN of ", "$(\\Phi, \\widehat{\\Phi})$")), 
                          TeX(paste0("FN of ", "$(\\Sigma_{\\Phi}, \\widehat{\\Sigma_{\\Phi}})$")), 
                          TeX(paste0("RV of ", "$(\\Sigma_{\\Phi}, \\widehat{\\Sigma_{\\Phi}})$")), 
                          TeX(paste0("FN of ", "$(\\Lambda_s, \\widehat{\\Lambda_s})$")), 
                          TeX(paste0("FN of ", "$(\\Sigma_{\\Lambda_s}, \\widehat{\\Sigma_{\\Lambda_s}})$")), 
                          TeX(paste0("RV of ", "$(\\Sigma_{\\Lambda_s}, \\widehat{\\Sigma_{\\Lambda_s}})$")), 
                          TeX(paste0("FN of ", "$(\\Sigma_{s}, \\widehat{\\Sigma_{s}})$"))
                          )

# ###############################Scenario 1################################
# Metrics in the main text
RV_Phi_sc1 <- extract_metric(metric = "RV_Phi", scenario = 1)
RV_Lambda_sc1 <- create_null_df("RV_Lambda_mean")
RV_SigmaMarginal_sc1 <- extract_metric(metric = "RV_SigmaMarginal_mean", scenario = 1)
Time_sc1 <- extract_metric(metric = "Elapsed_Time_sec", scenario = 1)
RAM_sc1 <- extract_metric(metric = "Peak_RAM_Used_MiB", scenario = 1)

# Metrics in the supplementary material
FN_Phi_sc1 <- extract_metric(metric = "FN_Phi", scenario = 1)
FN_SigmaPhi_sc1 <- extract_metric(metric = "FN_SigmaPhi", scenario = 1)
RV_SigmaPhi_sc1 <- extract_metric(metric = "RV_SigmaPhi", scenario = 1)
FN_Lambda_sc1 <- create_null_df("FN_Lambda_mean")
FN_SigmaLambda_mean_sc1 <- create_null_df("FN_SigmaLambda_mean")
RV_SigmaLambda_mean_sc1 <- create_null_df("RV_SigmaLambda_mean")
FN_SigmaMarginal_mean_sc1 <- extract_metric(metric = "FN_SigmaMarginal_mean", scenario = 1)

# merge
df_1_supp <- merge(FN_Phi_sc1, FN_SigmaPhi_sc1) %>%
  merge(RV_SigmaPhi_sc1) %>%
  merge(FN_Lambda_sc1) %>%
  merge(FN_SigmaLambda_mean_sc1) %>%
  merge(RV_SigmaLambda_mean_sc1) %>%
  merge(FN_SigmaMarginal_mean_sc1) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = supp_metrics,
                       labels = labels_supp),
         method = factor(method, levels = method_all, labels = method_lab))
df_1_eff <- merge(Time_sc1, RAM_sc1) %>% 
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_1_acc <- merge(RV_Phi_sc1, RV_Lambda_sc1) %>% 
  merge(RV_SigmaMarginal_sc1) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))

library(scales)

# Scenario setting
sc1_settings <- data.frame(
  "x" = rep(0, 7),
  "y" = 6:0/6,
  "text" = c(
    TeX("Scenario 2", output="character"),
    TeX("Based on MOM-SS", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 40", output="character"),
    TeX("$N_s = (100, 100, 100, 100)$", output = "character"),
    TeX("K = 4, Known covariates Q=2", output="character"),
    TeX("plus study-specific intercepts", output="character")
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
    axis.line = element_line(colour = "black", size=1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# ---------------------Scenario 2 -------------------------------------
# Metrics in the main text
RV_Phi_sc2 <- extract_metric(metric = "RV_Phi", scenario = 2)
RV_Lambda_sc2 <- extract_metric(metric = "RV_Lambda_mean", scenario = 2)
RV_SigmaMarginal_sc2 <- extract_metric(metric = "RV_SigmaMarginal_mean", scenario = 2)
Time_sc2 <- extract_metric(metric = "Elapsed_Time_sec", scenario = 2) 
RAM_sc2 <- extract_metric(metric = "Peak_RAM_Used_MiB", scenario = 2)

# Metrics in the supplementary material
FN_Phi_sc2 <- extract_metric(metric = "FN_Phi", scenario = 2)
FN_SigmaPhi_sc2 <- extract_metric(metric = "FN_SigmaPhi", scenario = 2)
RV_SigmaPhi_sc2 <- extract_metric(metric = "RV_SigmaPhi", scenario = 2)
FN_Lambda_sc2 <- extract_metric(metric = "FN_Lambda_mean", scenario = 2)
FN_SigmaLambda_mean_sc2 <- extract_metric(metric = "FN_SigmaLambda_mean", scenario = 2)
RV_SigmaLambda_mean_sc2 <- extract_metric(metric = "RV_SigmaLambda_mean", scenario = 2)
FN_SigmaMarginal_mean_sc2 <- extract_metric(metric = "FN_SigmaMarginal_mean", scenario = 2)

# merge
df_2_acc <- merge(RV_Phi_sc2, RV_Lambda_sc2) %>% 
  merge(RV_SigmaMarginal_sc2) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_2_eff <- merge(Time_sc2, RAM_sc2) %>% 
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))

df_2_supp <- merge(FN_Phi_sc2, FN_SigmaPhi_sc2) %>%
  merge(RV_SigmaPhi_sc2) %>%
  merge(FN_Lambda_sc2) %>%
  merge(FN_SigmaLambda_mean_sc2) %>%
  merge(RV_SigmaLambda_mean_sc2) %>%
  merge(FN_SigmaMarginal_mean_sc2) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = supp_metrics,
                       labels = labels_supp),
         method = factor(method, levels = method_all, labels = method_lab)
  )

sc2_settings <- data.frame(
  "x" = rep(0, 6),
  "y" = 5:0/5,
  "text" = c(
    TeX("Scenario 2", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 40", output="character"),
    TeX("$N_s = (100, 100, 100, 100)$", output = "character"),
    TeX("follwing BMSFA", output="character"),
    TeX("K = 4 and J_s=(1, 1, 1, 1)", output="character")
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
    axis.line = element_line(colour = "black", size=1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# ---------------------Scenario 3 -------------------------------------
# Metrics in the main text
RV_Phi <- extract_metric(metric = "RV_Phi", scenario = 3)
RV_Lambda <- extract_metric(metric = "RV_Lambda_mean", scenario = 3)
RV_SigmaMarginal <- extract_metric(metric = "RV_SigmaMarginal_mean", scenario = 3)
Time <- extract_metric(metric = "Elapsed_Time_sec", scenario = 3) 
RAM <- extract_metric(metric = "Peak_RAM_Used_MiB", scenario = 3)

# Metrics in the supplementary material
FN_Phi <- extract_metric(metric = "FN_Phi", scenario = 3) 
FN_SigmaPhi <- extract_metric(metric = "FN_SigmaPhi", scenario = 3)
RV_SigmaPhi <- extract_metric(metric = "RV_SigmaPhi", scenario = 3) 
FN_Lambda <- create_null_df("FN_Lambda_mean")
FN_SigmaLambda_mean <- extract_metric(metric = "FN_SigmaLambda_mean", scenario = 3)
RV_SigmaLambda_mean <- extract_metric(metric = "RV_SigmaLambda_mean", scenario = 3)
FN_SigmaMarginal_mean <- extract_metric(metric = "FN_SigmaMarginal_mean", scenario = 3)

# Merge
df_3_acc <- merge(RV_Phi, RV_Lambda) %>% 
  merge(RV_SigmaMarginal) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_3_eff <- merge(Time, RAM) %>% 
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))

df_3_supp <- merge(FN_Phi, FN_SigmaPhi) %>%
  merge(RV_SigmaPhi) %>%
  merge(FN_Lambda) %>%
  merge(FN_SigmaLambda_mean) %>%
  merge(RV_SigmaLambda_mean) %>%
  merge(FN_SigmaMarginal_mean) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = supp_metrics,
                       labels = labels_supp),
         method = factor(method, levels = method_all, labels = method_lab))

sc3_settings <- data.frame(
  "x" = rep(0, 6),
  "y" = 5:0/5,
  "text" = c(
    TeX("Scenario 1", output="character"),
    TeX("Based on PFA", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 40", output="character"),
    TeX("$N_s = (100, 100, 100, 100)$", output = "character"),
    TeX("K = 4 and Perturbations", output="character")
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
    axis.line = element_line(colour = "black", size=1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


# ---------------------Scenario 4 -------------------------------------
# Metrics in the main text
RV_Phi_sc4 <- extract_metric(metric = "RV_Phi", scenario = 4)
RV_Lambda_sc4 <- extract_metric(metric = "RV_Lambda_mean", scenario = 4)
RV_SigmaMarginal_sc4 <- extract_metric(metric = "RV_SigmaMarginal_mean", scenario = 4)
Time_sc4 <- extract_metric(metric = "Elapsed_Time_sec", scenario = 4) 
RAM_sc4 <- extract_metric(metric = "Peak_RAM_Used_MiB", scenario = 4)

# Metrics in the supplementary material
FN_Phi_sc4 <- extract_metric(metric = "FN_Phi", scenario = 4)
FN_SigmaPhi_sc4 <- extract_metric(metric = "FN_SigmaPhi", scenario = 4)
RV_SigmaPhi_sc4 <- extract_metric(metric = "RV_SigmaPhi", scenario = 4)
FN_Lambda_sc4 <- extract_metric(metric = "FN_Lambda_mean", scenario = 4)
FN_SigmaLambda_mean_sc4 <- extract_metric(metric = "FN_SigmaLambda_mean", scenario = 4)
RV_SigmaLambda_mean_sc4 <- extract_metric(metric = "RV_SigmaLambda_mean", scenario = 4)
FN_SigmaMarginal_mean_sc4 <- extract_metric(metric = "FN_SigmaMarginal_mean", scenario = 4)

# Merge
df_4_acc <- merge(RV_Phi_sc4, RV_Lambda_sc4) %>% 
  merge(RV_SigmaMarginal_sc4) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_4_eff <- merge(Time_sc4, RAM_sc4) %>% 
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_4_supp <- merge(FN_Phi_sc4, FN_SigmaPhi_sc4) %>%
  merge(RV_SigmaPhi_sc4) %>%
  merge(FN_Lambda_sc4) %>%
  merge(FN_SigmaLambda_mean_sc4) %>%
  merge(RV_SigmaLambda_mean_sc4) %>%
  merge(FN_SigmaMarginal_mean_sc4) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = supp_metrics,
                       labels = labels_supp),
         method = factor(method, levels = method_all, labels = method_lab))

sc4_settings <- data.frame(
  "x" = rep(0, 7),
  "y" = 6:0/6,
  "text" = c(
    TeX("Scenario 4 (mimic nutrition)", output="character"),
    TeX("Based on Tetris", output="character"),
    TeX("12 knwon covariates", output="character"),
    TeX("S = 12", output="character"),
    TeX("P = 42", output="character"),
    TeX("$ 373 <= N_s <= 3775$", output = "character"),
    TeX("K = 4, J_s= (1, ... 1), and 7 partial", output="character")
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
    axis.line = element_line(colour = "black", size=1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# ------------------------- Scenario 6 -------------------------
# Metrics in the main text
RV_Phi_sc6 <- extract_metric(metric = "RV_Phi", scenario = 6)
RV_Lambda_sc6 <- extract_metric(metric = "RV_Lambda_mean", scenario = 6)
RV_SigmaMarginal_sc6 <- extract_metric(metric = "RV_SigmaMarginal_mean", scenario = 6)
Time_sc6 <- extract_metric(metric = "Elapsed_Time_sec", scenario = 6) 
RAM_sc6 <- extract_metric(metric = "Peak_RAM_Used_MiB", scenario = 6)

# Metrics in the supplementary material
FN_Phi_sc6 <- extract_metric(metric = "FN_Phi", scenario = 6)
FN_SigmaPhi_sc6 <- extract_metric(metric = "FN_SigmaPhi", scenario = 6)
RV_SigmaPhi_sc6 <- extract_metric(metric = "RV_SigmaPhi", scenario = 6)
FN_Lambda_sc6 <- extract_metric(metric = "FN_Lambda_mean", scenario = 6)
FN_SigmaLambda_mean_sc6 <- extract_metric(metric = "FN_SigmaLambda_mean", scenario = 6)
RV_SigmaLambda_mean_sc6 <- extract_metric(metric = "RV_SigmaLambda_mean", scenario = 6)
FN_SigmaMarginal_mean_sc6 <- extract_metric(metric = "FN_SigmaMarginal_mean", scenario = 6)

# Merge
df_6_acc <- merge(RV_Phi_sc6, RV_Lambda_sc6) %>% 
  merge(RV_SigmaMarginal_sc6) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_6_eff <- merge(Time_sc6, RAM_sc6) %>% 
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_6_supp <- merge(FN_Phi_sc6, FN_SigmaPhi_sc6) %>%
  merge(RV_SigmaPhi_sc6) %>%
  merge(FN_Lambda_sc6) %>%
  merge(FN_SigmaLambda_mean_sc6) %>%
  merge(RV_SigmaLambda_mean_sc6) %>%
  merge(FN_SigmaMarginal_mean_sc6) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = supp_metrics,
                       labels = labels_supp),
         method = factor(method, levels = method_all, labels = method_lab))

sc6_settings <- data.frame(
  "x" = rep(0, 6),
  "y" = 5:0/5,
  "text" = c(
    TeX("Scenario 3", output="character"),
    TeX("Based on SUFA", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 40", output="character"),
    TeX("$N_s=(100, 100, 100, 100)$", output = "character"),
    TeX("K = 4 and J_s=(1, 1, 1, 1)", output="character")
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
    axis.line = element_line(colour = "black", size=1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


### ------------------------- Main plot - accuracy -------------------------
acc_6 <-
  df_6_acc %>% 
  ggplot(aes(x = method, y = value, colour = method)) +
  geom_boxplot() +
  facet_wrap(~metric, nrow = 1, labeller = label_parsed) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey95"),  # Keep major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines (optional)
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "grey95", color = NA), 
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),  # Rotate x-axis labels
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"),  # Smaller legend keys
        legend.text = element_text(size = 8),  # Smaller legend text
        plot.margin = margin(t = 2, r = 2, b = 0, l = 2),
        legend.title = element_text(size = 9)) +  # Place the legend on the right-hand side
  labs(title = "",
    x = "",
    y = "Value",
    fill = "Method") +
  scale_color_manual(values = method_colors)  # Use consistent color mapping

layout_acc <- "
ABBB
CDDD
EFFF
GHHH
IJJJ
"

(sc3_settings + acc_3 + 
  sc1_settings + acc_1 + 
  sc6_settings + acc_6 +
  sc4_settings + acc_4 +
  settings_sc5 + acc_5) +
  plot_layout(design=layout_acc)

ggplot2::ggsave("./Figs/box_acc_ylimFixed.pdf", width = 10, height = 14, units = "in", dpi = 300)


##------------------------------- Efficiency plots ---------------------------
eff_3<- df_3_eff %>%
  ggplot(aes(x = method, y = value, colour = method)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(
    ~metric,
    scales = "free",
    nrow = 1,
    labeller = label_parsed,
    strip.position = "top"
  ) +
  # scale_y_continuous(
  #   labels = scales::trans_format("identity", scales::math_format()),
  #   breaks = scales::trans_breaks("identity", identity)
  # ) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
    legend.position = "none",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  ) +
  labs(title = "Scenario 1",
    x = "",
    y = "Value",
    fill = "Method", color = "Method"
  ) +
  scale_color_manual(values = method_colors) +
  guides(colour = guide_legend(override.aes = list(size = 0.5))) 
# +
#   scale_y_log10()


layout_eff <- "
AADD
BBEE
CCFF
"
eff_3 + eff_1 + eff_6 + eff_4 + eff_5 + 
  plot_layout(design=layout_eff)
ggplot2::ggsave("./Figs/box_eff_originScale.pdf", width = 10, height = 9, units = "in", dpi = 300)

legend_plot <- cowplot::get_legend(eff_6)

# ---------------------------Supplement plots --------------------------------
supp_3 <- 
  df_3_supp %>%
  ggplot(aes(x = method, y = value, colour = method)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free", nrow = 1, labeller = label_parsed) +
  theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey95"),  # Keep major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines (optional)
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "grey95", color = NA), 
        axis.text.x = element_text(angle = 30, hjust = 1, size = 5),  # Rotate x-axis labels
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"),  # Smaller legend keys
        legend.text = element_text(size = 8),  # Smaller legend text
        legend.title = element_text(size = 9)) +  # Place the legend on the right-hand side
  labs(title = "Scenario 1",
       x = "",
       y = "Value",
       fill = "Method") +
  scale_color_manual(values = method_colors) +  # Use consistent color mapping
  guides(colour = guide_legend(override.aes = list(size = 0.5)))

layout_supp <- "
A
B
C
D
E
"
supp_3 + supp_1 + supp_6 + supp_4 + supp_5 + 
  plot_layout(design=layout_supp)
ggplot2::ggsave("./Figs/box_supp.pdf", width = 11, height = 12, units = "in", dpi = 300)

