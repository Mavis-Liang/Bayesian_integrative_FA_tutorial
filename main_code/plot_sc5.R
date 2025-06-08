# Scenario 5 (mimicing genomic data) has a large size. Therefore, we separately handle this scenario.
extract_metric_sc6 <- function(models, metric, seedRange = 1:50){
  file_paths <- paste0("./RDS/sc5", "/sc5_", models, "/sc5_",models, "_", seedRange, ".rds")
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
  metric_long <- do.call(rbind, metric_list) %>% as.data.frame() %>% 
    mutate(seed = (seedRange)*3) %>%
    pivot_longer(cols = -seed, names_to = "method", values_to = "metric_value") %>% 
    mutate(metric_value = sapply(metric_value, function(x) {
      # Convert logical values (e.g., TRUE/FALSE or NULL) to NA
      if (is.logical(x)) return(NA)
      # Return the numeric value if it's a list with one element
      return(as.numeric(x))
    })) %>% 
    dplyr::rename(!!sym(metric) := metric_value)
  # If there is only one method, than rename the method column
  if(length(unique(metric_long$method)) == 1){
    metric_long$method <- models
  }
  return(metric_long)
}

# This function create an NULL dataframe
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

# This is different from the plot_box
method_all <- c("StackFA", "IndFA", "PFA", "MOMSS", "SUFAfixJ", "SUFA", 
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

# Stack FA
RV_Phi_stackFA <- extract_metric_sc6(metric = "RV_Phi", models = "StackFA")
RV_Lambda_stackFA <- create_null_df("RV_Lambda_mean")
RV_SigmaMarginal_stackFA <- extract_metric_sc6(metric = "RV_SigmaMarginal_mean", models = "StackFA")
Time_stackFA <- extract_metric_sc6(metric = "Elapsed_Time_sec", models = "StackFA")
RAM_stackFA <- extract_metric_sc6(metric = "Peak_RAM_Used_MiB", models = "StackFA")

FN_Phi_stackFA <- extract_metric_sc6(metric = "FN_Phi", models = "StackFA")
FN_SigmaPhi_stackFA <- extract_metric_sc6(metric = "FN_SigmaPhi", models = "StackFA")
RV_SigmaPhi_stackFA <- extract_metric_sc6(metric = "RV_SigmaPhi", models = "StackFA")
FN_Lambda_stackFA <- create_null_df("FN_Lambda_mean")
FN_SigmaLambda_mean_stackFA <- create_null_df("FN_SigmaLambda_mean")
RV_SigmaLambda_mean_stackFA <- create_null_df("RV_SigmaLambda_mean")
FN_SigmaMarginal_mean_stackFA <- extract_metric_sc6(metric = "FN_SigmaMarginal_mean", models = "StackFA")


# MOM-SS and BMSFA
RV_Phi_MOMSSBMSFA <- extract_metric_sc6(metric = "RV_Phi", models = "MOMSS_BMSFA")
RV_Lambda_MOMSSBMSFA <- extract_metric_sc6(metric = "RV_Lambda_mean", models = "MOMSS_BMSFA")
RV_SigmaMarginal_MOMSSBMSFA <- extract_metric_sc6(metric = "RV_SigmaMarginal_mean", models = "MOMSS_BMSFA")
Time_MOMSSBMSFA <- extract_metric_sc6(metric = "Elapsed_Time_sec", models = "MOMSS_BMSFA")
RAM_MOMSSBMSFA <- extract_metric_sc6(metric = "Peak_RAM_Used_MiB", models = "MOMSS_BMSFA")

FN_Phi_MOMSSBMSFA <- extract_metric_sc6(metric = "FN_Phi", models = "MOMSS_BMSFA")
FN_SigmaPhi_MOMSSBMSFA <- extract_metric_sc6(metric = "FN_SigmaPhi", models = "MOMSS_BMSFA")
RV_SigmaPhi_MOMSSBMSFA <- extract_metric_sc6(metric = "RV_SigmaPhi", models = "MOMSS_BMSFA")
FN_Lambda_MOMSSBMSFA <- extract_metric_sc6(metric = "FN_Lambda_mean", models = "MOMSS_BMSFA")
FN_SigmaLambda_mean_MOMSSBMSFA <- extract_metric_sc6(metric = "FN_SigmaLambda_mean", models = "MOMSS_BMSFA")
RV_SigmaLambda_mean_MOMSSBMSFA <- extract_metric_sc6(metric = "RV_SigmaLambda_mean", models = "MOMSS_BMSFA")
FN_SigmaMarginal_mean_MOMSSBMSFA <- extract_metric_sc6(metric = "FN_SigmaMarginal_mean", models = "MOMSS_BMSFA")

# SUFA and Tetris_fixT (because SUFA takes so long so we only run 30 times)
RV_Phi_SUFA_Tetris <- extract_metric_sc6(metric = "RV_Phi", models = "SUFA_Tetris", seedRange = 1:29)
RV_Lambda_SUFA_Tetris <- extract_metric_sc6(metric = "RV_Lambda_mean", models = "SUFA_Tetris", seedRange = 1:29)
RV_SigmaMarginal_SUFA_Tetris <- extract_metric_sc6(metric = "RV_SigmaMarginal_mean", models = "SUFA_Tetris", seedRange = 1:29)
Time_SUFA_Tetris <- extract_metric_sc6(metric = "Elapsed_Time_sec", models = "SUFA_Tetris", seedRange = 1:29)
RAM_SUFA_Tetris <- extract_metric_sc6(metric = "Peak_RAM_Used_MiB", models = "SUFA_Tetris", seedRange = 1:29)

RV_SigmaPhi_SUFA_Tetris <- extract_metric_sc6(metric = "RV_SigmaPhi", models = "SUFA_Tetris", seedRange = 1:29)
RV_SigmaLambda_mean_SUFA_Tetris <- extract_metric_sc6(metric = "RV_SigmaLambda_mean", models = "SUFA_Tetris", seedRange = 1:29)
FN_Phi_SUFA_Tetris <- extract_metric_sc6(metric = "FN_Phi", models = "SUFA_Tetris", seedRange = 1:29)
FN_SigmaPhi_SUFA_Tetris <- extract_metric_sc6(metric = "FN_SigmaPhi", models = "SUFA_Tetris", seedRange = 1:29)
FN_Lambda_SUFA_Tetris <- extract_metric_sc6(metric = "FN_Lambda_mean", models = "SUFA_Tetris", seedRange = 1:29)
FN_SigmaLambda_mean_SUFA_Tetris <- extract_metric_sc6(metric = "FN_SigmaLambda_mean", models = "SUFA_Tetris", seedRange = 1:29)
FN_SigmaMarginal_mean_SUFA_Tetris <- extract_metric_sc6(metric = "FN_SigmaMarginal_mean", models = "SUFA_Tetris", seedRange = 1:29)

RV_Phi_TetrixfixT <- extract_metric_sc6(metric = "RV_Phi", 
                                        models = "TetrisfixT", seedRange = 31:50) %>% 
  mutate(method = "Tetris_fixT")
RV_Lambda_TetrixfixT <- extract_metric_sc6(metric = "RV_Lambda_mean", 
                                           models = "TetrisfixT",
                                           seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
RV_SigmaMarginal_TetrixfixT <- extract_metric_sc6(metric = "RV_SigmaMarginal_mean",
                                                  models = "TetrisfixT",
                                                  seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
Time_TetrixfixT <- extract_metric_sc6(metric = "Elapsed_Time_sec", 
                                      models = "TetrisfixT",
                                      seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
RAM_TetrixfixT <- extract_metric_sc6(metric = "Peak_RAM_Used_MiB", 
                                     models = "TetrisfixT",
                                     seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")

FN_Phi_TetrixfixT <- extract_metric_sc6(metric = "FN_Phi", 
                                         models = "TetrisfixT", seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
FN_SigmaPhi_TetrixfixT <- extract_metric_sc6(metric = "FN_SigmaPhi", 
                                               models = "TetrisfixT", seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
FN_Lambda_TetrixfixT <- extract_metric_sc6(metric = "FN_Lambda_mean", 
                                              models = "TetrisfixT", seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
FN_SigmaLambda_mean_TetrixfixT <- extract_metric_sc6(metric = "FN_SigmaLambda_mean", 
                                                        models = "TetrisfixT", seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
FN_SigmaMarginal_mean_TetrixfixT <- extract_metric_sc6(metric = "FN_SigmaMarginal_mean", 
                                                           models = "TetrisfixT", seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
RV_SigmaPhi_TetrixfixT <- extract_metric_sc6(metric = "RV_SigmaPhi", 
                                                models = "TetrisfixT", seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")
RV_SigmaLambda_mean_TetrixfixT <- extract_metric_sc6(metric = "RV_SigmaLambda_mean", 
                                                          models = "TetrisfixT", seedRange = 31:50) %>%
  mutate(method = "Tetris_fixT")

# IndFA
RV_Phi_IndFA <- create_null_df("RV_Phi")
RV_Lambda_IndFA <- extract_metric_sc6(metric = "RV_Lambda_mean", models = "IndFA")
RV_SigmaMarginal_IndFA <- extract_metric_sc6(metric = "RV_SigmaMarginal_mean", models = "IndFA")
RAM_IndFA <- extract_metric_sc6(metric = "Peak_RAM_Used_MiB", models = "IndFA")
Time_IndFA <- extract_metric_sc6(metric = "Elapsed_Time_sec", models = "IndFA")

RV_SigmaPhi_IndFA <- create_null_df("RV_SigmaPhi")
RV_SigmaLambda_mean_IndFA <- extract_metric_sc6(metric = "RV_SigmaLambda_mean", models = "IndFA")
FN_Phi_IndFA <- create_null_df("FN_Phi")
FN_SigmaPhi_IndFA <- create_null_df("FN_SigmaPhi")
FN_Lambda_IndFA <- extract_metric_sc6(metric = "FN_Lambda_mean", models = "IndFA")
FN_SigmaLambda_mean_IndFA <- extract_metric_sc6(metric = "FN_SigmaLambda_mean", models = "IndFA")
FN_SigmaMarginal_mean_IndFA <- extract_metric_sc6(metric = "FN_SigmaMarginal_mean", models = "IndFA")

# SUFA_fixJs
RV_Phi_SUFAfixJ <- extract_metric_sc6(metric = "RV_Phi", models = "SUFAfixJ", seedRange = 1:30)
RV_Lambda_SUFAfixJ <- extract_metric_sc6(metric = "RV_Lambda_mean", models = "SUFAfixJ", seedRange = 1:30)
RV_SigmaMarginal_SUFAfixJ <- extract_metric_sc6(metric = "RV_SigmaMarginal_mean", models = "SUFAfixJ", seedRange = 1:30)
Time_SUFAfixJ <- extract_metric_sc6(metric = "Elapsed_Time_sec", models = "SUFAfixJ", seedRange = 1:30)
RAM_SUFAfixJ <- extract_metric_sc6(metric = "Peak_RAM_Used_MiB", models = "SUFAfixJ", seedRange = 1:30)

RV_SigmaPhi_SUFAfixJ <- extract_metric_sc6(metric = "RV_SigmaPhi", models = "SUFAfixJ", seedRange = 1:30)
RV_SigmaLambda_mean_SUFAfixJ <- extract_metric_sc6(metric = "RV_SigmaLambda_mean", models = "SUFAfixJ", seedRange = 1:30)
FN_Phi_SUFAfixJ <- extract_metric_sc6(metric = "FN_Phi", models = "SUFAfixJ", seedRange = 1:30)
FN_SigmaPhi_SUFAfixJ <- extract_metric_sc6(metric = "FN_SigmaPhi", models = "SUFAfixJ", seedRange = 1:30)
FN_Lambda_SUFAfixJ <- extract_metric_sc6(metric = "FN_Lambda_mean", models = "SUFAfixJ", seedRange = 1:30)
FN_SigmaLambda_mean_SUFAfixJ <- extract_metric_sc6(metric = "FN_SigmaLambda_mean", models = "SUFAfixJ", seedRange = 1:30)
FN_SigmaMarginal_mean_SUFAfixJ <- extract_metric_sc6(metric = "FN_SigmaMarginal_mean", models = "SUFAfixJ", seedRange = 1:30)

# Method exceed 24 hours
RV_Phi_PFA <- create_null_df("RV_Phi")
RV_Lambda_PFA <- create_null_df("RV_Lambda_mean")
RV_SigmaMarginal_PFA <- create_null_df("RV_SigmaMarginal_mean")
Time_PFA <- create_null_df("Elapsed_Time_sec")
RAM_PFA <- create_null_df("Peak_RAM_Used_MiB")
FN_Phi_PFA <- create_null_df("FN_Phi")
FN_SigmaPhi_PFA <- create_null_df("FN_SigmaPhi")
FN_Lambda_PFA <- create_null_df("FN_Lambda_mean")
FN_SigmaLambda_mean_PFA <- create_null_df("FN_SigmaLambda_mean")
FN_SigmaMarginal_mean_PFA <- create_null_df("FN_SigmaMarginal_mean")
RV_SigmaPhi_PFA <- create_null_df("RV_SigmaPhi")
RV_SigmaLambda_mean_PFA <- create_null_df("RV_SigmaLambda_mean")

RV_Phi_Tetris <- create_null_df("RV_Phi")
RV_Lambda_Tetris <- create_null_df("RV_Lambda_mean")
RV_SigmaMarginal_Tetris <- create_null_df("RV_SigmaMarginal_mean")
Time_Tetris <- create_null_df("Elapsed_Time_sec")
RAM_Tetris <- create_null_df("Peak_RAM_Used_MiB")
FN_Phi_Tetris <- create_null_df("FN_Phi")
FN_SigmaPhi_Tetris <- create_null_df("FN_SigmaPhi")
FN_Lambda_Tetris <- create_null_df("FN_Lambda_mean")
FN_SigmaLambda_mean_Tetris <- create_null_df("FN_SigmaLambda_mean")
FN_SigmaMarginal_mean_Tetris <- create_null_df("FN_SigmaMarginal_mean")
RV_SigmaPhi_Tetris <- create_null_df("RV_SigmaPhi")
RV_SigmaLambda_mean_Tetris <- create_null_df("RV_SigmaLambda_mean")

# bind row
RV_Phi <- bind_rows(RV_Phi_stackFA, RV_Phi_IndFA, RV_Phi_SUFA_Tetris, 
                    RV_Phi_TetrixfixT, RV_Phi_MOMSSBMSFA, RV_Phi_SUFAfixJ,
                    RV_Phi_PFA, RV_Phi_Tetris)
RV_Lambda <- bind_rows(RV_Lambda_stackFA, RV_Lambda_IndFA, RV_Lambda_SUFA_Tetris, 
                       RV_Lambda_TetrixfixT, RV_Lambda_MOMSSBMSFA, RV_Lambda_SUFAfixJ,
                       RV_Lambda_PFA, RV_Lambda_Tetris)
RV_SigmaMarginal <- bind_rows(RV_SigmaMarginal_stackFA, RV_SigmaMarginal_IndFA, RV_SigmaMarginal_SUFA_Tetris, 
                              RV_SigmaMarginal_TetrixfixT, RV_SigmaMarginal_MOMSSBMSFA, RV_SigmaMarginal_SUFAfixJ,
                              RV_SigmaMarginal_PFA, RV_SigmaMarginal_Tetris)
Time <- bind_rows(Time_stackFA, Time_IndFA, Time_SUFA_Tetris,
                  Time_TetrixfixT, Time_MOMSSBMSFA, Time_SUFAfixJ,
                  RAM_PFA, Time_Tetris)
RAM <- bind_rows(RAM_stackFA, RAM_IndFA, RAM_SUFA_Tetris,
                 RAM_TetrixfixT, RAM_MOMSSBMSFA, RAM_SUFAfixJ,
                 RAM_PFA, RAM_Tetris)
FN_Phi_sc5 <- bind_rows(FN_Phi_stackFA, FN_Phi_IndFA, FN_Phi_SUFA_Tetris, 
                    FN_Phi_TetrixfixT, FN_Phi_MOMSSBMSFA, FN_Phi_SUFAfixJ,
                    FN_Phi_PFA, FN_Phi_Tetris)
FN_SigmaPhi_sc5 <- bind_rows(FN_SigmaPhi_stackFA, FN_SigmaPhi_IndFA, FN_SigmaPhi_SUFA_Tetris, 
                          FN_SigmaPhi_TetrixfixT, FN_SigmaPhi_MOMSSBMSFA, FN_SigmaPhi_SUFAfixJ,
                          FN_SigmaPhi_PFA, FN_SigmaPhi_Tetris)
FN_Lambda_sc5 <- bind_rows(FN_Lambda_stackFA, FN_Lambda_IndFA, FN_Lambda_SUFA_Tetris,
                       FN_Lambda_TetrixfixT, FN_Lambda_MOMSSBMSFA, FN_Lambda_SUFAfixJ,
                       FN_Lambda_PFA, FN_Lambda_Tetris)
FN_SigmaLambda_mean_sc5 <- bind_rows(FN_SigmaLambda_mean_stackFA, FN_SigmaLambda_mean_IndFA, FN_SigmaLambda_mean_SUFA_Tetris,
                                  FN_SigmaLambda_mean_TetrixfixT, FN_SigmaLambda_mean_MOMSSBMSFA, FN_SigmaLambda_mean_SUFAfixJ,
                                  FN_SigmaLambda_mean_PFA, FN_SigmaLambda_mean_Tetris)
FN_SigmaMarginal_mean_sc5 <- bind_rows(FN_SigmaMarginal_mean_stackFA, FN_SigmaMarginal_mean_IndFA, FN_SigmaMarginal_mean_SUFA_Tetris,
                                    FN_SigmaMarginal_mean_TetrixfixT, FN_SigmaMarginal_mean_MOMSSBMSFA, FN_SigmaMarginal_mean_SUFAfixJ,
                                    FN_SigmaMarginal_mean_PFA, FN_SigmaMarginal_mean_Tetris)
RV_SigmaPhi_sc5 <- bind_rows(RV_SigmaPhi_stackFA, RV_SigmaPhi_IndFA, RV_SigmaPhi_SUFA_Tetris,
                         RV_SigmaPhi_TetrixfixT, RV_SigmaPhi_MOMSSBMSFA, RV_SigmaPhi_SUFAfixJ,
                         RV_SigmaPhi_PFA, RV_SigmaPhi_Tetris)
RV_SigmaLambda_mean_sc5 <- bind_rows(RV_SigmaLambda_mean_stackFA, RV_SigmaLambda_mean_IndFA, RV_SigmaLambda_mean_SUFA_Tetris,
                                  RV_SigmaLambda_mean_TetrixfixT, RV_SigmaLambda_mean_MOMSSBMSFA, RV_SigmaLambda_mean_SUFAfixJ,
                                  RV_SigmaLambda_mean_PFA, RV_SigmaLambda_mean_Tetris)


# Further combining data
df_5_acc <- merge(RV_Phi, RV_Lambda, all = TRUE) %>% 
  merge(RV_SigmaMarginal) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))
df_5_eff <- merge(Time, RAM, all = TRUE) %>% 
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = main_metrics,
                       labels = labels_main),
         method = factor(method, levels = method_all, labels = method_lab))

df_5_supp <- merge(FN_Phi_sc5, FN_SigmaPhi_sc5, all = TRUE) %>%
  merge(RV_SigmaPhi_sc5, all = TRUE) %>%
  merge(FN_Lambda_sc5, all = TRUE) %>%
  merge(FN_SigmaLambda_mean_sc5, all = TRUE) %>%
  merge(RV_SigmaLambda_mean_sc5, all = TRUE) %>%
  merge(FN_SigmaMarginal_mean_sc5, all = TRUE) %>%
  pivot_longer(cols = -c(seed, method), names_to = "metric", values_to = "value") %>% 
  mutate(metric=factor(metric, levels = supp_metrics,
                       labels = labels_supp),
         method = factor(method, levels = method_all, labels = method_lab))


acc_5 <-
  df_5_acc %>%
  ggplot(aes(x = method, y = value, colour = method)) +
  geom_boxplot() +
  geom_point(aes(x = method, y = Inf, colour = method), size = 0, alpha = 0) +
  facet_wrap(~metric, nrow = 1, labeller = label_parsed) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey95"),  # Keep major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines (optional)
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "grey95", color = NA), 
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),  # Rotate x-axis labels
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"),  # Smaller legend keys
        legend.text = element_text(size = 8),  # Smaller legend text
        plot.margin = margin(t = 2, r = 2, b = 0, l = 2),
        legend.title = element_text(size = 9)) +  
  lims(y = c(0, 1)) +
  labs(title = "",
       x = "",
       y = "Value",
       fill = "method",
       colour = "method") +
  scale_color_manual(values = method_colors,
                     drop = FALSE          )  # Use consistent color mapping


eff_5 <- df_5_eff %>%
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
    panel.border = element_rect(color = "black", fill = NA, size = 0.4),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
    legend.position = "none",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  ) +
  labs(title = "Scenario 5 (mimic genomics)",
       x = "",
       y = "Value",
       fill = "method"
  ) +
  scale_color_manual(values = method_colors) +
  guides(colour = guide_legend(override.aes = list(size = 0.5))) 
# +
  # scale_y_log10()



settings_sc5 <- data.frame(
  "x" = rep(0, 7),
  "y" = 6:0/6,
  "text" = c(
    TeX("Scenario 5 (mimic genomics)", output="character"),
    TeX("Based on Tetris", output="character"),
    TeX("S = 4", output="character"),
    TeX("P = 1060", output="character"),
    TeX("$N_s = (157, 195, 285, 117)$", output = "character"),
    TeX("K = 15, J_s=(2, 2, 2, 2)", output="character"),
    TeX("3 partial factors", output="character")
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

#----------------- Supplementary plots----------------------------
levels(df_5_supp$method)
supp_5 <- df_5_supp %>%
  mutate(method = factor(method)) %>% 
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
  labs(title = "Scenario 5",
       x = "",
       y = "Value",
       fill = "Method") +
  scale_color_manual(values = method_colors) +  # Use consistent color mapping
  guides(colour = guide_legend(override.aes = list(size = 0.5)))
  