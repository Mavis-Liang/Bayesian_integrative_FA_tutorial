# This file generate the Table 2 (summarizing model's performance in estimating K and J_s)
# Based on the model results in RDS folder (names with XX_mis)

# num of factors for Stack FA, Ind FA, and BMSFA
fun_eigen <- function(Sig_mean) {
  val_eigen <- eigen(Sig_mean)$values
  prop_var <- val_eigen/sum(val_eigen)
  choose_K <- length(which(prop_var > 0.05))
  return(choose_K)
}

# numbers of factors for MOM-SS (MOMSS has already post-processed for the number of factors)
count_col <- function(Phi){
  return(ncol(Phi))
}

# numbers of factors for SUFA (use same function as MOM-SS)
# SUFA

# numbers of factors for Tetris (use same function as MOM-SS)

# Counting num of f for all the models
count_f <- function(point_est_list, models){
  df <- data.frame(matrix(ncol = length(models), nrow = 1))
  colnames(df) <- models
  for (model in models){
    if (model == "BMSFA" |  model == "stackFA"){
      num_f <- fun_eigen(point_est_list[[model]]$SigmaPhi) }
      else if (model == "MOMSS" | model == "SUFA" | model == "Tetris" | model== "PFA"){
      num_f <- count_col(point_est_list[[model]]$Phi)
    }
    
    df[,model] <- num_f
  }
  return(df)
}

# scenario 1
# For seed=1-50
models_sc1 <- c("stackFA","PFA", "MOMSS", "SUFA", "BMSFA", "Tetris")
df_sc1 <- NULL
for (i in 1:50){
  sc1_mis <- readRDS(paste0("./RDS/sc1_mis/sc1_mis_", i, ".rds"))$point_est
  sc1_mis_Tetris <- readRDS(paste0("./RDS/sc1/Tetris/sc1_Tetris_", i, ".rds"))$point_est
  all <- c(sc1_mis, sc1_mis_Tetris)
  curr <- count_f(all, models= models_sc1)
  df_sc1 <- rbind(df_sc1, curr)
}

# scenario 3
# For seed=1-50
models_sc3 <- c("stackFA","PFA", "MOMSS", "SUFA", "BMSFA", "Tetris")
df_sc3 <- NULL
for (i in 1:50){
  sc3_mis <- readRDS(paste0("./RDS/sc3_mis/sc3_mis_", i, ".rds"))$point_est
  sc3_mis_Tetris <- readRDS(paste0("./RDS/sc3/Tetris/sc3_Tetris_", i, ".rds"))$point_est
  all <- c(sc3_mis, sc3_mis_Tetris)
  curr <- count_f(all, models = models_sc3)
  df_sc3 <- rbind(df_sc3, curr)
}

# scenario 4
# For seed=1-50
sc4_mis_1 <- readRDS("./RDS/sc4_mis/sc4_mis_1.rds")
models_sc4 <- c("stackFA","PFA", "MOMSS", "SUFA", "BMSFA")
df_sc4 <- NULL
for (i in 1:50){
  sc4_mis <- readRDS(paste0("./RDS/sc4_mis/sc4_mis_", i, ".rds"))$point_est
  curr <- count_f(sc4_mis, models = models_sc4)
  df_sc4 <- rbind(df_sc4, curr)
}

# scenario 5
MOMSSBMSFA <- readRDS("./RDS/sc5_mis/sc5_MOMSS_BMSFA_mis_1.rds")
models_sc5_MOMSSBMSFA <- names(MOMSSBMSFA$point_est)
df_sc5_MOMSSBMSFA <- NULL
for (i in 1:50){
  sc5_mis <- readRDS(paste0("./RDS/sc5_mis/sc5_MOMSS_BMSFA_mis_", i, ".rds"))$point_est
  curr <- count_f(sc5_mis, models = models_sc5_MOMSSBMSFA)
  df_sc5_MOMSSBMSFA <- rbind(df_sc5_MOMSSBMSFA, curr)
}

df_sc5_StackFA <- NULL
for (i in 1:50){
  sc5_mis <- readRDS(paste0("./RDS/sc5_mis/sc5_StackFA_mis_", i, ".rds"))$point_est
  curr <- count_f(sc5_mis, models = "stackFA")
  df_sc5_StackFA <- rbind(df_sc5_StackFA, curr)
}

df_sc5_SUFA <- NULL
for (i in 1:30){
  sc5_mis <- readRDS(paste0("./RDS/sc5_mis/sc5_SUFA_mis_", i, ".rds"))$point_est
  curr <- count_f(sc5_mis, models = "SUFA")
  df_sc5_SUFA <- rbind(df_sc5_SUFA, curr)
}



# scenario 6
models_sc6 <- c("stackFA","PFA", "MOMSS", "SUFA", "BMSFA", "Tetris")
df_sc6 <- NULL
for (i in c(1:27,29:50)){
  sc6_mis <- readRDS(paste0("./RDS/sc6_mis/sc6_mis_", i, ".rds"))$point_est
  sc6_mis_Tetris <- readRDS(paste0("./RDS/sc6/Tetris/sc6_Tetris_", i, ".rds"))$point_est
  all <- c(sc6_mis, sc6_mis_Tetris)
  curr <- count_f(all, models = models_sc6)
  df_sc6 <- rbind(df_sc6, curr)
}


# ------------------------------------------------------------------------------
# Estimations of J_s
count_J_s <- function(point_est_list, models){
  S <- length(point_est_list[[models[1]]]$SigmaMarginal)
  df_Js <- data.frame(matrix(ncol = S+1, 
                             nrow = length(models))) %>% 
    # add a column for model names
    rename(model_name = paste0("X", S+1)) %>% 
    mutate(model_name = models)
  for (model in models){
    if (model == "BMSFA" |  model == "IndFA"){
      num_J_s <- rbind(
        lapply(point_est_list[[model]]$SigmaLambdaList, function(x) fun_eigen(x))
      )
    } else if (model == "SUFA" | model == "Tetris"){
      num_J_s <- rbind(
        lapply(point_est_list[[model]]$LambdaList, function(x) count_col(x))
      )
    }
    
    df_Js[df_Js$model_name == model,1:S] <- num_J_s
  }
  return(df_Js)
}

# scenario 1
models_Js_sc1 <- c("BMSFA", "IndFA", "SUFA", "Tetris")
df_sc1_J_s <- NULL
for (i in 1:50){
  sc1_mis <- readRDS(paste0("./RDS/sc1_mis/sc1_mis_", i, ".rds"))$point_est
  sc1_mis_Tetris <- readRDS(paste0("./RDS/sc1/Tetris/sc1_Tetris_", i, ".rds"))$point_est
  all <- c(sc1_mis, sc1_mis_Tetris)
  curr <- count_J_s(all, models_Js_sc1)
  df_sc1_J_s <- rbind(df_sc1_J_s, curr)
}


# scenario 3
models_Js_sc3 <- c("BMSFA", "IndFA", "SUFA", "Tetris")
df_sc3_J_s <- NULL
for (i in 1:50){
  sc3_mis <- readRDS(paste0("./RDS/sc3_mis/sc3_mis_", i, ".rds"))$point_est
  sc3_mis_Tetris <- readRDS(paste0("./RDS/sc3/Tetris/sc3_Tetris_", i, ".rds"))$point_est
  all <- c(sc3_mis, sc3_mis_Tetris)
  curr <- count_J_s(all, models_Js_sc3)
  df_sc3_J_s <- rbind(df_sc3_J_s, curr)
}

# scenario 4
models_Js_sc4 <- c("BMSFA", "IndFA", "SUFA")
df_sc4_J_s <- NULL
for (i in 1:50){
  sc4_mis <- readRDS(paste0("./RDS/sc4_mis/sc4_mis_", i, ".rds"))$point_est
  curr <- count_J_s(sc4_mis, models_Js_sc4)
  df_sc4_J_s <- rbind(df_sc4_J_s, curr)
}

# scenario 5
df_sc5_Js_BMSFA <- NULL
for (i in 1:50){
  sc5_mis <- readRDS(paste0("./RDS/sc5_mis/sc5_MOMSS_BMSFA_mis_", i, ".rds"))$point_est
  curr <- count_J_s(sc5_mis, "BMSFA")
  df_sc5_Js_BMSFA <- rbind(df_sc5_Js_BMSFA, curr)
}

df_sc5_Js_IndFA <- NULL
for (i in 1:50){
  sc5_mis <- readRDS(paste0("./RDS/sc5_mis/sc5_IndFA_mis", i, ".rds"))$point_est
  curr <- count_J_s(sc5_mis, "IndFA")
  df_sc5_Js_IndFA <- rbind(df_sc5_Js_IndFA, curr)
}

df_sc5_Js_SUFA <- NULL
for (i in 1:30){
  sc5_mis <- readRDS(paste0("./RDS/sc5_mis/sc5_SUFA_mis_", i, ".rds"))$point_est
  curr <- count_J_s(sc5_mis, "SUFA")
  df_sc5_Js_SUFA <- rbind(df_sc5_Js_SUFA, curr)
}

# scenario 6
models_Js_sc6 <- c("BMSFA", "IndFA", "SUFA", "Tetris")
df_sc6_J_s <- NULL
for (i in 1:50){
  sc6_mis <- readRDS(paste0("./RDS/sc6_mis/sc6_mis_", i, ".rds"))$point_est
  sc6_mis_Tetris <- readRDS(paste0("./RDS/sc6/Tetris/sc6_Tetris_", i, ".rds"))$point_est
  all <- c(sc6_mis, sc6_mis_Tetris)
  curr <- count_J_s(all, models_Js_sc6)
  df_sc6_J_s <- rbind(df_sc6_J_s, curr)
}


######## ------ calculate mean and sd and organize the table ------ ########
calculate_mean_sd <- function(df) {
  df %>%
    summarise(across(everything(), 
                     list(mean = ~ mean(.), sd = ~ sd(.)), 
                     .names = "{.col}_{.fn}")) %>%
    pivot_longer(cols = everything(), 
                 names_to = c("Model", ".value"), 
                 names_sep = "_") %>%
    mutate(Mean_SD = sprintf("%.2f(%.2f)", mean, sd)) %>%
    select(Model, Mean_SD)
}
calculate_mean_sd_J <- function(df_js) {
  df_js %>%
    pivot_longer(cols = starts_with("X"),
                 names_to = "Metric", values_to = "Value") %>%
    group_by(model_name, Metric) %>%
    summarise(
      Mean = mean(Value),
      SD = sd(Value),
      .groups = "drop"
    ) %>%
    mutate(Mean_SD = sprintf("%.2f(%.2f)", Mean, SD)) %>%
    group_by(model_name) %>%
    summarise(
      J_s = paste(Mean_SD, collapse = ", "),
      .groups = "drop"
    )
}


df_sc1_mean_sd <- calculate_mean_sd(df_sc1)
df_sc1_Js_mean_sd <- calculate_mean_sd_J(df_sc1_J_s) %>% 
  mutate(Scenario = "Scenario 1")


df_sc3_mean_sd <- calculate_mean_sd(df_sc3)
df_sc3_J_mean_sd <- calculate_mean_sd_J(df_sc3_J_s) %>% 
  mutate(Scenario = "Scenario 3")


df_sc4_mean_sd <- calculate_mean_sd(df_sc4)
df_sc4_J_mean_sd <- calculate_mean_sd_J(df_sc4_J_s) %>% 
  mutate(Scenario = "Scenario 4")


df_sc5_SUFA_mean_sd <- calculate_mean_sd(df_sc5_SUFA)
df_sc5_StackFA_mean_sd <- calculate_mean_sd(df_sc5_StackFA)
df_sc5_MOMSSBMSFA_mean_sd <- calculate_mean_sd(df_sc5_MOMSSBMSFA)
df_sc5_combined <- bind_rows(df_sc5_MOMSSBMSFA_mean_sd, df_sc5_SUFA_mean_sd, df_sc5_StackFA_mean_sd)

df_sc5_Js_SUFA_mean_sd <- calculate_mean_sd_J(df_sc5_Js_SUFA)
df_sc5_Js_BMSFA_mean_sd <- calculate_mean_sd_J(df_sc5_Js_BMSFA)
df_sc5_Js_IndFA_mean_sd <- calculate_mean_sd_J(df_sc5_Js_IndFA)

df_sc5_Js_combined <- bind_rows(df_sc5_Js_BMSFA_mean_sd, df_sc5_Js_SUFA_mean_sd, df_sc5_Js_IndFA_mean_sd)%>%
  mutate(Scenario = "Scenario 5")

df_sc6_mean_sd <- calculate_mean_sd(df_sc6)
df_sc6_J_mean_sd <- calculate_mean_sd_J(df_sc6_J_s) %>% 
  mutate(Scenario = "Scenario 6")


# ------------------------------------------------------------------------------
# Sample R code for combining K and J_s tables into a unified format
library(dplyr)
library(tidyr)
library(xtable)

# Ensure Scenario column exists in both data frames before joining
prepare_data <- function(df, scenario_name, is_js = FALSE) {
  if (!"Scenario" %in% colnames(df)) {
    df <- df %>% mutate(Scenario = scenario_name)
  }
  if (is_js && !"model_name" %in% colnames(df)) {
    stop("The J_s data frame must have a 'model_name' column.")
  }
  return(df)
}

# Combine K and J_s tables for a scenario using a full join
combine_k_js <- function(df_k, df_js, scenario_name) {
  df_k <- prepare_data(df_k, scenario_name)
  df_js <- prepare_data(df_js, scenario_name, is_js = TRUE)
  
  # Adjust column names if necessary
  if (!"Mean_SD" %in% colnames(df_k)) {
    df_k <- df_k %>% rename(Mean_SD = K_Value)
  }
  
  df_combined <- full_join(
    df_k %>% select(Scenario, Model, Mean_SD),
    df_js %>% select(Scenario, model_name, J_s) %>% rename(Model = model_name),
    by = c("Scenario", "Model")
  )
  return(df_combined)
}


# Repeat for other scenarios as needed
# Example for combining multiple scenarios
final_combined <- bind_rows(
  combine_k_js(df_sc1_mean_sd, df_sc1_Js_mean_sd, "Scenario 1"),  # Rename Scenario 1 to Scenario 2
  combine_k_js(df_sc3_mean_sd, df_sc3_J_mean_sd, "Scenario 3"),  # Rename Scenario 3 to Scenario 1
  combine_k_js(df_sc6_mean_sd, df_sc6_J_mean_sd, "Scenario 6"),  # Rename Scenario 6 to Scenario 3
  combine_k_js(df_sc4_mean_sd, df_sc4_J_mean_sd, "Scenario 4"),
  combine_k_js(df_sc5_combined, df_sc5_Js_combined, "Scenario 5")
) %>%
  mutate(Scenario = case_when(
    Scenario == "Scenario 1" ~ "Scenario 2",
    Scenario == "Scenario 3" ~ "Scenario 1",
    Scenario == "Scenario 6" ~ "Scenario 3",
    TRUE ~ Scenario
  )) %>%
  mutate(Model = factor(Model, levels = c("stackFA", "IndFA", "PFA", "MOMSS", "SUFA", "BMSFA", "Tetris"))) %>%
  arrange(Scenario, Model) %>%
  mutate(
    Model = recode(Model,
                   "stackFA" = "Stack FA",
                   "IndFA" = "Ind FA",
                   "MOMSS" = "MOM-SS"
    )
  )%>% 
  rename("Estimated K" = Mean_SD, "Estimated J_s" = J_s)%>%
  mutate(
    `True K` = case_when(
      Scenario == "Scenario 2" ~ 4,
      Scenario == "Scenario 1" ~ 4,
      Scenario == "Scenario 3" ~ 4,
      Scenario == "Scenario 4" ~ 4,
      Scenario == "Scenario 5" ~ 15
    ),
    `True J_s` = case_when(
      Scenario == "Scenario 2" ~ "0, 0, 0, 0",
      Scenario == "Scenario 1" ~ "0, 0, 0, 0",
      Scenario == "Scenario 3" ~ "1, 1, 1, 1",
      Scenario == "Scenario 4" ~ "1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1",
      Scenario == "Scenario 5" ~ "2, 2, 2, 2"
    )
  )


