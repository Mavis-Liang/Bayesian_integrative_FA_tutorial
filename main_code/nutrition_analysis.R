load("./RDS/dataLAT_projale2.rda")

count_na_and_negatives <- function(df) {
  # Count NA values
  na_count <- sum(is.na(df))
  # Count negative values (ignoring NAs)
  negative_count <- sum(df < 0, na.rm = TRUE)
  
  # Print counts
  cat("Number of NAs:", na_count, "\n")
  cat("Number of negative values:", negative_count, "\n")
}

# Apply the function to each data frame in the list
invisible(lapply(X_s2, count_na_and_negatives))

process_study_data <- function(df) {
  # Remove rows where all values are NA
  cleaned_df <- df[!apply(df, 1, function(row) all(is.na(row))), , drop = FALSE]
  # Count remaining rows
  remaining_rows <- nrow(cleaned_df)
  
  # Print results for the study
  cat("Remaining rows:", remaining_rows, "\n")
  return(cleaned_df)
}

Y_list <- lapply(X_s2, process_study_data)

# Replace negative values with 0, then log(x+0.01) + 0.01
replace_negatives <- function(df) {
  # Replace negative values with 0
  df[df < 0] <- 0
  # Apply log transformation
  transformed_df <- log(df + 0.01)
  return(transformed_df)
}

Z_list <- lapply(Y_list, replace_negatives)
saveRDS(Z_list, "./RDS/nutrition_processed.rds")

# Running each methods with over-specified numbers of factors


# Further Post-processing
# Numbers of factors - StackFA, IndFA, BMSFA
fun_eigen <- function(Sig_mean) {
  val_eigen <- eigen(Sig_mean)$values
  prop_var <- val_eigen/sum(val_eigen)
  choose_K <- length(which(prop_var > 0.05))
  return(choose_K)
}
SigmaPhi_StackFA <- readRDS("./RDS/results_nutrition/Rnutrition_StackFA.rds")$SigmaPhi
K_StackFA <- fun_eigen(SigmaPhi_StackFA)
K_StackFA

SigmaLambda_IndFA <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA.rds")$SigmaLambda
Js_IndFA <- lapply(SigmaLambda_IndFA, fun_eigen)
Js_IndFA

SigmaPhi_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA.rds")$SigmaPhi
K_BMSFA <- fun_eigen(SigmaPhi_BMSFA)
SigmaLambda_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA.rds")$SigmaLambda
Js_BMSFA <- lapply(SigmaLambda_BMSFA, fun_eigen)
K_BMSFA
Js_BMSFA

# Numbers of factors - PFA
# columns have all loadings less than 10^-3
fun_neibour <- function(Phi, threshold = 1e-3) {
  return(
    sum(apply(Phi, 2, function(x) {
      # if this column has all loadings less than threshold
      sum(abs(x) <= threshold) < length(x)}
    ))
  )
}
Phi_PFA <- readRDS("./RDS/results_nutrition/Rnutrition_PFA.rds")$Phi
K_PFA <- fun_neibour(Phi_PFA)
K_PFA

# ------------------------- Visualizations -------------------------
source("./functions/plot_single_loading.R")
library(patchwork)
reorder_loading <- function(loading) {
  variance_explained <- colSums(loading^2)
  order_indices <- order(variance_explained, decreasing = TRUE)
  loading_reordered <- loading[, order_indices]
  return(loading_reordered)
}

# Make the majority sign positive of each column
switch_sign <- function(loading){
  loading <- apply(loading, 2, function(x) x * sign(sum(sign(x))))
  return(loading)
}

# StackFA
Phi_StackFA <- readRDS("./RDS/results_nutrition/Rnutrition_StackFA_2.rds")$Phi %>% 
  reorder_loading() %>% 
  switch_sign()
Phi_StackFA[, c(2, 3)] <- - Phi_StackFA[, c(2, 3)]  
P_stackFA <- plot_single(Phi_StackFA, y_label = colnames(X_s2[[1]]), show_colorbar = TRUE)+ ggtitle(TeX("$\\Phi$"))
heatStackFA <- P_stackFA + plot_annotation("(a) Stack FA")
ggplot2::ggsave("./Figs/nutrition/heatStackFA.pdf", width = 4, height = 5, units = "in", dpi = 300)

# IndFA
Lambda1_IndFA <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA_2.rds")$LambdaList[[1]] %>% reorder_loading() %>% switch_sign()
L1_IndFA <- plot_single(Lambda1_IndFA, y_label = colnames(X_s2[[1]])) + ggtitle(TeX("$\\Lambda_1$"))
Lambda2_IndFA <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA_2.rds")$LambdaList[[2]] %>% reorder_loading() %>% switch_sign()
Lambda2_IndFA[, 3] <- - Lambda2_IndFA[, 3]
L2_IndFA <- plot_single(Lambda2_IndFA)+ ggtitle(TeX("$\\Lambda_2$"))
Lambda3_IndFA <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA_2.rds")$LambdaList[[3]] %>% reorder_loading() %>% switch_sign()
L3_IndFA <- plot_single(Lambda3_IndFA)+ ggtitle(TeX("$\\Lambda_3$"))
Lambda4_IndFA <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA_2.rds")$LambdaList[[4]] %>% reorder_loading() %>% switch_sign()
Lambda4_IndFA[, 3] <- - Lambda4_IndFA[, 3]
L4_IndFA <- plot_single(Lambda4_IndFA)+ ggtitle(TeX("$\\Lambda_4$"))
Lambda5_IndFA <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA_2.rds")$LambdaList[[5]] %>% reorder_loading() %>% switch_sign()
L5_IndFA <- plot_single(Lambda5_IndFA)+ ggtitle(TeX("$\\Lambda_5$"))
Lambda6_IndFA <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA_2.rds")$LambdaList[[6]] %>% reorder_loading() %>% switch_sign()
Lambda6_IndFA[, 3] <- - Lambda6_IndFA[, 3]
L6_IndFA <- plot_single(Lambda6_IndFA, show_colorbar = TRUE)+ ggtitle(TeX("$\\Lambda_6$"))

layout_IndFA <- "
ABCDEF
"
heatIndFA <- L1_IndFA + L2_IndFA + L3_IndFA + L4_IndFA + L5_IndFA + L6_IndFA +
  plot_layout(design=layout_IndFA) + plot_annotation("(b) Ind FA")

ggplot2::ggsave("./Figs/nutrition/heatIndFA.pdf", width = 10, height = 5, units = "in", dpi = 300)

# PFA
# This is actually Phi for the first study
Phi_PFA <- readRDS("./RDS/results_nutrition/Rnutrition_PFA.rds")$Phi
Phi_PFA <- reorder_loading(Phi_PFA) %>% switch_sign()
P_PFA <- plot_single(Phi_PFA, 
                     y_label = colnames(X_s2[[1]]), 
                     show_colorbar = TRUE)+ 
  ggtitle(TeX("$\\Phi$"))
heatPFA <- P_PFA + plot_annotation("(c) PFA")
ggplot2::ggsave("./Figs/nutrition/heatPFA.pdf", width = 4.2, height = 5, units = "in", dpi = 300)
# MOM-SS
P_MOMSS <- readRDS("./RDS/results_nutrition/Rnutrition_MOMSS.rds")$Phi %>% reorder_loading() %>% switch_sign()
P_MOMSS <- plot_single(Phi_MOMSS, y_label = colnames(X_s2[[1]]), fill_limits = c(-8, 8), show_colorbar = TRUE) + ggtitle(TeX("$\\Phi$"))
heatMOMSS <- P_MOMSS + plot_annotation("(d) MOM-SS")
ggplot2::ggsave("./Figs/nutrition/heatMOMSS.pdf", width = 4.2, height = 5, units = "in", dpi = 300)


# SUFA
Phi_SUFA <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA.rds")$Phi %>% reorder_loading() %>% switch_sign()
P_SUFA <- plot_single(Phi_SUFA, y_label = colnames(X_s2[[1]])) + ggtitle(TeX("$\\Phi$"))

# No need to re-order since only one factor is estimated
Lambda1_SUFA <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA.rds")$LambdaList[[1]] %>% switch_sign()
L1_SUFA <- plot_single(Lambda1_SUFA) + ggtitle(TeX("$\\Lambda_1$"))
Lambda2_SUFA <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA.rds")$LambdaList[[2]] %>% switch_sign()
L2_SUFA <- plot_single(Lambda2_SUFA) + ggtitle(TeX("$\\Lambda_2$"))
Lambda3_SUFA <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA.rds")$LambdaList[[3]] %>% switch_sign()
L3_SUFA <- plot_single(Lambda3_SUFA) + ggtitle(TeX("$\\Lambda_3$"))
Lambda4_SUFA <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA.rds")$LambdaList[[4]] %>% switch_sign()
L4_SUFA <- plot_single(Lambda4_SUFA) + ggtitle(TeX("$\\Lambda_4$"))
Lambda5_SUFA <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA.rds")$LambdaList[[5]] %>% switch_sign()
L5_SUFA <- plot_single(Lambda5_SUFA) + ggtitle(TeX("$\\Lambda_5$"))
Lambda6_SUFA <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA.rds")$LambdaList[[6]] %>% switch_sign()
L6_SUFA <- plot_single(Lambda6_SUFA, show_colorbar = TRUE) + ggtitle(TeX("$\\Lambda_6$"))

layout_SUFA <- "AAAAAABCDEFG"
heatSUFA <- P_SUFA + L1_SUFA + L2_SUFA + L3_SUFA + L4_SUFA + L5_SUFA + L6_SUFA +
  plot_layout(design=layout_SUFA) + plot_annotation("(e) SUFA")
ggplot2::ggsave("./Figs/nutrition/heatSUFA.pdf", width = 6.8, height = 5, units = "in", dpi = 300)


# BMSFA
Phi_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")$Phi %>% 
  reorder_loading() %>% 
  switch_sign()
Phi_BMSFA[, c(2,4)] <- - Phi_BMSFA[, c(2,4)]  
P_BMSFA <- plot_single(Phi_BMSFA, y_label = colnames(X_s2[[1]]))+ ggtitle(TeX("$\\Phi$"))

Lambda1_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")$LambdaList[[1]] %>% reorder_loading() %>% switch_sign()
Lambda1_BMSFA[, 1] <- - Lambda1_BMSFA[, 1]
L1_BMSFA <- plot_single(Lambda1_BMSFA) + ggtitle(TeX("$\\Lambda_1$"))
Lambda2_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")$LambdaList[[2]] %>% reorder_loading() %>% switch_sign()
Lambda2_BMSFA[, 1] <- - Lambda2_BMSFA[, 1]
L2_BMSFA <- plot_single(Lambda2_BMSFA) + ggtitle(TeX("$\\Lambda_2$"))
Lambda3_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")$LambdaList[[3]] %>% reorder_loading() %>% switch_sign()
Lambda3_BMSFA[, 1] <- - Lambda3_BMSFA[, 1]
L3_BMSFA <- plot_single(Lambda3_BMSFA) + ggtitle(TeX("$\\Lambda_3$"))
Lambda4_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")$LambdaList[[4]] %>% reorder_loading() %>% switch_sign()
Lambda4_BMSFA[, 1] <- - Lambda4_BMSFA[, 1]
L4_BMSFA <- plot_single(Lambda4_BMSFA) + ggtitle(TeX("$\\Lambda_4$"))
Lambda5_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")$LambdaList[[5]] %>% reorder_loading() %>% switch_sign()
Lambda5_BMSFA[, 1] <- - Lambda5_BMSFA[, 1]
L5_BMSFA <- plot_single(Lambda5_BMSFA) + ggtitle(TeX("$\\Lambda_5$"))
Lambda6_BMSFA <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_2.rds")$LambdaList[[6]] %>% reorder_loading() %>% switch_sign()
Lambda6_BMSFA[, 1] <- - Lambda6_BMSFA[, 1]
L6_BMSFA <- plot_single(Lambda6_BMSFA, show_colorbar = TRUE) + ggtitle(TeX("$\\Lambda_6$"))

layout_BMSFA <- "AAABCDEFG"
heatBMSFA <- P_BMSFA + L1_BMSFA + L2_BMSFA + L3_BMSFA + L4_BMSFA + L5_BMSFA + L6_BMSFA +
  plot_layout(design=layout_BMSFA) + plot_annotation("(f) BMSFA")
ggplot2::ggsave("./Figs/nutrition/heatBMSFA.pdf", width = 7, height = 5, units = "in", dpi = 300)

# Tetris
Phi_Tetris <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris.rds")$Phi %>% reorder_loading() %>% switch_sign()
P_Tetris <- plot_single(Phi_Tetris, y_label = colnames(X_s2[[1]])) + ggtitle(TeX("$\\Phi$"))

Lambda1_Tetris <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris.rds")$LambdaList[[1]] # J_1=0
L1_Tetris <- plot_single(matrix(data = 0, ncol=1, nrow = 42)) + ggtitle(TeX("$\\Lambda_1$"))
Lambda2_Tetris <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris.rds")$LambdaList[[2]] %>% switch_sign()# J_2=1 
L2_Tetris <- plot_single(Lambda2_Tetris) + ggtitle(TeX("$\\Lambda_2$"))
Lambda3_Tetris <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris.rds")$LambdaList[[3]] # J_3=0
L3_Tetris <- plot_single(matrix(data = 0, ncol=1, nrow = 42)) + ggtitle(TeX("$\\Lambda_3$"))
Lambda4_Tetris <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris.rds")$LambdaList[[4]]%>% switch_sign() # J_4=1
L4_Tetris <- plot_single(Lambda4_Tetris) + ggtitle(TeX("$\\Lambda_4$"))
Lambda5_Tetris <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris.rds")$LambdaList[[5]] %>% switch_sign()# J_5=1
L5_Tetris <- plot_single(Lambda5_Tetris) + ggtitle(TeX("$\\Lambda_5$"))
Lambda6_Tetris <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris.rds")$LambdaList[[6]]%>% switch_sign() # J_6=1
L6_Tetris <- plot_single(Lambda6_Tetris, show_colorbar = TRUE) + ggtitle(TeX("$\\Lambda_6$"))

layout_Tetris <- "AAAAABCDEFG"
heatTetris <- P_Tetris + L1_Tetris + L2_Tetris + L3_Tetris + L4_Tetris + L5_Tetris + L6_Tetris +
  plot_layout(design=layout_Tetris) + plot_annotation("(g) Tetris")

ggplot2::ggsave("./Figs/nutrition/heatTetris.pdf", width = 7, height = 5, units = "in", dpi = 300)

# -------------------------- Predictions and MSE --------------------------
train_ratio <- 0.7
train_list <- list()
test_list <- list()
set.seed(6)
for (s in seq_along(Z_list)) {
  N_s <- nrow(Z_list[[s]])  # Number of rows in the study
  train_indices <- sample(1:N_s, size = floor(train_ratio * N_s), replace = FALSE)
  
  train_list[[s]] <- Z_list[[s]][train_indices, ]
  test_list[[s]] <- Z_list[[s]][-train_indices, ]
}
saveRDS(train_list, "./RDS/nutrition_train.rds")
saveRDS(test_list, "./RDS/nutrition_test.rds")

# Fit the models on the training data again

# Test data has to be centered
test_list <- lapply(test_list, as.matrix)
test_list_scaled <- lapply(
  test_list, function(x) scale(x, center = TRUE, scale = FALSE)
)

# Stack FA
fit_StackFA_train <- readRDS("./RDS/results_nutrition/Rnutrition_StackFA_train.rds")
Phi <- fit_StackFA_train$Phi
Psi <- fit_StackFA_train$Psi

mse_stackFA <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    scores <- test_list_scaled[[s]] %*% solve(Psi) %*% Phi %*% mnormt::pd.solve(signif(t(Phi) %*% solve(Psi) %*% Phi))
    norm(test_list_scaled[[s]] - scores %*% t(Phi), "F")^2
  })
)

# Ind FA
fit_IndFA_train <- readRDS("./RDS/results_nutrition/Rnutrition_IndFA_train.rds")
LambdaList <- fit_IndFA_train$LambdaList
Psi <- fit_IndFA_train$Psi
mse_IndFA <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    scores <- test_list_scaled[[s]] %*% solve(Psi[[s]]) %*% LambdaList[[s]] %*% mnormt::pd.solve(signif(t(LambdaList[[s]]) %*% solve(Psi[[s]]) %*% LambdaList[[s]]))
    
    norm(test_list_scaled[[s]] - scores %*% t(LambdaList[[s]]), "F")^2
  })
)

# PFA
fit_PFA_train <- readRDS("./RDS/results_nutrition/Rnutrition_PFA_train.rds")
Phi <- fit_PFA_train$Phi
Psi <- fit_PFA_train$Psi
Q_list <- fit_PFA_train$Q
mse_PFA <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    scores <- test_list_scaled[[s]] %*% t(Q_list[[s]]) %*% Phi %*% mnormt::pd.solve(signif(t(Phi) %*% solve(Psi) %*% Phi))
    Y_est <- scores %*% t(Phi) %*% t(solve(Q_list[[s]]))
    norm(test_list_scaled[[s]] - Y_est, "F")^2
  })
)

# BMSFA
fit_BMSFA_train <- readRDS("./RDS/results_nutrition/Rnutrition_BMSFA_train.rds")
Phi <- fit_BMSFA_train$Phi
LambdaList <- fit_BMSFA_train$LambdaList
Psi <- fit_BMSFA_train$PsiList %>%  lapply(as.vector) %>% lapply(diag)
mse_BMSFA <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    Omega <- cbind(Phi, LambdaList[[s]])
    scores <- test_list_scaled[[s]] %*% solve(Psi[[s]]) %*% Omega %*% mnormt::pd.solve(signif(t(Omega) %*% solve(Psi[[s]]) %*% Omega))
    
    norm(test_list_scaled[[s]] - scores %*% t(Omega), "F")^2
  })
)

# SUFA
fit_SUFA_train <- readRDS("./RDS/results_nutrition/Rnutrition_SUFA_train.rds")
Phi <- fit_SUFA_train$Phi
LambdaList <- fit_SUFA_train$LambdaList
Psi <- fit_SUFA_train$Psi
mse_SUFA <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    Omega <- cbind(Phi, LambdaList[[s]])
    scores <- test_list_scaled[[s]] %*% solve(Psi) %*% Omega %*% mnormt::pd.solve(signif(t(Omega) %*% solve(Psi) %*% Omega))
    
    norm(test_list_scaled[[s]] - scores %*% t(Omega), "F")^2
  })
)

# MOM-SS
fit_MOMSS_train <- readRDS("./RDS/results_nutrition/Rnutrition_MOMSS_train.rds")
Phi <- fit_MOMSS_train$Phi
Psi <- fit_MOMSS_train$Psi %>% lapply(diag)
alpha <- lapply(1:6, function(s) {
  fit_MOMSS_train$alpha[,s]
})

# Results here not so correct
mse_MOMSS <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    scores <- t(apply(test_list[[s]], 1, function(row) {row - alpha[[s]]})) %*% solve(Psi[[s]]) %*% Phi %*% mnormt::pd.solve(signif(t(Phi) %*% solve(Psi[[s]]) %*% Phi))
    Y_est <- t(apply(scores %*% t(Phi), 1, function(row) {row + alpha[[s]]}))
    norm(test_list[[s]] - Y_est, "F")^2
  })
)

mse_MOMSS <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    scores <- test_list_scaled[[s]] %*% solve(Psi[[s]]) %*% Phi %*% mnormt::pd.solve(signif(t(Phi) %*% solve(Psi[[s]]) %*% Phi))
    norm(test_list_scaled[[s]] - scores %*% t(Phi), "F")^2
  })
)


# Tetris
fit_Tetris_train <- readRDS("./RDS/results_nutrition/Rnutrition_Tetris_train.rds")
Phi <- fit_Tetris_train$Phi
LambdaList <- fit_Tetris_train$LambdaList
Marginal <- fit_Tetris_train$SigmaMarginal
SigmaPhi <- fit_Tetris_train$SigmaPhi
SigmaLambdaList <- fit_Tetris_train$SigmaLambdaList
Psi <- lapply(1:6, function(s) {
  Marginal[[s]] - SigmaPhi - SigmaLambdaList[[s]]
})
T_s_list <- fit_Tetris_train$T_s
mse_Tetris <- 1/(42 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:6, function(s){
    Omega <- cbind(Phi, LambdaList[[s]])
    scores <- test_list_scaled[[s]] %*% solve(Psi[[s]]) %*% Omega %*% mnormt::pd.solve(signif(t(Omega) %*% solve(Psi[[s]]) %*% Omega))
    
    norm(test_list_scaled[[s]] - scores %*% t(Omega), "F")^2
  })
)
  
  
