# Figures
# Heatmap for cross-product of matrix Phi (true vs estimated)

# Load libraries
library(ggplot2)
library(reshape2)
library(latex2exp)
library(gridExtra)
library(tidyverse)


################# Senerio 1: Dense ###################

# Load true data
data_sen1_dense <- readRDS("./RDS/data_sen1_dense.rds")
Phi_sen1_dense <- data_sen1_dense$Phi
covPhi_sen1_dense <- tcrossprod(Phi_sen1_dense)

# Load estimated data
data_sen1_dense_MOMSS <- readRDS("./RDS/result_MOMSS_sen1_dense.rds")
Phi_sen1_dense_MOMSS <- data_sen1_dense_MOMSS$M## Un-scaled, thus wrong
covPhi_sen1_dense_MOMSS <- tcrossprod(Phi_sen1_dense_MOMSS)
data_sen1_dense_BMSFA <- readRDS("./RDS/result_BMSFA_sen1_dense.rds")
covPhi_sen1_dense_BMSFA <- data_sen1_dense_BMSFA$SigmaPhi

# Heatmap
p1 <- ggplot(data = melt(covPhi_sen1_dense), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  ggtitle(TeX("True $\\Phi\\Phi^T$ (dense)")) +
  labs(x="", y="") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))

p2 <- ggplot(data = melt(covPhi_sen1_dense_MOMSS), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)  +
  ggtitle(TeX("Estimated $\\Phi\\Phi^T$ (dense) by MOM-SS")) +
  labs(x="", y="") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))

p3 <- ggplot(data = melt(covPhi_sen1_dense_BMSFA), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)  +
  ggtitle(TeX("Estimated $\\Phi\\Phi^T$ (dense) by BMSFA")) +
  labs(x="", y="") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))

grid.arrange(p1, p2, p3, ncol=3)
# ggsave("./figs/sen1_dense_heatmap.png", width = 12, height = 3)


################# Senerio 2: Sparse ###################

# Load true data
data_sen1_sparse <- readRDS("./RDS/data_sen1_sparse.rds")
Phi_sen1_sparse <- data_sen1_sparse$Phi
covPhi_sen1_sparse <- tcrossprod(Phi_sen1_sparse)

# Load estimated data
data_sen1_sparse_MOMSS <- readRDS("./RDS/result_MOMSS_sen1_sparse.rds")
Phi_sen1_sparse_MOMSS <- scale(data_sen1_sparse_MOMSS$M)
covPhi_sen1_sparse_MOMSS <- tcrossprod(Phi_sen1_sparse_MOMSS)
data_sen1_sparse_BMSFA <- readRDS("./RDS/result_BMSFA_sen1_sparse.rds")
covPhi_sen1_sparse_BMSFA <- data_sen1_sparse_BMSFA$SigmaPhi

# Heatmap
p1 <- ggplot(data = melt(covPhi_sen1_sparse), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  ggtitle(TeX("True $\\Phi\\Phi^T$ (sparse)")) +
  labs(x="", y="") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))

p2 <- ggplot(data = melt(covPhi_sen1_sparse_MOMSS), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)  +
  ggtitle(TeX("Estimated $\\Phi\\Phi^T$ (sparse) by MOM-SS")) +
  labs(x="", y="") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))

p3 <- ggplot(data = melt(covPhi_sen1_sparse_BMSFA), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)  +
  ggtitle(TeX("Estimated $\\Phi\\Phi^T$ (sparse) by BMSFA")) +
  labs(x="", y="") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))

grid.arrange(p1, p2, p3, ncol=3)


##################### Summaries (box plots) for simulations ####################
##################### Senerio 1: Dense ###################
# Load data
sen1_dense <- readRDS("./RDS/sen1_dense.rds")
# boxplot - runtime
p1 <- sen1_dense %>% as.data.frame() %>% 
  select(run_time_MOMSS, run_time_BMSFA) %>% 
  rename_with(~ gsub("run_time_", "", .x), everything()) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "value") %>%
  unnest(value) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  labs(title = "Runtime comparison (dense)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))+
  labs(x = "Method", y = "Runtime (seconds)")

# boxplot - peak_RAM
p2 <- sen1_dense %>% as.data.frame() %>% 
  select(peak_RAM_MOMSS, peak_RAM_BMSFA) %>% 
  rename_with(~ gsub("peak_RAM_", "", .x), everything()) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "value") %>%
  unnest(value) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  labs(title = "Peak RAM comparison (dense)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))+
  labs(x = "Method", y = "Peak RAM (MB)")

# boxplot - RV of Phi
p3 <- sen1_dense %>% as.data.frame() %>% 
  dplyr::select(RV_MOMSS_Phi, RV_BMSFA_Phi) %>% 
  rename_with(~ gsub("RV_(.*)_Phi", "\\1", .x), everything()) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "value") %>%
  unnest(value) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  ggtitle(TeX("RV of $\\Phi\\Phi^T$ (dense)\\ comparison")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))+
  labs(x = "Method", y = "RV")

grid.arrange(p1, p2, p3, ncol=3)


##################### Senerio 2: Sparse ###################
# Load data
sen1_sparse <- readRDS("./RDS/sen1_sparse.rds")
# boxplot - runtime
p1 <- sen1_sparse %>% as.data.frame() %>% 
  select(run_time_MOMSS, run_time_BMSFA) %>% 
  rename_with(~ gsub("run_time_", "", .x), everything()) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "value") %>%
  unnest(value) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  labs(title = "Runtime comparison (sparse)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))+
  labs(x = "Method", y = "Runtime (seconds)")

# boxplot - peak_RAM
p2 <- sen1_sparse %>% as.data.frame() %>% 
  select(peak_RAM_MOMSS, peak_RAM_BMSFA) %>% 
  rename_with(~ gsub("peak_RAM_", "", .x), everything()) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "value") %>%
  unnest(value) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  labs(title = "Peak RAM comparison (sparse)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))+
  labs(x = "Method", y = "Peak RAM (MB)")

# boxplot - RV of Phi
p3 <- sen1_sparse %>% as.data.frame() %>% 
  select(RV_MOMSS_Phi, RV_BMSFA_Phi) %>% 
  rename_with(~ gsub("RV_(.*)_Phi", "\\1", .x), everything()) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "value") %>%
  unnest(value) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  ggtitle(TeX("RV of $\\Phi\\Phi^T$ (sparse) comparison")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=10))+
  labs(x = "Method", y = "RV")

grid.arrange(p1, p2, p3, ncol=3)


