# Load libraries
library(ggplot2)
library(reshape2)
library(latex2exp)
library(gridExtra)
library(tidyverse)
source("./functions/plot_metrics.R")


##################### Summaries (box plots) for simulations ####################
##################### Senerio 1: Dense ###################
sen1_dense <- readRDS("./RDS/sen1_dense_7.9.rds")
plot_metrics(sen1_dense)
##################### Senerio 1: Sparse ###################
sen1_sparse <- readRDS("./RDS/sen1_sparse_7.9.rds")
plot_metrics(sen1_sparse)




# Heatmap for cross-product of matrix Phi (true vs estimated)
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
