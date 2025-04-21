library(curatedOvarianData)
library(tidyverse)
# Load datasets
data(package="curatedOvarianData")
data(GSE13876_eset)
data(GSE26712_eset)
data(GSE9891_eset)
data(PMID17290060_eset)

# find intersection of genes
inter_genes <- Reduce(intersect, list(featureNames(GSE13876_eset),
                                      featureNames(GSE26712_eset), 
                                      featureNames(GSE9891_eset),
                                      featureNames(PMID17290060_eset)))

GSE13876_eset <- GSE13876_eset[inter_genes,]
GSE26712_eset <- GSE26712_eset[inter_genes,]
GSE9891_eset <- GSE9891_eset[inter_genes,]
PMID17290060_eset <- PMID17290060_eset[inter_genes,]

study1 <- t(exprs(GSE13876_eset))
study2 <- t(exprs(GSE26712_eset))
study3 <- t(exprs(GSE9891_eset))
study4 <- t(exprs(PMID17290060_eset))


# calculate Coefficient of Variation of each gene
cv1 <- apply(study1, 2, sd) / apply(study1, 2, mean)
cv2 <- apply(study2, 2, sd) / apply(study2, 2, mean)
cv3 <- apply(study3, 2, sd) / apply(study3, 2, mean)
cv4 <- apply(study4, 2, sd) / apply(study4, 2, mean)
cv_matrix <- rbind(cv1, cv2, cv3, cv4)

# Find genes with CV >= threshold in at least one study
threshold <- 0.16
genes_to_keep <- apply(cv_matrix, 2, function(cv) any(cv >= threshold))
sum(genes_to_keep)
# Filtered
study1 <- GSE13876_eset[genes_to_keep,]
study2 <- GSE26712_eset[genes_to_keep,]
study3 <- GSE9891_eset[genes_to_keep,]
study4 <- PMID17290060_eset[genes_to_keep,]

df1 <- study1 %>% exprs() %>% t() %>% 
  log() %>% as.data.frame()
df2 <- study2 %>% exprs() %>% t() %>% 
  log() %>% as.data.frame()
df3 <- study3 %>% exprs() %>% t() %>% 
  log()  %>% as.data.frame()
df4 <- study4 %>% exprs() %>% t() %>% 
  log()  %>% as.data.frame()
list_gene <- list(df1, df2, df3, df4)
saveRDS(list(df1, df2, df3, df4), "CuratedOvarian_processed.rds")
sapply(list_gene, dim)

# Train set
train_ratio <- 0.7
train_list <- list()
test_list <- list()
Z_list <- list(df1, df2, df3, df4)
set.seed(6)
for (s in seq_along(Z_list)) {
  N_s <- nrow(Z_list[[s]])  # Number of rows in the study
  train_indices <- sample(1:N_s, size = floor(train_ratio * N_s), replace = FALSE)
  
  train_list[[s]] <- Z_list[[s]][train_indices, ]
  test_list[[s]] <- Z_list[[s]][-train_indices, ]
}
saveRDS(train_list, "./RDS/curated_train.rds")
saveRDS(test_list, "./RDS/curated_test.rds")

# Further Post-processing
# Numbers of factors - StackFA, IndFA, BMSFA
fun_eigen <- function(Sig_mean) {
  val_eigen <- eigen(Sig_mean)$values
  prop_var <- val_eigen/sum(val_eigen)
  choose_K <- length(which(prop_var > 0.05))
  return(choose_K)
}
SigmaPhi_StackFA <- readRDS("./RDS/results_curated/RCuratedOvarian_stackFA_scaled.rds")$SigmaPhi
K_StackFA <- fun_eigen(SigmaPhi_StackFA)
K_StackFA

SigmaLambda_IndFA <- readRDS("./RDS/results_curated/RCuratedOvarian_IndFA_scaled.rds")$SigmaLambda
Js_IndFA <- lapply(SigmaLambda_IndFA, fun_eigen)
Js_IndFA

SigmaPhi_BMSFA <- readRDS("./RDS/results_curated/RCuratedOvarian_BMSFA_scaled.rds")$SigmaPhi
K_BMSFA <- fun_eigen(SigmaPhi_BMSFA)
SigmaLambda_BMSFA <- readRDS("./RDS/results_curated/RCuratedOvarian_BMSFA_scaled.rds")$SigmaLambda
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
Phi_PFA <- readRDS("./RDS/results_curated/RCuratedOvarian_PFA.rds")$Phi
K_PFA <- fun_neibour(Phi_PFA)
K_PFA

# Numbers of factors - SUFA
Phi_SUFA <- readRDS("./RDS/results_curated/RCuratedOvarian_SUFA.rds")$Phi
Lambda_SUFA <- readRDS("./RDS/results_curated/RCuratedOvarian_SUFA.rds")$LambdaList

# number of factors - MOM-SS
Phi_MOMSS <- readRDS("./RDS/results_curated/RCuratedOvarian_MOMSS.rds")$Phi

# number of factors - Tetris
Phi_Tetris <- readRDS("./RDS/results_curated/RCuratedOvarian_Tetris.rds")$Phi
Lambda_Tetris <- readRDS("./RDS/results_curated/RCuratedOvarian_Tetris.rds")$LambdaList

# ----------------- MSE -----------------

# Test data has to be centered
test_list <- lapply(test_list, as.matrix)
test_list_scaled <- lapply(
  test_list, function(x) scale(x, center = TRUE, scale = FALSE)
)

# Stack FA
fit_StackFA_train <- readRDS("./RDS/results_curated/RCuratedOvarian_stackFA_train.rds")
Phi <- fit_StackFA_train$Phi
Psi <- fit_StackFA_train$Psi

mse_stackFA <- 1/(1060 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:4, function(s){
    scores <- test_list_scaled[[s]] %*% solve(Psi) %*% Phi %*% mnormt::pd.solve(signif(t(Phi) %*% solve(Psi) %*% Phi))
    norm(test_list_scaled[[s]] - scores %*% t(Phi), "F")^2
  })
)

# Ind FA (sp_fa() does not work with K=1)
fit_IndFA_train <- readRDS("./RDS/results_curated/RCuratedOvarian_IndFA2_train.rds")
LambdaList <- fit_IndFA_train$LambdaList
Psi <- fit_IndFA_train$Psi
mse_IndFA <- 1/(1060 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:4, function(s){
    scores <- test_list_scaled[[s]] %*% solve(Psi[[s]]) %*% LambdaList[[s]] %*% mnormt::pd.solve(signif(t(LambdaList[[s]]) %*% solve(Psi[[s]]) %*% LambdaList[[s]]))

    norm(test_list_scaled[[s]] - scores %*% t(LambdaList[[s]]), "F")^2
  })
)

# PFA
fit_PFA_train <- readRDS("./RDS/results_curated/RCuratedOvarian_PFA_train.rds")
Phi <- fit_PFA_train$Phi
Psi <- fit_PFA_train$Psi
Q_list <- fit_PFA_train$Q
mse_PFA <- 1/(1060 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:4, function(s){
    scores <- test_list_scaled[[s]] %*% t(Q_list[[s]]) %*% Phi %*% mnormt::pd.solve(signif(t(Phi) %*% solve(Psi) %*% Phi))
    Y_est <- scores %*% t(Phi) %*% t(solve(Q_list[[s]]))
    norm(test_list_scaled[[s]] - Y_est, "F")^2
  })
)

# BMSFA
fit_BMSFA_train <- readRDS("./RDS/results_curated/RCuratedOvarian_BMSFA_train.rds")
Phi <- fit_BMSFA_train$Phi
LambdaList <- fit_BMSFA_train$LambdaList
Psi <- fit_BMSFA_train$PsiList %>%  lapply(as.vector) %>% lapply(diag)
mse_BMSFA <- 1/(1060 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:4, function(s){
    Omega <- cbind(Phi, LambdaList[[s]])
    scores <- test_list_scaled[[s]] %*% solve(Psi[[s]]) %*% Omega %*% mnormt::pd.solve(signif(t(Omega) %*% solve(Psi[[s]]) %*% Omega))
    
    norm(test_list_scaled[[s]] - scores %*% t(Omega), "F")^2
  })
)

# SUFA
fit_SUFA_train <- readRDS("./RDS/results_curated/RCuratedOvarian_SUFA_train.rds")
Phi <- fit_SUFA_train$Phi
LambdaList <- fit_SUFA_train$LambdaList
Psi <- fit_SUFA_train$Psi
mse_SUFA <- 1/(1060 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:4, function(s){
    Omega <- cbind(Phi, LambdaList[[s]])
    scores <- test_list_scaled[[s]] %*% solve(Psi) %*% Omega %*% MASS::ginv(signif(t(Omega) %*% solve(Psi) %*% Omega))
    
    norm(test_list_scaled[[s]] - scores %*% t(Omega), "F")^2
  })
)

# MOM-SS
fit_MOMSS_train <- readRDS("./RDS/results_curated/RCuratedOvarian_MOMSS_train.rds")
Phi <- fit_MOMSS_train$Phi
Psi <- fit_MOMSS_train$Psi %>% lapply(diag)
alpha <- lapply(1:4, function(s) {
  fit_MOMSS_train$alpha[,s]
})

mse_MOMSS <- 1/(1060 * sum(sapply(test_list, nrow)))*sum(
  sapply(1:4, function(s){
    scores <- test_list_scaled[[s]] %*% solve(Psi[[s]]) %*% Phi %*% mnormt::pd.solve(signif(t(Phi) %*% solve(Psi[[s]]) %*% Phi))
    norm(test_list_scaled[[s]] - scores %*% t(Phi), "F")^2
  })
)



# Tetris
fixed <- readRDS("./RDS/results_curated/Tetris_train_run_fixed.rds") # A matrix is all 0
fit_Tetris_train <- readRDS("./RDS/results_curated/Tetris_train_run_fixed.rds")
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


# Filtering genes for visualization in Gephi
genenames <- colnames(Z_list[[1]])

# StackFA
SigmaPhi_StackFA_curated <- readRDS("./RDS/results_curated/RCuratedOvarian_stackFA2_scaled.rds")$SigmaPhi
colnames(SigmaPhi_StackFA_curated) <- rownames(SigmaPhi_StackFA_curated) <- genenames
diag(SigmaPhi_StackFA_curated) <- NA
above_thresh <- SigmaPhi_StackFA_curated > 0.85
keep_genes <- apply(above_thresh, 1, function(x) any(x, na.rm = TRUE))
sum(keep_genes)
SigmaPhi_StackFA_curated[abs(SigmaPhi_StackFA_curated) < 0.85] <- 0
SigmaPhi_StackFA_curated_filtered <- SigmaPhi_StackFA_curated[keep_genes, keep_genes]
write.csv(SigmaPhi_StackFA_curated_filtered, "SigmaPhi_StackFA_curated_filtered.csv")

# MOM-SS
SigmaPhiMOMSS <- readRDS("./RDS/results_curated/RCuratedOvarian_MOMSS_scaled.rds")$SigmaPhi
colnames(SigmaPhiMOMSS) <- rownames(SigmaPhiMOMSS) <- genenames
diag(SigmaPhiMOMSS) <- NA
above_thresh <- abs(SigmaPhiMOMSS) > 0.95
keep_genes <- apply(above_thresh, 1, function(x) any(x, na.rm = TRUE))
sum(keep_genes)
SigmaPhiMOMSS[abs(SigmaPhiMOMSS) < 0.95] <- 0
SigmaPhiMOMSS_filtered <- SigmaPhiMOMSS[keep_genes, keep_genes]
write.csv(SigmaPhiMOMSS_filtered, "SigmaPhiMOMSS_filtered.csv")

# BMSFA
SigmaPhiBMSFA <- readRDS("./RDS/results_curated/RCuratedOvarian_BMSFA2_scaled.rds")$SigmaPhi
colnames(SigmaPhiBMSFA) <- rownames(SigmaPhiBMSFA) <- genenames
diag(SigmaPhiBMSFA) <- NA
above_thresh <- abs(SigmaPhiBMSFA)  > 0.5
keep_genes <- apply(above_thresh, 1, function(x) any(x, na.rm = TRUE))
sum(keep_genes)
SigmaPhiBMSFA[abs(SigmaPhiBMSFA) < 0.5] <- 0
SigmaPhiBMSFA_filtered <- SigmaPhiBMSFA[keep_genes, keep_genes]
write.csv(SigmaPhiBMSFA_filtered, "SigmaPhiBMSFA_filtered.csv")

# SUFA
SigmaPhiSUFA <- readRDS("./RDS/results_curated/RCuratedOvarian_SUFA_scaled.rds")$SigmaPhi
colnames(SigmaPhiSUFA) <- rownames(SigmaPhiSUFA) <- genenames
diag(SigmaPhiSUFA) <- NA
above_thresh <- abs(SigmaPhiSUFA) > 0.28
keep_genes <- apply(above_thresh, 1, function(x) any(x, na.rm = TRUE))
sum(keep_genes)
SigmaPhiSUFA[abs(SigmaPhiSUFA) < 0.28] <- 0
SigmaPhiSUFA_filtered <- SigmaPhiSUFA[keep_genes, keep_genes]
write.csv(SigmaPhiSUFA_filtered, "SigmaPhiSUFA_filtered.csv")

# PFA
SigmaPhiPFA <- readRDS("./RDS/results_curated/RCuratedOvarian_PFA_scaled.rds")$SigmaPhi
colnames(SigmaPhiPFA) <- rownames(SigmaPhiPFA) <- genenames
diag(SigmaPhiPFA) <- NA
above_thresh <- abs(SigmaPhiPFA) > 0.55
keep_genes <- apply(above_thresh, 1, function(x) any(x, na.rm = TRUE))
sum(keep_genes)
SigmaPhiPFA[abs(SigmaPhiPFA) < 0.55] <- 0
SigmaPhiPFA_filtered <- SigmaPhiPFA[keep_genes, keep_genes]
write.csv(SigmaPhiPFA_filtered, "SigmaPhiPFA_filtered.csv")


# Exploratory
# Melt data for plotting
df1$study <- "GSE13876"
df2$study <- "GSE26712"
df3$study <- "GSE9891"
df4$study <- "PMID17290060"

df_combine <- 
  rbind(df1, df2, df3, df4) %>%
  reshape2::melt(id.vars = "study") 

# Define the custom color palette
custom_colors <- c(
  "#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", 
  "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30"
)

# Add the custom colors and alternate them
df_combine %>%
  ggplot(aes(x = value, color = variable)) +
  geom_density(stat = "density", alpha = 0.6) +
  facet_wrap(~study, scales = "free", ncol = 1) +
  scale_color_manual(values = rep(custom_colors, 
                                  length.out = length(unique(df_combine$variable)))) +
  labs(title = "Density Plots of Gene Expression Across Studies", 
       x = "Log Expression Level (centered)", y = "Density") +
  theme_minimal() +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(hjust = 1),
    strip.text = element_text(size = 10, face = "bold")
  )
ggsave("Density_log_centered.svg", width = 10, height = 10, dpi = 300)

# Boxplot
ggplot(df_combine, aes(x = value , y =variable, color = "darkgrey")) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 21, 
               outlier.fill = "darkslategray", 
               fill = "gray60", color = "gray40") +
  facet_wrap(~study, ncol = 1) +
  labs(
    title = "Boxplot of Expression Levels Across Studies",
    y = "Genes",
    x = "Log Expression Levels"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    legend.position = "none"
  )


# Boxplot -alternitive
# Boxplot with study-specific colors, including outliers
boxplot(value ~ variable * study, data = df_combine,
        col = rep(study_colors, each = length(unique(df_combine$variable))),
        border = rep(study_colors, each = length(unique(df_combine$variable))),
        xaxt = "n", las = 2, xlab = "", ylab = "Log Expression Levels",
        main = "Boxplot of Expression Levels Across Studies",
        outline = TRUE, cex.axis = 0.7, pch = 20, cex = 0.5)

# Customize X-axis labels (e.g., show every nth gene)
unique_genes <- unique(df_combine$variable)
n_genes <- length(unique_genes)

# Define tick positions (aligning with the boxplot groups)
tick_positions <- seq(1, n_genes * length(unique(df_combine$study)), 
                      by = length(unique(df_combine$study)))

# Add angled X-axis labels
text(x = tick_positions, y = par("usr")[3] - 0.5,  # Adjust 'y' for label placement
     labels = unique_labels, srt = 45, adj = 1, xpd = TRUE, cex = 0.8) 


# Explore sample phenotype data
pD <- phenoData(study2)
varLabels(pD)
table(pData(pD)$summarygrade)





