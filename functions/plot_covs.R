library(ggplot2)
library(reshape2)
library(gridExtra)
library(latex2exp)
library(dplyr)
library(patchwork)
all_methods <- c("stackFA", "IndFA", "PFA", "MOMSS", "SUFA_fixJs", 
                 "SUFA", "BMSFA", "Tetris_fixT", "Tetris")
# This function plots the heatmaps for the estimated covariance matrices (both common and 
# study-specific) of all the fitted models and the true model, within a single plot.
plot_covs <- function(fitted_lists, sim_data) {
  # Arrange data
  phi_data <- arrange_SigmaPhi(fitted_lists, sim_data)
  lambda_data <- arrange_SigmaLambda(fitted_lists, sim_data)
  marginal_data <- arrange_SigmaMarginal(fitted_lists, sim_data)
  
  # Convert the CovType and Method factors to use LaTeX in facet titles
  phi_data$CovType <- factor(phi_data$CovType, labels = c(TeX("$\\Sigma_{\\Phi}$")))
  lambda_data$CovType <- factor(lambda_data$CovType, labels = sapply(1:length(sim_data$SigmaMarginal),
                                                                     function(i) TeX(paste0("$\\Sigma_{\\Lambda_{", i, "}}$"))))
  marginal_data$CovType <- factor(marginal_data$CovType, labels = sapply(1:length(sim_data$SigmaMarginal),
                                                                         function(i) TeX(paste0("$\\Sigma_{", i, "}$"))))
  
  phi_data$Method <- factor(phi_data$Method, levels = c("True", all_methods))
  lambda_data$Method <- factor(lambda_data$Method, levels = c("True", all_methods))
  marginal_data$Method <- factor(marginal_data$Method, levels = c("True", all_methods))
  
  # Plot for SigmaPhi (Common Covariance)
  p_phi <- ggplot(phi_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(CovType ~ Method, labeller = label_parsed, switch = "y") +
    scale_x_continuous(breaks = NULL) +  # Remove x ticks
    scale_y_continuous(breaks = NULL) +  # Remove y ticks
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, size = 12, margin = margin(r = 8, l = 10)),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),  # Remove axis titles
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          plot.title = element_text(hjust = 0.5)) +  # Center title
    scale_fill_gradient2(low = "#a1d76a", mid = "#f7f7f7", high = "#e9a3c9", midpoint = 0, na.value = "grey50", name = TeX("$\\Sigma_{\\Phi}$\\ values")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4)) +
    ggtitle("Common Covariance")
  
  # Plot for SigmaLambda (Study-Specific Covariance)
  p_lambda <- ggplot(lambda_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(CovType ~ Method, labeller = label_parsed, switch = "y") +
    scale_x_continuous(breaks = NULL) +  # Remove x ticks
    scale_y_continuous(breaks = NULL) +  # Remove y ticks
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, size = 12, margin = margin(r = 8, l = 10)),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),  # Remove axis titles
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          plot.title = element_text(hjust = 0.5)) +  # Center title
    scale_fill_gradient2(low = "#998ec3", mid = "#f7f7f7", high = "#f1a340", 
                         midpoint = 0, na.value = "grey50", name = TeX("$\\Sigma_{\\Lambda}$\\ values")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4)) +
    ggtitle("Study-Specific Covariance")
  
  # Plot for SigmaMarginal (Marginal Covariance)
  p_marginal <- ggplot(marginal_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(CovType ~ Method, labeller = label_parsed, switch = "y") +
    scale_x_continuous(breaks = NULL) +  # Remove x ticks
    scale_y_continuous(breaks = NULL) +  # Remove y ticks
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, size = 12, margin = margin(r = 8, l = 10)),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),  # Remove axis titles
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          plot.title = element_text(hjust = 0.5)) +  # Center title
    scale_fill_gradient2(low = "#67a9cf", mid = "#f7f7f7", high = "#ef8a62", midpoint = 0, 
                         na.value = "grey50", name = TeX("$\\Sigma_{s}$\\ values")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4)) +
    ggtitle("Marginal Covariance")
  
  # Combine the plots using patchwork
  p_combined <- p_phi / p_lambda / p_marginal + plot_layout(ncol = 1, heights = c(0.15, 0.425, 0.425))
  
  return(p_combined)
}

arrange_SigmaPhi <- function(fitted_lists, sim_data) {
  nrow_true <- ncol(sim_data$SigmaPhi)
  
  phi_melted_list <- lapply(c("True", all_methods), function(method) {
    if (method == "True") {
      SigmaPhi_matrix <- sim_data$SigmaPhi
    } else if (!is.null(fitted_lists[[method]]$SigmaPhi)) {
      SigmaPhi_matrix <- fitted_lists[[method]]$SigmaPhi
    } else {
      SigmaPhi_matrix <- matrix(NA, nrow = nrow_true, ncol = nrow_true)
    }
    
    df <- melt(SigmaPhi_matrix)
    df$CovType <- "SigmaPhi"
    df$Method <- method
    return(df)
  })
  
  phi_melted <- do.call(rbind, phi_melted_list)
  return(phi_melted)
}

arrange_SigmaLambda <- function(fitted_lists, sim_data) {
  num_studies <- length(sim_data$SigmaMarginal)
  nrow_true <- ncol(sim_data$SigmaMarginal[[1]])
  
  lambda_melted_list <- lapply(c("True", all_methods), function(method) {
    if (method == "True") {
      if (is.null(sim_data$SigmaLambdaList)) {
        SigmaLambdaList <- replicate(num_studies, matrix(NA, nrow = nrow_true, ncol = nrow_true), simplify = FALSE)
      } else {
        SigmaLambdaList <- sim_data$SigmaLambdaList
      }
    } else if (!is.null(fitted_lists[[method]]$SigmaLambdaList)) {
      SigmaLambdaList <- fitted_lists[[method]]$SigmaLambdaList
    } else {
      SigmaLambdaList <- replicate(num_studies, matrix(NA, nrow = nrow_true, ncol = nrow_true), simplify = FALSE)
    }
    
    melted_list <- lapply(seq_along(SigmaLambdaList), function(i) {
      # if (is.null(SigmaLambdaList[[i]]) || all(is.na(SigmaLambdaList[[i]]))) {
      #   return(NULL)
      # }
      df <- melt(SigmaLambdaList[[i]])
      if (nrow(df) == 0) {
        SigmaLambdaList[[i]] <- matrix(0, nrow = K, ncol = max_ncol)
        df <- melt(SigmaLambdaList[[i]])
      }
      df$CovType <- paste0("SigmaLambda_", i)
      df$Method <- method
      return(df)
    })
    
    return(do.call(rbind, melted_list))
  })
  

  lambda_melted <- do.call(rbind, lambda_melted_list)
  return(lambda_melted)
}

arrange_SigmaMarginal <- function(fitted_lists, sim_data) {
  num_studies <- length(sim_data$SigmaMarginal)
  nrow_true <- ncol(sim_data$SigmaMarginal[[1]])
  
  marginal_melted_list <- lapply(c("True", all_methods), function(method) {
    if (method == "True") {
      SigmaMarginalList <- sim_data$SigmaMarginal
    } else if (!is.null(fitted_lists[[method]]$SigmaMarginal)) {
      SigmaMarginalList <- fitted_lists[[method]]$SigmaMarginal
    } else {
      SigmaMarginalList <- replicate(num_studies, matrix(NA, nrow = nrow_true, ncol = nrow_true), simplify = FALSE)
    }
    
    melted_list <- lapply(seq_along(SigmaMarginalList), function(i) {
      if (is.null(SigmaMarginalList[[i]]) || all(is.na(SigmaMarginalList[[i]]))) {
        return(NULL)
      }
      df <- melt(SigmaMarginalList[[i]])
      df$CovType <- paste0("SigmaMarginal_", i)
      df$Method <- method
      return(df)
    })
    
    return(do.call(rbind, melted_list))
  })
  
  marginal_melted <- do.call(rbind, marginal_melted_list)
  return(marginal_melted)
}

# Example usage
#plot_covs(postfit_scPFA_list, sim_PFA1)
