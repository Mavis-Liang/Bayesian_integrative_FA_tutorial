library(ggplot2)
library(reshape2)
library(gridExtra)
library(latex2exp)
library(dplyr)
library(patchwork)

# This function plots the heatmaps for the loading matrices (both common and 
# study-specific) of all the fitted models and the true model, within a single plot.

# Main function for plotting the loading matrices
plot_loadings <- function(fitted_lists, sim_data) {
  # Prepare Phi and Lambda data
  phi_data <- arrange_Phi(fitted_lists, sim_data)
  lambda_data <- arrange_Lambda(fitted_lists, sim_data)
  
  # Convert the Lambda and Method factors to use LaTeX in facet titles
  phi_data$Lambda <- factor(phi_data$Lambda, labels = c(TeX("$\\Phi$")))
  lambda_data$Lambda <- factor(lambda_data$Lambda, labels = sapply(unique(lambda_data$Lambda), function(i) TeX(paste0("$\\Lambda_", gsub("Lambda_", "", i), "$"))))
  
  phi_data$Method <- factor(phi_data$Method, levels = c("True", names(fitted_lists)))
  lambda_data$Method <- factor(lambda_data$Method, levels = c("True", names(fitted_lists)))
  
  # Plot for Phi matrices
  # Determine the range of the data
  x_range <- range(phi_data$Var2, na.rm = TRUE)
  # Set intervals to avoid overlap 
  x_breaks <- seq(x_range[1], x_range[2], by = ceiling((x_range[2] - x_range[1])/5))

  p_phi <- ggplot(phi_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(Lambda ~ Method, labeller = label_parsed, switch = "y") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_reverse(breaks = seq(min(phi_data$Var1, na.rm = TRUE), max(phi_data$Var1, na.rm = TRUE), by = 1), labels = NULL) +
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, size = 12, margin = margin(r = 8, l = 10)),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),  # Remove axis titles
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          plot.title = element_text(hjust = 0.5)) +  # Center title
    scale_fill_gradient2(low = "navy", mid = "white", high = "darkgoldenrod1", midpoint = 0, na.value = "grey80", name = TeX("$\\phi$\\ values")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4))+
    ggtitle("Common Loading matrix")
  
  
  # Plot for Lambda matrices
  x_range <- range(lambda_data$Var2, na.rm = TRUE)
  if(x_range[2] - x_range[1] > 15) {
    x_breaks = NULL
  } else {
    x_breaks <- seq(x_range[1], x_range[2], by = ceiling((x_range[2] - x_range[1])/5))
  }
  p_lambda <- ggplot(lambda_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(Lambda ~ Method, labeller = label_parsed, switch = "y") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_reverse(breaks = seq(min(lambda_data$Var1, na.rm = TRUE), max(lambda_data$Var1, na.rm = TRUE), by = 1), labels = NULL) +
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, size = 12, margin = margin(r = 8, l = 10)),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +  # Center title
    scale_fill_gradient2(low = "darkslateblue", mid = "white", high = "chocolate", midpoint = 0, na.value = "grey80", name = TeX("$\\lambda$\\ values")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4))+
    ggtitle("Study-specific Loading matrix")
  
  # Combine the plots using patchwork
  p_combined <- p_phi / p_lambda + plot_layout(ncol = 1, heights = c(0.33, 1))  # Adjust heights to compress Phi plot
  
  return(p_combined)
}

arrange_Phi <- function(fitted_lists, sim_data) {
  # Number of rows and columns in Phi (based on the true data)
  nrow_true <- nrow(sim_data$Phi)
  max_ncol <- max(ncol(sim_data$Phi), sapply(fitted_lists, function(fit) {
    if (!is.null(fit$Phi)) {
      ncol(fit$Phi)
    } else {
      0
    }
  }), na.rm = TRUE)
  
  # Prepare the Phi data for true and fitted methods
  phi_melted_list <- lapply(c("True", names(fitted_lists)), function(method) {
    if (method == "True") {
      Phi_matrix <- sim_data$Phi
    } else if (!is.null(fitted_lists[[method]]$Phi)) {
      Phi_matrix <- fitted_lists[[method]]$Phi
    } else {
      # Create an NA matrix if Phi does not exist for this method
      Phi_matrix <- matrix(NA, nrow = nrow_true, ncol = max_ncol)
    }
    
    df <- melt(Phi_matrix)
    df$Lambda <- "Phi"
    df$Method <- method
    return(df)
  })
  
  # Combine all Phi data
  phi_melted <- do.call(rbind, phi_melted_list)
  
  return(phi_melted)
}



arrange_Lambda <- function(fitted_lists, sim_data) {
  # Number of studies in sim_data
  S <- length(sim_data$SigmaMarginal)
  
  # Determine the number of rows (from sim_data$Phi) and the maximum number of columns
  K <- nrow(sim_data$Phi)
  max_ncol <- max(sapply(fitted_lists, function(fit) {
    if (!is.null(fit$LambdaList)) {
      max(sapply(fit$LambdaList, ncol))
    } else {
      0
    }
  }), na.rm = TRUE)
  
  # Handle the case where true LambdaList is NULL
  if (is.null(sim_data$LambdaList)) {
    true_LambdaList <- replicate(S, matrix(NA, nrow = K, ncol = max_ncol), simplify = FALSE)
  } else {
    true_LambdaList <- sim_data$LambdaList
  }
  
  # Prepare the true data
  true_melted <- lapply(seq_len(S), function(i) {
    df <- melt(true_LambdaList[[i]])
    df$Lambda <- paste0("Lambda_", i)
    df$Method <- "True"
    return(df)
  })
  
  true_melted <- do.call(rbind, true_melted)
  
  # Prepare the fitted data for all methods
  fitted_melted_list <- lapply(names(fitted_lists), function(method) {
    if (!is.null(fitted_lists[[method]]$LambdaList)) {
      LambdaList <- fitted_lists[[method]]$LambdaList
    } else {
      # Create a list of NA matrices if LambdaList doesn't exist
      LambdaList <- replicate(S, matrix(NA, nrow = K, ncol = max_ncol), simplify = FALSE)
    }
    
    melted_list <- lapply(seq_along(LambdaList), function(i) {
      df <- melt(LambdaList[[i]])
      if (nrow(df) == 0) {
        LambdaList[[i]] <- matrix(0, nrow = K, ncol = max_ncol)
        df <- melt(LambdaList[[i]])
      }
      df$Lambda <- paste0("Lambda_", i)
      df$Method <- method
      return(df)
    })
    
    do.call(rbind, melted_list)
  })
  
  fitted_melted <- do.call(rbind, fitted_melted_list)
  
  # Combine true and fitted Lambda data
  combined_melted <- rbind(true_melted, fitted_melted)
  
  return(combined_melted)
}


# 
# # Example usage:
# # sc3
# fit_sc3 <- readRDS("./RDS/sc3/sc3_1.rds")
# plot_loadings(fit_sc3$point_est, fit_sc3$true)
# ggplot2::ggsave("./Figs/sc3_loadings.png", width = 6, height = 7, units = "in")
# # sc2
# fit_sc2 <- readRDS("./RDS/sc2/sc2_1.rds")
# plot_loadings(fit_sc2$point_est, fit_sc2$true)
# ggplot2::ggsave("./Figs/sc2_loadings.png", width = 8, height = 7, units = "in")
# # sc1
# fit_sc1 <- readRDS("./RDS/sc1/sc1_1.rds")
# plot_loadings(fit_sc1$point_est, fit_sc1$true)
# ggplot2::ggsave("./Figs/sc1_loadings.png", width = 6, height = 7, units = "in")
