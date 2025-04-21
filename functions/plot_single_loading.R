library(ggplot2)
library(reshape2)
library(gridExtra)
library(latex2exp)
library(dplyr)
# Main function for plotting the loading matrices
plot_single <- function(matrix, y_label = NULL, fill_limits = c(-1.3, 1.5), show_colorbar = FALSE) {
  df <- melt(matrix)
  # Determine the range of the data
  x_range <- range(df$Var2, na.rm = TRUE)
  y_range <- range(df$Var1, na.rm = TRUE)
  x_breaks <- seq(x_range[1], x_range[2], by = ceiling((x_range[2] - x_range[1]) / 5))
  # Set intervals to avoid overlap 
  if (is.null(y_label)) {
    y_breaks <- seq(y_range[2], y_range[1], by = -5)  # Use breaks every 5
    y_labels <- y_breaks  # Labels same as breaks
  } else {
    y_breaks <- seq_along(y_label)  # Use indices of the labels as breaks
    y_labels <- y_label  # Use the provided custom labels
  }

  # Plot for Phi matrices
  p <- ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_reverse(
      breaks = y_breaks,  # Define where the ticks appear
      labels = y_label  # Define custom labels
    )  +  # Reverse Y-axis with intervals
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, size = 15, margin = margin(r = 8, l = 10)),
          axis.text.x = element_text(size = 8, margin = margin(t = -15)),  # Make X-axis numbers smaller
          axis.text.y = element_text(size = 7),  # Make Y-axis numbers smaller
          axis.ticks = element_blank(),  # Remove axis ticks
          axis.title = element_blank(),  # Remove axis titles
          #legend.position = "none",  # Move legend to the right
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          panel.grid = element_blank(),                     # Remove all grid lines
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 13, margin = margin(b = -10))) +  # Center title
    scale_fill_gradient2(
      low = "#998ec3", mid = "#f7f7f7", high = "#f1a340", 
      midpoint = 0,
      na.value = "grey50",
      limits = fill_limits,                  # Set unified color bar limits
      name = TeX("Loadings")
    ) +
    guides(fill = guide_colorbar(barwidth = 0.3, barheight = 3, title.theme = element_text(size = 7),
                                 label.theme = element_text(size = 6)))
  
  if (!show_colorbar) {
    p <- p + guides(fill = "none")  # Turn off the color bar
  }
  
  return(p)
}


# Usage
# 
# grid.arrange(plot_single(sim_B_small$LambdaList[[1]])+
#                ggtitle(TeX("$\\Lambda_1$"), size = 12),
#              plot_single(sim_B_small$LambdaList[[2]])+
#                ggtitle(TeX("$\\Lambda_2$")),
#              plot_single(sim_B_small$LambdaList[[3]])+
#                ggtitle(TeX("$\\Lambda_3$")),
#              nrow = 1)
# 


