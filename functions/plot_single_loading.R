library(ggplot2)
library(reshape2)
library(gridExtra)
library(latex2exp)
library(dplyr)
# Main function for plotting the loading matrices
plot_single <- function(matrix) {
  df <- melt(matrix)
  
  # Determine the range of the data
  x_range <- range(df$Var2, na.rm = TRUE)
  y_range <- range(df$Var1, na.rm = TRUE)
  
  # Set intervals to avoid overlap 
  x_breaks <- seq(x_range[1], x_range[2], by = ceiling((x_range[2] - x_range[1])/5))
  y_breaks <- seq(y_range[1], y_range[2], by = ceiling((y_range[2] - y_range[1])/10))
  
  # Plot for Phi matrices
  p <- ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_reverse(breaks = y_breaks) +  # Reverse Y-axis with intervals
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, size = 15, margin = margin(r = 8, l = 10)),
          axis.text.x = element_text(size = 10),  # Make X-axis numbers smaller
          axis.text.y = element_text(size = 10),  # Make Y-axis numbers smaller
          axis.ticks = element_blank(),  # Remove axis ticks
          axis.title = element_blank(),  # Remove axis titles
          #legend.position = "none",  # Move legend to the right
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          plot.title = element_text(hjust = 0.5, size = 30)) +  # Center title
    scale_fill_gradient2(low = "royalblue", mid = "white", high = "darkorange", midpoint = 0, 
                         na.value = "grey50", name = TeX("Values")) +
    guides(fill = guide_colorbar(barwidth = 0.6, barheight = 4, title.theme = element_text(size = 10),
                                 label.theme = element_text(size = 8)))
  
  return(p)
}


# Usage
plot_single(sim_PFA1$Phi)+
  ggtitle(TeX("\\Phi"))
# 
# grid.arrange(plot_single(sim_B_small$LambdaList[[1]])+
#                ggtitle(TeX("$\\Lambda_1$"), size = 12),
#              plot_single(sim_B_small$LambdaList[[2]])+
#                ggtitle(TeX("$\\Lambda_2$")),
#              plot_single(sim_B_small$LambdaList[[3]])+
#                ggtitle(TeX("$\\Lambda_3$")),
#              nrow = 1)
# 
grid.arrange(plot_single(fit_sc1$true$Psi_list[[1]])+
               ggtitle(TeX("$\\Psi_1$")),
             plot_single(fit_sc1$true$Psi_list[[2]])+
               ggtitle(TeX("$\\Psi_2$")),
             plot_single(fit_sc1$true$Psi_list[[3]])+
               ggtitle(TeX("$\\Psi_3$")),
             plot_single(fit_sc1$true$Psi_list[[4]])+
               ggtitle(TeX("$\\Psi_4$")),
             nrow = 2)
# sc1
plot_single(fit_sc1$true$Phi) +
  ggtitle(TeX("$\\Phi$"))

