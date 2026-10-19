# Load required libraries
library(ggplot2)
library(reshape2)

# Title and axis labels
title <- "Cell type proportions of BiSC mothers DNAm - mealtiming (N=389)"
axis_x <- "Samples"
axis_y <- "Proportions"

# Load metadata
metadata_path <- "db/metadata_plaEWAS_MT_complete_cases_with_scores_110012025.rds"
metadata <- readRDS(metadata_path)

# Define cell types
cell_type_vars <- c('Basename', 'Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC', 'Syncytiotrophoblast')
cell_type_data <- metadata[, cell_type_vars]

# Order samples by Syncytiotrophoblast proportion
ordered_samples <- cell_type_data[order(cell_type_data$Syncytiotrophoblast, decreasing = FALSE), "Basename"]
cell_type_data$Basename <- factor(cell_type_data$Basename, levels = ordered_samples)

# Melt data
cell_type_long <- melt(cell_type_data, id.vars = "Basename")
colnames(cell_type_long) <- c("Basename", "CellType", "Proportion")

# Reverse the stacking order (Syncytiotrophoblast on top)
reversed_order <- c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC', 'Syncytiotrophoblast')
cell_type_long$CellType <- factor(cell_type_long$CellType, levels = reversed_order)

# Define matching colors (in same reverse order)
cell_type_colors <- c("#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4")

# Plot
cell_type_plot <- ggplot(cell_type_long, aes(x = Basename, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = cell_type_colors) +
  labs(title = title, x = axis_x, y = axis_y, fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  ) +
  ylim(0, 1)

# Save the plot
save_to <- "results/2_quality_control/hist_CTprop_syncytio_top.png"
ggsave(save_to, plot = cell_type_plot, width = 7, height = 10, dpi = 300)

# Print plot
print(cell_type_plot)
