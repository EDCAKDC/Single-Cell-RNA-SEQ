# Clear environment and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(clustree)
  library(celldex)
  library(singleseqgset)
  library(devtools)
  library(tidyr)
  library(ggthemes)
  library(gplots)
})

# Create output folder and set working directory
output_dir <- "5-prop"
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

# Load annotated Seurat object with celltype and sample metadata
seurat_obj <- readRDS("../3-Celltype/sce_celltype.rds")

# Optional: Set factor level order for samples
seurat_obj$sample <- factor(seurat_obj$sample, levels = c("A", "B", "C"))

# -------------------------------
# Proportion visualization section
# -------------------------------

# Tabulate number of cells by (sample × celltype)
cell_counts <- table(seurat_obj$sample, seurat_obj$celltype)

# Balloon plot for quick inspection
balloonplot(cell_counts)

# Prepare data for bar plot
bar_df <- as.data.frame(cell_counts)
colnames(bar_df) <- c("Sample", "CellType", "Count")

# Compute relative proportions
bar_df <- bar_df %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count),
         Proportion = Count / Total)

# Custom color palette for cell types
cell_colors <- c(
  "#3176B7", "#F78000", "#3FA116", "#CE2820", "#9265C1",
  "#885649", "#DD76C5", "#BBBE00", "#41BED1"
)

# Plot relative proportion (percentage-based bar)
plot_relative <- ggplot(bar_df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = cell_colors) +
  labs(x = NULL, y = "Relative proportion", fill = NULL) +
  theme_few() +
  theme(
    legend.position = "top",
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Plot absolute cell counts
plot_absolute <- ggplot(bar_df, aes(x = Sample, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = cell_colors) +
  labs(x = NULL, y = "Cell count", fill = NULL) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.ticks = element_blank()
  )

# Combine the two plots
library(patchwork)
combined_plot <- plot_relative + plot_absolute
print(combined_plot)

# Export combined plot to PDF
ggsave("prop.pdf", plot = combined_plot, width = 12, height = 7)

# Restore working directory
setwd("..")
