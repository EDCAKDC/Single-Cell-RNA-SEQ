# Clear the workspace
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(clustree)
  library(cowplot)
  library(dplyr)
  library(SingleR)
  library(celldex)
  library(singleseqgset)
  library(devtools)
  library(grid)
  library(gridExtra)
})

# Set working directory and create output folder
output_dir <- "4-plot"
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

# Load Seurat object with cell type annotations
seurat_obj <- readRDS("../3-Celltype/sce_celltype.rds")

# Define marker genes for visualization
marker_genes <- c('EPCAM', 'NKG7', 'LYZ', 'CD79A', 'CLDN5', 'DCN')

# Basic expression plots for each marker
fp_basic <- FeaturePlot(seurat_obj,
                        features = marker_genes,
                        cols = c("lightgrey", "#DE1F1F"),
                        ncol = 3,
                        raster = FALSE)
ggsave("FeaturePlot_marker.pdf", plot = fp_basic, width = 12, height = 8)

# FeaturePlot with fixed scale range for all genes
plot_list <- FeaturePlot(seurat_obj, features = marker_genes, combine = FALSE, raster = FALSE)
fixed_scale <- scale_color_gradientn(colors = c("lightgrey", "#DE1F1F"), limits = c(0, 6))
plot_list_scaled <- lapply(plot_list, \(p) p + fixed_scale)
CombinePlots(plot_list_scaled)

# Clean visualization: no legend/axes, black panel border
FeaturePlot(seurat_obj,
            features = marker_genes,
            cols = c("lightgrey", "red"),
            ncol = 3) &
  NoLegend() & NoAxes() &
  theme(panel.border = element_rect(color = "black", size = 1))

# Blend expression of two genes to highlight co-expression
FeaturePlot(seurat_obj,
            features = c("S100A9", "S100A8"),
            cols = c("lightgrey", "green", "orange"),
            blend = TRUE,
            blend.threshold = 0)

# UMAP colored by multiple genes in sequence using ggnewscale
library(ggnewscale)

umap_df <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
umap_df$celltype <- seurat_obj$celltype

gene_data <- GetAssayData(seurat_obj, slot = "data")[c("S100A9", "S100A8", "CXCL8"), ]
plot_df <- cbind(t(as.matrix(gene_data)), umap_df)

ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = S100A9)) +
  geom_point(size = 0.3, alpha = 1) +
  scale_color_gradientn(colors = c("lightgrey", "green"), limits = c(0, 0.3), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(color = S100A8), size = 0.3, alpha = 0.7) +
  scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.1, 0.2), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(color = CXCL8), size = 0.3, alpha = 0.1) +
  scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, 0.3), oob = scales::squish) +
  theme_classic()

# -------------------------
# Heatmap of top markers
# -------------------------

# Set identity to cell type
seurat_obj <- SetIdent(seurat_obj, value = "celltype")
print(table(Idents(seurat_obj)))

# Keep specific cell types for visualization
selected_types <- c("B", "Endothelial", "Epithelial", "Fibro", "Myeloid", "T&NK")
subset_obj <- subset(seurat_obj, idents = selected_types)

# Perform differential expression or load existing results
if (!file.exists("sce.markers.csv")) {
  marker_results <- FindAllMarkers(subset_obj,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25)
  write.csv(marker_results, "sce.markers.csv")
} else {
  marker_results <- read.csv("sce.markers.csv", row.names = 1)
}

# Select top 5 genes per cluster
top_markers <- marker_results %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Downsample to 100 cells per cluster for better heatmap readability
downsampled_obj <- subset(subset_obj, downsample = 100)
scaled_obj <- ScaleData(downsampled_obj, features = top_markers$gene)

# Plot heatmap
heatmap_plot <- DoHeatmap(scaled_obj,
                          features = top_markers$gene,
                          assay = "RNA",
                          label = TRUE) +
  scale_fill_gradientn(colors = c("white", "grey", "firebrick3"))
ggsave("markers_heatmap.pdf", plot = heatmap_plot, width = 10, height = 7)

# DotPlot for same markers
dot_plot <- DotPlot(seurat_obj, features = unique(top_markers$gene)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
ggsave("markers_top5_dotplot.pdf", plot = dot_plot, width = 10, height = 7)

# Return to parent directory
setwd("..")
