# Clear session and load required libraries
rm(list = ls())
dir.create("9-monocle", showWarnings = FALSE)
setwd("9-monocle")

library(Seurat)
library(monocle)
library(dplyr)
library(tidyverse)
library(Matrix)
library(cowplot)
library(clustree)
library(ggplot2)

# Load Seurat object with epithelial subclusters
seu_obj <- readRDS("../7-epi/epi_sce_celltype.rds")
Idents(seu_obj) <- seu_obj$celltype
table(Idents(seu_obj))

# Optionally downsample cells per group (commented out for now)
# sampled_cells <- unlist(lapply(levels(Idents(seu_obj)), function(ct) {
#   sample(WhichCells(seu_obj, idents = ct), size = 10)
# }))
# seu_obj <- subset(seu_obj, cells = sampled_cells)

# -----------------------------
# Convert Seurat to Monocle CellDataSet
# -----------------------------
expr_mat <- as(GetAssayData(seu_obj, slot = "counts"), "sparseMatrix")

gene_annot <- data.frame(
  gene_id = rownames(expr_mat),
  gene_short_name = rownames(expr_mat),
  stringsAsFactors = FALSE
)
rownames(gene_annot) <- gene_annot$gene_id
fData <- new("AnnotatedDataFrame", data = gene_annot)

cell_annot <- seu_obj@meta.data
pData <- new("AnnotatedDataFrame", data = cell_annot)

cds <- newCellDataSet(
  expr_mat,
  phenoData = pData,
  featureData = fData,
  expressionFamily = negbinomial.size()
)

# -----------------------------
# Estimate size factors and dispersions
# -----------------------------
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# -----------------------------
# Select ordering genes based on expression
# -----------------------------
disp_table <- dispersionTable(cds)
order_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, order_genes$gene_id)

# -----------------------------
# Dimensionality reduction (DDRTree) and trajectory ordering
# -----------------------------
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)

# -----------------------------
# Visualize pseudotime trajectory
# -----------------------------
p1 <- plot_cell_trajectory(cds, color_by = "celltype", cell_size = 1)
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1)
p1 + p2
ggsave("monocle.pdf", width = 12, height = 8)

# Facet by sample
plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~sample", nrow = 1)
ggsave("monocle_sample.pdf", width = 14, height = 8)

# -----------------------------
# Manually assign root (State 3 here as example)
# -----------------------------
cds_rooted <- orderCells(cds, root_state = 3)

# Add pseudotime gene expression (e.g., TGFBR2) to metadata
pData(cds_rooted)$TGFBR2 <- log2(exprs(cds_rooted)["TGFBR2", ] + 1)

# Plot trajectory again with defined root
p1 <- plot_cell_trajectory(cds_rooted, color_by = "celltype", cell_size = 1)
p2 <- plot_cell_trajectory(cds_rooted, color_by = "Pseudotime", cell_size = 1)
p1 + p2

# Return to project root directory
setwd("..")
