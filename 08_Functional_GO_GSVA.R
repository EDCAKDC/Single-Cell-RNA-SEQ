# Clear workspace and load packages
rm(list = ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(clustree)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSVA)
library(msigdbr)
library(limma)
library(patchwork)
library(reshape2)
library(pheatmap)
library(ggthemes)

source("scRNA_scripts/lib.R")
source("scRNA_scripts/mycolors.R")

# Set working directory
setwd("7-epi/")
seurat_obj <- readRDS("epi_sce_celltype.rds")
Idents(seurat_obj) <- seurat_obj$celltype

# ------------------------
# GO Enrichment Analysis
# ------------------------

# Load markers or calculate if not available
if (!file.exists("epimarkers.csv")) {
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  write.csv(markers, "epimarkers.csv")
} else {
  markers <- read.csv("epimarkers.csv", row.names = 1)
}

# Filter genes for enrichment (p < 0.01 & log2FC > 2)
selected_genes <- markers %>%
  filter(p_val < 0.01 & avg_log2FC > 2) %>%
  pull(gene)

# Convert gene symbols to Entrez IDs
gene_map <- bitr(selected_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
filtered_markers <- markers[markers$gene %in% gene_map$SYMBOL, ]
filtered_markers <- filtered_markers[!duplicated(filtered_markers$gene), ]
filtered_markers$ENTREZID <- gene_map$ENTREZID[match(filtered_markers$gene, gene_map$SYMBOL)]

# Group genes by cluster
gene_sets <- split(filtered_markers$ENTREZID, filtered_markers$cluster)

# Run GO Biological Process enrichment
go_result <- compareCluster(
  geneCluster = gene_sets,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Visualize GO enrichment
dotplot(go_result) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# ------------------------
# GSVA Pathway Enrichment
# ------------------------

# Retrieve MSigDB C2 gene sets
msig_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame() %>%
  split(.$gs_name)

# Get average expression for each cell type
avg_expr <- AverageExpression(seurat_obj, assays = "RNA", slot = "data")[[1]]
avg_expr <- avg_expr[rowSums(avg_expr) > 0, ]  # filter 0-sum genes
avg_expr <- as.matrix(avg_expr)

# Run GSVA
gsva_scores <- gsva(expr = avg_expr, gset.idx.list = msig_gene_sets, method = "gsva")
saveRDS(gsva_scores, "gsva.res.rds")
write.csv(data.frame(Geneset = rownames(gsva_scores), gsva_scores), "gsva_res.csv", row.names = FALSE)

# Subset and visualize random 30 gene sets
set.seed(123)
selected_pathways <- sample(rownames(gsva_scores), 30)
heat_data <- gsva_scores[selected_pathways, ]

pheatmap(
  mat = heat_data,
  scale = "row",
  show_colnames = TRUE,
  angle_col = 45,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)

# Convert GSVA results to long format
gsva_melted <- melt(heat_data)
colnames(gsva_melted) <- c("Pathway", "Celltype", "Score")

# Bubble plot of pathway activity
ggplot(gsva_melted, aes(x = Celltype, y = Pathway, size = Score, color = Score)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(1, 6)) +
  scale_color_gradient(low = "#427183", high = "#D2D470") +
  theme_minimal() +
  labs(x = "Cell Type", y = "Pathway", size = "GSVA Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# Return to parent directory
setwd("..")
