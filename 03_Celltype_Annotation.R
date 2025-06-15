# Clear workspace
rm(list = ls())

# Load helper scripts
source("scRNA_scripts/lib.R")
source("scRNA_scripts/mycolors.R")

# Load integrated Seurat object (output from Harmony integration)
seurat_obj <- readRDS("2-harmony/sce.all_int.rds")

# Use clustering resolution for downstream labeling
cluster_res <- "RNA_snn_res.0.5"
seurat_obj <- SetIdent(seurat_obj, value = cluster_res)
print(table(Idents(seurat_obj)))

# Prepare working directory
output_dir <- "./3-Celltype"
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

# Rename for clarity
scRNA <- seurat_obj

# Define markers for different lineages
marker_genes <- c(
  # Epithelial
  "EPCAM", "KRT19", "CLDN4", "SCGB1A1",
  # Endothelial
  "PECAM1", "COL1A2", "VWF", "CDH5", "CLDN5",
  # Fibroblasts
  "LUM", "FGF7", "MME",
  # T cells
  "CD3D", "CD3E", "CD8A", "CD4", "CD2",
  # Myeloid
  "AIF1", "C1QC", "C1QB", "LYZ", "CD68", "CSF1R", "CSF3R",
  # Cycling cells
  "MKI67", "STMN1", "PCNA",
  # Mast cells
  "CPA3", "CST3", "KIT", "TPSAB1", "TPSB2",
  # Neutrophils
  "GOS2", "S100A9", "S100A8", "CXCL8",
  # NK cells
  "KLRD1", "GNLY", "KLRF1", "AREG", "XCL2", "HSPA6",
  # B/Plasma
  "MS4A1", "CD19", "CD79A", "IGHG1", "MZB1", "SDC1", "IGHD"
) %>% toupper()

# Plot markers to guide cell type annotation
dotplot <- DotPlot(scRNA, features = unique(marker_genes)) + coord_flip()
print(dotplot)

# Optional: annotate UMAP axis direction
source("../scRNA_scripts/Bottom_left_axis.R")
axis_info <- left_axes(scRNA)

# UMAP by cluster
umap_cluster <- DimPlot(scRNA, reduction = "umap", group.by = cluster_res,
                        cols = my36colors, label = TRUE, label.box = TRUE, pt.size = 0.8) +
  NoAxes() +
  geom_line(data = axis_info$axes, aes(x = x, y = y, group = group),
            arrow = arrow(length = unit(0.1, "inches"), ends = "last")) +
  geom_text(data = axis_info$label, aes(x = x, y = y, label = lab, angle = angle), fontface = "italic") +
  theme(aspect.ratio = 1, plot.title = element_blank())
ggsave("RNA_snn_res.0.5_umap.pdf", umap_cluster, width = 9, height = 7)

# UMAP by sample
umap_sample <- DimPlot(scRNA, reduction = "umap", group.by = "sample",
                       cols = my36colors, pt.size = 0.8) +
  NoAxes() +
  geom_line(data = axis_info$axes, aes(x = x, y = y, group = group),
            arrow = arrow(length = unit(0.1, "inches"), ends = "last")) +
  geom_text(data = axis_info$label, aes(x = x, y = y, label = lab, angle = angle), fontface = "italic") +
  theme(aspect.ratio = 1, plot.title = element_blank())
ggsave("sample_umap.pdf", umap_sample, width = 9, height = 7)

# Define mapping from cluster IDs to cell types
cluster_to_type <- tibble(
  ClusterID = 0:22,
  CellType = case_when(
    ClusterID %in% c(3,6,9,13,14,15,16,20,22) ~ "Myeloid",
    ClusterID %in% c(5,11,12,17,19) ~ "Epithelial",
    ClusterID %in% c(0,1,21) ~ "T/NK",
    ClusterID == 8 ~ "Fibroblast",
    ClusterID == 18 ~ "Proliferating",
    ClusterID == 10 ~ "Plasma",
    ClusterID == 2 ~ "B",
    ClusterID == 7 ~ "Endothelial",
    ClusterID == 4 ~ "Mast",
    TRUE ~ "Unknown"
  )
)

# Attach cell type annotation to metadata
scRNA <- scRNA %>%
  AddMetaData(metadata = left_join(
    tibble(barcode = colnames(.), cluster = pull(.[[cluster_res]])),
    cluster_to_type, by = c("cluster" = "ClusterID")
  ) %>% column_to_rownames("barcode") %>% select(CellType), col.name = "celltype")

# Check cell type distribution
table(scRNA$celltype)

# UMAP colored by cell type
umap_celltype <- DimPlot(scRNA, reduction = "umap", group.by = "celltype",
                         cols = my36colors, label = TRUE, pt.size = 0.8) +
  NoAxes() +
  geom_line(data = axis_info$axes, aes(x = x, y = y, group = group),
            arrow = arrow(length = unit(0.1, "inches"), ends = "last")) +
  geom_text(data = axis_info$label, aes(x = x, y = y, label = lab, angle = angle), fontface = "italic") +
  theme(aspect.ratio = 1, plot.title = element_blank())
ggsave("umap_by_celltype.pdf", umap_celltype, width = 9, height = 7)

# Combine UMAPs
library(patchwork)
combined_umap <- umap_cluster + umap_sample + umap_celltype
ggsave("combine_umap.pdf", combined_umap, width = 15, height = 7)

# Save final annotated Seurat object
saveRDS(scRNA, "sce_celltype.rds")

# Return to parent directory
setwd("../")
