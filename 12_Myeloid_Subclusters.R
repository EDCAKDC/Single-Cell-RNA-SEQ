# Myeloid and T/NK Cell Subtype Analysis and Proportional Comparison

```r
rm(list=ls())
options(stringsAsFactors = FALSE)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
source('scRNA_scripts/lib.R')
source('scRNA_scripts/mycolors.R')

dir.create("12-T")
setwd("12-T/")
set.seed(12345)

sce.all <- readRDS("../3-Celltype/sce_celltype.rds")

# Select Myeloid and Mast cells
sce1 <- sce.all[, sce.all$celltype %in% c('Myeloid', 'Mast')]
LayerData(sce1, assay = "RNA", layer = "counts")
sce1 <- JoinLayers(sce1)

sce <- sce1
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(sce))

ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.1)
sce <- RunTSNE(sce, dims = 1:15)
sce <- RunUMAP(sce, dims = 1:5)

DimPlot(sce, label = TRUE, cols = mycolors)

# Define known marker genes
marker_genes <- c('CCL3L1', 'FABP4', 'CPA3', 'CST3', 'KIT', 'TPSAB1', 'TPSB2',
                  'G0S2', 'S100A9', 'S100A8', 'CXCL8', 'S100B', 'CD68', 'CD163',
                  'CD14', 'MKI67', 'PCNA', 'VCAN', 'FCN1', 'CD300E', 'TXN')
marker_genes <- toupper(marker_genes)

DotPlot(sce, features = unique(marker_genes), assay = 'RNA') + coord_flip()
ggsave('myeloid_check_markers.pdf', height = 10, width = 7)

# Add UMAP axis annotation
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(sce)
umap <- DimPlot(sce, reduction = "umap", cols = my36colors, pt.size = 0.8,
                group.by = "RNA_snn_res.0.1", label = TRUE, label.box = TRUE) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = result$axes, aes(x = x, y = y, group = group),
            arrow = arrow(length = unit(0.1, "inches"), ends = "last")) +
  geom_text(data = result$label, aes(x = x, y = y, angle = angle, label = lab),
            fontface = 'italic')
ggsave('myeloid_RNA_snn_res.0.1_umap.pdf', width = 9, height = 7)

# Annotate clusters manually
annotation <- data.frame(ClusterID = 0:9,
                         celltype = c('M2','Mast','M1','proliferating','Granulocyte',
                                      'S100B_DC','Mast','M1','TXN_DC','TXN_DC'))

sce$celltype <- annotation$celltype[match(sce$RNA_snn_res.0.1, annotation$ClusterID)]

# UMAP colored by annotated cell types
result <- left_axes(sce)
DimPlot(sce, reduction = "umap", cols = my36colors, pt.size = 0.5,
        group.by = "celltype", label = TRUE) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = result$axes, aes(x = x, y = y, group = group),
            arrow = arrow(length = unit(0.1, "inches"), ends = "last")) +
  geom_text(data = result$label, aes(x = x, y = y, angle = angle, label = lab),
            fontface = 'italic')
ggsave('myeloid_umap_by_celltype.pdf', width = 9, height = 7)

saveRDS(sce, "myeloidsce_celltype.rds")

# Merge metadata with T/NK data for compositional comparison
myeloidsce <- readRDS("myeloidsce_celltype.rds")
myeloidsce$group <- 'myeloid'
myeloidmeta <- myeloidsce@meta.data

TNKsce <- readRDS("Tsce_celltype.rds")
TNKsce$group <- 'T&NK'
TNKmeta <- TNKsce@meta.data

meta <- rbind(myeloidmeta, TNKmeta)
tb <- table(meta$sample, meta$celltype, meta$group)
balloonplot(tb)
bar_data <- as.data.frame(tb)

bar_per <- bar_data %>%
  group_by(Var2) %>%
  mutate(percent = Freq / sum(Freq))

bar_per <- bar_per[bar_per$Freq != 0, ]
bar_per$Var1 <- factor(bar_per$Var1, levels = c('AIS', 'MIA', 'IAC'))
bar_per$Var2 <- factor(bar_per$Var2, levels = c("TXN_DC", "S100B_DC", "M1", "M2",
                                                "Mast", "proliferating", "Granulocyte",
                                                "CCL4L2_CD8", "GZMK_CD8", "memory_CD4",
                                                "GNLY_NK", "Treg", "Fibro-like T",
                                                "exhausted_CD4", "FGFBP2_NK"))

# Generate faceted barplot
p1 <- ggplot(bar_per, aes(x = percent, y = Var2)) +
  geom_bar(aes(fill = Var1), stat = "identity") +
  coord_flip() +
  labs(y = " ", fill = NULL, x = 'Relative proportion(%)') +
  scale_fill_manual(values = c("#3176B7", "#F78000", "#3FA116", "#CE2820", "#9265C1",
                               "#885649", "#DD76C5", "#BBBE00", "#41BED1")) +
  facet_grid(~Var3, scales = "free_x", space = "free_x") +
  theme_pubr(base_size = 10) +
  theme_few() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        legend.key.size = unit(3, "pt"))

ggsave(filename = "prop.pdf", width = 12, height = 7)
setwd('../')
```
