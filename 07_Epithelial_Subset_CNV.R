# Clear environment and load packages
rm(list = ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(clustree)
source("scRNA_scripts/lib.R")
source("scRNA_scripts/mycolors.R")

# Create working directory
dir.create("7-epi", showWarnings = FALSE)
setwd("7-epi/")
set.seed(12345)

# Load annotated Seurat object
full_obj <- readRDS("../3-Celltype/sce_celltype.rds")

# Subset epithelial cells
epi_obj <- subset(full_obj, celltype == "Epithelial")
epi_obj <- JoinLayers(epi_obj)
table(epi_obj$orig.ident)

# Normalization and feature selection
epi_obj <- NormalizeData(epi_obj)
epi_obj <- FindVariableFeatures(epi_obj, selection.method = "vst", nfeatures = 2000)
epi_obj <- ScaleData(epi_obj)
epi_obj <- RunPCA(epi_obj)
ElbowPlot(epi_obj)

# Clustering
epi_obj <- FindNeighbors(epi_obj, dims = 1:15)
epi_obj <- FindClusters(epi_obj, resolution = 0.1)
epi_obj <- RunUMAP(epi_obj, dims = 1:5)
epi_obj <- RunTSNE(epi_obj, dims = 1:15)

# Marker inspection (from reference paper)
marker_genes <- toupper(c("TM4SF1", "CRABP2", "UBE2C", "TOP2A", "MKI67",
                          "CAV1", "CLDN18", "CAPS", "SCGB1A1"))
DotPlot(epi_obj, features = marker_genes) + coord_flip()
ggsave("check_markers.pdf", width = 7, height = 9)

# Annotated UMAP
source("../scRNA_scripts/Bottom_left_axis.R")
axis_info <- left_axes(epi_obj)
DimPlot(epi_obj, reduction = "umap", label = TRUE, cols = my36colors) +
  geom_line(data = axis_info$axes, aes(x, y, group = group),
            arrow = arrow(length = unit(0.1, "inches"), ends = "last")) +
  geom_text(data = axis_info$label, aes(x, y, angle = angle, label = lab), fontface = "italic") +
  NoAxes() +
  ggsave("RNA_snn_res.0.1_umap.pdf", width = 9, height = 7)

# Violin plot of marker expression
VlnPlot(epi_obj, features = marker_genes, pt.size = 0, ncol = 4, cols = mycolors)

# ----------------------------
# inferCNV CNV analysis
# ----------------------------
dir.create("CNV", showWarnings = FALSE)
setwd("CNV")

# Assign cluster names
epi_obj$cnv_cluster <- paste0("C", Idents(epi_obj))
Idents(epi_obj) <- epi_obj$cnv_cluster

# Downsample epithelial cells
cnv_epi <- subset(epi_obj, downsample = 200)

# Reference cells: endothelial and T
endo_ref <- subset(full_obj, celltype == "Endothelial")
t_ref <- subset(full_obj, celltype == "T&NK", downsample = 500)

# Combine raw count matrix
ref1 <- as.data.frame(endo_ref[["RNA"]]$counts)
ref2 <- as.data.frame(t_ref[["RNA"]]$counts)
target <- as.data.frame(cnv_epi[["RNA"]]$counts)
gene_ids <- intersect(rownames(target), rownames(ref1))
expr_mat <- cbind(target[gene_ids, ], ref1[gene_ids, ], ref2[gene_ids, ])

# Write inferCNV files
anno <- data.frame(
  v1 = colnames(expr_mat),
  v2 = c(cnv_epi$cnv_cluster, rep("ref-1", ncol(ref1)), rep("ref-2", ncol(ref2)))
)
write.table(anno, "groupFiles.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

library(AnnoProbe)
gene_info <- annoGene(rownames(expr_mat), "SYMBOL", "human") %>%
  arrange(chr, start) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(SYMBOL, chr, start, end)
expr_mat <- expr_mat[rownames(expr_mat) %in% gene_info$SYMBOL, ]
expr_mat <- expr_mat[match(gene_info$SYMBOL, rownames(expr_mat)), ]

write.table(expr_mat, "expFile.txt", sep = "\t", quote = FALSE)
write.table(gene_info, "geneFile.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Run inferCNV
library(infercnv)
cnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = "expFile.txt",
  annotations_file = "groupFiles.txt",
  delim = "\t",
  gene_order_file = "geneFile.txt",
  ref_group_names = c("ref-1", "ref-2")
)

cnv_result <- run(cnv_obj,
                  cutoff = 0.1,
                  out_dir = "infercnv_output",
                  cluster_by_groups = TRUE,
                  hclust_method = "ward.D2",
                  plot_steps = TRUE)
save(expr_mat, gene_info, anno, cnv_epi, file = "infercnv.Rdata")

# Load inferCNV output and compute CNV scores
cnv_final <- readRDS("infercnv_output/run.final.infercnv_obj")
expr_data <- cnv_final@expr.data

ref_data <- expr_data[, unlist(cnv_final@reference_grouped_cell_indices)]
ref_mean <- rowMeans(ref_data)
ref_sd <- apply(ref_data, 1, sd)
cutoff_low <- mean(ref_mean) - 2 * mean(ref_sd)
cutoff_high <- mean(ref_mean) + 2 * mean(ref_sd)
cnv_range <- cutoff_high - cutoff_low

thresholds <- c(
  a1 = cutoff_low - 2 * cnv_range,
  a2 = cutoff_low - 1 * cnv_range,
  a3 = cutoff_high + 1 * cnv_range,
  a4 = cutoff_high + 2 * cnv_range
)

# Discretize CNV values
label_mat <- expr_data
label_mat[expr_data > 0 & expr_data < thresholds["a2"]] <- "A"
label_mat[expr_data >= thresholds["a2"] & expr_data < cutoff_low] <- "B"
label_mat[expr_data >= cutoff_low & expr_data < cutoff_high] <- "C"
label_mat[expr_data >= cutoff_high & expr_data <= thresholds["a3"]] <- "D"
label_mat[expr_data > thresholds["a3"] & expr_data <= thresholds["a4"]] <- "E"
label_mat[expr_data > thresholds["a4"]] <- "F"

score_mat <- matrix(0, nrow = nrow(label_mat), ncol = ncol(label_mat))
score_mat[label_mat == "A"] <- 2
score_mat[label_mat == "B"] <- 1
score_mat[label_mat == "D"] <- 1
score_mat[label_mat == "E"] <- 2
score_mat[label_mat == "F"] <- 2

cell_scores <- data.frame(
  cnv_score = colSums(score_mat),
  group = gsub("ref-1|spike-1", "Endo", anno$v2)
)
cell_scores$group <- gsub("ref-2|spike-2", "T", cell_scores$group)

ggplot(cell_scores, aes(x = group, y = cnv_score, fill = group)) +
  geom_boxplot() +
  theme_base() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = NULL, y = "CNV Score")

setwd("..")

# ----------------------------
# Manual subcluster annotation
# ----------------------------

cluster_map <- data.frame(
  ClusterID = 0:8,
  celltype = c(
    "Cancer(TM4SF1+)", "AT1", "Ciliated cells", "AT2", "Cancer(UBE2C+)",
    "Cancer(TM4SF1+)", "Clara-like cancer", "Clara cells", "Cancer(CRABP2+)"
  )
)

epi_obj$celltype <- "NA"
for (i in 1:nrow(cluster_map)) {
  epi_obj$celltype[which(epi_obj$RNA_snn_res.0.1 == cluster_map$ClusterID[i])] <- cluster_map$celltype[i]
}

# UMAP colored by new annotation
axis_info <- left_axes(epi_obj)
p <- DimPlot(epi_obj, reduction = "umap", group.by = "celltype", cols = my36colors, label = TRUE, pt.size = 0.8) +
  geom_line(data = axis_info$axes, aes(x, y, group = group),
            arrow = arrow(length = unit(0.1, "inches"), ends = "last")) +
  geom_text(data = axis_info$label, aes(x, y, angle = angle, label = lab), fontface = "italic") +
  NoAxes() + theme(aspect.ratio = 1)
ggsave("umap_by_celltype.pdf", plot = p, width = 9, height = 7)

# Save object
saveRDS(epi_obj, "epi_sce_celltype.rds")

# Identify DEGs per subcluster
Idents(epi_obj) <- epi_obj$celltype
if (!file.exists("epimarkers.csv")) {
  markers <- FindAllMarkers(epi_obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  write.csv(markers, "epimarkers.csv")
}

# ----------------------------
# Pie chart of proportions
# ----------------------------
source("prop_pie.R")
prop_pie(epi_obj, seurat_by = "sample", pheno_by = "celltype", palette = "nejm")

# Relative proportions + trend lines
library(ggsci)
library(ggrepel)
plot_data <- epi_obj@meta.data %>%
  dplyr::count(sample, celltype) %>%
  group_by(celltype) %>%
  mutate(Proportion = n / sum(n) * 100)

ggplot(plot_data, aes(x = sample, y = Proportion, color = sample)) +
  geom_line(aes(group = 1), size = 1.5, color = "#64A36C") +
  geom_point(size = 5) +
  facet_wrap(~ celltype) +
  theme_bw() +
  scale_color_nejm()
ggsave("pie_prob.pdf", width = 12, height = 7)

# Return to parent directory
setwd("..")
