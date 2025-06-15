# Clear environment
rm(list = ls())

# Load helper functions and set up working directory
source("scRNA_scripts/lib.R")
dir.create("6-DEG", showWarnings = FALSE)
setwd("6-DEG")

# Load Seurat object with annotations
seurat_obj <- readRDS("../3-Celltype/sce_celltype.rds")
library(tidyverse)

# Set sample group as identity for DEG
Idents(seurat_obj) <- seurat_obj$sample
sample_levels <- c("C", "B", "A")

# Run DEG analysis between all pairwise sample combinations
if (!file.exists("markers_list.Rdata")) {
  deg_results <- list()
  
  combn(sample_levels, 2, simplify = FALSE) %>%
    walk(~{
      res <- FindMarkers(seurat_obj, group.by = "sample",
                         ident.1 = .x[1], ident.2 = .x[2],
                         logfc.threshold = 0.1)
      key <- paste0(.x[1], "_vs_", .x[2])
      deg_results[[key]] <- res
    })
  
  save(deg_results, file = "markers_list.Rdata")
} else {
  load("markers_list.Rdata")
}

# Filter significant genes
sig_deg <- map(deg_results, ~ filter(.x, p_val_adj < 0.01))

# DEG summary stats (up/down-regulated counts)
deg_stats <- map_df(sig_deg, function(df) {
  tibble(
    Up = sum(df$avg_log2FC > 1),
    Down = sum(df$avg_log2FC < -1),
    Total = n()
  )
}, .id = "Comparison")

print(deg_stats)

# Annotate up/down for each contrast
annotate_group <- function(df) {
  df %>%
    mutate(group = case_when(
      avg_log2FC > 1 & p_val_adj < 0.01 ~ "up",
      avg_log2FC < -1 & p_val_adj < 0.01 ~ "down",
      TRUE ~ NA_character_
    )) %>%
    drop_na(group)
}

deg_C_vs_A <- annotate_group(sig_deg[[1]])
deg_C_vs_B <- annotate_group(sig_deg[[2]])
deg_A_vs_B <- annotate_group(sig_deg[[3]])

# Get upregulated gene sets
genes_up_C_A <- rownames(deg_C_vs_A[deg_C_vs_A$group == "up", ])
genes_up_C_B <- rownames(deg_C_vs_B[deg_C_vs_B$group == "up", ])
genes_up_A_B <- rownames(deg_A_vs_B[deg_A_vs_B$group == "up", ])

# Venn diagram (upregulated)
library(VennDetail)
venn_up <- venndetail(list(
  C_vs_B_up = genes_up_C_B,
  C_vs_A_up = genes_up_C_A,
  A_vs_B_up = genes_up_A_B
))

plot(venn_up, type = "vennpie")
dplot(venn_up, order = TRUE, textsize = 4)
plot(venn_up, type = "upset")

# Repeat for downregulated genes
genes_down_C_A <- rownames(deg_C_vs_A[deg_C_vs_A$group == "down", ])
genes_down_C_B <- rownames(deg_C_vs_B[deg_C_vs_B$group == "down", ])
genes_down_A_B <- rownames(deg_A_vs_B[deg_A_vs_B$group == "down", ])

venn_down <- venndetail(list(
  C_vs_B_down = genes_down_C_B,
  C_vs_A_down = genes_down_C_A,
  A_vs_B_down = genes_down_A_B
))

plot(venn_down, type = "vennpie")
dplot(venn_down, order = TRUE, textsize = 4)
plot(venn_down, type = "upset")

# PCA on average expression of DEGs across sample-patient groups
seurat_obj$sample_patient <- paste(seurat_obj$sample, seurat_obj$patient, sep = "_")
avg_expr <- AverageExpression(seurat_obj, group.by = "sample_patient")$RNA
all_deg_genes <- unique(c(
  rownames(deg_C_vs_A),
  rownames(deg_C_vs_B),
  rownames(deg_A_vs_B)
))
expr_mat <- avg_expr[rownames(avg_expr) %in% all_deg_genes, ]

pca <- prcomp(t(expr_mat), scale. = TRUE)
pca_df <- as.data.frame(pca$x)

library(scatterplot3d)
group_colors <- ifelse(str_detect(rownames(pca_df), "B"), "purple",
                       ifelse(str_detect(rownames(pca_df), "C"), "orange", "green"))
scatterplot3d(pca_df[, 1:3], color = group_colors, pch = 16, angle = 30)
legend("topleft", legend = c("B", "C", "A"),
       fill = c("purple", "orange", "green"), box.col = NA, cex = 0.8)

# Compare DEGs across two contrasts
df1 <- deg_C_vs_A
df2 <- deg_A_vs_B
common_genes <- intersect(rownames(df1), rownames(df2))
genes_to_plot <- sample(common_genes, 200)

# Load marker genes and filter
marker_df <- read.csv("../4-plot/sce.markers.csv", row.names = 1)
marker_df <- marker_df[marker_df$gene %in% genes_to_plot, ]
marker_df <- marker_df[!duplicated(marker_df$gene), ]

# Create dataframe for correlation plot
df_corr <- tibble(
  gene = marker_df$gene,
  celltype = marker_df$cluster,
  logFC_CA = df1[marker_df$gene, "avg_log2FC"],
  logFC_AB = df2[marker_df$gene, "avg_log2FC"]
)

library(ggrepel)
ggplot(df_corr, aes(x = logFC_AB, y = logFC_CA, color = celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = gene), size = 2, color = "black", fontface = "italic") +
  xlim(-2, 3) + ylim(-3, 2) +
  labs(x = "Log2FC B vs A", y = "Log2FC C vs B") +
  theme_few()

# Summarize top DEGs across comparisons
filtered_list <- map(deg_results, ~ filter(.x, abs(avg_log2FC) > 0.5))
filtered_list <- map2(filtered_list, names(filtered_list), ~ mutate(.x, group = .y, gene = rownames(.x)))
all_df <- bind_rows(filtered_list)
all_df$label <- ifelse(all_df$p_val_adj < 0.01, "adj p < 0.01", "adj p ≥ 0.01")

top_genes <- all_df %>%
  group_by(group) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 15)

all_df$size <- ifelse(all_df$gene %in% top_genes$gene, 2, 1)

# Background grid bars
bg1 <- tibble(x = 1:3, y = c(10, 10, 11))
bg2 <- tibble(x = 1:3, y = c(-7.5, -7, -5))

# Label positions
labels_df <- tibble(
  x = 1:3,
  y = 0,
  label = c("B_vs_A", "C_vs_A", "C_vs_B"),
  color = c("#00A0877F", "#3C54887F", "#F39B7F7F")
)

# Volcano dot-like plot
ggplot() +
  geom_col(data = bg1, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = bg2, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_jitter(data = filter(all_df, size == 1), aes(x = group, y = avg_log2FC, color = label), size = 0.85, width = 0.4) +
  geom_jitter(data = filter(all_df, size == 2), aes(x = group, y = avg_log2FC, color = label), size = 1.2, width = 0.4) +
  geom_tile(data = labels_df, aes(x = x, y = y), height = 0.8, fill = labels_df$color, color = "black") +
  geom_text(data = labels_df, aes(x = x, y = y, label = label), color = "white", size = 4) +
  labs(x = "Comparison", y = "Average log2FC") +
  scale_color_manual(values = c("red", "black")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.line.y = element_line(size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top"
  )

# Return to upper directory
setwd("..")
