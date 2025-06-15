# -----------------------------------------
# SCENIC: Regulatory Network Inference from scRNA-seq
# -----------------------------------------

rm(list = ls())
set.seed(123)

# Create and set working directory
dir.create("10-SCENIC", showWarnings = FALSE)
setwd("10-SCENIC")
dir.create("int")

# Load required packages
packages <- c("GenomicFeatures","AUCell","RcisTarget","GENIE3","zoo","mixtools",
              "rbokeh","DT","NMF","R2HTML","Rtsne","doMC","doRNG","SCENIC")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# Optional GitHub packages
if (!requireNamespace("SCopeLoomR", quietly = TRUE)) {
  devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
}
library(SCopeLoomR)

# Load expression data
library(Seurat)
library(tidyverse)
library(patchwork)

seu <- readRDS("../7-epi/epi_sce_celltype.rds")
table(seu$celltype)

# -------------------------------
# Subset target cells (Cancer)
# -------------------------------
target_cells <- c("Cancer(CRABP2+)", "Cancer(TM4SF1+)", "Cancer(UBE2C+)", "Clara-like cancer")
seu_sub <- subset(seu, subset = celltype %in% target_cells)
Idents(seu_sub) <- seu_sub$celltype

# Downsample to 50 cells per group to reduce computational cost
cell_ids <- unlist(lapply(levels(Idents(seu_sub)), function(cl) {
  sample(WhichCells(seu_sub, idents = cl), 50)
}))
seu_sub <- subset(seu_sub, cells = cell_ids)
saveRDS(seu_sub, "scRNAsub.rds")

# Save metadata
meta_info <- seu_sub@meta.data[, c("sample", "celltype")]
saveRDS(meta_info, file = "int/cellInfo.Rds")

# Extract raw expression matrix
expr_mat <- as.matrix(GetAssayData(seu_sub, slot = "counts"))
dim(expr_mat)

# -------------------------------
# Set database location & SCENIC options
# -------------------------------
db_dir <- "/home/data/t050453/paper_plot_redraw/scRNA/GSE189357-LUAD-scRNA-ST/10-SCENIC"
db_files <- c("hg19-500bp-upstream-7species.mc9nr.feather",
              "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(db_files) <- c("500bp", "10kb")

scenic_opts <- initializeScenic(org = "hgnc",
                                nCores = 60,
                                dbDir = db_dir,
                                dbs = db_files,
                                datasetTitle = "SCENIC")
saveRDS(scenic_opts, "int/scenicOptions.rds")

# -------------------------------
# Step 1: Gene Filtering
# -------------------------------
genes_pass <- geneFiltering(expr_mat, scenic_opts,
                            minCountsPerGene = 0.03 * ncol(expr_mat),
                            minSamples = 0.01 * ncol(expr_mat))
expr_mat_filtered <- expr_mat[genes_pass, ]

# -------------------------------
# Step 2: TF Network Inference
# -------------------------------
runCorrelation(expr_mat_filtered, scenic_opts)
runGenie3(log2(expr_mat_filtered + 1), scenic_opts, nParts = 20)

# -------------------------------
# Step 3: Regulon Identification
# -------------------------------
runSCENIC_1_coexNetwork2modules(scenic_opts)
runSCENIC_2_createRegulons(scenic_opts)

# -------------------------------
# Step 4: Regulon Activity Scoring
# -------------------------------
expr_all <- as.matrix(GetAssayData(seu, slot = "counts"))
expr_all <- log2(expr_all + 1)
saveRDS(expr_all, "exprMat_all.rds")

runSCENIC_3_scoreCells(scenic_opts, exprMat = expr_all)
runSCENIC_4_aucell_binarize(scenic_opts, exprMat = expr_all)

# -------------------------------
# Step 5: Regulon Analysis
# -------------------------------
regulons <- loadInt(scenic_opts, "aucell_regulons")
regulon_auc <- loadInt(scenic_opts, "aucell_regulonAUC")
regulon_auc <- regulon_auc[onlyNonDuplicatedExtended(rownames(regulon_auc)), ]

# Aggregated regulon activity by cell type
meta_info <- readRDS("int/cellInfo.Rds")
activity_matrix <- sapply(split(rownames(meta_info), meta_info$celltype), function(cells) {
  rowMeans(getAUC(regulon_auc)[, cells, drop = FALSE])
})

# Select top 20 regulons for heatmap visualization
activity_df <- as.data.frame(activity_matrix)
top_20 <- activity_df[sample(rownames(activity_df), 20), ]
scaled_mat <- t(scale(t(top_20)))

library(ComplexHeatmap)
dev.new()
Heatmap(scaled_mat, name = "Regulon activity")
dev.off()

# -------------------------------
# Step 6: Regulon Specificity Score (RSS)
# -------------------------------
celltypes <- meta_info[colnames(regulon_auc), "celltype"]
rss <- calcRSS(getAUC(regulon_auc), celltypes)
rss_plot <- plotRSS(rss)
print(rss_plot$plot)

# -------------------------------
# Step 7: Inspect specific TF regulon (e.g., Stat6)
# -------------------------------
target_info <- loadInt(scenic_opts, "regulonTargetsInfo")
stat6_targets <- target_info[TF == "Stat6" & highConfAnnot == TRUE]
head(stat6_targets)

# Save code history
savehistory("scenic_code.txt")
