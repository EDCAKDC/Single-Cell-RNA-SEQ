# Define marker gene lists for different cell types and datasets
GASTRIC_CANCER_MARKERS <- c("PTPRC", "MUC2", "ITLN1", "FABP1", "APOA1", "CEACAM5", "CEACAM6",
                            "EPCAM", "KRT18", "MUC1", "MUC6", "TFF2", "PGA4", "PGA3",
                            "MUC5AC", "TFF1", "CHGA", "CHGB")

MYO_MARKERS <- c("Krt17", "Krt14", "Krt5", "Acta2", "Myl9", "Mylk", "Myh11")
LUM_MARKERS <- c("Krt19", "Krt18", "Krt8")
HS_MARKERS <- c("Prlr", "Cited1", "Pgr", "Prom1", "Esr1")
AV_MARKERS <- c("Mfge8", "Trf", "Csn3", "Wfdc18", "Elf5", "Ltf")
LP_MARKERS <- c("Kit", "Aldh1a3", "Cd14")
FIB_MARKERS <- c("Col1a1", "Col1a2", "Col3a1", "Fn1")

GSE150580_BREAST_CANCER_LIST <- list(Myo = MYO_MARKERS, Lum = LUM_MARKERS, Hs = HS_MARKERS, 
                                     AV = AV_MARKERS, Lp = LP_MARKERS, Fib = FIB_MARKERS)

# Myeloid marker definitions from SCP1661 and additional sources
SCP1661_MYLOID_LIST <- list(
  macrophages = c("Adgre1", "Cd14", "Fcgr3"),
  cDCs = c("Xcr1", "Flt3", "Ccr7"),
  pDCs = c("Siglech", "Clec10a", "Clec12a"),
  monocytes = c("Ly6c2", "Spn"),
  neutrophils = c("Csf3r", "S100a8", "Cxcl3")
)

# Further immune cell marker lists...
# (continued from your long code block, split into smaller chunks if necessary)

# Main logic for checking marker expression and generating DotPlots
if (species == "human") {
  marker_format <- str_to_upper
} else if (species == "mouse") {
  marker_format <- str_to_title
} else {
  stop("Only 'human' or 'mouse' species are supported.")
}

plot_markers <- function(marker_group, obj, group_name, coord_flip_plot = TRUE) {
  formatted_genes <- marker_format(marker_group)
  dotplot <- DotPlot(obj, features = formatted_genes) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (coord_flip_plot) {
    dotplot <- dotplot + coord_flip()
  }
  return(dotplot)
}

# Example usage for all groups
for (marker_name in c("GASTRIC_CANCER_MARKERS")) {
  marker_genes <- get(marker_name)
  p <- plot_markers(marker_genes, sce.all.int, marker_name)
  ggsave(paste0("check_for_", marker_name, ".pdf"), p)
}

# Example: Composite plot with UMAP
final_plot <- plot_markers(last_markers_to_check, sce.all.int, "last") + p_umap
ggsave("last_markers_and_umap.pdf", final_plot)

# Quality control visualization
if ("percent_mito" %in% colnames(sce.all.int@meta.data)) {
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
  p1 <- VlnPlot(sce.all.int, features = feats[1:2], pt.size = 0, ncol = 2) + NoLegend()
  ggsave("qc-Vlnplot1.pdf", p1)
  p2 <- VlnPlot(sce.all.int, features = feats[3:5], pt.size = 0, ncol = 3) +
    scale_y_continuous(breaks = seq(0, 100, 5)) + NoLegend()
  ggsave("qc-Vlnplot2.pdf", p2)
  p3 <- FeatureScatter(sce.all.int, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5)
  ggsave("qc-Scatterplot.pdf", p3)
}

# Cluster-specific marker analysis using COSG
if (requireNamespace("cosg", quietly = TRUE)) {
  marker_cosg <- cosg::cosg(sce.all.int, groups = 'all', assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 100)
  save(marker_cosg, file = "qc-_marker_cosg.Rdata")
  top10 <- unique(unlist(apply(marker_cosg$names, 2, head, 10)))
  top3 <- unique(unlist(apply(marker_cosg$names, 2, head, 3)))
  
  for (top_genes in list(top10 = top10, top3 = top3)) {
    scaled_obj <- ScaleData(sce.all.int, features = top_genes)
    heatmap <- DoHeatmap(scaled_obj, features = top_genes, size = 3)
    dotplot <- DotPlot(sce.all.int, features = top_genes, assay = 'RNA') + coord_flip()
    suffix <- ifelse(length(top_genes) > 5, "top10", "top3")
    ggsave(paste0("qc-DoHeatmap_check_", suffix, ".pdf"), heatmap)
    ggsave(paste0("qc-DotPlot_check_", suffix, ".pdf"), dotplot)
  }
}
