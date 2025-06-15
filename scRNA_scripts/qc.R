basic_qc <- function(input_sce) {
  # Calculate mitochondrial gene percentage
  mito_genes <- rownames(input_sce)[grep("^MT-", rownames(input_sce), ignore.case = TRUE)]
  print(mito_genes)
  input_sce <- PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
  fivenum(input_sce@meta.data$percent_mito)
  
  # Calculate ribosomal gene percentage
  ribo_genes <- rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce), ignore.case = TRUE)]
  print(ribo_genes)
  input_sce <- PercentageFeatureSet(input_sce, features = ribo_genes, col.name = "percent_ribo")
  fivenum(input_sce@meta.data$percent_ribo)
  
  # Calculate hemoglobin gene percentage
  Hb_genes <- rownames(input_sce)[grep("^Hb[^p]", rownames(input_sce), ignore.case = TRUE)]
  print(Hb_genes)
  input_sce <- PercentageFeatureSet(input_sce, features = Hb_genes, col.name = "percent_hb")
  fivenum(input_sce@meta.data$percent_hb)
  
  # Visualize metrics
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1 <- VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + NoLegend()
  w <- length(unique(input_sce$orig.ident)) / 3 + 5
  ggsave("Vlnplot1.pdf", p1, width = w, height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2 <- VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims = TRUE) +
    scale_y_continuous(breaks = seq(0, 100, 5)) + NoLegend()
  w <- length(unique(input_sce$orig.ident)) / 2 + 5
  ggsave("Vlnplot2.pdf", p2, width = w, height = 5)
  
  # Feature scatter plot
  p3 <- FeatureScatter(input_sce, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
  ggsave("Scatterplot.pdf", p3)
  
  # Optional: filter low-quality cells/genes
  if (FALSE) {
    selected_c <- WhichCells(input_sce, expression = nFeature_RNA > 500)
    selected_f <- rownames(input_sce)[Matrix::rowSums(input_sce@assays$RNA$counts > 0) > 3]
    input_sce.filt <- subset(input_sce, features = selected_f, cells = selected_c)
    dim(input_sce)
    dim(input_sce.filt)
  }
  
  input_sce.filt <- input_sce
  
  # Boxplot of top 50 highly expressed genes
  C <- subset(input_sce.filt, downsample = 100)@assays$RNA$counts
  C <- Matrix::t(Matrix::t(C) / Matrix::colSums(C)) * 100
  most_expressed <- order(apply(C, 1, median), decreasing = TRUE)[50:1]
  
  pdf("TOP50_most_expressed_gene.pdf", width = 14)
  boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
          cex = 0.1, las = 1,
          xlab = "% total count per cell",
          col = (scales::hue_pal())(50)[50:1],
          horizontal = TRUE)
  dev.off()
  rm(C)
  
  # Filter by mitochondrial/ribosomal/hemoglobin gene expression
  selected_mito <- WhichCells(input_sce.filt, expression = percent_mito < 25)
  selected_ribo <- WhichCells(input_sce.filt, expression = percent_ribo > 3)
  selected_hb <- WhichCells(input_sce.filt, expression = percent_hb < 1)
  
  input_sce.filt <- subset(input_sce.filt, cells = selected_mito)
  input_sce.filt <- subset(input_sce.filt, cells = selected_hb)
  dim(input_sce.filt)
  table(input_sce.filt$orig.ident)
  
  # Post-filtering visualization
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1_filtered <- VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + NoLegend()
  w <- length(unique(input_sce.filt$orig.ident)) / 3 + 5
  ggsave("Vlnplot1_filtered.pdf", p1_filtered, width = w, height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2_filtered <- VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + NoLegend()
  w <- length(unique(input_sce.filt$orig.ident)) / 2 + 5
  ggsave("Vlnplot2_filtered.pdf", p2_filtered, width = w, height = 5)
  
  return(input_sce.filt)
}
