# Single-Cell RNA-seq Analysis Workflows

A modular R-based collection of workflows for **single-cell transcriptomic analysis in cancer**, using a human lung adenocarcinoma example dataset.

This repository is maintained as a reproducible methods reference covering common steps from quality control to cell-state interpretation and downstream translational analysis.

## Analysis Modules

```text
LungTumor-scRNAseq/
├── 01_QC_02_Integration_Harmony.R
├── 03_Celltype_Annotation.R
├── 04_Feature_Heatmap_Plot.R
├── 05_Cluster_Composition.R
├── 06_DEG_Bulk_Comparison.R
├── 07_Epithelial_Subset_CNV.R
├── 08_Functional_GO_GSVA.R
├── 09_Trajectory_Monocle.R
├── 10_Regulon_SCENIC.R
├── 11_TCGA_Prognosis_Validation.R
├── 12_Myeloid_Subclusters.R
├── 13_Tcell_Subclusters.R
└── 14_CellCell_Communication.R
```

Utility functions are organized under `scRNA_scripts/`.

## Covered Methods

- quality control and filtering
- Seurat preprocessing
- Harmony batch correction
- PCA, UMAP, and t-SNE visualization
- marker-based cell-type annotation
- differential expression analysis
- cluster-composition analysis
- epithelial-cell CNV inference
- GO / KEGG / GSVA pathway analysis
- Monocle pseudotime analysis
- SCENIC regulon inference
- TCGA survival validation
- myeloid and T-cell subclustering
- CellChat cell-cell communication analysis

## Main Tools

```text
Seurat
Harmony
Monocle
inferCNV
SCENIC
GSVA
clusterProfiler
CellChat
ggplot2
ggpubr
cowplot
patchwork
pheatmap
```

## Notes

The lung adenocarcinoma dataset is used as an example implementation. The repository is intended as a methods-oriented workflow, so sample metadata, thresholds, annotations, and comparison groups should be adapted before use with a new dataset.

## Research Context

The methods represented here are directly relevant to tumor microenvironment analysis, immune-cell state characterization, cancer genomics, and translational single-cell studies.
