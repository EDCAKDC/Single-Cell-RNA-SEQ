# Single-Cell-RNA-SEQ
👋 Hi there, I'm EDCAKDC
🎓 I'm currently learning bioinformatics and sharing my learning journey.

🔬 Featured Project: Human Lung scRNA-seq Analysis (LUAD)
This repository contains the LUAD single-cell RNA-seq analysis pipeline I had learned before.

🧾 Sample Group Info
Actually, the specific cell types are not the focus here — the main purpose is to learn the underlying principles of analysis. 
You can replace the sample groups with your own data as needed.
Each group contains multiple biological replicates.

## 📁 Project Structure

```text
LungTumor-scRNAseq/
├── 01_QC_02_Integration_Harmony.R       # QC, filtering, and Harmony batch correction
├── 03_Celltype_Annotation.R             # Manual annotation with known markers
├── 04_Feature_Heatmap_Plot.R            # Visualization: DotPlot, FeaturePlot, Heatmap
├── 05_Cluster_Composition.R             # Cell proportion analysis across tumor stages
├── 06_DEG_Bulk_Comparison.R             # DEG analysis between groups
├── 07_Epithelial_Subset_CNV.R          # Subset epithelial cells & inferCNV analysis
├── 08_Functional_GO_GSVA.R             # GO, KEGG, and GSVA enrichment
├── 09_Trajectory_Monocle.R             # Pseudotime inference (Monocle2)
├── 10_Regulon_SCENIC.R                 # Regulon inference with SCENIC
├── 11_TCGA_Prognosis_Validation.R      # Survival analysis using TCGA-LUAD
├── 12_Myeloid_Subclusters.R            # Subclustering myeloid cells
├── 13_Tcell_Subclusters.R              # Subclustering T cells
├── 14_CellCell_Communication.R         # CellChat-based interaction analysis

scRNA_scripts/
├── Bottom_left_axis.R                  # Draw bottom-left axis on UMAPs
├── check-all-markers.R                 # Canonical marker list and DotPlot generation
├── harmony.R                           # Harmony integration functions
├── lib.R                               # Shared utility functions
├── mycolors.R                          # Custom ggplot2 color palettes
├── qc.R                                # QC helper: mito/ribo/hb filtering, violin plots

main.Rproj                              # RStudio project file

📊 Output Highlights
UMAP & t-SNE: For cell visualization post-Harmony

DotPlots & ViolinPlots: Marker-based annotation

CNV Heatmaps: inferCNV to distinguish malignant epithelial cells

Trajectory Plots: Monocle2 for lineage relationships

Regulon Activity: SCENIC heatmaps & RSS scores

Survival Curves: TCGA-LUAD Kaplan–Meier plots

CellChat Networks: Stage-specific cell–cell communication diagrams

🧪 Tools & Packages
Seurat, harmony, monocle, infercnv, SCENIC, GSVA, clusterProfiler, CellChat

ggplot2, ggpubr, cowplot, patchwork, pheatmap
