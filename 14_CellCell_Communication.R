### Interaction Between Cancer and TME Using CellChat

```r
rm(list=ls())
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(CellChat)

options(stringsAsFactors = FALSE)
dir.create("13-cellchat")
setwd("13-cellchat/")

# Load pre-annotated Seurat objects
Tsce <- readRDS("../12-T/Tsce_celltype.rds")
myeloidsce <- readRDS("../12-T/myeloidsce_celltype.rds")
epi <- readRDS("../7-epi/epi_sce_celltype.rds")

# Subset cancer cells from epithelial data
target_cancer_types <- c('Cancer(CRABP2+)', 'Cancer(TM4SF1+)', 'Cancer(UBE2C+)', 'Clara-like cancer')
cancer <- epi[, epi$celltype %in% target_cancer_types]

# Combine T, myeloid, and cancer cells
sce <- merge(Tsce, list(myeloidsce, cancer), add.cell.ids = c('T', 'myeloid', 'cancer'))
Idents(sce) <- sce$celltype
sce <- JoinLayers(sce)

# Subset by sample type
AIS <- sce[, sce$sample == 'AIS']
MIA <- sce[, sce$sample == 'MIA']
IAC <- sce[, sce$sample == 'IAC']

save(sce, file = 'input.Rdata')

# CellChat wrapper for individual group
cellchatA <- function(sce_obj) {
  cellchat <- createCellChat(sce_obj@assays$RNA$data, meta = sce_obj@meta.data, group.by = "celltype")
  cellchat@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

# Run CellChat on each subtype
AIScellchat <- cellchatA(AIS)
MIAcellchat <- cellchatA(MIA)
IACcellchat <- cellchatA(IAC)

# Merge results for comparison
object.list <- list(AIS = AIScellchat, MIA = MIAcellchat, IAC = IACcellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Compare communication count and weight
gg1 <- compareInteractions(cellchat, group = c(1,2,3), measure = "count", show.legend = FALSE)
gg2 <- compareInteractions(cellchat, group = c(1,2,3), measure = "weight", show.legend = FALSE)
ggsave("Overview_number_strength.pdf", gg1 + gg2, width = 10, height = 7)

# Differential network visualization
par(mfrow = c(1,3))
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")

# Circle plots for each sample group
par(mfrow = c(1,3))
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
for (i in seq_along(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Interactions - ", names(object.list)[i]))
}

# Rank and compare pathway activity
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
ggsave("Compare_pathway_strengh.pdf", gg1 + gg2, width = 10, height = 6)

# Visualization of TGFb pathway
pathways.show <- "TGFb"
weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathways.show)
par(mfrow = c(1,3), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Heatmap comparison for TGFb
par(mfrow = c(1,3), xpd = TRUE)
ht <- lapply(seq_along(object.list), function(i) {
  netVisual_heatmap(object.list[[i]], signaling = pathways.show,
                    color.heatmap = "Reds",
                    title.name = paste(pathways.show, "signaling", names(object.list)[i]))
})
ComplexHeatmap::draw(Reduce(`+`, ht), ht_gap = unit(0.5, "cm"))
```
