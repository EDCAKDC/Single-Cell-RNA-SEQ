# Clear environment
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(Matrix)
  library(data.table)
})

# Set the data directory
data_dir <- "xxxx/"
sample_dirs <- list.files(data_dir, full.names = TRUE)
sample_names <- basename(sample_dirs)

# Print all detected sample folder names
print(sample_names)

# Read each sample and store as Seurat objects in a list
seurat_list <- vector("list", length(sample_dirs))
names(seurat_list) <- sample_names

for (i in seq_along(sample_dirs)) {
  message("Processing sample: ", sample_names[i])
  count_data <- Read10X(data.dir = sample_dirs[i])
  
  # If count_data is multimodal (e.g. RNA and protein), use RNA counts only
  if (is.list(count_data)) {
    count_data <- count_data[[1]]
  }
  
  # Create Seurat object with filtering
  seurat_obj <- CreateSeuratObject(
    counts = count_data,
    project = sample_names[i],
    min.cells = 5,
    min.features = 300
  )
  seurat_list[[i]] <- seurat_obj
}

# Show dimensions (genes x cells) for each sample
lapply(seurat_list, dim)

# Merge all Seurat objects into one combined object
combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

# Print the merged object
print(combined)

# Join layers (required for some Seurat versions with multi-layer support)
combined <- JoinLayers(combined)

# Print the dimension of the count matrix
print(dim(combined[["RNA"]]$counts))

# Display a preview of the count matrix and metadata
head(as.data.frame(combined[["RNA"]]$counts[1:10, 1:2]))
head(combined@meta.data)

# Show sample distribution
table(combined$orig.ident)
