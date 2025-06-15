# TCGA-LUAD Survival Analysis using scRNA-seq Derived Gene Signatures

rm(list = ls())

# ---------------------------------------------
# Step 1: Setup & Load Data
# ---------------------------------------------
dir.create("11-TCGA_LUAD", showWarnings = FALSE)
setwd("11-TCGA_LUAD")
source("../scRNA_scripts/mycolors.R")

library(tidyverse)
library(data.table)

# Define project
proj <- "TCGA-LUAD"
dir.create("input", showWarnings = FALSE)

# Download data if not already done
if (FALSE) {
  download.file(paste0("https://gdc.xenahubs.net/download/", proj, ".htseq_counts.tsv.gz"),
                destfile = paste0("input/", proj, ".htseq_counts.tsv.gz"))
  download.file(paste0("https://gdc.xenahubs.net/download/", proj, ".GDC_phenotype.tsv.gz"),
                destfile = paste0("input/", proj, ".GDC_phenotype.tsv.gz"))
  download.file(paste0("https://gdc.xenahubs.net/download/", proj, ".survival.tsv"),
                destfile = paste0("input/", proj, ".survival.tsv"))
}

clinical <- fread(paste0("input/", proj, ".GDC_phenotype.tsv.gz"), sep = "\t")
surv <- fread(paste0("input/", proj, ".survival.tsv"))

# ---------------------------------------------
# Step 2: Process Expression Matrix
# ---------------------------------------------
dat <- fread(paste0("input/", proj, ".htseq_counts.tsv.gz"))
dat <- dat[1:(nrow(dat) - 5), ]
rownames(dat) <- dat$Ensembl_ID
a <- as.matrix(2^dat[, -1] - 1)
exp <- apply(a, 2, as.integer)
rownames(exp) <- substr(rownames(dat), 1, 15)

# Annotate genes
library(AnnoProbe)
library(tinyarray)
re <- annoGene(rownames(exp), ID_type = "ENSEMBL")
exp <- trans_array(exp, ids = re, from = "ENSEMBL", to = "SYMBOL")
proj <- "tcga-luad"
save(exp, file = paste0(proj, ".htseq_counts.rdata"))

# ---------------------------------------------
# Step 3: Filter Tumor Samples with Survival Info
# ---------------------------------------------
rm(list = ls())
proj <- "tcga-luad"
load(paste0(proj, ".htseq_counts.rdata"))
Group <- ifelse(as.numeric(str_sub(colnames(exp), 14, 15)) < 10, "tumor", "normal")
exprSet <- exp[, Group == "tumor"]

clinical <- fread("input/TCGA-LUAD.GDC_phenotype.tsv.gz", sep = "\t")
surv <- fread("input/TCGA-LUAD.survival.tsv")
meta <- left_join(surv, clinical, by = c("sample" = "submitter_id.samples"))

meta <- meta[meta$OS.time >= 30 & !is.na(meta$OS.time) & !is.na(meta$OS), ]
meta <- meta[, c("sample", "OS", "OS.time")]
colnames(meta) <- c("ID", "event", "time")
meta$time <- meta$time / 30
rownames(meta) <- meta$ID

s <- intersect(rownames(meta), colnames(exprSet))
exprSet <- exprSet[, s]
meta <- meta[s, ]
save(exprSet, meta, file = paste0(proj, ".for_survival.rdata"))

# ---------------------------------------------
# Step 4: GSVA & Survival Analysis
# ---------------------------------------------
library(survival)
library(survminer)
library(GSEABase)
library(GSVA)

proj <- "tcga-luad"
load(paste0(proj, ".for_survival.rdata"))
phe <- meta
mySurv <- with(phe, Surv(time, event))

# Gene sets from scRNA-seq markers
markers <- read.csv("../7-epi/epimarkers.csv", row.names = 1)
markers$cluster <- as.factor(markers$cluster)
deg_list <- split(markers$gene, markers$cluster)
gs <- lapply(deg_list, toupper)

geneset <- GeneSetCollection(mapply(function(geneIds, clusterId) {
  GeneSet(geneIds, geneIdType = EntrezIdentifier(),
          collectionType = KEGGCollection(clusterId),
          setName = clusterId)
}, gs, names(gs)))

es.max <- gsva(as.matrix(exprSet), geneset, mx.diff = FALSE, verbose = FALSE, parallel.sz = 4)
pheatmap::pheatmap(es.max)

# Median-based survival analysis
splots <- lapply(names(deg_list), function(i) {
  v <- as.numeric(es.max[i, ])
  phe$sub_group <- ifelse(v < 0, "low", "high")
  fit <- survfit(Surv(time, event) ~ sub_group, data = phe)
  ggsurvplot(fit, data = phe, pval = TRUE, risk.table = TRUE,
             conf.int = TRUE, palette = "jco", ggtheme = theme_bw(),
             title = paste0("cluster_", i))
})

layout_dim <- ceiling(sqrt(length(splots)))
all_plot <- arrange_ggsurvplots(splots, print = FALSE, ncol = layout_dim, nrow = layout_dim)
ggsave(all_plot, filename = "all_survival_plot.pdf", width = layout_dim * 5, height = layout_dim * 5)

# Optimal cutpoint survival analysis
csplots <- lapply(names(deg_list), function(i) {
  v <- as.numeric(es.max[i, ])
  phe$v <- v
  cut <- surv_cutpoint(phe, time = "time", event = "event", variables = "v")
  cat <- surv_categorize(cut)
  sfit <- survfit(Surv(time, event) ~ v, data = cat)
  ggsurvplot(sfit, data = phe, pval = TRUE, risk.table = TRUE,
             conf.int = TRUE, palette = "jco", ggtheme = theme_bw(),
             title = paste0("cluster_", i))
})

layout_dim <- ceiling(sqrt(length(csplots)))
all_cut_plot <- arrange_ggsurvplots(csplots, print = FALSE, ncol = layout_dim, nrow = layout_dim)
ggsave(all_cut_plot, filename = "all_cut_point_survival_plot.pdf", width = layout_dim * 5, height = layout_dim * 5)

setwd("../")
