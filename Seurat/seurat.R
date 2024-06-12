#!/usr/bin/env Rscript

# Load libraries
# install.packages("hdf5r")
# BiocManager::install("glmGamPoi")
# devtools::install_github("immunogenomics/presto")

# remotes::install_github("paodan/studySeu")

library(Seurat)
library(tidyverse)
library(sctransform)
library(ggplot2)
library(gridExtra)

# install.packages("remotes")

args <- commandArgs(trailingOnly=TRUE)
h5_path <- args[1]
annot_path <- args[2]
tag <- args[3]
 
# h5_path = 'X:/scratch/delete90/paguirigan_a/smoorthi/cromwell-executions/SCG_10x_AML_Genotyping/8b5ef016-66a4-48c9-8cec-feeca490ce7b/call-cellRangerGEX//shard-2/cacheCopy/execution/glob-ca744b32743fde4269f708a8451f389e/filtered_feature_bc_matrix.h5'
# annot_path = 'X:/fast/paguirigan_a/sanaz/Research/GeneExpression/genotype_npm1_unfiltered/10k_mix_5050_npm1_genotype.txt'
# tag='mol50'  

# Initialize Seurat object with the raw (non-normalized data)
Sobj <- CreateSeuratObject(counts = Read10X_h5(filename=h5_path), project='test1', min.cells = 3, min.features = 200)

# Load annotation file
ann <- read.csv(annot_path, header=TRUE, sep = '\t')
ann <- ann[!duplicated(ann$cellBarcode),]
ann <- ann[c('cellBarcode', 'genotype')]
ann <- ann %>%  mutate(bar=paste0(cellBarcode, '-1'))
annot_dict <- c('0/0' = 'WT', '0/1' = 'Heterozygous mutation', '1/1' = 'Homozygous mutation')
ann <- ann %>%  mutate(genotype=recode(genotype, !!!annot_dict))

# Merge annotation with Sobj
vec <- setNames(ann$genotype, ann$bar)
barcodes <- colnames(Sobj)
Sobj@meta.data$genotype <- vec[barcodes] 
Sobj@meta.data$genotype[is.na(Sobj@meta.data$genotype)] <- "Unknown"


## Quality Control
Sobj[["percent.mt"]] <- PercentageFeatureSet(Sobj, pattern = "^MT-")

# View(Sobj@meta.data)
# VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
# FeatureScatter(Sobj, feature1='nCount_RNA', feature2='nFeature_RNA') + geom_smooth(method='lm')

# Use doublet finder

Sobj <- subset(Sobj, subset = nFeature_RNA > 200 & percent.mt < 5)   # Filtering

################
# Sobj <- NormalizeData(Sobj, normalization.method="LogNormalize", scale.factor=10000)   # Normalize data
# Sobj <- FindVariableFeatures(Sobj, selection.method = "vst", nfeatures=2000)      # Identify HVG
# Sobj <- ScaleData(Sobj, features = rownames(Sobj), ) #vars.to.regress=c("percent.mt"))       # scale - so that HVG do not dominate, regress out artifcats - biological (cell cycle), technical (batch effect)

Sobj <- SCTransform(Sobj, vars.to.regress = "percent.mt", verbose = TRUE, return.only.var.genes = TRUE)


################

Sobj <- RunPCA(Sobj, features = VariableFeatures(object=Sobj))                   # PCA - Linear Dimensionality Reduction

# Visualize PCA Results
# ElbowPlot(Sobj)

Sobj <- FindNeighbors(Sobj, dims=1:20)     # Find Neighbours
Sobj <- FindClusters(Sobj, resolution=0.05)  # Find Cluster (Louvain)
Sobj <- RunUMAP(Sobj, dims=1:20)          # Non-linear dim. reduction

plot1 <- DimPlot(Sobj, reduction = "umap") + ggtitle(tag)
plot2 <- DimPlot(Sobj, reduction = "umap", group.by = 'genotype')



ggsave(paste("umap_", tag, "_seurat_cluster.png", sep = ""), plot = arrangeGrob(plot1, plot2, nrow = 1), width=16, height=8)
## for DE Analysis

markers <- FindAllMarkers(Sobj, test.use='wilcox', only.pos=FALSE, logfc.threshold=0.5, min.pct=0.25, return.thresh=1e-5)
write.csv(markers, paste('seurat_DEG_',tag,'.csv', sep=""), row.names = FALSE)



