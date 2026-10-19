####################################################
### Step 5:: scRNA-seq Data Integration:: ....... ###
## Perform SCT Integration ..................
BiocManager::install("glmGamPoi")

library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(Matrix)
library(SeuratWrappers)
library(limma)
library(monocle3)
library(glmGamPoi)
######################################
# load dataset
BLCA_M5 <- readRDS("BLCA_M5.RDS")

###################################
###################################
# split the RNA measurements into layers
BLCA_M5[["RNA"]] <- split(BLCA_M5[["RNA"]], f = BLCA_M5$Batch)
BLCA_M5

# Increase the future.globals.maxSize limit
options(future.globals.maxSize = 25*1024^3)  # 25 GB

BLCA_M5 <- SCTransform(BLCA_M5)
BLCA_M5 <- RunPCA(BLCA_M5)
BLCA_M5 <- RunUMAP(BLCA_M5, dims = 1:30)
DimPlot(BLCA_M5, reduction = "umap", group.by = c("Batch", "seurat_clusters"))

# integrate datasets
BLCA_M5 <- IntegrateLayers(object = BLCA_M5, method = CCAIntegration, 
                           normalization.method = "SCT", verbose = F)
###################################
# Specify that we will perform downstream analysis on the integrated data 
DefaultAssay(BLCA_M5) <- "SCT"

# re-join layers after integration
BLCA_M5[["RNA"]] <- JoinLayers(BLCA_M5[["RNA"]])

# Re-run the pipeline on the integrated data
BLCA_M5 <- BLCA_M5 |> 
  ScaleData() |> 
  RunPCA() |> 
  FindNeighbors(reduction = "integrated.dr", dims = 1:30) |> 
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = "integrated.dr") |> 
  RunTSNE(dims = 1:30, reduction = "integrated.dr")
###################################
Idents(BLCA_M4) <- "seurat_clusters"
plot12 <- DimPlot(BLCA_M4, reduction = "umap")
plot13 <- DimPlot(BLCA_M5, reduction = "umap")
plot14 <- (plot12 | plot13) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC60.png", plot = plot14, width = 12, height = 5)
####################################
plot12 <- DimPlot(BLCA_M5, reduction = "umap", 
                  group.by = c("Batch", "seurat_clusters"))
plot12 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC61.png", plot = plot12, width = 12, height = 5)
###################################
plot12 <- FeaturePlot(BLCA_M5, features = "TIMP1", 
                      label = TRUE, repel = TRUE,
                      split.by = "Batch",
                      cols = c('lightgrey', '#943bcc'), 
                      reduction = "umap")
plot12 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC62.png", plot = plot12, width = 18, height = 6)
##########################################
# Dim Plots
DimPlot(BLCA_M5, group.by = "orig.ident") + DimPlot(BLCA_M5, group.by = "seurat_clusters")
##########################################
# Visualization
DimPlot(BLCA_M5, reduction = "umap", group.by = c("Batch", "seurat_clusters"))
DimPlot(BLCA_M5, reduction = "umap", split.by = "Batch")
##########################################
plot12 <- DimPlot(BLCA_M5, reduction = "umap", split.by = "Batch")
plot12 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC63.png", plot = plot12, width = 12, height = 5)
##########################################
plot12 <- DimPlot(BLCA_M5, reduction = "umap", split.by = "Batch", label = TRUE)
plot12 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2") & NoLegend()
ggsave(filename = "QC64.png", plot = plot12, width = 12, height = 6)
##########################################
saveRDS(BLCA_M5, file = "BLCA_M6.RDS")
BLCA_M6 <- readRDS("BLCA_M6.RDS")
##########################################