### Step 4:: Principal Component Analysis ###
##############################################
# Run PCA and identify strong expressed genes
# Install and load required packages:

library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(limma)
#############################################
BLCA_M3 <- readRDS("BLCA_M2.RDS")
#############################################
#############################################
### Step 4:: Principal Component Analysis ###
#############################################
# Run PCA
# BLCA_M3 <- RunPCA(BLCA_M3)
BLCA_M3 <- RunPCA(BLCA_M3, npcs = 50,
                  assay = "RNA", 
                  features = VariableFeatures(object = BLCA_M3),
                  ndims.print = 1:10, verbose = TRUE, 
                  nfeatures.print = 30)
##################################################################
##################################################################
## View Top Genes Contributing to PCs
# Display the most important genes contributing to the top 10 principal components
print("Top genes contributing to the first 6 principal components:")
print(BLCA_M3[["pca"]], dims = 1:10, nfeatures = 10)

#####################################################################
#####################################################################
# What are the genes that determine PC1 and PC2?
## Visualize PCA Loadings
# Plot the top 15 genes contributing to each of the first 5 PCs
VizDimLoadings(BLCA_M3, dims = 1:4)  # dims = 1:4

plot9 <- VizDimLoadings(BLCA_M3, dims = 1:4, nfeatures = 30,
                        col = "blue", reduction = "pca", combine = TRUE)
ggsave(filename = "QC17.png", plot = plot9, dpi = 600,
       width = 18, height = 10)
#####################################################################
#####################################################################
# Plot PCs 1 and 5 
# DimPlot(), and DimHeatmap()
DimPlot(BLCA_M3)
DimPlot(BLCA_M3, reduction = "pca", dims = c(1,3))
DimPlot(BLCA_M3, reduction = "pca", dims = c(1,5))

Idents(BLCA_M3) <- "orig.ident"
plot6 <- DimPlot(BLCA_M3, reduction = "pca", dims = c(1,3))
ggsave(filename = "QC18.png", plot = plot6, dpi = 600,
       width = 14, height = 10)

plot7 <- DimPlot(BLCA_M3, reduction = "pca", dims = c(1,3),
                 group.by = "Batch", raster = FALSE)
ggsave(filename = "QC19.png", plot = plot7, dpi = 600,
       width = 14, height = 10)

plot8 <- plot6+plot7
ggsave(filename = "QC20.png", plot = plot8, dpi = 600,
       width = 25, height = 10)
#####################################################################
#####################################################################
## Dimensional Heatmap
# Dimensional heatmap for the first 3 PCs with a subset of 500 cells for clarity
DimHeatmap(BLCA_M3, dims = 1)
DimHeatmap(BLCA_M3, dims = 3)

DimHeatmap(BLCA_M3, dims = 1:3, cells = 500, balanced = TRUE)

png("QC21.png", width = 8000, height = 2800, 
    res = 600)
DimHeatmap(BLCA_M3, dims = c(1:3), cells = 500, balanced = TRUE)
dev.off() 
#####################################################################
#####################################################################
## Determine Optimal Number of PCs ...... very important
# Elbow plot to determine the optimal number of principal components
ElbowPlot(BLCA_M3)
ElbowPlot(BLCA_M3, ndims = 50)

plot7 <- ElbowPlot(BLCA_M3, ndims = 50, reduction = "pca")
ggsave(filename = "QC22.png", plot = plot7, width = 10, height = 8)
######################################################################
#####################################################################
BLCA_M3

# Save our Seurat object
saveRDS(BLCA_M3, file = "BLCA_M3.RDS")
BLCA_M4 <- readRDS("BLCA_M3.RDS")
#####################################
#####################################
# FeaturePlot of PCA top features 
FeaturePlot(BLCA_M3, features = c("ST14", "KLF5"), dims = c(1,3))

FeaturePlot(BLCA_M3, features = c("TIMP1", "MMP9"), dims = c(1,3))

plot8 <- FeaturePlot(BLCA_M3, features = c("TIMP1", "CD63"), blend = TRUE,
                     dims = c(1,3), raster = TRUE)
ggsave(filename = "QC23.png", plot = plot8, width = 15, height = 5, dpi = 600)

plot8 <- FeaturePlot(BLCA_M3, features = c("TIMP1", "CD63", "MMP9",
                                           "MMP14", "IL6", "FN1"), 
                     cols = c("lightgrey", "#943bcc"), #f2a350
                     dims = c(1,3), raster = TRUE)
ggsave(filename = "QC24.png", plot = plot8, width = 10, height = 15, dpi = 600)

plot8 <- FeaturePlot(BLCA_M3, features = c("TIMP1", "CD63", "MMP9",
                                           "MMP14", "IL6", "FN1"), 
                     cols = c("lightgrey", "#943bcc"), 
                     dims = c(1,3), split.by = "Batch")
ggsave(filename = "QC25.png", plot = plot8, width = 10, height = 20, dpi = 600)
#######################################################################
#######################################################################