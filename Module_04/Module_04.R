### Step 4:: Clustering the Cells ###
#####################################
# Install and load required packages:

library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(limma)
library(SingleR)
library(pheatmap)
#####################################################
BLCA_M4 <- readRDS("BLCA_M3.RDS")

BLCA_M4[["RNA"]] <- JoinLayers(BLCA_M4[["RNA"]])
DefaultAssay(BLCA_M4) <- "RNA"
BLCA_M4

BLCA_M4@assays$RNA$data
#################################
#################################
## Clustering and Dimensional Reduction
## Nearest neighbor graph, clustering, and embedding ###

#### Step 1: Find Neighbors and Identify Clusters
# Identify neighbors using the first 30 principal components
BLCA_M4 <- FindNeighbors(BLCA_M4, dims = 1:30)

# Find the clusters
# BLCA_M4 <- FindClusters(BLCA_M4, resolution = c(0.1,0.3, 0.5, 0.7, 1))

# Perform clustering with a resolution of 0.5
BLCA_M4 <- FindClusters(BLCA_M4, resolution = 0.5)
################################
DimPlot(BLCA_M4, group.by = "RNA_snn_res.0.5", dims = c(1,3),
        label = FALSE)

# setting identity of clusters
Idents(BLCA_M4)
Idents(BLCA_M4) <- "RNA_snn_res.0.5"
Idents(BLCA_M4)

# Display cluster IDs for the first 5 cells
head(Idents(BLCA_M4), 5)
########################################################
########################################################
#### Step 2: Run Dimensional Reduction (UMAP and t-SNE)
########################################################
# Run UMAP for non-linear dimensional reduction
BLCA_M4 <- RunUMAP(BLCA_M4, dims = 1:30, verbose = TRUE)

# Plot the UMAP with clustering
DimPlot(BLCA_M4, reduction = "umap", dims = c(1,2))
DimPlot(BLCA_M4, reduction = "umap", label = FALSE, 
        repel = TRUE, group.by = "Batch")

#### Explore Cluster Metadata
# Display the number of cells per cluster
table(BLCA_M4@meta.data$seurat_clusters)
######################################
plot9 <- DimPlot(BLCA_M4, reduction = "umap", 
                 label = TRUE,
                 repel = TRUE, 
                 group.by = c("GSE_ID", "Batch"),
                 label.size = 4) & NoLegend() & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC26.png", plot = plot9, width = 15, height = 8)

plot9 <- DimPlot(BLCA_M4, reduction = "umap", 
                 label = FALSE, 
                 repel = TRUE,
                 group.by = c("GSE_ID", "Batch"),
                 label.size = 4) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC27.png", plot = plot9, width = 18, height = 8)
######################################
plot10 <- DimPlot(BLCA_M4, reduction = "umap", 
                  label = TRUE, 
                  repel = TRUE, 
                  group.by = c("seurat_clusters", "Batch"),
                  label.size = 4) & NoLegend() & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC28.png", plot = plot10, width = 15, height = 8)

plot10 <- DimPlot(BLCA_M4, reduction = "umap", 
                  label = FALSE, 
                  repel = TRUE, 
                  group.by = c("seurat_clusters", "Batch"),
                  label.size = 4) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC29.png", plot = plot10, width = 18, height = 8)
######################################
# setting identity of clusters
Idents(BLCA_M4)
Idents(BLCA_M4) <- "GSE_ID"
Idents(BLCA_M4)

features = c("TIMP1", "CD63", "MMP9", 
             "MMP14", "IL6", "FN1")

plot12 <- FeaturePlot(BLCA_M4, features = features, 
                      label = TRUE, repel = TRUE,
                      split.by = "Batch", label.size = 3,
                      cols = c('lightgrey', '#943bcc'), 
                      reduction = "umap") & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC30.png", plot = plot12, width = 10, height = 20)

plot12 <- FeaturePlot(BLCA_M4, features = features, 
                      label = FALSE, repel = TRUE,
                      split.by = "Batch", label.size = 3,
                      cols = c('lightgrey', '#943bcc'), 
                      reduction = "umap") & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC31.png", plot = plot12, width = 10, height = 20)
####################################################################
####################################################################
# Get the TSNE embedding .............
## check_duplicates = FALSE
BLCA_M4 <- RunTSNE(BLCA_M4, dims = 1:30, 
                   verbose = TRUE)
#########################################
#########################################
# Plot the TSNE with clustering
DimPlot(BLCA_M4, reduction = "tsne")

plot10 <- DimPlot(object = BLCA_M4, reduction = "tsne", 
                  group.by = c("seurat_clusters", "Batch"),
                  label = TRUE, repel = TRUE, label.size = 4) & NoLegend()
ggsave(filename = "QC32.png", plot = plot10, width = 20, height = 10)

plot10 <- DimPlot(object = BLCA_M4, reduction = "tsne", 
                  group.by = c("GSE_ID", "Batch"),
                  label = TRUE, repel = TRUE, label.size = 4) & NoLegend()
ggsave(filename = "QC33.png", plot = plot10, width = 20, height = 10)

plot10 <- DimPlot(object = BLCA_M4, reduction = "tsne", 
                  group.by = c("seurat_clusters", "Batch"),
                  repel = TRUE, label.size = 4)
ggsave(filename = "QC34.png", plot = plot10, width = 22, height = 10)

plot10 <- DimPlot(object = BLCA_M4, reduction = "tsne", 
                  group.by = c("GSE_ID", "Batch"),
                  repel = TRUE, label.size = 4)
ggsave(filename = "QC35.png", plot = plot10, width = 22, height = 10)

plot10 <- DimPlot(object = BLCA_M4, reduction = "tsne",
                  repel = TRUE, label.size = 4)
ggsave(filename = "QC36.png", plot = plot10, width = 10, height = 8)

plot10 <- DimPlot(object = BLCA_M4, reduction = "umap",
                  repel = TRUE, label.size = 4)
ggsave(filename = "QC37.png", plot = plot10, width = 10, height = 8)
########################################
########################################   ****** here ------------------------
# Define color palettes
library(RColorBrewer)
library(pheatmap)
# WtOrRd_pal <- colorRampPalette(c("#a742f5", brewer.pal(9, "#EFE9C3")))(100)
WtOrRd_pal <- colorRampPalette(c("white", "#943bcc"))(100)

# Compare Batch and Cluster ID
compTable <- table(BLCA_M4$orig.ident, BLCA_M4$seurat_clusters)
compTable <- (compTable / rowSums(compTable)) * 100
plot11 <- pheatmap(compTable, color = WtOrRd_pal)
ggsave(filename = "QC38.png", plot = plot11, width = 15, height = 10)
#########################################
##########################################
plot12 <- FeaturePlot(BLCA_M4, features = features, 
                      label = TRUE, repel = TRUE,
                      split.by = "Batch", label.size = 3,
                      cols = c('lightgrey', '#943bcc'), 
                      reduction = "tsne") 
ggsave(filename = "QC39.png", plot = plot12, width = 10, height = 20)

plot12 <- FeaturePlot(BLCA_M4, features = features, 
                      label = FALSE, repel = TRUE,
                      split.by = "Batch",
                      cols = c('lightgrey', '#943bcc'), 
                      reduction = "tsne") 
ggsave(filename = "QC40.png", plot = plot12, width = 10, height = 20)

plot12 <- FeaturePlot(BLCA_M4, features = features, 
                      label = FALSE, repel = TRUE,
                      cols = c('lightgrey', '#943bcc'), 
                      reduction = "tsne") 
ggsave(filename = "QC41.png", plot = plot12, width = 14, height = 20)

View(BLCA_M4@meta.data)
#########################################
#########################################
Idents(BLCA_M4) <- "orig.ident"

rami1 <- DimPlot(BLCA_M4, reduction = "umap", label = FALSE)
plot13 <- (rami1) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC42.png", plot = plot13, width = 11, height = 8)

rami2 <- DimPlot(BLCA_M4, reduction = "umap", label = FALSE, 
                 group.by = "Batch")
plot13 <- (rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC43.png", plot = plot13, width = 10, height = 8)
########################
########################
plot12 <- FeaturePlot(BLCA_M4, features=c("TIMP1", "CD63", "MMP9", 
                                          "MMP14", "IL6", "FN1"), 
                      repel = TRUE, 
                      cols = c('lightgrey', '#943bcc'), reduction = "umap")
plot12 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2")  
ggsave(filename = "QC44.png", plot = plot12, width = 10, height = 14)
#######################
#######################
plot13 <- FeaturePlot(BLCA_M4, features=c("TIMP1", "CD63", "MMP9", 
                                          "MMP14", "IL6", "FN1"),
                      cols = c('lightgrey', '#943bcc'), reduction = "tsne")
ggsave(filename = "QC45.png", plot = plot13, width = 10, height = 14)
###########################################################################
plot10 <- DimPlot(object = BLCA_M4, reduction = "tsne", 
                  group.by = "seurat_clusters", 
                  label = TRUE, repel = TRUE, label.size = 6) + NoLegend()

plot11 <- DimPlot(object = BLCA_M4, reduction = "tsne", 
                  group.by = "Batch", 
                  label = TRUE, repel = TRUE, label.size = 6) + NoLegend()


plot14 <- plot10 + plot13
plot15 <- plot11 + plot13

ggsave(filename = "QC46.png", plot = plot14, width = 22, height = 12)
ggsave(filename = "QC47.png", plot = plot15, width = 22, height = 12)
##########################################################################
Idents(BLCA_M4) <- "orig.ident"

rami1 <- DimPlot(BLCA_M4, reduction = "tsne", label = FALSE)
ggsave(filename = "QC48.png", plot = rami1, width = 11, height = 8)

rami2 <- DimPlot(BLCA_M4, reduction = "tsne", label = FALSE, 
                 group.by = "Batch")
ggsave(filename = "QC49.png", plot = rami2, width = 10, height = 8)

Idents(BLCA_M4) <- "GSE_ID"

rami1 <- DimPlot(BLCA_M4, reduction = "tsne", label = FALSE)
ggsave(filename = "QC50.png", plot = rami1, width = 11, height = 8)
##########################################################################
Idents(BLCA_M4) <- "orig.ident"

rami1 <- DimPlot(BLCA_M4, reduction = "tsne", repel = TRUE,
                 label.size = 3,
                 label = TRUE) & NoLegend()
ggsave(filename = "QC51.png", plot = rami1, width = 9, height = 8) 

rami2 <- DimPlot(BLCA_M4, reduction = "tsne", label = TRUE,
                 label.size = 3, repel = TRUE,
                 group.by = "Batch") & NoLegend()
ggsave(filename = "QC52.png", plot = rami2, width = 9, height = 8) 

rami2 <- DimPlot(BLCA_M4, reduction = "tsne", label = FALSE, 
                 group.by = "seurat_clusters")
ggsave(filename = "QC53.png", plot = rami2, width = 10, height = 8)
###########################################################################
plot6 <- DimPlot(BLCA_M4, reduction = "pca", dims = c(1,3), 
                 group.by = "Batch",
                 label = FALSE, 
                 repel = TRUE, label.size = 6) + NoLegend()

plot8 <- FeaturePlot(BLCA_M4, features = features, 
                     raster = TRUE,
                     cols = c("lightgrey", "#943bcc"), 
                     reduction = "pca", dims = c(1,3))

plot16 <- plot6 + plot8
ggsave(filename = "QC54.png", plot = plot16, width = 18, height = 10)
############################################################################
plot9 <- DimPlot(BLCA_M4, reduction = "umap", label = TRUE,
                 repel = TRUE, 
                 group.by = "seurat_clusters",
                 label.size = 6) + NoLegend()

plot10 <- (plot9) & labs(x = "UMAP_1", y = "UMAP_2")
plot16 <- FeaturePlot(BLCA_M4, features = features, 
                      cols = c("lightgrey", "#943bcc"), reduction = "umap")

plot18 <- (plot16) & labs(x = "UMAP_1", y = "UMAP_2")
plot17 <- plot10 + plot18
ggsave(filename = "QC55.png", plot = plot17, width = 20, height = 12)
################################################
### Advanced Data visualization methods in Seurat
################################################
# Ridge plots - Visualize single cell expression distributions in each cluster
RidgePlot(BLCA_M4, features = features, ncol = 2)

rami1 <- RidgePlot(BLCA_M4, features = "TIMP1", ncol = 1) & NoLegend()
ggsave(filename = "QC56.png", plot = rami1, width = 16, height = 12)

rami1 <- RidgePlot(BLCA_M4, features = "TIMP1", group.by = "Batch",
                   ncol = 1) & NoLegend()
ggsave(filename = "QC57.png", plot = rami1, width = 10, height = 8)
#################################################
#################################################
# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(BLCA_M4, features = features)

rami1 <- VlnPlot(BLCA_M4, features = features, group.by = "Batch")
ggsave(filename = "QC58.png", plot = rami1, width = 16, height = 12)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(BLCA_M4, features = features) + RotatedAxis()

rami1 <- DotPlot(BLCA_M4, features = features, 
                 group.by = "Batch") + RotatedAxis()
ggsave(filename = "QC59.png", plot = rami1, width = 10, height = 8)
ggsave(filename = "QC59.jpg", plot = rami1, width = 10, height = 8)
##################################################
##################################################
# Adjust the contrast in the plot
FeaturePlot(BLCA_M4, features = "TIMP1", cols = c("lightgrey", "#943bcc"),
            min.cutoff = 1, max.cutoff = 3)

# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
FeaturePlot(BLCA_M4, features = features, cols = c("lightgrey", "#943bcc"),
            min.cutoff = "q10", max.cutoff = "q90")
################################################
################################################
saveRDS(BLCA_M4, file = "BLCA_M5.RDS")
BLCA_M5 <- readRDS("BLCA_M5.RDS")
################################################