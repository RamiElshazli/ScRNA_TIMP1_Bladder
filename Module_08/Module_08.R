#########################################
### Get the cell cycle scores ..........
########################################
## Load libraries ##
library(Seurat)
library(VennDiagram)
library(tidyverse)
library(ggpubr)
library(phateR)
library(plotly)
library(Rmagic)
library(UpSetR)
library(cowplot)
library(SingleCellExperiment)
library(patchwork)
library(pheatmap)
library(enrichR)
library(RColorBrewer)
library(reticulate)
library(phateR)

# LOad your data
BLCA_M9 <- readRDS("BLCA_M8.RDS")
features <- c("TIMP1", "CD63", "MMP9", "MMP14", "IL6", "FN1")

# Get the cell cycle scores
BLCA_M9 <- CellCycleScoring(BLCA_M9, search = TRUE,
                            s.features = cc.genes$s.genes,
                            g2m.features = cc.genes$g2m.genes)
################################################################
################################################################
plot53 <- DimPlot(BLCA_M9, group.by = "Phase") + DimPlot(BLCA_M9, group.by = "cell_type")

plot54 <- (plot53) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC121.png", plot = plot54, width = 22, height = 10)
#########################
#########################
plot52 <- DimPlot(BLCA_M9, group.by = "Phase")
plot53 <- DimPlot(BLCA_M9, group.by = "cell_type", 
                  label = TRUE, repel = TRUE, label.size = 5) & NoLegend()

plot54 <- (plot52 | plot53) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC122.png", plot = plot54, width = 22, height = 10)
########################
########################
plot53 <- DimPlot(BLCA_M9, group.by = "Phase") & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC123.png", plot = plot53, width = 9, height = 8)

plot53 <- DimPlot(BLCA_M9, group.by = "Phase", reduction = "tsne")
ggsave(filename = "QC124.png", plot = plot53, width = 9, height = 8)
########################
########################
FeaturePlot(BLCA_M9, features = features, split.by = "Phase")

rami2 <- FeaturePlot(BLCA_M9, features = features,
                     split.by = "Phase", cols = c("lightgrey", "#943bcc"))
rami3 <- (rami2) & labs(x = "UMAP_1", y = "UMAP_2")

ggsave(filename = "QC125.png", plot = rami3, width = 14, height = 25)
#############################
rami2 <- FeaturePlot(BLCA_M9, features = features,
                     label = TRUE, repel = TRUE, label.size = 3,
                     split.by = "Phase", cols = c("lightgrey", "#943bcc"))
rami3 <- (rami2) & labs(x = "UMAP_1", y = "UMAP_2")

ggsave(filename = "QC126.png", plot = rami3, width = 14, height = 25)
##############################
##############################
rami4 <- FeaturePlot(BLCA_M9, features = features, reduction = "tsne",
                     split.by = "Phase", cols = c("lightgrey", "#943bcc"))

ggsave(filename = "QC127.png", plot = rami4, width = 14, height = 25)
##############################
##############################
rami4 <- FeaturePlot(BLCA_M9, features = features, reduction = "tsne",
                     label = TRUE, repel = TRUE, label.size = 3,
                     split.by = "Phase", cols = c("lightgrey", "#943bcc"))

ggsave(filename = "QC128.png", plot = rami4, width = 14, height = 25)
#####################################################
#####################################################
## https://www.scembroider.org/visualization/
#####################################################
# Plot with clusters
compTable.int <- table(BLCA_M9$cell_type, BLCA_M9$Phase)
compTable.int <- (compTable.int / rowSums(compTable.int)) * 100
plot29 <- pheatmap(compTable.int, color = colorRampPalette(c('white','#943bcc'))(100), breaks = 1:20)
ggsave(filename = "QC129.png", plot = plot29, width = 20, height = 10)

# FeaturePlot for showing mod score
rami5 <- FeaturePlot(BLCA_M9, features = c("S.Score", "G2M.Score"), 
                     reduction = "umap")
rami6 <- FeaturePlot(BLCA_M9, features = c("S.Score", "G2M.Score"), 
                     reduction = "tsne")

rami7 <- (rami5) & labs(x = "UMAP_1", y = "UMAP_2")

ggsave(filename = "QC130.png", plot = rami7, width = 20, height = 10)
ggsave(filename = "QC131.png", plot = rami6, width = 20, height = 10)
###############################
###############################
rami5 <- FeaturePlot(BLCA_M9, features = c("S.Score", "G2M.Score"),
                     label = TRUE, repel = TRUE, label.size = 4.5,
                     reduction = "umap")
rami7 <- (rami5) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC132.png", plot = rami7, width = 20, height = 10)

rami6 <- FeaturePlot(BLCA_M9, features = c("S.Score", "G2M.Score"), 
                     label = TRUE, repel = TRUE, label.size = 4.5,
                     reduction = "tsne")
ggsave(filename = "QC133.png", plot = rami6, width = 20, height = 10)
######################################################
rami5 <- FeaturePlot(BLCA_M9, features = c("S.Score", "G2M.Score"), 
                     reduction = "umap", label = TRUE, repel = TRUE)
rami7 <- (rami5) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC134.png", plot = rami7, width = 20, height = 10)

rami5 <- FeaturePlot(BLCA_M9, features = c("S.Score", "G2M.Score"), 
                     reduction = "umap", cols = c("lightgrey", "#943bcc"),
                     label = TRUE, repel = TRUE)
rami7 <- (rami5) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC135.png", plot = rami7, width = 20, height = 10)

rami5 <- FeaturePlot(BLCA_M9, features = c("S.Score", "G2M.Score"), 
                     reduction = "tsne", cols = c("lightgrey", "#943bcc"),
                     label = TRUE, repel = TRUE)
ggsave(filename = "QC136.png", plot = rami5, width = 20, height = 10)
#######################################################
#######################################################
rami1 <- DotPlot(BLCA_M9, features = features, group.by = "cell_type") + RotatedAxis()
rami2 <- rami1 + scale_colour_gradient2(low = "red", mid = "#89CFF0", high = "#943bcc")
ggsave(filename = "QC137.jpg", plot = rami2, width = 8, height = 6)
ggsave(filename = "QC137.png", plot = rami2, width = 8, height = 6)

rami1 <- DotPlot(BLCA_M9, features = features, group.by = "Batch") + RotatedAxis() 
rami2 <- rami1 + scale_colour_gradient2(low = "red", mid = "#89CFF0", high = "#943bcc")
ggsave(filename = "QC138.png", plot = rami2, width = 8, height = 6)
ggsave(filename = "QC138.jpg", plot = rami2, width = 8, height = 6)

rami1 <- DotPlot(BLCA_M9, features = features, group.by = "Batch") + RotatedAxis() 
rami2 <- rami1 + scale_colour_gradient2()
ggsave(filename = "QC139.png", plot = rami2, width = 8, height = 6)
ggsave(filename = "QC139.jpg", plot = rami2, width = 8, height = 6)

rami1 <- DotPlot(BLCA_M9, features = features, group.by = "cell_type") + RotatedAxis()
rami2 <- rami1 + scale_colour_gradient2()
ggsave(filename = "QC140.jpg", plot = rami2, width = 8, height = 6)
ggsave(filename = "QC140.png", plot = rami2, width = 8, height = 6)
#################################
#################################
# Ridge plots - from ggridges. 
# Visualize single cell expression distributions in each cluster
rami1 <- RidgePlot(BLCA_M9, features = features, ncol = 2)
ggsave(filename = "QC141.png", plot = rami1, width = 14, height = 18)

rami1 <- RidgePlot(BLCA_M9, features = "TIMP1", ncol = 1) & NoLegend()
ggsave(filename = "QC142.png", plot = rami1, width = 6, height = 6)

rami1 <- RidgePlot(BLCA_M9, features = c("TIMP1", "CD63", "MMP14", "FN1"), 
                   ncol = 2) & NoLegend()
ggsave(filename = "QC143.png", plot = rami1, width = 10, height = 10)
########################################
# Violin plot - Visualize single cell expression distributions in each cluster
rami1 <- VlnPlot(BLCA_M9, features = features)
ggsave(filename = "QC144.png", plot = rami1, width = 24, height = 12)

# Feature plot - visualize feature expression in low-dimensional space
rami1 <- FeaturePlot(BLCA_M9, features = features,
                     label = TRUE, repel = TRUE,
                     cols = c("lightgrey", "#943bcc"))
rami1 <- (rami1) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC145.png", plot = rami1, width = 18, height = 25)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
rami1 <- DotPlot(BLCA_M9, features = features,
                 cols = c("lightgrey", "#943bcc")) + RotatedAxis()
ggsave(filename = "QC146.png", plot = rami1, width = 10, height = 8)
ggsave(filename = "QC146.jpg", plot = rami1, width = 10, height = 8)

# metadata and pass it to the split.by argument
rami1 <- VlnPlot(BLCA_M9, features = "percent.mt", split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC147.png", plot = rami1, width = 14, height = 8)

rami1 <- VlnPlot(BLCA_M9, features = "TIMP1", split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC148.png", plot = rami1, width = 14, height = 8)
########################################
########################################
rami1 <- VlnPlot(BLCA_M9, features = features,split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC149.png", plot = rami1, width = 14, height = 8)
#######################################
#######################################
rami1 <- VlnPlot(BLCA_M9, features = "CD63", split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC150.png", plot = rami1, width = 14, height = 8)

rami1 <- VlnPlot(BLCA_M9, features = "MMP9", split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC151.png", plot = rami1, width = 14, height = 8)

rami1 <- VlnPlot(BLCA_M9, features = "MMP14", split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC152.png", plot = rami1, width = 14, height = 8)

rami1 <- VlnPlot(BLCA_M9, features = "IL6", split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC153.png", plot = rami1, width = 14, height = 8)

rami1 <- VlnPlot(BLCA_M9, features = "FN1", split.plot = TRUE,
                 split.by = "Batch")
ggsave(filename = "QC154.png", plot = rami1, width = 14, height = 8)
######################################
######################################
saveRDS(BLCA_M9, file = "BLCA_M10.RDS")
BLCA_M10 <- readRDS("BLCA_M10.RDS")
#######################################
#######################################