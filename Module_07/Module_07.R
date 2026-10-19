##################################################################
### Run Magic ...........
#########################
#########################
library(Rmagic)
library(phateR)
library(reticulate)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(Matrix)
library(RColorBrewer)
library(SeuratWrappers)

# pip list
# pip install magic-impute
# pip install phate
# pip install phate --break-system-packages

features = c("TIMP1", "CD63", "MMP9","MMP14", "IL6", "FN1")
BLCA_M7 <- readRDS("BLCA_M7.RDS")

# Why MAGIC?
FeatureScatter(BLCA_M7, 
               group.by = "Batch",
               feature1 = "TIMP1", 
               feature2 = "CD68")

DefaultAssay(BLCA_M7) <- "RNA"
BLCA_M7@assays$RNA$data

# Run MAGIC
data_MAGIC <- magic(t(BLCA_M7@assays$RNA$data), knn=15, seed = 42)

# Create Seurat with MAGIC data
BLCA_MAGIC <- SetAssayData(BLCA_M7, 
                           layer = "data", 
                           new.data = t(as.matrix(data_MAGIC$result)))

# Recalculate with MAGIC
BLCA_MAGIC <- BLCA_MAGIC %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.7) %>%
  RunUMAP(dims = 1:20, n.neighbors = 60) %>% 
  RunTSNE()
##############################################################################
# DimPlot MAGIC vs Normal
DimPlot(BLCA_MAGIC, group.by = "cell_type", 
        label = TRUE) + DimPlot(BLCA_M7, group.by = "Batch")

rami1 <- DimPlot(BLCA_MAGIC, group.by = "cell_type", 
                 label = TRUE, repel = TRUE,
                 label.size = 5) + labs(x = "UMAP_1",
                                        y = "UMAP_2", 
                                        title = "Cell_type with MAGIC") + NoLegend()
rami2 <- DimPlot(BLCA_M7, group.by = "cell_type", 
                 label = TRUE, repel = TRUE,
                 label.size = 6) + labs(x = "UMAP_1",
                                        y = "UMAP_2",
                                        title = "Cell_type without MAGIC") + NoLegend()
rami3 <- (rami1 + rami2) 
ggsave(filename = "QC97.png", plot = rami3, width = 17, height = 8)
#################################################
#################################################
rami1 <- DimPlot(BLCA_MAGIC, group.by = "Batch", 
                 label = FALSE, repel = TRUE,
                 label.size = 5) + labs(x = "UMAP_1",
                                        y = "UMAP_2", 
                                        title = "Batch with MAGIC")
rami2 <- DimPlot(BLCA_M7, group.by = "Batch", 
                 label = FALSE, repel = TRUE,
                 label.size = 6) + labs(x = "UMAP_1",
                                        y = "UMAP_2",
                                        title = "Batch without MAGIC")
rami3 <- (rami1 + rami2) 
ggsave(filename = "QC98.png", plot = rami3, width = 24, height = 10)
###################################################################
###################################################################
rami4 <- DimPlot(BLCA_MAGIC, group.by = "Batch", 
                 label = FALSE, reduction = "tsne",
                 label.size = 6) + labs(title = "Batch with MAGIC")
rami5 <- DimPlot(BLCA_M7, group.by = "Batch", 
                 label = FALSE, reduction = "tsne",
                 label.size = 6) + labs(title = "Batch without MAGIC")
rami6 <- (rami4 + rami5) 
ggsave(filename = "QC99.png", plot = rami6, width = 24, height = 10)
#################################################
#################################################
rami4 <- DimPlot(BLCA_MAGIC, group.by = "cell_type", 
                 label = TRUE, reduction = "tsne",
                 label.size = 4.5) + labs(title = "Cell_type with MAGIC") + NoLegend()
rami5 <- DimPlot(BLCA_M7, group.by = "cell_type", 
                 label = TRUE,reduction = "tsne",
                 label.size = 4.5) + labs(title = "Cell_type without MAGIC") + NoLegend()
rami6 <- (rami4 + rami5) 
ggsave(filename = "QC100.png", plot = rami6, width = 20, height = 10)
###############################################
# FeaturePlot of gene expression compared
features <- c("TIMP1", "CD63", "MMP9", "MMP14", "IL6", "FN1")

cols <- c("#afafed", "#803b58", "#2f6942")

fs1 <- FeatureScatter(BLCA_MAGIC, 
                      cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "CD63") + NoLegend()

fs2 <- FeatureScatter(BLCA_MAGIC, 
                      cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "MMP9") + NoLegend()

fs3 <- FeatureScatter(BLCA_MAGIC, 
                      cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "MMP14") + NoLegend()

fs4 <- FeatureScatter(BLCA_MAGIC,
                      cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "IL6") + NoLegend()

fs5 <- FeatureScatter(BLCA_MAGIC, 
                      cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "FN1")

fs7 <- fs1 | fs2 | fs3 | fs4 | fs5 
ggsave(filename = "QC101.png", plot = fs7, width = 30, height = 6)
#######################################
#######################################
fs1 <- FeatureScatter(BLCA_MAGIC, 
                      group.by = "cell_type",
                      feature1 = "TIMP1", 
                      feature2 = "CD63") + NoLegend()

fs2 <- FeatureScatter(BLCA_MAGIC, 
                      group.by = "cell_type",
                      feature1 = "TIMP1", 
                      feature2 = "MMP9") + NoLegend()

fs3 <- FeatureScatter(BLCA_MAGIC, 
                      group.by = "cell_type",
                      feature1 = "TIMP1", 
                      feature2 = "MMP14") + NoLegend()

fs4 <- FeatureScatter(BLCA_MAGIC,
                      group.by = "cell_type",
                      feature1 = "TIMP1", 
                      feature2 = "IL6") + NoLegend()

fs5 <- FeatureScatter(BLCA_MAGIC, 
                      group.by = "cell_type",
                      feature1 = "TIMP1", 
                      feature2 = "FN1")

fs7 <- fs1 | fs2 | fs3 | fs4 | fs5 
ggsave(filename = "QC102.png", plot = fs7, width = 30, height = 6)

Idents(BLCA_MAGIC) <- "cell_type"
#######################################################
#######################################################
fx1 <- FeatureScatter(BLCA_M7, 
                      #cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "CD63") + NoLegend()

fx2 <- FeatureScatter(BLCA_M7, 
                      #cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "MMP9") + NoLegend()

fx3 <- FeatureScatter(BLCA_M7, 
                      #cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "MMP14") + NoLegend()

fx4 <- FeatureScatter(BLCA_M7,
                      #cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "IL6") + NoLegend()

fx5 <- FeatureScatter(BLCA_M7, 
                      #cols = cols,
                      group.by = "Batch",
                      feature1 = "TIMP1", 
                      feature2 = "FN1") 

fs6 <- fx1 | fx2 | fx3 | fx4 | fx5
ggsave(filename = "QC103.png", plot = fs6, width = 25, height = 5)
#############################################
saveRDS(BLCA_MAGIC, file = "BLCA_MAGIC.RDS")

BLCA_M7 <- readRDS("BLCA_M7.RDS")
BLCA_MAGIC <- readRDS("BLCA_MAGIC.RDS")
#######################################
#######################################
## Running phate ..........
###########################
###########################
library(phateR)

# Run PHATE
BLCA_phate <- phate(t(BLCA_M7@assays$RNA$data))

# Get Embedding
emd <- BLCA_phate$embedding

# Add into SRT
BLCA_M7$PHATE1 <- emd[,c(1)]
BLCA_M7$PHATE2 <- emd[,c(2)]

# Plot the PHATE embedding
p1 <- DimPlot(BLCA_M7, group.by = "cell_type") + NoLegend() 
p2 <- FeatureScatter(BLCA_M7, feature1 = "PHATE1",
                     group.by = "cell_type",
                     plot.cor = FALSE, 
                     feature2 = "PHATE2") 

p3 <- DimPlot(BLCA_M7, group.by = "Batch") + NoLegend() 
p4 <- FeatureScatter(BLCA_M7, feature1 = "PHATE1", plot.cor = FALSE, 
                     group.by = "Batch",
                     feature2 = "PHATE2")

p5 <- p1+p2+p3+p4
ggsave(filename = "QC104.png", plot = p5, width = 15, height = 10)
#############################
#############################
saveRDS(BLCA_M7, file = "BLCA_M8.RDS")
BLCA_M8 <- readRDS("BLCA_M8.RDS")
######################################
######################################
p1 <- DimPlot(BLCA_MAGIC, group.by = "cell_type",
              label = FALSE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2",
                                        title = "MAGIC clustered by cell_type")

p2 <- FeatureScatter(BLCA_M7, feature1 = "PHATE1",
                     plot.cor = FALSE, group.by = "cell_type",
                     feature2 = "PHATE2") & labs(x = "PHATE1", 
                                                 y = "PHATE2",
                                                 title = "PHATE clustered by cell_type")

p3 <- DimPlot(BLCA_MAGIC, group.by = "Batch", 
              label = FALSE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2", 
                                          title = "MAGIC clustered by Batch")

p4 <- FeatureScatter(BLCA_M7, feature1 = "PHATE1", plot.cor = FALSE, 
                     group.by = "Batch",
                     feature2 = "PHATE2") & labs(x = "PHATE1", y = "PHATE2",
                                                 title = "PHATE clustered by Batch")

p5 <- p1+p2+p3+p4
ggsave(filename = "QC105.png", plot = p5, width = 14, height = 10)
###########################################
###########################################
p1 <- DimPlot(BLCA_M7, group.by = "cell_type", repel = TRUE,
              label = FALSE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2",
                                                 title = "UMAP clustered by cell_type")

p2 <- FeatureScatter(BLCA_M7, feature1 = "PHATE1",
                     plot.cor = FALSE, group.by = "cell_type",
                     feature2 = "PHATE2") & labs(x = "PHATE1", 
                                                 y = "PHATE2",
                                                 title = "PHATE clustered by cell_type")

p3 <- DimPlot(BLCA_M7, group.by = "Batch", 
              label = FALSE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2", 
                                                 title = "UMAP clustered by Sample")

p4 <- FeatureScatter(BLCA_M7, feature1 = "PHATE1", plot.cor = FALSE, 
                     group.by = "Batch",
                     feature2 = "PHATE2") & labs(x = "PHATE1", y = "PHATE2",
                                                 title = "PHATE clustered by Sample")

p5 <- p1+p2+p3+p4
ggsave(filename = "QC106.png", plot = p5, width = 14, height = 10)
#####################################################
p1 <- DimPlot(BLCA_M7, group.by = "cell_type", repel = TRUE, label.size = 4,
              label = TRUE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2",
                                                 title = "UMAP clustered by cell_type")

p2 <- DimPlot(BLCA_MAGIC, group.by = "cell_type", repel = TRUE,label.size = 4,
              label = TRUE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2",
                                                 title = "MAGIC clustered by cell_type")

p3 <- DimPlot(BLCA_M7, group.by = "Batch", label.size = 3.5,
              label = FALSE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2", 
                                                 title = "UMAP clustered by Batch")

p4 <- DimPlot(BLCA_MAGIC, group.by = "Batch", label.size = 3.5,
              label = FALSE) + NoLegend() & labs(x = "UMAP_1", y = "UMAP_2", 
                                                 title = "MAGIC clustered by Batch")

p5 <- p1+p2+p3+p4
ggsave(filename = "QC107.png", plot = p5, width = 15, height = 12)
############################################
features = c("TIMP1", "CD63", "MMP9", 
             "MMP14", "IL6", "FN1")

DimPlot(BLCA_M7, reduction = "umap", label = TRUE, group.by = "orig.ident")

plot12 <- FeaturePlot(BLCA_M7, features = features, 
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "umap")

plot14 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC108.png", plot = plot14, width = 10, height = 14)

plot13 <- FeaturePlot(BLCA_M7, features = features,
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "tsne")
ggsave(filename = "QC109.png", plot = plot13, width = 12, height = 16)
###########################################################################
plot12 <- FeaturePlot(BLCA_MAGIC, features = features,
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "umap")

plot14 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC110.png", plot = plot14, width = 10, height = 14)

plot13 <- FeaturePlot(BLCA_MAGIC, features = features,
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "tsne")
ggsave(filename = "QC111.png", plot = plot13, width = 12, height = 16)
############################################################################
Idents(BLCA_MAGIC) <- "cell_type"
plot12 <- FeaturePlot(BLCA_MAGIC, features = features,
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "umap")

plot14 <- (plot12) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC112.png", plot = plot14, width = 10, height = 14)

plot13 <- FeaturePlot(BLCA_MAGIC, features = features,
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "tsne")
ggsave(filename = "QC113.png", plot = plot13, width = 12, height = 16)
############################################################################
plot10 <- DimPlot(object = BLCA_MAGIC, reduction = "tsne", 
                  group.by = "cell_type", 
                  label = TRUE, repel = TRUE, label.size = 4) + NoLegend()

plot11 <- DimPlot(object = BLCA_MAGIC, reduction = "tsne", 
                  group.by = "Batch", 
                  label = TRUE, repel = TRUE, label.size = 4) + NoLegend()


plot14 <- plot10 + plot13
plot15 <- plot11 + plot13

ggsave(filename = "QC114.png", plot = plot14, width = 18, height = 10)
####################################
####################################
plot10 <- DimPlot(object = BLCA_MAGIC, reduction = "umap", 
                  group.by = "cell_type", 
                  label = TRUE, repel = TRUE, label.size = 4) + NoLegend()

plot13 <- FeaturePlot(BLCA_MAGIC, features = features,
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "umap")

plot14 <- (plot10 + plot13) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC115.png", plot = plot14, width = 20, height = 12)
######################################################################
plot9 <- DimPlot(BLCA_M7, reduction = "umap", label = TRUE, 
                 repel = TRUE, group.by = "cell_type", label.size = 6) + NoLegend()

plot10 <- (plot9) & labs(x = "UMAP_1", y = "UMAP_2")
plot16 <- FeaturePlot(BLCA_M7, features = features, 
                      label = TRUE, repel = TRUE,
                      cols = c("lightgrey", "#943bcc"), reduction = "umap")

plot18 <- (plot16) & labs(x = "UMAP_1", y = "UMAP_2")
plot17 <- plot10 + plot18
ggsave(filename = "QC116.png", plot = plot17, width = 24, height = 14)
#####################################################################
plot10 <- DimPlot(object = BLCA_MAGIC, reduction = "umap", 
                  group.by = "cell_type", 
                  label = TRUE, repel = TRUE, label.size = 4) + NoLegend()

plot11 <- DimPlot(object = BLCA_MAGIC, reduction = "umap", 
                  group.by = "Batch", 
                  label = TRUE, repel = TRUE, label.size = 4) + NoLegend()

plot13 <- FeaturePlot(BLCA_MAGIC, features = features,
                      label = TRUE, repel = TRUE, label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "umap")

plot14 <- (plot10 + plot13) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC117.png", plot = plot14, width = 18, height = 10)
#####################################################################
plot9 <- DimPlot(BLCA_M7, reduction = "umap", label = FALSE, 
                 repel = TRUE, group.by = "Batch", label.size = 6)

plot10 <- (plot9) & labs(x = "UMAP_1", y = "UMAP_2")
plot16 <- FeaturePlot(BLCA_M7, features = features, 
                      label = TRUE, repel = TRUE,
                      cols = c("lightgrey", "#943bcc"), reduction = "umap")

plot18 <- (plot16) & labs(x = "UMAP_1", y = "UMAP_2")
plot17 <- plot10 + plot18
ggsave(filename = "QC118.png", plot = plot17, width = 24, height = 14)
#######################################################################
plot10 <- DimPlot(object = BLCA_MAGIC, reduction = "umap", 
                  group.by = "Batch", 
                  label = FALSE, repel = TRUE, label.size = 4) 

plot13 <- FeaturePlot(BLCA_MAGIC, features = features,
                      label = TRUE, repel = TRUE,label.size = 3,
                      cols = c('lightgrey', '#943bcc'), reduction = "umap")

plot14 <- (plot10 + plot13) & labs(x = "UMAP_1", y = "UMAP_2")

ggsave(filename = "QC119.png", plot = plot14, width = 18, height = 10)
########################################################################
plot9 <- DimPlot(BLCA_M7, reduction = "tsne", label = TRUE, 
                 repel = TRUE, group.by = "cell_type", label.size = 5) + NoLegend()

plot10 <- (plot9) & labs(x = "UMAP_1", y = "UMAP_2")
plot16 <- FeaturePlot(BLCA_M7, features = features, 
                      label = TRUE, repel = TRUE, label.size = 3, 
                      cols = c("lightgrey", "#943bcc"), reduction = "tsne")

plot18 <- (plot16) & labs(x = "UMAP_1", y = "UMAP_2")
plot17 <- plot10 + plot18
ggsave(filename = "QC120.png", plot = plot17, width = 24, height = 14)
########################################################################
########################################################################
saveRDS(BLCA_M7, file = "BLCA_M8.RDS")
BLCA_M8 <- readRDS("BLCA_M8.RDS")
######################################
######################################