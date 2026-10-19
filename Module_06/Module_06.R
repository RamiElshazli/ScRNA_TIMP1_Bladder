# Install and load required packages:

library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(limma)
library(SingleR)
library(pheatmap)
library(monocle3)
library(celldex)
##################################
# saveRDS(BLCA_M4, file = "BLCA_M5.RDS")
BLCA_M6 <- readRDS("BLCA_M6.RDS")

##################################
### Step 6: Automated classification using celldex and singleR ###
# Predict using the MonacoImmuneData from celldex
ref <- celldex::HumanPrimaryCellAtlasData()

# ref <- celldex::DatabaseImmuneCellExpressionData()
# ref <- celldex::HumanPrimaryCellAtlasData()
# ref <- celldex::MonacoImmuneData()

# What is the structure of ref?
ref
View(ref)
View(as.data.frame(ref@colData))
View(as.data.frame(colData(ref)))
colData(ref)
##################################
##################################
### Method 01 ........
######################
# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the data-set ...
BLCA_M5_counts <- GetAssayData(BLCA_M5, layer = 'counts')
BLCA_M5_counts <- as.SingleCellExperiment(BLCA_M5)

pred <- SingleR(test = BLCA_M5_counts,
                ref = ref,
                labels = ref$label.main)

pred
OC_M5$singleR.labels <- pred$labels[match(rownames(OC_M5@meta.data), 
                                            rownames(pred))]
DimPlot(OC_M5, reduction = 'umap', 
        group.by = 'singleR.labels')
#########################################################################
#########################################################################
### Method 02 ........
######################
# run SingleR (default mode) ---------
DefaultAssay(BLCA_M6) <- "RNA"

###################################
###################################
# 1. Predict using the Broad labels
BLCA_M6_counts <- as.SingleCellExperiment(BLCA_M6)
results1 <- SingleR(test = BLCA_M6_counts, 
                    ref = ref, 
                    labels = ref$label.main)

plotScoreHeatmap(results1,
                 max.labels = 25,
                 clusters = BLCA_M6_counts$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)
######################################################
######################################################
rami1 <- plotScoreHeatmap(results1, max.labels = 25,
                          clusters = BLCA_M6_counts$seurat_clusters,
                          order.by = "clusters", show_colnames = FALSE)

ggsave(filename = "QC65.png", plot = rami1, width = 20, height = 15)
##########################################################################
##########################################################################
# 2. Predict using the Fine-grained labels
results2 <- SingleR(test = BLCA_M6_counts, 
                    ref = ref, 
                    labels = ref$label.fine)

plotScoreHeatmap(results2,
                 max.labels = 25,
                 clusters = BLCA_M6_counts$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)
#####################################################
#####################################################
rami2 <- plotScoreHeatmap(results2, max.labels = 25,
                          clusters = BLCA_M6_counts$seurat_clusters,
                          order.by = "clusters", show_colnames = FALSE)

ggsave(filename = "QC66.png", plot = rami2, width = 30, height = 45)
############################################################################
# 3. Combine the labels for improved accuracy
### combineRecomputedResults .....
results3 <- combineCommonResults(list("Broad" = results1,
                                      "Fine" = results2))
plotScoreHeatmap(results3,
                 max.labels = 25,
                 clusters = BLCA_M6_counts$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)

rami3 <- plotScoreHeatmap(results3,
                          max.labels = 25, 
                          clusters = BLCA_M6_counts$seurat_clusters,
                          order.by = "clusters", show_colnames = FALSE)

ggsave(filename = "QC67.png", plot = rami3, width = 20, height = 15)
#################################################################
#################################################################
# Add the predicted labels and compare to clusters
# You are looking for clusters which have unambiguous labels
library(pheatmap)
BLCA_M6$predicted_id <- results3$pruned.labels

tab1 <- table(BLCA_M6$predicted_id, BLCA_M6$seurat_clusters)
rami33 <- pheatmap(log10(tab1+10), color = colorRampPalette(c('grey','darkred'))(10))
ggsave(filename = "QC68.png", plot = rami33, width = 18, height = 20)

table(tab1)
View(as.data.frame(tab1))
write.csv(tab1,"filename.csv", row.names = TRUE)

# Select final IDs
Idents(BLCA_M6) <- "seurat_clusters"

new.cluster.ids <- c("Epithelial_cells","Epithelial_cells","Endothelial_cells",
                     "T_cells","Epithelial_cells","Macrophage","Tissue_stem_cells",
                     "Epithelial_cells","T_cells","Neutrophils","Neutrophils",
                     "Fibroblasts","Neurons","Epithelial_cells","Epithelial_cells",
                     "Epithelial_cells","Epithelial_cells","Neutrophils","Fibroblasts",
                     "Monocytes","NK_cells","Epithelial_cells","Monocytes",
                     "Epithelial_cells","Fibroblasts","CMP","B_cells","Epithelial_cells")

names(new.cluster.ids) <- levels(BLCA_M6)

BLCA_M6 <- RenameIdents(BLCA_M6, new.cluster.ids)
BLCA_M6$cell_type <- Idents(BLCA_M6)
DimPlot(BLCA_M6)

rami1 <- DimPlot(BLCA_M6, label.size = 4.5, 
                 label = TRUE, repel = TRUE) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC69.png", plot = rami1, width = 10, height = 8)

plot13 <- FeaturePlot(BLCA_M6, features = "TIMP1", label = TRUE, 
                      repel = TRUE, label.size = 4.5, label.color = "black", 
                      alpha = 0.7, cols = c('lightgrey', 
                                            '#943bcc')) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC70.png", plot = plot13, width = 10, height = 8)
################################################################
################################################################
saveRDS(BLCA_M6, file = "BLCA_M7.RDS")
BLCA_M6 <- readRDS("BLCA_M7.RDS")
#######################################
#######################################
DimPlot(BLCA_M6, label.size = 4.5, label = TRUE) + NoLegend()
FeaturePlot(BLCA_M6, features = "TIMP1", label = TRUE,
            repel = TRUE,label.size = 4.5,
            cols = c('lightgrey', '#943bcc'))

BLCA_M6$singlr_labels <- BLCA_M6$cell_type
##########################################
plot6 <- DimPlot(BLCA_M6, reduction = "umap", 
                 split.by = "Batch",
                 label = TRUE, 
                 repel = TRUE, 
                 label.size = 4) & labs(x = "UMAP_1", y = "UMAP_2") & NoLegend()
ggsave(filename = "QC71.png", plot = plot6, width = 16, height = 7)

plot6 <- DimPlot(BLCA_M6, reduction = "umap", 
                 label = TRUE, 
                 repel = TRUE, 
                 label.size = 4) & labs(x = "UMAP_1", y = "UMAP_2") & NoLegend()
ggsave(filename = "QC72.png", plot = plot6, width = 9, height = 8)

plot6 <- DimPlot(BLCA_M6, reduction = "umap", 
                 label.size = 4) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC73.png", plot = plot6, width = 10, height = 8)
#############################################
features = c("TIMP1", "CD63", "MMP9", 
             "MMP14", "IL6", "FN1")

plot6 <- FeaturePlot(BLCA_M6, 
                     features = "TIMP1", label = TRUE,
                     cols = c('lightgrey', '#943bcc')) & labs(x = "UMAP_1", y = "UMAP_2") & NoLegend()
ggsave(filename = "QC74.png", plot = plot6, width = 8, height = 8)

plot6 <- FeaturePlot(BLCA_M6, 
                     features = features, label = TRUE, repel = TRUE,
                     label.size = 3.5,
                     cols = c('lightgrey', '#943bcc')) & labs(x = "UMAP_1", y = "UMAP_2") & NoLegend()
ggsave(filename = "QC75.png", plot = plot6, width = 14, height = 20)


plot6 <- FeaturePlot(BLCA_M6, 
                     features = "TIMP1", label = TRUE,
                     min.cutoff = "q1", max.cutoff = "q10",
                     label.size = 3.5, 
                     repel = TRUE,
                     split.by = "Batch",
                     cols = c('lightgrey', '#943bcc')) & labs(x = "UMAP_1", y = "UMAP_2") & NoLegend()
ggsave(filename = "QC76.png", plot = plot6, width = 20, height = 8)

plot6 <- FeaturePlot(BLCA_M6, 
                     features = "TIMP1", label = TRUE,
                     label.size = 3.5, 
                     repel = TRUE,
                     split.by = "Batch",
                     cols = c('lightgrey', '#943bcc')) & labs(x = "UMAP_1", y = "UMAP_2") & NoLegend()
ggsave(filename = "QC77.png", plot = plot6, width = 20, height = 8)
##########################################
plot6 <- DimPlot(BLCA_M6, reduction = "tsne", 
                 split.by = "Batch",
                 label = TRUE, 
                 repel = TRUE, 
                 label.size = 4) & NoLegend()
ggsave(filename = "QC78.png", plot = plot6, width = 18, height = 7)

plot6 <- DimPlot(BLCA_M6, reduction = "tsne", 
                 label = TRUE, 
                 repel = TRUE, 
                 label.size = 5) & NoLegend()
ggsave(filename = "QC79.png", plot = plot6, width = 9, height = 8)

plot6 <- DimPlot(BLCA_M6, reduction = "tsne", 
                 label.size = 4)
ggsave(filename = "QC80.png", plot = plot6, width = 10, height = 8)
##########################################
plot6 <- FeaturePlot(BLCA_M6, reduction = "tsne",
                     features = "TIMP1", label = TRUE,
                     cols = c('lightgrey', '#943bcc')) & NoLegend()
ggsave(filename = "QC81.png", plot = plot6, width = 8, height = 8)

plot6 <- FeaturePlot(BLCA_M6, reduction = "tsne",
                     features = features, label = TRUE, repel = TRUE,
                     label.size = 3.5,
                     cols = c('lightgrey', '#943bcc')) & NoLegend()
ggsave(filename = "QC82.png", plot = plot6, width = 14, height = 20)


plot6 <- FeaturePlot(BLCA_M6, reduction = "tsne",
                     features = "TIMP1", label = TRUE,
                     label.size = 3.5,
                     repel = TRUE,
                     split.by = "Batch",
                     cols = c('lightgrey', '#943bcc')) & NoLegend()
ggsave(filename = "QC83.png", plot = plot6, width = 20, height = 7)
##########################################
##########################################
rami1 <- VlnPlot(BLCA_M6, features = "TIMP1", 
                 layer = "counts", log = TRUE, group.by = "Batch")
ggsave(filename = "QC84.png", plot = rami1, width = 15, height = 9)

rami1 <- VlnPlot(BLCA_M6, features = features, 
                 layer = "counts", log = TRUE, group.by = "Batch")
ggsave(filename = "QC85.png", plot = rami1, width = 15, height = 9)
############################################
rami2 <- VlnPlot(BLCA_M6, features = "TIMP1", 
                 layer = "counts", log = TRUE)
ggsave(filename = "QC86.png", plot = rami2, width = 15, height = 8)

rami2 <- VlnPlot(BLCA_M6, features = features, 
                 layer = "counts", log = TRUE)
ggsave(filename = "QC87.png", plot = rami2, width = 18, height = 8)
#################################################
rami1 <- RidgePlot(BLCA_M6, features = "TIMP1", ncol = 1) & NoLegend()
ggsave(filename = "QC88.png", plot = rami1, width = 9, height = 8)

rami1 <- DotPlot(BLCA_M6, features = features) + RotatedAxis()
ggsave(filename = "QC89.png", plot = rami1, width = 8, height = 6)
ggsave(filename = "QC89.jpg", plot = rami1, width = 8, height = 6)
###################################################
# Single cell heatmap of feature expression
DoHeatmap(subset(BLCA_M6, downsample = 100), features = features, size = 3, 
          group.by = "singlr_labels")

features = c("TIMP1", "CD63", "MMP9", 
             "MMP14", "IL6", "FN1")

p1 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "CD63") + NoLegend()
p2 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "MMP9") + NoLegend()
p3 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "MMP14") + NoLegend()
p4 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "IL6") + NoLegend()
p5 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "FN1")
p6 <- p1|p2|p3|p4|p5

ggsave(filename = "QC90.png", plot = p6, width = 30, height = 6)
############################################
############################################
p0 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "CD63") + NoLegend()
p1 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "CD63") 
p2 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "MMP9") 
p3 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "MMP14") 
p4 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "IL6") 
p5 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "FN1")

ggsave(filename = "QC91.png", plot = p1, width = 10, height = 8)
ggsave(filename = "QC92.png", plot = p2, width = 10, height = 8)
ggsave(filename = "QC93.png", plot = p3, width = 10, height = 8)
ggsave(filename = "QC94.png", plot = p4, width = 10, height = 8)
ggsave(filename = "QC95.png", plot = p5, width = 10, height = 8)
##########################################
p1 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "CD63", 
                     group.by = "Batch") + NoLegend()
p2 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "MMP9",
                     group.by = "Batch") + NoLegend()
p3 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "MMP14",
                     group.by = "Batch") + NoLegend()
p4 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "IL6",
                     group.by = "Batch") + NoLegend()
p5 <- FeatureScatter(BLCA_M6, feature1 = "TIMP1", feature2 = "FN1", group.by = "Batch")
p6 <- p1|p2|p3|p4|p5

ggsave(filename = "QC96.png", plot = p6, width = 30, height = 6)
##########################################

saveRDS(BLCA_M6, file = "BLCA_M7.RDS")
BLCA_M7 <- readRDS("BLCA_M7.RDS")
##########################################
## convert seurat into H5ad python
library(SeuratData)
library(SeuratDisk)
#SaveH5Seurat(BLCA_M7, filename = "BLCA_M7.h5Seurat")
#Convert("BLCA_M7.h5Seurat", dest = "h5ad")
##########################################
##########################################