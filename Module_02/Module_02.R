# Load required packages:

library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(limma)
library(cowplot)
############################################
# Read our Seurat object 
BLCA_M1 <- readRDS("BLCA_M1.RDS")
BLCA_M1
####################################
## 2. Normalization and Identification of Highly Variable Genes
####################################
View(BLCA_M1@meta.data)
range(BLCA_M1$nFeature_RNA)
range(BLCA_M1$nCount_RNA)
range(BLCA_M1$Batch)
##########################
##########################
table(BLCA_M1$orig.ident, BLCA_M1$Batch)
##########################
Idents(BLCA_M1) <- "orig.ident"

## free unused memory
gc()

##
class(BLCA_M1[["RNA"]])
#######################################
######## Data normalization (Normalize.Data)
# Log-transform the counts = Log-normalize the data
BLCA_M2 <- NormalizeData(object = BLCA_M1,
                         normalization.method = "LogNormalize", # Standard log normalization
                         scale.factor = 10000,                  # Default (10,000 transcripts)
                         verbose = TRUE)                        # Print progress messages

# Find Variable Features
?FindVariableFeatures

# BLCA_M2 <- FindVariableFeatures(BLCA_M2)
# Identify highly variable features using Variance-stabilizing transformation 'vst' method
BLCA_M2 <- FindVariableFeatures(BLCA_M2, 
                                selection.method = "vst", nfeatures = 3000,
                                verbose = TRUE)

# Extract the most top 10 highly variable genes
top10 <- head(VariableFeatures(BLCA_M2), 60)
top10

# plot variable features with and without labels
# Visualize Highly Variable Features
plot1 <- VariableFeaturePlot(BLCA_M2)
plot1
ggsave(filename = "QC9.png", plot = plot1, width = 12, height = 8.5)

# plot variable features with and without labels
# Label the top 10 variable genes
plot1 <- VariableFeaturePlot(BLCA_M2, cols = c("grey", "#6a5cd1")) + NoLegend()
plot1
ggsave(filename = "QC10.png", plot = plot1, width = 12, height = 8.5)
###########################################
###########################################
plot2 <- LabelPoints(plot = plot1, points = top10, 
                     #color = "brown", 
                     repel = TRUE, xnudge = 0, ynudge = 0)
plot2
ggsave(filename = "QC11.png", plot = plot2, width = 12, height = 8.5)

plot3 <- plot1 + plot2
ggsave(filename = "QC12.png", plot = plot3, width = 18, height = 8)
###########################################
###########################################
### Scale the Data
# Scale only the variable features to reduce noise and improve computation time
#BLCA_M3 <- ScaleData(BLCA_M3)
BLCA_M2 <- ScaleData(BLCA_M2, 
                     features = VariableFeatures(object = BLCA_M2))
###########################################
plot5 <- FeatureScatter(BLCA_M2, feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA", raster = FALSE)

ggsave(filename = "QC13.png", plot = plot5, width = 12, height = 8.5)

plot6 <- FeatureScatter(BLCA_M2, feature1 = "nCount_RNA",
                        feature2 = "percent.mt", raster = FALSE)

ggsave(filename = "QC14.png", plot = plot6, width = 12, height = 8.5)

plot7 <- FeatureScatter(BLCA_M2, feature1 = "nCount_RNA",
                        feature2 = "percent.rb", raster = FALSE)

ggsave(filename = "QC15.png", plot = plot7, width = 12, height = 8.5)

plot8 <- plot5 + plot6 + plot7
ggsave(filename = "QC16.png", plot = plot8, width = 40, height = 10)
#############################################

# Save our Seurat object
saveRDS(BLCA_M2, file = "BLCA_M2.RDS")
BLCA_M3 <- readRDS("BLCA_M2.RDS")
#######################################
#######################################