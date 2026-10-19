# Install and load required packages:
######################################
install.packages("Seurat")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"), force = TRUE)
BiocManager::install("limma")

install.packages("ggpubr")
install.packages("tidyverse")
############################################
# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)

options(timeout=9999999)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)


BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'), force = TRUE)

remotes::install_github('cole-trapnell-lab/monocle3')

library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(Matrix)
library(SeuratWrappers)
library(limma)
library(monocle3)
#######################################
getwd()
setwd()
#################################################
## 1. Loading data from different sources .....
#################################################
# Load the BLCA dataset [Male, Age late 50s] Bladder BC001 .............
P1_X10Genomics <- Read10X_h5("datasets/4plex_Bladder_matrix.h5",
                                 use.names = TRUE,
                                 unique.features = TRUE)

## GSE130001 [two sample] https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130001
P2_GSE130001 <- Read10X(data.dir = "datasets/GSE130001_RAW/sample1",
                              unique.features = TRUE)
P3_GSE130001 <- Read10X(data.dir = "datasets/GSE130001_RAW/sample2",
                              unique.features = TRUE)

## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129845 [three]
N1_GSE129845 <- Read10X(data.dir = "datasets/GSE129845_RAW/sample1",
                              unique.features = TRUE)
N2_GSE129845 <- Read10X(data.dir = "datasets/GSE129845_RAW/sample2",
                              unique.features = TRUE)
N3_GSE129845 <- Read10X(data.dir = "datasets/GSE129845_RAW/sample3",
                              unique.features = TRUE)
#############################################################################################
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267718
U1_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample1",
                            unique.features = TRUE)
U2_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample2",
                            unique.features = TRUE)
U3_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample3",
                            unique.features = TRUE)
U4_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample4",
                            unique.features = TRUE)
U5_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample5",
                            unique.features = TRUE)
U6_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample6",
                            unique.features = TRUE)
U7_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample9",
                            unique.features = TRUE)
U8_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/urine-derived cells/sample10",
                            unique.features = TRUE)
###############################################################################################
P4_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/tumor/sample1",
                              unique.features = TRUE)
P5_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/tumor/sample2",
                              unique.features = TRUE)
P6_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/tumor/sample3",
                              unique.features = TRUE)
P7_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/tumor/sample4",
                              unique.features = TRUE)
P8_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/tumor/sample5",
                              unique.features = TRUE)
P9_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/tumor/sample8",
                              unique.features = TRUE)
P10_GSE267718 <- Read10X(data.dir = "datasets/GSE267718_RAW2/tumor/sample9",
                               unique.features = TRUE)
###################################################################################
## cleanup environment to free up memory .............
## rm(list = ls()); gc()
##########################################
# Replace underscores with dashes in feature names
rownames(P1_X10Genomics) <- gsub("_", "-", rownames(P1_X10Genomics)) 
rownames(P2_GSE130001) <- gsub("_", "-", rownames(P2_GSE130001))
rownames(P3_GSE130001) <- gsub("_", "-", rownames(P3_GSE130001))
rownames(P4_GSE267718) <- gsub("_", "-", rownames(P4_GSE267718))
rownames(P5_GSE267718) <- gsub("_", "-", rownames(P5_GSE267718))
rownames(P6_GSE267718) <- gsub("_", "-", rownames(P6_GSE267718))
rownames(P7_GSE267718) <- gsub("_", "-", rownames(P7_GSE267718))
rownames(P8_GSE267718) <- gsub("_", "-", rownames(P8_GSE267718))
rownames(P9_GSE267718) <- gsub("_", "-", rownames(P9_GSE267718))
rownames(P10_GSE267718) <- gsub("_", "-", rownames(P10_GSE267718))
##########################################################
rownames(N1_GSE129845) <- gsub("_", "-", rownames(N1_GSE129845))
rownames(N2_GSE129845) <- gsub("_", "-", rownames(N2_GSE129845))
rownames(N3_GSE129845) <- gsub("_", "-", rownames(N3_GSE129845))
###########################################################
rownames(U1_GSE267718) <- gsub("_", "-", rownames(U1_GSE267718))
rownames(U2_GSE267718) <- gsub("_", "-", rownames(U2_GSE267718))
rownames(U3_GSE267718) <- gsub("_", "-", rownames(U3_GSE267718))
rownames(U4_GSE267718) <- gsub("_", "-", rownames(U4_GSE267718))
rownames(U5_GSE267718) <- gsub("_", "-", rownames(U5_GSE267718))
rownames(U6_GSE267718) <- gsub("_", "-", rownames(U6_GSE267718))
rownames(U7_GSE267718) <- gsub("_", "-", rownames(U7_GSE267718))
rownames(U8_GSE267718) <- gsub("_", "-", rownames(U8_GSE267718))
###########################################################
# Initialize the Seurat object with the raw (non-normalized data).
# Create a Seurat object with basic filtering for low-quality cells
## Make a function .........

seurat_M <- function(seurat, project_name = "my_seurat_project") {
  CreateSeuratObject(counts = seurat, 
                     project = project_name,
                     min.cells = 3, min.features = 200)
}
#######################################################################
P1_X10Genomics <- seurat_M(P1_X10Genomics, project_name = "P1-X10Genomics")
P2_GSE130001 <- seurat_M(P2_GSE130001, project_name = "P2-GSE130001")
P3_GSE130001 <- seurat_M(P3_GSE130001, project_name = "P3-GSE130001")
P4_GSE267718 <- seurat_M(P4_GSE267718, project_name = "P4-GSE267718")
P5_GSE267718 <- seurat_M(P5_GSE267718, project_name = "P5-GSE267718")
P6_GSE267718 <- seurat_M(P6_GSE267718, project_name = "P6-GSE267718")
P7_GSE267718 <- seurat_M(P7_GSE267718, project_name = "P7-GSE267718")
P8_GSE267718 <- seurat_M(P8_GSE267718, project_name = "P8-GSE267718")
P9_GSE267718 <- seurat_M(P9_GSE267718, project_name = "P9-GSE267718")
P10_GSE267718 <- seurat_M(P10_GSE267718, project_name = "P10-GSE267718")
######################################
######################################
N1_GSE129845 <- seurat_M(N1_GSE129845, project_name = "N1-GSE129845")
N2_GSE129845 <- seurat_M(N2_GSE129845, project_name = "N2-GSE129845")
N3_GSE129845 <- seurat_M(N3_GSE129845, project_name = "N3-GSE129845")
#######################################
#######################################
U1_GSE267718 <- seurat_M(U1_GSE267718, project_name = "U1-GSE267718")
U2_GSE267718 <- seurat_M(U2_GSE267718, project_name = "U2-GSE267718")
U3_GSE267718 <- seurat_M(U3_GSE267718, project_name = "U3-GSE267718")
U4_GSE267718 <- seurat_M(U4_GSE267718, project_name = "U4-GSE267718")
U5_GSE267718 <- seurat_M(U5_GSE267718, project_name = "U5-GSE267718")
U6_GSE267718 <- seurat_M(U6_GSE267718, project_name = "U6-GSE267718")
U7_GSE267718 <- seurat_M(U7_GSE267718, project_name = "U7-GSE267718")
U8_GSE267718 <- seurat_M(U8_GSE267718, project_name = "U8-GSE267718")
########################################
########################################
## Add new column called GSE_ID
P1_X10Genomics$GSE_ID <- "BLCA-P1"
P2_GSE130001$GSE_ID <- "BLCA-P2"
P3_GSE130001$GSE_ID <- "BLCA-P3"
P4_GSE267718$GSE_ID <- "BLCA-P4"
P5_GSE267718$GSE_ID <- "BLCA-P5" 
##############################################
P6_GSE267718$GSE_ID <- "BLCA-P6"
P7_GSE267718$GSE_ID <- "BLCA-P7"
P8_GSE267718$GSE_ID <- "BLCA-P8"
P9_GSE267718$GSE_ID <- "BLCA-P9"
P10_GSE267718$GSE_ID <- "BLCA-P10"
##################################################
N1_GSE129845$GSE_ID <- "BLCA-N1"
N2_GSE129845$GSE_ID <- "BLCA-N2"
N3_GSE129845$GSE_ID <- "BLCA-N3"
##################################################
U1_GSE267718$GSE_ID <- "BLCA-U1"
U2_GSE267718$GSE_ID <- "BLCA-U2"
U3_GSE267718$GSE_ID <- "BLCA-U3"
U4_GSE267718$GSE_ID <- "BLCA-U4"
U5_GSE267718$GSE_ID <- "BLCA-U5"
U6_GSE267718$GSE_ID <- "BLCA-U6"
U7_GSE267718$GSE_ID <- "BLCA-U7"
U8_GSE267718$GSE_ID <- "BLCA-U8"
##################################################
#################################
## Add new column called Batch
P1_X10Genomics$Batch <- "Bladder-carcinoma"
P2_GSE130001$Batch <- "Bladder-carcinoma"
P3_GSE130001$Batch <- "Bladder-carcinoma"
P4_GSE267718$Batch <- "Bladder-carcinoma"
P5_GSE267718$Batch <- "Bladder-carcinoma" 
##############################################
P6_GSE267718$Batch <- "Bladder-carcinoma"
P7_GSE267718$Batch <- "Bladder-carcinoma"
P8_GSE267718$Batch <- "Bladder-carcinoma"
P9_GSE267718$Batch <- "Bladder-carcinoma"
P10_GSE267718$Batch <- "Bladder-carcinoma"
##################################################
N1_GSE129845$Batch <- "Normal-Bladder tissues"
N2_GSE129845$Batch <- "Normal-Bladder tissues"
N3_GSE129845$Batch <- "Normal-Bladder tissues"
##################################################
U1_GSE267718$Batch <- "Urine-derived cells"
U2_GSE267718$Batch <- "Urine-derived cells"
U3_GSE267718$Batch <- "Urine-derived cells"
U4_GSE267718$Batch <- "Urine-derived cells"
U5_GSE267718$Batch <- "Urine-derived cells"
U6_GSE267718$Batch <- "Urine-derived cells"
U7_GSE267718$Batch <- "Urine-derived cells"
U8_GSE267718$Batch <- "Urine-derived cells"
##################################################
########## optional ....... ######################
##################################################
# Merge Seurat Objects
BLCA_M <- merge(x = P1_X10Genomics, 
                y = list(P2_GSE130001, P3_GSE130001, P4_GSE267718, 
                         P5_GSE267718, P6_GSE267718, P7_GSE267718, 
                         P8_GSE267718, P9_GSE267718, P10_GSE267718,
                         N1_GSE129845, N2_GSE129845, N3_GSE129845, 
                         U1_GSE267718, U2_GSE267718, U3_GSE267718, U4_GSE267718, 
                         U5_GSE267718, U6_GSE267718, U7_GSE267718, U8_GSE267718),
                add.cell.ids = c(1:21),
                #add.cell.ids = NULL,
                merge.data = TRUE,
                project = "BLCA_M")

table(BLCA_M$GSE_ID, BLCA_M$Batch)
ls()
BLCA_M
##############################
### Step 3: Calculate Quality Control (QC) Metrics

# store mitochondrial gene percentage in object meta data:
BLCA_M <- PercentageFeatureSet(BLCA_M, pattern = "^MT-", col.name = "percent.mt")

# store ribosomal gene percentage in object meta data:
BLCA_M <- PercentageFeatureSet(BLCA_M, pattern = "^RP[SL]", col.name = "percent.rb")

View(BLCA_M@meta.data)
range(BLCA_M$percent.mt)
range(BLCA_M$percent.rb)
##############################
##############################
### Step 4: Visualize Initial QC Metrics
# Visualize QC metrics as a violin plot to inspect distributions
### Initial Quality Control Metrics ###
PLT1 <- VlnPlot(BLCA_M, features = c("nFeature_RNA", "nCount_RNA", 
                                     "percent.mt", "percent.rb"), 
                ncol = 4, group.by = "Batch", raster = FALSE, 
                assay = "RNA", layer = "counts")

ggsave(filename = "QC1.png", plot = PLT1, width = 15, height = 8)
##########################
##########################
PLT2 <- VlnPlot(BLCA_M, features = c("nFeature_RNA", "nCount_RNA", 
                                     "percent.mt", "percent.rb"), 
                ncol = 4, raster = FALSE, assay = "RNA", layer = "counts")

ggsave(filename = "QC2.png", plot = PLT2, width = 40, height = 10,
       limitsize = TRUE)
##########################
Idents(BLCA_M) <- "GSE_ID"
levels(BLCA_M)

PLT3 <- VlnPlot(BLCA_M, features = c("nFeature_RNA", "nCount_RNA", 
                                     "percent.mt", "percent.rb"), 
                ncol = 4, 
                raster = FALSE, 
                assay = "RNA", layer = "counts")

ggsave(filename = "QC3.png", plot = PLT3, width = 40, height = 10,
       limitsize = TRUE)
###########################
# Scatter plots for QC metric relationships
Idents(BLCA_M) <- "orig.ident"
levels(BLCA_M)

plot1 <- FeatureScatter(BLCA_M, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt", raster=FALSE) + 
  ggtitle("RNA Count vs Mitochondrial Percentage") + theme_minimal()

plot2 <- FeatureScatter(BLCA_M, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA", raster=FALSE) + 
  ggtitle("RNA Count vs Feature Count") + theme_minimal()

plot3 <- FeatureScatter(BLCA_M, feature1 = "nCount_RNA", 
                        feature2 = "percent.rb", raster=FALSE) + 
  ggtitle("RNA Count vs Ribosomal Percentage") + theme_minimal()

PLT4 <- plot1 + plot2 + plot3

ggsave(filename = "QC4.png", plot = PLT4, width = 40, height = 10,
       limitsize = TRUE)
##############################################
##############################################
### Step 5: Filter Based on QC Metrics
seurat_QC <- function(seurat) {
  seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
  seurat <- PercentageFeatureSet(seurat, pattern = "^RP[SL]", col.name = "percent.rb")
  seurat <- subset(seurat, 
                   subset = (subset = (nFeature_RNA < 4000 &
                                         nCount_RNA < 20000 &
                                         percent.mt < 15 & 
                                         percent.rb < 30)))
   return(seurat)
}
#####################################################
#####################################################
## Perform Quality control ..............
P1_X10Genomics <- seurat_QC(P1_X10Genomics)
P2_GSE130001 <- seurat_QC(P2_GSE130001)
P3_GSE130001 <- seurat_QC(P3_GSE130001)
P4_GSE267718 <- seurat_QC(P4_GSE267718)
P5_GSE267718 <- seurat_QC(P5_GSE267718)
P6_GSE267718 <- seurat_QC(P6_GSE267718)
P7_GSE267718 <- seurat_QC(P7_GSE267718)
P8_GSE267718 <- seurat_QC(P8_GSE267718)
P9_GSE267718 <- seurat_QC(P9_GSE267718)
P10_GSE267718 <- seurat_QC(P10_GSE267718)
######################################
######################################
N1_GSE129845 <- seurat_QC(N1_GSE129845)
N2_GSE129845 <- seurat_QC(N2_GSE129845)
N3_GSE129845 <- seurat_QC(N3_GSE129845)
#######################################
#######################################
U1_GSE267718 <- seurat_QC(U1_GSE267718)
U2_GSE267718 <- seurat_QC(U2_GSE267718)
U3_GSE267718 <- seurat_QC(U3_GSE267718)
U4_GSE267718 <- seurat_QC(U4_GSE267718)
U5_GSE267718 <- seurat_QC(U5_GSE267718)
U6_GSE267718 <- seurat_QC(U6_GSE267718)
U7_GSE267718 <- seurat_QC(U7_GSE267718)
U8_GSE267718 <- seurat_QC(U8_GSE267718)
######################################## 
# Merge Seurat Objects ***************
BLCA_M <- merge(x = P1_X10Genomics, 
                y = list(P2_GSE130001, P3_GSE130001, P4_GSE267718, 
                         P5_GSE267718, P6_GSE267718, P7_GSE267718, 
                         P8_GSE267718, P9_GSE267718, P10_GSE267718,
                         N1_GSE129845, N2_GSE129845, N3_GSE129845, 
                         U1_GSE267718, U2_GSE267718, U3_GSE267718, U4_GSE267718, 
                         U5_GSE267718, U6_GSE267718, U7_GSE267718, U8_GSE267718),
                add.cell.ids = c("P1", "P2", "P3", "P4", "P5", 
                                 "P6", "P7", "P8", "P9", "P10",
                                 "N1", "N2", "N3", 
                                 "U1", "U2", "U3", "U4", 
                                 "U5", "U6", "U7", "U8"),
                merge.data = TRUE,
                project = "BLCA_M")

table(BLCA_M$GSE_ID, BLCA_M$Batch)
ls()
BLCA_M
##############################
#############################
Idents(BLCA_M)
levels(BLCA_M)

### Quality Control Metrics after Filtering ###
PLT1 <- VlnPlot(BLCA_M, features = c("nFeature_RNA", "nCount_RNA", 
                                     "percent.mt", "percent.rb"), 
                ncol = 4, group.by = "Batch", raster = FALSE, 
                assay = "RNA", layer = "counts")

ggsave(filename = "QC5.png", plot = PLT1, width = 15, height = 8)
##########################
##########################
PLT2 <- VlnPlot(BLCA_M, features = c("nFeature_RNA", "nCount_RNA", 
                                     "percent.mt", "percent.rb"), 
                ncol = 4, raster = FALSE, assay = "RNA", layer = "counts")

ggsave(filename = "QC6.png", plot = PLT2, width = 40, height = 10,
       limitsize = TRUE)
##########################
Idents(BLCA_M) <- "GSE_ID"
levels(BLCA_M)

PLT3 <- VlnPlot(BLCA_M, features = c("nFeature_RNA", "nCount_RNA", 
                                     "percent.mt", "percent.rb"), 
                ncol = 4, raster = FALSE, 
                assay = "RNA", layer = "counts")

ggsave(filename = "QC7.png", plot = PLT3, width = 40, height = 10,
       limitsize = TRUE)
##########################
##########################
# Scatter plots for QC metric relationships
Idents(BLCA_M) <- "orig.ident"
levels(BLCA_M)

plot1 <- FeatureScatter(BLCA_M, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt", raster=FALSE) + 
  ggtitle("RNA Count vs Mitochondrial Percentage") + theme_minimal()

plot2 <- FeatureScatter(BLCA_M, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA", raster=FALSE) + 
  ggtitle("RNA Count vs Feature Count") + theme_minimal()

plot3 <- FeatureScatter(BLCA_M, feature1 = "nCount_RNA", 
                        feature2 = "percent.rb", raster=FALSE) + 
  ggtitle("RNA Count vs Ribosomal Percentage") + theme_minimal()

PLT4 <- plot1 + plot2 + plot3

ggsave(filename = "QC8.png", plot = PLT4, width = 40, height = 10,
       limitsize = TRUE)
##########################
##########################
view(BLCA_M@meta.data)
class(BLCA_M[["RNA"]])
#########################
#########################
# Save our Seurat object 
saveRDS(BLCA_M, file = "BLCA_M1.RDS")

rm(list = ls()); gc()

# Read our Seurat object 
BLCA_M1 <- readRDS("BLCA_M1.RDS")
#########################################
#########################################