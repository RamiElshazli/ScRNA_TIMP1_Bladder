################################
## Human MSigDB Collections GSEA

## https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html
## https://igordot.github.io/msigdbr/articles/msigdbr-intro.html
################################
## Load libraries ##
BiocManager::install("msigdb")
install.packages("msigdbdf", repos = "https://igordot.r-universe.dev")

library(Seurat)
library(VennDiagram)
library(tidyverse)
library(ggpubr)
library(plotly)
library(UpSetR)
library(cowplot)
library(SingleCellExperiment)
library(patchwork)
library(pheatmap)
library(enrichR)
library(RColorBrewer)
library(reticulate)
library(msigdbr)
library(msigdbdf)

BLCA_M10 <- readRDS("BLCA_M10.RDS")

features <- c("TIMP1", "CD63", "MMP9", 
              "MMP14", "IL6", "FN1")
########################################
## https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
# Module Scoring and Comparison:
# 1. Get the "HALLMARK_ANGIOGENESIS" gene set 
## HALLMARK_ANGIOGENESIS .....

h_gene_sets <- msigdbr(category = c("H")) # collection
angio_genes <- h_gene_sets %>%
  filter(gs_name == "HALLMARK_ANGIOGENESIS") %>%
  pull(gene_symbol)

# 2. Calculate the module score
BLCA_M10 <- AddModuleScore(BLCA_M10, name = "angiogenesis_score", assay = "SCT", 
                           features = list(angio_genes),
                           search = TRUE)
rami1 <- FeaturePlot(BLCA_M10, cols = c("#DDE5B2", "#943bcc"),
                     features = "angiogenesis_score1") 
rami7 <- (rami1) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC155.png", plot = rami7, width = 10, height = 8)
##############################
##############################
rami1 <- FeaturePlot(BLCA_M10, cols = c("#DDE5B2", "#943bcc"),
                     label = TRUE, repel = TRUE,
                     features = "angiogenesis_score1") 
rami7 <- (rami1) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC156.png", plot = rami7, width = 10, height = 8)
####################################################################
####################################################################
# 3. Compare to clusters using VlnPlot and DimPlot
DimPlot(BLCA_M10, group.by = "cell_type") + 
  VlnPlot(BLCA_M10, features = "angiogenesis_score1",
          group.by = "cell_type")

rami1 <- FeaturePlot(BLCA_M10, cols = c("#DDE5B2", "#943bcc"),
                     label = TRUE, repel = TRUE,
                     features = "angiogenesis_score1") 
rami1 <- (rami1) & NoLegend() & labs(x = "UMAP_1", y = "UMAP_2")
rami2 <- VlnPlot(BLCA_M10, features = "angiogenesis_score1",
                 group.by = "cell_type") & NoLegend()

rami3 <- rami1 | rami2
ggsave(filename = "QC157.png", plot = rami3, width = 20, height = 10)
#####################################################################
#####################################################################
rami1 <- FeaturePlot(BLCA_M10, features = "TIMP1", label = TRUE, 
                     repel = TRUE, cols = c("lightgrey", "#943bcc")) 
rami2 <- FeaturePlot(BLCA_M10, features = "angiogenesis_score1",
                     label = TRUE, repel = TRUE, cols = c("#DDE5B2", "#943bcc")) 
rami3 <- (rami1 | rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC158.png", plot = rami3, width = 20, height = 10)
################################
################################
rami1 <- FeaturePlot(BLCA_M10, features = "angiogenesis_score1", 
                     label = TRUE, repel = TRUE, cols = c("#DDE5B2", "#943bcc"))
rami2 <- VlnPlot(BLCA_M10, features = "angiogenesis_score1",
                 group.by = "cell_type") & NoLegend()

rami3 <- (rami1 | rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC159.png", plot = rami3, width = 20, height = 10)
############################################
############################################
# 4. Scatter plot comparing G2M module score with the G2M score found by Cell.CycleScoring
rami4 <- FeatureScatter(BLCA_M10, feature1 = "angiogenesis_score1", 
                        feature2 = "G2M.Score",
                        group.by = "Phase")

ggsave(filename = "QC160.png", plot = rami4, width = 12, height = 10)
ggsave(filename = "QC160.jpg", plot = rami4, width = 12, height = 10)
######################################################################
#####################################################################
# Module Scoring and Comparison:
# 1. Get the "HALLMARK_APOPTOSIS" gene set
h_gene_sets <- msigdbr(category = c("H")) # collection
apopt_genes <- h_gene_sets %>%
  filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol)
#########################################
#########################################
## HALLMARK_APOPTOSIS // HALLMARK_HYPOXIA
# 2. Calculate the Hypoxia module score
BLCA_M10 <- AddModuleScore(BLCA_M10, name = "apoptosis_score", assay = "SCT", 
                           features = list(apopt_genes),
                           search = TRUE)
FeaturePlot(BLCA_M10, features = "apoptosis_score1") 
#########################################
#########################################
# 3. Compare to clusters using VlnPlot and DimPlot
DimPlot(BLCA_M10, group.by = "cell_type") + 
  VlnPlot(BLCA_M10, features = "apoptosis_score1",
          group.by = "cell_type")
#########################################
#########################################
rami1 <- FeaturePlot(BLCA_M10, cols = c("#DDE5B2", "#943bcc"),
                     label = TRUE, repel = TRUE,
                     features = "apoptosis_score1") 
rami1 <- (rami1) & NoLegend() & labs(x = "UMAP_1", y = "UMAP_2")

rami2 <- VlnPlot(BLCA_M10, features = "apoptosis_score1",
                 group.by = "cell_type") & NoLegend()

rami3 <- rami1 | rami2
ggsave(filename = "QC161.png", plot = rami3, width = 20, height = 10)
##########################################
##########################################
rami1 <- FeaturePlot(BLCA_M10, features = "TIMP1", label = TRUE, 
                     repel = TRUE, cols = c("lightgrey", "#943bcc")) 
rami2 <- FeaturePlot(BLCA_M10, features = "apoptosis_score1",
                     label = TRUE, repel = TRUE, cols = c("#DDE5B2", "#943bcc")) 
rami3 <- (rami1 | rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC162.png", plot = rami3, width = 20, height = 10)
##########################################
##########################################
rami1 <- FeaturePlot(BLCA_M10, features = "apoptosis_score1", 
                     label = TRUE, repel = TRUE, cols = c("#DDE5B2", "#943bcc"))
rami2 <- VlnPlot(BLCA_M10, features = "apoptosis_score1",
                 group.by = "cell_type") & NoLegend()

rami3 <- (rami1 | rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC163.png", plot = rami3, width = 20, height = 10)
###############################
###############################
# 4. Scatter plot comparing G2M module score with the G2M score found by Cell.CycleScoring
rami4 <- FeatureScatter(BLCA_M10, feature1 = "apoptosis_score1", 
                        feature2 = "G2M.Score",
                        group.by = "Phase")

ggsave(filename = "QC164.png", plot = rami4, width = 12, height = 10)
ggsave(filename = "QC164.jpg", plot = rami4, width = 12, height = 10)
######################################################################
######################################################################
# Module Scoring and Comparison:
# 1. Get the "HALLMARK_INFLAMMATORY_RESPONSE" gene set
h_gene_sets <- msigdbr(category = c("H")) # collection
inflammatory_genes <- h_gene_sets %>%
  filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>%
  pull(gene_symbol)

## HALLMARK_APOPTOSIS // HALLMARK_HYPOXIA
# 2. Calculate the Hypoxia module score
BLCA_M10 <- AddModuleScore(BLCA_M10, name = "inflammatory_score", assay = "SCT", 
                           features = list(inflammatory_genes),
                           search = TRUE)
FeaturePlot(BLCA_M10, features = "inflammatory_score1") 
###########################################
###########################################
# 3. Compare to clusters using VlnPlot and DimPlot
DimPlot(BLCA_M10, group.by = "cell_type") + 
  VlnPlot(BLCA_M10, features = "inflammatory_score1",
          group.by = "cell_type")
###########################################
###########################################
rami1 <- FeaturePlot(BLCA_M10, cols = c("#DDE5B2", "#943bcc"),
                     label = TRUE, repel = TRUE,
                     features = "inflammatory_score1") 
rami1 <- (rami1) & NoLegend() & labs(x = "UMAP_1", y = "UMAP_2")

rami2 <- VlnPlot(BLCA_M10, features = "inflammatory_score1",
                 group.by = "cell_type") & NoLegend()

rami3 <- rami1 | rami2
ggsave(filename = "QC165.png", plot = rami3, width = 20, height = 10)
###########################################
###########################################
rami1 <- FeaturePlot(BLCA_M10, features = "TIMP1", label = TRUE, 
                     repel = TRUE, cols = c("lightgrey", "#943bcc")) 
rami2 <- FeaturePlot(BLCA_M10, features = "inflammatory_score1",
                     label = TRUE, repel = TRUE, cols = c("#DDE5B2", "#943bcc")) 
rami3 <- (rami1 | rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC166.png", plot = rami3, width = 20, height = 10)
############################################
############################################
rami1 <- FeaturePlot(BLCA_M10, features = "inflammatory_score1", 
                     label = TRUE, repel = TRUE, cols = c("#DDE5B2", "#943bcc"))
rami2 <- VlnPlot(BLCA_M10, features = "inflammatory_score1",
                 group.by = "cell_type") & NoLegend()

rami3 <- (rami1 | rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC167.png", plot = rami3, width = 18, height = 10)
##############################################
##############################################
# 4. Scatter plot comparing G2M module score with the G2M score found by Cell.CycleScoring
rami4 <- FeatureScatter(BLCA_M10, feature1 = "inflammatory_score1", 
                        feature2 = "G2M.Score",
                        group.by = "Phase")

ggsave(filename = "QC168.png", plot = rami4, width = 12, height = 10)
ggsave(filename = "QC168.jpg", plot = rami4, width = 12, height = 10)
#####################################################################
VlnPlot(BLCA_M10, features = c("angiogenesis_score1", "TIMP1"),
        group.by = "cell_type")

ggsave(filename = "QC169.png", width = 20, height = 8)
dev.off()
##############################
##############################
VlnPlot(BLCA_M10, features = c("apoptosis_score1", "TIMP1"),
        group.by = "cell_type")

ggsave(filename = "QC170.png", width = 20, height = 8)
dev.off()
##############################
##############################
VlnPlot(BLCA_M10, features = c("inflammatory_score1", "TIMP1"),
        group.by = "cell_type")

ggsave(filename = "QC171.png", width = 20, height = 8)
dev.off()
########################################
########################################
# cp_gene_sets <- msigdbr(category = c("C8")) ## H C1-C8
# emt_genes <- cp_gene_sets %>%
#    filter(gs_name == "TAVAZOIE_METASTASIS") %>%
#    pull(gene_symbol)
#########################################
saveRDS(BLCA_M10, file = "BLCA_M11.RDS")
BLCA_M11 <- readRDS("BLCA_M11.RDS")
#########################################